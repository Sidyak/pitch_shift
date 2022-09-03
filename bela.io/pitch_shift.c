/*

fft based pitch shifting
Algorithm is based on dafx - U. Zoelzer page 279 ff, block-by-block pitch shifting apporach w/ resampling
Project is deployed on 
 ____  _____ _        _    
| __ )| ____| |      / \   
|  _ \|  _| | |     / _ \  
| |_) | |___| |___ / ___ \ 
|____/|_____|_____/_/   \_\

The platform for ultra-low latency audio and sensor processing

http://bela.io

The Bela software is distributed under the GNU Lesser General Public License
(LGPL 3.0), available here: https://www.gnu.org/licenses/lgpl-3.0.txt
*/

#include <Bela.h>
#include <ne10/NE10.h>                    // NEON FFT library
#include <Midi.h>
#include <math.h>
#include <stdio.h>

#define BUFFER_SIZE (16384)

// set the frequency of the oscillators
float gInputBuffer[BUFFER_SIZE];
int gInputBufferPointer = 0;
float gOutputBuffer[BUFFER_SIZE];
float inL, inR, outL, outR;
int gOutputBufferWritePointer = 0;
int gOutputBufferReadPointer = 0;
int gSampleCount = 0;

// These variables used internally in the example:
int gFFTSize = 2048;
// phase parameters
int Ha = gFFTSize/4-0; /* analysis hopsize */
int Hs = gFFTSize/4;   /* synthisis hopsize */
int gHopSize = Hs;
int gPeriod = gHopSize;
float tstretch = 1.0f;
float resample = 1.0f;
int pitch = 0;
int prevPitch = 0;
int lx;

float gFFTScaleFactor = 0;
float *gInputAudio = NULL;
float *gWindowBuffer;

void process_pitch_shift_background(void *);

// FFT vars
ne10_fft_cpx_float32_t* timeDomainIn;
ne10_fft_cpx_float32_t* timeDomainOut;
ne10_fft_cpx_float32_t* frequencyDomain;
ne10_fft_cfg_float32_t cfg;

// phase processing vars
float *deltaPhi;
float *phi0;
float *phi;
float *psi;
float *omega;
float *amplitude;

// Set the analog channels to read from
const int gAnalogIn = 0;
int gAudioFramesPerAnalogFrame = 0;

int gReadPtr = 0;        // Position of last read sample from file
AuxiliaryTask gFFTTask;
int gFFTInputBufferPointer = 0;
int gFFTOutputBufferPointer = 0;

// instantiate the scope
//Scope scope;

// userData holds an opaque pointer to a data structure that was passed
// in from the call to initAudio().
//
// Return true on success; returning false halts the program.
bool setup(BelaContext* context, void* userData)
{
    printf("go setup\n");
    printf("context->audioFrames = %d\n", context->audioFrames);
    printf("context->audioSampleRate = %f\n", context->audioSampleRate);
    printf("context->audioInChannels = %d\n", context->audioInChannels);
    printf("context->audioOutChannels = %d\n", context->audioOutChannels);
    // Check that we have the same number of inputs and outputs.
    if(context->audioInChannels != context->audioOutChannels ||
            context->analogInChannels != context-> analogOutChannels){
        printf("Error: for this project, you need the same number of input and output channels.\n");
        return false;
    }

    gFFTScaleFactor = 1.0f / (float)gFFTSize;
    gOutputBufferWritePointer = Hs;
    gOutputBufferReadPointer = 0;

    timeDomainIn = (ne10_fft_cpx_float32_t*) NE10_MALLOC (gFFTSize * sizeof (ne10_fft_cpx_float32_t));
    timeDomainOut = (ne10_fft_cpx_float32_t*) NE10_MALLOC (gFFTSize * sizeof (ne10_fft_cpx_float32_t));
    frequencyDomain = (ne10_fft_cpx_float32_t*) NE10_MALLOC (gFFTSize * sizeof (ne10_fft_cpx_float32_t));
    cfg = ne10_fft_alloc_c2c_float32_neon (gFFTSize);
    
    memset(timeDomainIn, 0, gFFTSize * sizeof (ne10_fft_cpx_float32_t));
    memset(timeDomainOut, 0, gFFTSize * sizeof (ne10_fft_cpx_float32_t));
    memset(gOutputBuffer, 0, BUFFER_SIZE * sizeof(float));
    memset(gInputBuffer, 0, BUFFER_SIZE * sizeof(float));

    // Allocate phase processing buffer and init vars
    psi = (float *)malloc(gFFTSize * sizeof(float));
    if(psi == 0)
        return false;
    phi = (float *)malloc(gFFTSize * sizeof(float));
    if(phi == 0)
        return false;
    amplitude = (float *)malloc(gFFTSize * sizeof(float));
    if(amplitude == 0)
        return false;
    phi0 = (float *)malloc(gFFTSize * sizeof(float));
    if(phi0 == 0)
        return false;
    deltaPhi = (float *)malloc(gFFTSize * sizeof(float));
    if(deltaPhi == 0)
        return false;
    omega = (float *)malloc(gFFTSize * sizeof(float));
    if(omega == 0)
        return false;
        
    // Allocate buffer to mirror and modify the input
    gInputAudio = (float *)malloc(context->audioFrames * context->audioOutChannels * sizeof(float));
    if(gInputAudio == 0)
        return false;

    // Allocate the window buffer based on the FFT size
    gWindowBuffer = (float *)malloc(gFFTSize * sizeof(float));
    if(gWindowBuffer == 0)
        return false;

    // Calculate a Hann window
    for(int n = 0; n < gFFTSize; n++) {
        gWindowBuffer[n] = 0.5f * (1.0f - cosf(2.0f * M_PI * (float)n / (float)(gFFTSize)));
        omega[n] = 2.0f*(float)M_PI*Hs*(float)n/(float)gFFTSize;
        phi0[n] = 0.0f;
        psi[n] = 0.0f;
    }

    // Initialise auxiliary tasks
    if((gFFTTask = Bela_createAuxiliaryTask(&process_pitch_shift_background, 90, "fft-calculation")) == 0)
        return false;


    // Check if analog channels are enabled
    if(context->analogFrames == 0 || context->analogFrames > context->audioFrames) {
        rt_printf("Error: this example needs analog enabled, with 4 or 8 channels\n");
        return false;
    }
    // Useful calculations
    if(context->analogFrames)
        gAudioFramesPerAnalogFrame = context->audioFrames / context->analogFrames;
        
    printf("bye setup\n");
    return true;
}

float princarg(float phase){
    return fmod(phase + M_PI, -2*M_PI) + M_PI;
}

// This function handles the FFT based pitch shifting processing
void process_pitch_shift(float *inBuffer, int inWritePointer, float *outBuffer, int outWritePointer)
{
    // Copy buffer into FFT input
    int pointer = (inWritePointer - gFFTSize + BUFFER_SIZE) % BUFFER_SIZE;
    for(int n = 0; n < gFFTSize; n++) {
        timeDomainIn[n].r = (ne10_float32_t) inBuffer[pointer] * gWindowBuffer[n];
        timeDomainIn[n].i = 0.0f;

        pointer++;
        if(pointer >= BUFFER_SIZE)
            pointer = 0;
    }

    // Run the FFT
    ne10_fft_c2c_1d_float32_neon (frequencyDomain, timeDomainIn, cfg, 0);

    for(int n = 0; n < gFFTSize; n++) {
        amplitude[n] = sqrtf(frequencyDomain[n].r * frequencyDomain[n].r + frequencyDomain[n].i * frequencyDomain[n].i);
        phi[n] = atan2(frequencyDomain[n].i, frequencyDomain[n].r); //rand()/(float)RAND_MAX * 2.f* M_PI;
    }
    
/*
  Algorithm is based on dafx - U. Zoelzer page 279 ff, block-by-block pitch shifting apporach w/ resampling
  phi = atan2(imag(Y1), real(Y1)); % arg(Y1); %
  a = sqrt(real(Y1).*real(Y1) + imag(Y1).*imag(Y1)); % abs(Y1);
  deltaPhi = omega + princarg(phi-phi0-omega);
  phi0 = phi;
  psi = princarg(psi+deltaPhi*tstretch); 
  Y1_ = a .* exp(1i*(psi));
  % done phase processing
  if(tstretch == 1)
    y1_ = fftshift(real(ifft(Y1_, N)));
  else
    y1_ = fftshift(real(ifft(Y1_, N))) .* w';
  end  
*/

    // In case the pitch would change during processing phi0 and psi have to be recalculated. Otherwise artifacts arise.
    if(prevPitch != pitch)
    {
        //rt_printf("psi reset cause pitch_offset changed from %d to %d\n", prevPitch, pitch);
        prevPitch = pitch;
        memset(phi0, 0, sizeof(phi0));
        memset(psi, 0, sizeof(psi));
    }
    
    tstretch = (float)(Ha+pitch)/(float)Hs;
    resample = 1.0f/tstretch;
    for(int n = 0; n < gFFTSize; n++) {
        deltaPhi[n] = omega[n] + princarg(phi[n]-phi0[n]-omega[n]);
          phi0[n] = phi[n];
          psi[n] = princarg(psi[n]+deltaPhi[n]*tstretch); 
          //Y1_[n] = amplitude[n] .* exp(1i*(psi[n]));
          frequencyDomain[n].r = cosf(psi[n]) * amplitude[n];
          frequencyDomain[n].i = sinf(psi[n]) * amplitude[n];
          //y1_[n] = fftshift(real(ifft(Y1_, N)));
    }
    ne10_fft_c2c_1d_float32_neon (timeDomainOut, frequencyDomain, cfg, 1);
    //y1_[n] = fftshift(real(ifft(Y1_, N))) .* w[n]';
    if(pitch){
        for(int n = 0; n < gFFTSize; n++) {
            timeDomainOut[n].r = timeDomainOut[n].r * gWindowBuffer[n];
        }
    }

/*
% for linear interpolation of a grain of length N
% resampling factor Hs/Ha
lx = floor(N*resample);
x = 1+(0:lx-1)'*N/lx;
ix = floor(x) ;
ix1 = ix+1;
dx = x-ix;
dx1 = 1-dx;
%----- interpolation (resampling)
grain2 = [grain'; 0];
grain3 = grain2(ix) .* dx1 + grain2(ix1) .* dx;
*/
    // Overlap-and-add timeDomainOut into the output buffer
    pointer = outWritePointer;
    lx = (int)floor((float)gFFTSize*resample);
    int n;
    for(n = 0; n < lx/*gFFTSize*/; n++) {
        float x = 0.0f + (float)n*((float)gFFTSize/(float)lx);
        int ix = (int)floor(x);
        int ix1 = ix+1;
        float dx = x-(float)ix;
        float dx1 = 1.0f-dx;
        outBuffer[pointer] += (timeDomainOut[ix].r * dx1 + timeDomainOut[ix1].r * dx);
        //outBuffer[pointer] += (timeDomainOut[n].r);// * gFFTScaleFactor;
        //if(timeDomainOut[n].i != 0)
        //    rt_printf("timeDomainOut[n].i not zero \n");        
        if(isnan(outBuffer[pointer]))
            rt_printf("outBuffer OLA\n");
        pointer++;
        if(pointer >= BUFFER_SIZE)
            pointer = 0;
    }
}

// Function to process the FFT in a thread at lower priority
void process_pitch_shift_background(void*) {
    process_pitch_shift(gInputBuffer, gFFTInputBufferPointer, gOutputBuffer, gFFTOutputBufferPointer);
}

void render(BelaContext *context, void *userData)
{
    // iterate over the audio frames and create three oscillators, seperated in phase by PI/2
    for(unsigned int n = 0; n < context->audioFrames; n++) {
        if(gAudioFramesPerAnalogFrame && !(n % gAudioFramesPerAnalogFrame)) {
            // read analog inputs and update pitch value
            pitch = (int)floor(map(analogRead(context, n/gAudioFramesPerAnalogFrame, gAnalogIn), 0, 1, -200, 200));
            //pitch = 0;
        }

        // Read audio inputs
        inL = audioRead(context,n,0);
        inR = audioRead(context,n,1);
#if 1 // apply pitch shifting if defined. otherwise it's bypassed
        gInputBuffer[gInputBufferPointer] = (inR+inL) * 0.5f;
            
        outL = gOutputBuffer[gOutputBufferReadPointer];
        outR = outL;
        
        // Clear the output sample in the buffer so it is ready for the next overlap-add
        gOutputBuffer[gOutputBufferReadPointer] = 0;
        
        gOutputBufferReadPointer++;
        if(gOutputBufferReadPointer >= (BUFFER_SIZE))
            gOutputBufferReadPointer = 0;

        gOutputBufferWritePointer++;
        if(gOutputBufferWritePointer >= (BUFFER_SIZE))
            gOutputBufferWritePointer = 0;

        gInputBufferPointer++;
        if(gInputBufferPointer >= (BUFFER_SIZE))
            gInputBufferPointer = 0;

        gSampleCount++;
        if(gSampleCount >= Hs) {
#if 0
            /* do not use scheduling */
            process_pitch_shift(gInputBuffer, gInputBufferPointer, gOutputBuffer, gOutputBufferWritePointer);
#else
            gFFTInputBufferPointer = gInputBufferPointer;
            gFFTOutputBufferPointer = gOutputBufferWritePointer;
            Bela_scheduleAuxiliaryTask(gFFTTask);
#endif
            gSampleCount = 0;
        }    
        
#else
        outL = inL;
        outR = inR;
#endif        
        audioWrite(context, n, 0, outL);
        audioWrite(context, n, 1, outR);
    }
}

// cleanup_render() is called once at the end, after the audio has stopped.
// Release any resources that were allocated in initialise_render().

void cleanup(BelaContext* context, void* userData)
{
    NE10_FREE(timeDomainIn);
    NE10_FREE(timeDomainOut);
    NE10_FREE(frequencyDomain);
    NE10_FREE(cfg);
    free(gInputAudio);
    free(gWindowBuffer);
    free(psi);
    free(phi);
    free(amplitude);
    free(phi0);
    free(deltaPhi);
    free(omega);
}
