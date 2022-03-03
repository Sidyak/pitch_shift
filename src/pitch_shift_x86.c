/*
fft based pitch shifting
Algorithm is based on dafx - U. Zoelzer page 279 ff, block-by-block pitch shifting apporach w/ resampling
Project is deployed on x86
*/

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#if defined(_MSC_VER)
#include <getopt.h>
#else
#include <unistd.h>
#endif

extern "C" {
    #include "wavreader.h"
    #include "wavwriter.h"
}

#include "kissfft/kiss_fft.h"

#define BUFFER_SIZE (16384)

void cleanup(void);

// set the frequency of the oscillators
int16_t gInputBuffer[BUFFER_SIZE];
int gInputBufferPointer = 0;
int16_t gOutputBuffer[BUFFER_SIZE];
int16_t inL, inR, outL, outR;
int gOutputBufferWritePointer = 0;
int gOutputBufferReadPointer = 0;
int gSampleCount = 0;

// These variables used internally in the example:
static const int gFFTSize = 9;//11;
// phase parameters
static const int Ha = (1<<gFFTSize)/4-0; /* analysis hopsize */
static const int Hs = (1<<gFFTSize)/4;   /* synthisis hopsize */
int gHopSize = Hs;
int gPeriod = gHopSize;
float tstretch = 1.0f;
float resample = 1.0f;
int pitch = 0;
int lx;

float gFFTScaleFactor = 0;
float *gWindowBuffer;

#if 0
void process_pitch_shift_background(void *);

// FFT vars
ne10_fft_cpx_float32_t* timeDomainIn;
ne10_fft_cpx_float32_t* timeDomainOut;
ne10_fft_cpx_float32_t* frequencyDomain;
ne10_fft_cfg_float32_t cfg;
#endif

kiss_fft_cfg fft_cfg;
kiss_fft_cfg ifft_cfg;

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
#if 0 // TODO
AuxiliaryTask gFFTTask;
#endif
int gFFTInputBufferPointer = 0;
int gFFTOutputBufferPointer = 0;

// instantiate the scope
//Scope scope;

// userData holds an opaque pointer to a data structure that was passed
// in from the call to initAudio().
//
// Return true on success; returning false halts the program.

bool setup(void)
{
#if 0
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
#endif

    fft_cfg = kiss_fft_alloc(1<<gFFTSize, 0, (void*)0, 0);
    ifft_cfg = kiss_fft_alloc(1<<gFFTSize, 1, (void*)0, 0);

    gFFTScaleFactor = 1.0f / (float)(1<<gFFTSize);
    gOutputBufferWritePointer = Hs;
    gOutputBufferReadPointer = 0;
#if 0
    timeDomainIn = (ne10_fft_cpx_float32_t*) NE10_MALLOC ((1<<gFFTSize) * sizeof (ne10_fft_cpx_float32_t));
    timeDomainOut = (ne10_fft_cpx_float32_t*) NE10_MALLOC ((1<<gFFTSize) * sizeof (ne10_fft_cpx_float32_t));
    frequencyDomain = (ne10_fft_cpx_float32_t*) NE10_MALLOC ((1<<gFFTSize) * sizeof (ne10_fft_cpx_float32_t));
    // https://community.vcvrack.com/t/complete-list-of-c-c-fft-libraries/9153
    cfg = ne10_fft_alloc_c2c_float32_neon ((1<<gFFTSize));

    memset(timeDomainIn, 0, (1<<gFFTSize) * sizeof (ne10_fft_cpx_float32_t));
    memset(timeDomainOut, 0, (1<<gFFTSize) * sizeof (ne10_fft_cpx_float32_t));
#endif

    memset(gOutputBuffer, 0, BUFFER_SIZE * sizeof(int16_t));
    memset(gInputBuffer, 0, BUFFER_SIZE * sizeof(int16_t));
    
#if 0        
    // Allocate buffer to mirror and modify the input
    gInputAudio = (float *)malloc(context->audioFrames * context->audioOutChannels * sizeof(float));
    if(gInputAudio == 0)
        return false;

#endif

    // Allocate phase processing buffer and init vars
    psi = (float *)malloc((1<<gFFTSize) * sizeof(float));
    if(psi == 0)
    {
        return false;
    }
    phi = (float *)malloc((1<<gFFTSize) * sizeof(float));
    if(phi == 0)
    {
        return false;
    }
    amplitude = (float *)malloc((1<<gFFTSize) * sizeof(float));
    if(amplitude == 0)
    {
        return false;
    }
    phi0 = (float *)malloc((1<<gFFTSize) * sizeof(float));
    if(phi0 == 0)
    {
        return false;
    }
    deltaPhi = (float *)malloc((1<<gFFTSize) * sizeof(float));
    if(deltaPhi == 0)
    {
        return false;
    }
    omega = (float *)malloc((1<<gFFTSize) * sizeof(float));
    if(omega == 0)
    {
        return false;
    }
    // Allocate the window buffer based on the FFT size
    gWindowBuffer = (float *)malloc((1<<gFFTSize) * sizeof(float));

    if(gWindowBuffer == NULL)
    {
        return false;
    }

    // Calculate a Hann window
    for(int n = 0; n < (1<<gFFTSize); n++)
    {
        gWindowBuffer[n] = 0.5f * (1.0f - cosf(2.0f * M_PI * (float)n / (float)((1<<gFFTSize))));
        omega[n] = 2.0f*(float)M_PI*Hs*(float)n/(float)(1<<gFFTSize);
        phi0[n] = 0.0f;
        psi[n] = 0.0f;
    }

#if 0    // Initialise auxiliary tasks
    if((gFFTTask = Bela_createAuxiliaryTask(&process_pitch_shift_background, 90, "fft-calculation")) == 0)
        return false;


    // Check if analog channels are enabled
    if(context->analogFrames == 0 || context->analogFrames > context->audioFrames) {
        printf("Error: this example needs analog enabled, with 4 or 8 channels\n");
        return false;
    }
    // Useful calculations
    if(context->analogFrames)
        gAudioFramesPerAnalogFrame = context->audioFrames / context->analogFrames;
        
    printf("bye setup\n");
#endif
    return true;

}

float princarg(float phase)
{
    return fmod(phase + M_PI, -2*M_PI) + M_PI;
}

// This function handles the FFT based pitch shifting processing
void process_pitch_shift(int16_t *inBuffer, int inWritePointer, int16_t *outBuffer, int outWritePointer)
{
    // Copy buffer into FFT input
    int pointer = (inWritePointer - (1<<gFFTSize) + BUFFER_SIZE) % BUFFER_SIZE;

#if 1 // TODO

#if 0
    // Run the FFT https://community.vcvrack.com/t/complete-list-of-c-c-fft-libraries/9153
    ne10_fft_c2c_1d_float32_neon (frequencyDomain, timeDomainIn, cfg, 0);
#else

    kiss_fft_cpx timeDomainIn[(1<<gFFTSize)];
    kiss_fft_cpx frequencyDomain[(1<<gFFTSize)];

    // interleave real imag parts
    for(int n = 0; n < (1<<gFFTSize); n++)
    {
        timeDomainIn[n].r = (int16_t)((float)inBuffer[pointer] * gWindowBuffer[n]);
        timeDomainIn[n].i = 0;

        pointer++;
        if(pointer >= BUFFER_SIZE)
        {
            pointer = 0;
        }
    }
    
    kiss_fft(fft_cfg, timeDomainIn , frequencyDomain);

#endif

    for(int n = 0; n < (1<<gFFTSize); n++)
    {
//       amplitude[n] = sqrtf(frequencyDomain[n].r * frequencyDomain[n].r + frequencyDomain[n].i * frequencyDomain[n].i);
        amplitude[n] = sqrtf(frequencyDomain[n].r * frequencyDomain[n].r + frequencyDomain[n].i * frequencyDomain[n].i);
 //       phi[n] = atan2(frequencyDomain[n].i, frequencyDomain[n].r); //rand()/(float)RAND_MAX * 2.f* M_PI;
        phi[n] = atan2(frequencyDomain[n].i, frequencyDomain[n].r); //rand()/(float)RAND_MAX * 2.f* M_PI;
    }
#endif
    
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
    tstretch = (float)(Ha+pitch)/(float)Hs;
    resample = 1.0f/tstretch;
    
    for(int n = 0; n < (1<<gFFTSize); n++)
    {
        deltaPhi[n] = omega[n] + princarg(phi[n]-phi0[n]-omega[n]);
        phi0[n] = phi[n];
        psi[n] = princarg(psi[n]+deltaPhi[n]*tstretch); 
        //Y1_[n] = amplitude[n] .* exp(1i*(psi[n]));
#if 1// TODO
        //frequencyDomain[n].r = cosf(psi[n]) * amplitude[n];
        frequencyDomain[n].r = (int16_t)(cosf(psi[n]) * amplitude[n]);
        //frequencyDomain[n].i = sinf(psi[n]) * amplitude[n];
        frequencyDomain[n].i = (int16_t)(sinf(psi[n]) * amplitude[n]);
#endif
        //y1_[n] = fftshift(real(ifft(Y1_, N)));
    }

    //ne10_fft_c2c_1d_float32_neon (timeDomainOut, frequencyDomain, cfg, 1);
    kiss_fft(ifft_cfg, frequencyDomain, timeDomainIn);

    //y1_[n] = fftshift(real(ifft(Y1_, N))) .* w[n]';

    //if(pitch)
    {
#if 1 // TODO
        for(int n = 0; n < (1<<gFFTSize); n++)
        {
            //timeDomainOut[n].r = timeDomainOut[n].r * gWindowBuffer[n];
            timeDomainIn[n].r = (int16_t)((float)timeDomainIn[n].r * gWindowBuffer[n]);
        }
#endif
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
#if 1 // TODO

    float absMax = 0.f;
    uint32_t ueqZero = 0;
    
    // Overlap-and-add timeDomainOut into the output buffer
    pointer = outWritePointer;
    lx = (int)floor((float)(1<<gFFTSize)*resample);
    int n;
    for(n = 0; n < lx/*(1<<gFFTSize)*/; n++)
    {
#if 1
        float x = 0.0f + (float)n*((float)(1<<gFFTSize)/(float)lx);
        int ix = (int)floor(x);
        int ix1 = ix+1;
        float dx = x-(float)ix;
        float dx1 = 1.0f-dx;
#else
        int x = 0 + n*((1<<gFFTSize)/lx);
        int ix = (int)(x);
        int ix1 = ix+1;
        int dx = x-ix;
        int dx1 = 1-dx;
#endif
        //outBuffer[pointer] += (timeDomainOut[ix].r * dx1 + timeDomainOut[ix1].r * dx);
        //outBuffer[pointer] += (timeDomainOut[n].r);// * gFFTScaleFactor;

        outBuffer[pointer] += (int16_t)(timeDomainIn[ix].r * dx1 + timeDomainIn[ix1].r * dx);
        //outBuffer[pointer] += (int16_t)(timeDomainIn[n].r);// * gFFTScaleFactor;

        //if(timeDomainOut[n].i != 0)
        //    printf("timeDomainOut[n].i not zero \n");

        if(timeDomainIn[n].i != 0)
        {
            absMax = (absMax < fabs(timeDomainIn[n].i)) ? fabs(timeDomainIn[n].i) : absMax;
            ueqZero++;
        }

        if(isnan(outBuffer[pointer]))
        {
            printf("outBuffer OLA\n");
        }

        pointer++;
        if(pointer >= BUFFER_SIZE)
        {
            pointer = 0;
        }
    }

    if(ueqZero)
    {
        printf("WARNING: timeDomainIn[N].i not zero %d times (max = %f)\n", ueqZero, absMax);
    }
#endif
    printf("end\n");
}

#if 0
// Function to process the FFT in a thread at lower priority
void process_pitch_shift_background(void*)
{
    process_pitch_shift(gInputBuffer, gFFTInputBufferPointer, gOutputBuffer, gFFTOutputBufferPointer);
}
#endif

void usage(const char* name)
{
    fprintf(stderr, "%s in.wav out.wav\n", name);
}

int main(int argc, char *argv[])
{
    const char *infile, *outfile;
    FILE *out;
    void *wavIn;
    void *wavOut;
    int format, sample_rate, channels, bits_per_sample;
    uint32_t data_length;
    int input_size;
    uint8_t* input_buf;
    int16_t* convert_buf;

    if(!setup())
    {
        fprintf(stderr, "setup failed\n");
    }
    
    if (argc - optind < 2)
    {
        fprintf(stderr, "Error: not enough parameter provided\n");
        usage(argv[0]);
        return 1;
    }
    
    infile = argv[optind];
    outfile = argv[optind + 1];

    wavIn = wav_read_open(infile);
    if (!wavIn)
    {
        fprintf(stderr, "Unable to open wav file %s\n", infile);
        return 1;
    }
    if (!wav_get_header(wavIn, &format, &channels, &sample_rate, &bits_per_sample, &data_length))
    {
        fprintf(stderr, "Bad wav file %s\n", infile);
        return 1;
    }
    if (format != 1)
    {
        fprintf(stderr, "Unsupported WAV format %d\n", format);
        return 1;
    }

    wavOut = wav_write_open(outfile, sample_rate, bits_per_sample, channels);

    if (!wavOut)
    {
        fprintf(stderr, "Unable to open wav file for writing %s\n", infile);
        return 1;
    }

    input_size = data_length;//channels*2*frameLength;//info.frameLength;
    input_buf = (uint8_t*) malloc(input_size);
    convert_buf = (int16_t*) malloc(input_size);

    if (input_buf == NULL || convert_buf == NULL)
    {
        fprintf(stderr, "Unable to allocate memory for buffer\n");
        return 1;
    }

    int read = wav_read_data(wavIn, input_buf, input_size);

    fprintf(stderr, "data_length = %d\tread = %d\tinput_size = %d \n", data_length, read, input_size);
    fprintf(stderr, "sample_rate = %d\tbits_per_sample = %d\tchannels = %d \n", sample_rate, bits_per_sample, channels);

    int numSamples = read/2;
    for(unsigned int n = 0; n < numSamples; n++)
    {
        const uint8_t* in = &input_buf[2*n];
        convert_buf[n] = in[0] | (in[1] << 8);
    }

    pitch = 0; //floor(map(analogRead(context, n/gAudioFramesPerAnalogFrame, gAnalogIn), 0, 1, -200, 200));

    // iterate over the audio frames and create three oscillators, seperated in phase by PI/2
    for(unsigned int n = 0; n < numSamples; n+=channels)
    {
        // Read audio inputs
        if(channels == 1)
        {
            inL = (float)convert_buf[n];///(1<<15);
            inR = inL;
        }
        if(channels == 2)
        {
            // interleaved left right channel
            inL = (float)convert_buf[n];///(1<<15);
            inR = (float)convert_buf[n+1];///(1<<15);
        }
        else
        {
            fprintf(stderr, "channel = %d\n", channels);
            return -1;
        }

#if 1 // apply pitch shifting if defined. otherwise it's bypassed
        gInputBuffer[gInputBufferPointer] = (inR+inL)/2;

        outL = gOutputBuffer[gOutputBufferReadPointer];
        outR = outL;
        
        // Clear the output sample in the buffer so it is ready for the next overlap-add
        gOutputBuffer[gOutputBufferReadPointer] = 0;
        
        gOutputBufferReadPointer++;
        if(gOutputBufferReadPointer >= (BUFFER_SIZE))
        {
            gOutputBufferReadPointer = 0;
        }   
        gOutputBufferWritePointer++;
        if(gOutputBufferWritePointer >= (BUFFER_SIZE))
        {
            gOutputBufferWritePointer = 0;
        }
        gInputBufferPointer++;
        if(gInputBufferPointer >= (BUFFER_SIZE))
        {
            gInputBufferPointer = 0;
        }    
        gSampleCount++;
        if(gSampleCount >= Hs)
        {
#if 1
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
    	//outL *= (1<<15);
        //outR *= (1<<15);
        int16_t oL = (int16_t)outL;
        int16_t oR = (int16_t)outR;
        wav_write_data(wavOut, (unsigned char*)&oL, 2);
        wav_write_data(wavOut, (unsigned char*)&oR, 2);
    }    

    free(convert_buf);
    free(input_buf);
    
    cleanup();

    wav_write_close(wavOut);
    wav_read_close(wavIn);

    return 0;
}
// cleanup_render() is called once at the end, after the audio has stopped.
// Release any resources that were allocated in initialise_render().

void cleanup()
{
#if 0 // TODO
    NE10_FREE(timeDomainIn);
    NE10_FREE(timeDomainOut);
    NE10_FREE(frequencyDomain);
    NE10_FREE(cfg);

    free(gInputAudio);
#endif

    kiss_fft_free(fft_cfg);

    free(gWindowBuffer);
    free(psi);
    free(phi);
    free(amplitude);
    free(phi0);
    free(deltaPhi);
    free(omega);
}
