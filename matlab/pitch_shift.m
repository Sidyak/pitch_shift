clear
clc
close all

[y, fs] = audioread('batman.wav');
%Y = fft(y);
Y = dct(y);
subplot(211), plot(y);
L = length(Y);
subplot(212), plot( abs(Y(1:L/2)) );
