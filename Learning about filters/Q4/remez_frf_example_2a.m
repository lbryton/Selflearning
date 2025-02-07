% Remez filter design example
% Design a low pass FIR filter with the following characteristics
% filter length is a multiple of 30
% sample rate 4.5 MHz
% passband edge 50 kHz
% passband ripple 0.1 dB
% stopband edge 100 kHz
% rejection 60 dB
% the weights are adjusted to optimize the inband vs. stopband performance
h1=remez(279-1,[0, 1, 2, 50]/50, [1 1 0 0], [1 10]);

figure(1)
plot((-0.5:1/4096:0.5-1/4096),20*log10(abs(fftshift(fft(h1,4096)))))
axis([-0.5 0.5 -100 5])
title('Equal Ripple Remez Filter Design Example')
grid on
zoom on
% now give the filter a 1/f rolloff in the stopband
% myfrf is a modified version of the remezfrf function to provide the 1/f rolloff
%   in myfrf, removed linse start with "%---"
%   and new lines end with "%+++ ..."
%   to make the changes visible
% the weights are adjusted to optimize the inband vs. stopband performance
h2=remez(279-1,[0, 1, 2, 50]/50, {'myfrf', [1 1 0 0]}, [1 10]);
figure(2)
plot((-0.5:1/4096:0.5-1/4096),20*log10(abs(fftshift(fft(h2,4096)))))
axis([-0.5 0.5 -100 5])
title('1/f Stopband Rolloff Remez Filter Design Example')
grid on
zoom on

figure(3)
subplot(2,1,1)
plot((-0.5:1/4096:0.5-1/4096),20*log10(abs(fftshift(fft(h1,4096)))),'b','linewidth',1.5)
axis([-0.25 0.25 -100 5])
title('Equal Ripple Remez Filter Design Example')
xlabel('Normalized Frequency')
ylabel('Log Magnitude (dB)')
grid on
axes('position',[0.6 0.75 0.2 0.15])
plot((-0.5:1/4096:0.5-1/4096),20*log10(abs(fftshift(fft(h1,4096)))),'b','linewidth',1.5)
axis([-0.02 0.02 -.2 .2])
title('Zoom to Pass Band')
ylabel('dB')
grid on

subplot(2,1,2)
plot((-0.5:1/4096:0.5-1/4096),20*log10(abs(fftshift(fft(h2,4096)))),'b','linewidth',1.5)
axis([-0.25 0.25 -100 5])
title('1/f Stopband Rolloff Remez Filter Design Example')
xlabel('Normalized Frequency')
ylabel('Log Magnitude (dB)')
grid on
axes('position',[0.6 0.27 0.2 0.15])
plot((-0.5:1/4096:0.5-1/4096),20*log10(abs(fftshift(fft(h2,4096)))),'b','linewidth',1.5)
axis([-0.02 0.02 -.2 .2])
title('Zoom to Pass Band')
ylabel('dB')
grid on

figure(4)
subplot(2,1,1)
plot(0:278,h1,'b','linewidth',2)
grid on
axis([-5 285 -0.01 0.035])
title('Impulse Response, Remez Algorithm','fontsize',12)
xlabel('Time Index','fontsize',12)
ylabel('Amplitude','fontsize',12)

axes('position',[0.17 0.7 0.20 0.15])
stem(0:10, h1(1:11),'b','linewidth',2)
grid on 
axis([-1 10 -0.0005 0.00])
title('End Samples Details','fontsize',12)
xlabel('Time Index','fontsize',12)
ylabel('Amplitude','fontsize',12)

subplot(2,1,2)
plot(0:278,h2,'b','linewidth',2)
grid on
axis([-5 285 -0.01 0.035])
title('Impulse Response, Remez Algorithm with myfrf','fontsize',12)
xlabel('Time Index','fontsize',12)
ylabel('Amplitude','fontsize',12)

axes('position',[0.17 0.2 0.20 0.15])
stem(0:10, h2(1:11),'b','linewidth',2)
grid on 
axis([-1 10 -0.0005 0.00])
title('End Samples Details','fontsize',12)
xlabel('Time Index','fontsize',12)
ylabel('Amplitude','fontsize',12)
