clc;
clf reset;
% xtickformat('%.0f');

sf = 80e3;
cycles = 20;
num_samples = 320;
samples_5k = cos_wav_call(5e3, sf, 0, num_samples);
% samples_10k = cos_wav_call(10e3, sf, 0, num_samples);
% const_t = 0.5 * ones(1, num_samples);
x = samples_5k;
% combined_signal = samples_10k + samples_5k + const_t;
% disp(combined_signal(1:100))


figure(5)
subplot(4,1,1);
plot(0:100,real(x(1:101)),'b','linewidth',2)
hold on;
plot(0:100,imag(x(1:101)),'r','linewidth',2)
hold off;
grid on
axis([0 100 -1.2 1.2])
title('Complex input Sinusoid at f_C= 5, and f_S= 80')

subplot(4,1,2);
ww=kaiser(320,8)';
ww=ww/sum(ww);
plot([-0.5:1/1000:0.5-1/1000]*80,fftshift(abs(fft(x.*ww,1000))),'b','linewidth',2)
grid on
axis([-40 40 0 1.2])
title('Magnitude of windowed FFT of Complex input Sinusoid at f_C= 5, and f_S= 80')

subplot(4,1,3);
plot([-0.5:1/1000:0.5-1/1000]*80,fftshift(real(fft(x.*ww,1000))),'b','linewidth',2)
grid on
axis([0 10 -1.2 1.2])
title('Real Part of windowed FFT of Complex input Sinusoid at f_C= 5, and f_S= 80')

subplot(4,1,4);
plot([-0.5:1/1000:0.5-1/1000]*80,fftshift(imag(fft(x.*ww,1000))),'b','linewidth',2)
grid on
axis([0 10 -1.2 1.2])
title('Imag Part of windowed FFT of Complex input Sinusoid at f_C= 5, and f_S= 80')


figure(6)
subplot(4,1,1);
plot(0:100,real(x(1:101)),'b','linewidth',2)
hold on;
plot(0:100,imag(x(1:101)),'r','linewidth',2)
hold off;
grid on
axis([0 100 -1.2 1.2])
title('Complex input Sinusoid at f_C= 5, and f_S= 80')

subplot(4,1,2);
ww=kaiser(320,8)';
ww=ww/sum(ww);
plot([-0.5:1/1000:0.5-1/1000]*80,fftshift(abs(fft(real(x).*ww,1000))),'b','linewidth',2)
grid on
axis([-40 40 0 1.2])
title('Magnitude of windowed FFT of Real input Sinusoid at f_C= 5, and f_S= 80')

subplot(4,1,3);
plot([-0.5:1/1000:0.5-1/1000]*80,fftshift(real(fft(real(x.*ww),1000))),'b','linewidth',2)
grid on
axis([-10 10 -1.2 1.2])
title('Real Part of windowed FFT of Real input Sinusoid at f_C= 5, and f_S= 80')

subplot(4,1,4);
plot([-0.5:1/1000:0.5-1/1000]*80,fftshift(imag(fft(real(x.*ww),1000))),'b','linewidth',2)
grid on
axis([-10 10 -1.2 1.2])
title('Imag Part of windowed FFT of Real input Sinusoid at f_C= 5, and f_S= 80')


figure(7)
subplot(4,1,1);
plot(0:100,real(x(1:101)),'b','linewidth',2)
hold on;
plot(0:100,imag(x(1:101)),'r','linewidth',2)
hold off;
grid on
axis([0 100 -1.2 1.2])
title('Complex input Sinusoid at f_C= 5, and f_S= 80')

subplot(4,1,2);
ww=kaiser(320,8)';
ww=ww/sum(ww);
plot([-0.5:1/1000:0.5-1/1000]*80,fftshift(abs(fft(imag(x).*ww,1000))),'b','linewidth',2)
grid on
axis([-40 40 0 1.2])
title('Magnitude of windowed FFT of imaginary input Sinusoid at f_C= 5, and f_S= 80')

subplot(4,1,3);
plot([-0.5:1/1000:0.5-1/1000]*80,fftshift(real(fft(imag(x.*ww),1000))),'b','linewidth',2)
grid on
axis([-10 10 -1.2 1.2])
title('Real Part of windowed FFT of imaginary input Sinusoid at f_C= 5, and f_S= 80')

subplot(4,1,4);
plot([-0.5:1/1000:0.5-1/1000]*80,fftshift(imag(fft(imag(x.*ww),1000))),'b','linewidth',2)
grid on
axis([-10 10 -1.2 1.2])
title('Imag Part of windowed FFT of Imaginary input Sinusoid at f_C= 5, and f_S= 80')