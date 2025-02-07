% Clearing data
clc;
close all;
clf reset;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up
% Math for stopband edge frequency
% W1 = (stopband frequency) / (sampling rate / 2)
f_s = 80e3;
f_pass = 5e3;
f_stop = 10e3;
f_ss = f_s/2;
stopband_ripple = 60;
Rs= 60;     % Stopband ripple (dB)

% Filter Order
% N = (f_s/(f_stop-f_pass)) * attenuation(db) / 22
N = ceil(f_s/(f_stop - f_pass)*Rs/22);


% We do freq * 2 / samp_freq to get ratio between passband frequency and
% nquist frequency

% Remez filter 1
f_s1 = 80e3;
f_pass1 = 5e3;
f_stop1 = 10e3;
Sr1 = 60;

% Filter Order
% N = (f_s/(f_stop-f_pass)) * attenuation(db) / 22
N1 = ceil(f_s1/(f_stop1 - f_pass1) * Sr1 / 22);
h1 = firpm(N1, [0 f_pass1 f_stop1 f_s1/2]/(f_s1/2), [1 1 0 0], [1 10]);

% Remez filter 2
f_s2 = 80e3;
f_pass2 = 5e3;
f_stop2 = 10e3;
Sr2 = 80;
N2 = ceil(f_s2/(f_stop2 - f_pass2) * Sr2 / 22);
h2 = firpm(N2, [0 f_pass2 f_stop2 f_s2/2]/(f_s2/2), [1 1 0 0], [10 100]);

% Quantization bits
bits1 = Sr1 / 5;
bits1 = 2 ^(bits1-1);
bits2 = Sr2 / 5;
bits2 = 2 ^(bits2-1);

% Nonscaled
quant_noscale1 = round(h1 * bits1) / bits1;
quant_noscale2 = round(h2 * bits2) / bits2;

% Scaled 
% Max value for scaling
scaled_coef1 = max(abs(h1));
scaled_coef2 = max(abs(h2));

h1_scale = h1/scaled_coef1;
h2_scale = h2/scaled_coef2;

quant_scale1 = round(h1_scale * bits1) / bits1 * scaled_coef1;
quant_scale2 = round(h2_scale * bits2) / bits2 * scaled_coef2;

%% Filter 1 no scaling plots
figure(1)
subplot(3,1,1)

% Impulse Response
plot(0:length(quant_noscale1)-1,quant_noscale1,'linewidth',2)
grid on
% axis([0 50 -0.05 0.20])
set(gca,'fontsize',14)
set(gca,'xtick',[0:5:50])
title('Impulse response, 60db Rs, non-scaled','fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

% Log magnitude frequency response
subplot(3,1,2)
fft_h1_ns = 20*log10(abs(fftshift(fft(quant_noscale1,4096))));
plot(-40:80/4096:40-80/4096,fft_h1_ns,'linewidth',2)
grid on
axis([-20 20 -100 10])
set(gca,'fontsize',14)
title('Log magnitude frequency response for 60db Rs, non-scaled')
xlabel('Normalized frequency')
ylabel('Log Magnitude (dB)')

% Zoom to passband ripple
subplot(3,2,5)
% axes('position',[0.596 0.25 0.231 0.128])
plot(-40:80/4096:40-80/4096,fft_h1_ns,'linewidth',2)
grid on
axis([-6 6 -0.2 0.2])
title('Zoom to Passband Ripple, 60db Rs, non-scaled','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('(dB)','fontsize',14)

% Pole Zero plot
set(gcf, 'WindowState', 'maximized');
subplot(3,2,6)
plot(0,0)
plot(exp(1i*2*pi*(0:0.01:1)),'b','linewidth',2)
hold on
plot(roots(quant_noscale1),'ro','linewidth',2)
hold off
grid on
axis('equal')
axis([-1.2 1.2 -1.2 1.2])
title('Zero-pole: 60db Rs, non-scaled')

%% Filter 2 no scaling plots
figure(2)
subplot(3,1,1)

% Impulse Response
plot(0:length(quant_noscale2)-1,quant_noscale2,'linewidth',2)
grid on
% axis([0 50 -0.05 0.20])
set(gca,'fontsize',14)
set(gca,'xtick',[0:5:50])
title('Impulse response, 80db Rs, non-scaled','fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

% Log magnitude frequency response
subplot(3,1,2)
fft_h2_ns = 20*log10(abs(fftshift(fft(quant_noscale2,4096))));
plot(-40:80/4096:40-80/4096,fft_h2_ns,'linewidth',2)
grid on
axis([-20 20 -100 10])
set(gca,'fontsize',14)
title('Log magnitude frequency response for 80db Rs, non-scaled')
xlabel('Normalized frequency')
ylabel('Log Magnitude (dB)')

% Zoom to passband ripple
subplot(3,2,5)
% axes('position',[0.596 0.25 0.231 0.128])
plot(-40:80/4096:40-80/4096,fft_h2_ns,'linewidth',2)
grid on
axis([-6 6 -0.2 0.2])
title('Zoom to Passband Ripple, 80db Rs, non-scaled','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('(dB)','fontsize',14)

% Pole Zero plot
set(gcf, 'WindowState', 'maximized');
subplot(3,2,6)
plot(0,0)
plot(exp(1i*2*pi*(0:0.01:1)),'b','linewidth',2)
hold on
plot(roots(quant_noscale2),'ro','linewidth',2)
hold off
grid on
axis('equal')
axis([-1.2 1.2 -1.2 1.2])
title('Zero-pole: 80db Rs, non-scaled')

%% Filter 1 scaling plots
figure(3)
subplot(3,1,1)

% Impulse Response
plot(0:length(quant_scale1)-1,quant_scale1,'linewidth',2)
grid on
% axis([0 50 -0.05 0.20])
set(gca,'fontsize',14)
set(gca,'xtick',[0:5:50])
title('Impulse response, 60db Rs, scaled','fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

% Log magnitude frequency response
subplot(3,1,2)
fft_h1_s = 20*log10(abs(fftshift(fft(quant_scale1,4096))));
plot(-40:80/4096:40-80/4096,fft_h1_s,'linewidth',2)
grid on
axis([-20 20 -100 10])
set(gca,'fontsize',14)
title('Log magnitude frequency response for 60db Rs, scaled')
xlabel('Normalized frequency')
ylabel('Log Magnitude (dB)')

% Zoom to passband ripple
subplot(3,2,5)
% axes('position',[0.596 0.25 0.231 0.128])
plot(-40:80/4096:40-80/4096,fft_h1_s,'linewidth',2)
grid on
axis([-6 6 -0.2 0.2])
title('Zoom to Passband Ripple, 60db Rs, scaled','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('(dB)','fontsize',14)

% Pole Zero plot
set(gcf, 'WindowState', 'maximized');
subplot(3,2,6)
plot(0,0)
plot(exp(1i*2*pi*(0:0.01:1)),'b','linewidth',2)
hold on
plot(roots(quant_scale1),'ro','linewidth',2)
hold off
grid on
axis('equal')
axis([-1.2 1.2 -1.2 1.2])
title('Zero-pole: 60db Rs, scaled')

%% Filter 2 scaling plots
figure(4)
subplot(3,1,1)

% Impulse Response
plot(0:length(quant_scale2)-1,quant_scale2,'linewidth',2)
grid on
% axis([0 50 -0.05 0.20])
set(gca,'fontsize',14)
set(gca,'xtick',[0:5:50])
title('Impulse response, 80db Rs, scaled','fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

% Log magnitude frequency response
subplot(3,1,2)
fft_h2_s = 20*log10(abs(fftshift(fft(quant_scale2,4096))));
plot(-40:80/4096:40-80/4096,fft_h2_s,'linewidth',2)
grid on
axis([-20 20 -100 10])
set(gca,'fontsize',14)
title('Log magnitude frequency response for 80db Rs, scaled')
xlabel('Normalized frequency')
ylabel('Log Magnitude (dB)')

% Zoom to passband ripple
subplot(3,2,5)
% axes('position',[0.596 0.25 0.231 0.128])
plot(-40:80/4096:40-80/4096,fft_h2_s,'linewidth',2)
grid on
axis([-6 6 -0.2 0.2])
title('Zoom to Passband Ripple, 80db Rs, scaled','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('(dB)','fontsize',14)

% Pole Zero plot
set(gcf, 'WindowState', 'maximized');
subplot(3,2,6)
plot(0,0)
plot(exp(1i*2*pi*(0:0.01:1)),'b','linewidth',2)
hold on
plot(roots(quant_scale2),'ro','linewidth',2)
hold off
grid on
axis('equal')
axis([-1.2 1.2 -1.2 1.2])
title('Zero-pole: 80db Rs, scaled')