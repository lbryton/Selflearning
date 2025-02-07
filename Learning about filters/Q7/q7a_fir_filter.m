% Clearing data
clc;
close all;
clf reset;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up
% Constants:
f_s = 80e3;
Rs= 60;     % Stopband ripple (dB)

% Filter 1
f_pass = 5e3;
f_stop = 10e3;
% Filter Order
% N = (f_s/(f_stop-f_pass)) * attenuation(db) / 22
N = ceil(f_s/(f_stop - f_pass)*Rs/22);
h1 = firpm(N, [0 f_pass f_stop f_s/2]/(f_s/2), [1 1 0 0], [1 10]);

% Filter 2
f_pass = 10e3;
f_stop = 15e3;
% Filter Order
% N = (f_s/(f_stop-f_pass)) * attenuation(db) / 22
N = ceil(f_s/(f_stop - f_pass)*Rs/22);
h2 = firpm(N, [0 f_pass f_stop f_s/2]/(f_s/2), [1 1 0 0], [1 10]);

% Filter 3
f_pass = 15e3;
f_stop = 20e3;
% Filter Order
% N = (f_s/(f_stop-f_pass)) * attenuation(db) / 22
N = ceil(f_s/(f_stop - f_pass)*Rs/22);
h3 = firpm(N, [0 f_pass f_stop f_s/2]/(f_s/2), [1 1 0 0], [1 10]);

% Filter 4
f_pass = 20e3;
f_stop = 25e3;
% Filter Order
% N = (f_s/(f_stop-f_pass)) * attenuation(db) / 22
N = ceil(f_s/(f_stop - f_pass)*Rs/22);
h4 = firpm(N, [0 f_pass f_stop f_s/2]/(f_s/2), [1 1 0 0], [1 10]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Impulse response
figure(1)
% Filter 1
subplot(4,1,1)
plot(0:length(h1)-1,h1,'linewidth',2)
grid on
% axis([0 50 -0.05 0.20])
set(gca,'fontsize',14)
set(gca,'xtick',[0:5:50])
title('Impulse response, filter 1','fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

% Filter 2
subplot(4,1,2)
plot(0:length(h2)-1,h2,'linewidth',2)
grid on
% axis([0 50 -0.05 0.20])
set(gca,'fontsize',14)
set(gca,'xtick',[0:5:50])
title('Impulse response, filter 2','fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

% Filter 3
subplot(4,1,3)
plot(0:length(h3)-1,h3,'linewidth',2)
grid on
% axis([0 50 -0.05 0.20])
set(gca,'fontsize',14)
set(gca,'xtick',[0:5:50])
title('Impulse response, filter 3','fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

% Filter 4
subplot(4,1,4)
plot(0:length(h4)-1,h4,'linewidth',2)
grid on
% axis([0 50 -0.05 0.20])
set(gca,'fontsize',14)
set(gca,'xtick',[0:5:50])
title('Impulse response, filter 4','fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Log magnitude frequency response
% Filter 1
figure(2)
subplot(4,1,1)
fft_h1 = 20*log10(abs(fftshift(fft(h1,4096))));
plot(-40:80/4096:40-80/4096,fft_h1,'linewidth',2)
grid on
axis([-40 40 -100 10])
set(gca,'fontsize',14)
title('Log magnitude frequency response for filter 1')
xlabel('Normalized frequency')
ylabel('Log Magnitude (dB)')

% Filter 2
subplot(4,1,2)
fft_h2 = 20*log10(abs(fftshift(fft(h2,4096))));
plot(-40:80/4096:40-80/4096,fft_h2,'linewidth',2)
grid on
axis([-40 40 -100 10])
set(gca,'fontsize',14)
title('Log magnitude frequency response for filter 2')
xlabel('Normalized frequency')
ylabel('Log Magnitude (dB)')

% Filter 3
subplot(4,1,3)
fft_h3 = 20*log10(abs(fftshift(fft(h3,4096))));
plot(-40:80/4096:40-80/4096,fft_h3,'linewidth',2)
grid on
axis([-40 40 -100 10])
set(gca,'fontsize',14)
title('Log magnitude frequency response for filter 3')
xlabel('Normalized frequency')
ylabel('Log Magnitude (dB)')

% Filter 4
subplot(4,1,4)
fft_h4 = 20*log10(abs(fftshift(fft(h4,4096))));
plot(-40:80/4096:40-80/4096,fft_h4,'linewidth',2)
grid on
axis([-40 40 -100 10])
set(gca,'fontsize',14)
title('Log magnitude frequency response for filter 1')
xlabel('Normalized frequency')
ylabel('Log Magnitude (dB)')




