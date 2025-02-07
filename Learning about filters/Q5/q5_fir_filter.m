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
h2 = firpm(N2, [0 f_pass2 f_stop2 f_s2/2]/(f_s2/2), [1 1 0 0], [1 10]);


% Remez filter 3
f_s3 = 80e3;
f_pass3 = 5e3;
f_stop3 = 10e3;
Sr3 = 60;
N3 = ceil(f_s3/(f_stop3 - f_pass3) * Sr3 / 22);
h3 = firpm(N3, [0 f_pass3 f_stop3 f_s3/2]/(f_s3/2), [1 1 0 0], [10 10]);

% Remez filter 4
f_s4 = 80e3;
f_pass4 = 5e3;
f_stop4 = 10e3;
Sr4 = 80;
N4 = ceil(f_s4/(f_stop4 - f_pass4) * Sr4 / 22);
h4 = firpm(N4, [0 f_pass4 f_stop4 f_s4/2]/(f_s4/2), [1 1 0 0], [10 100]);

%% Filter 1 plots
figure(1)
subplot(3,1,1)

% Impulse Response
plot(0:length(h1)-1,h1,'linewidth',2)
grid on
% axis([0 50 -0.05 0.20])
set(gca,'fontsize',14)
set(gca,'xtick',[0:5:50])
title('Impulse response, filter 1','fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

% Log magnitude frequency response
subplot(3,1,2)
fft_h1 = 20*log10(abs(fftshift(fft(h1,4096))));
plot(-40:80/4096:40-80/4096,fft_h1,'linewidth',2)
grid on
axis([-20 20 -100 10])
set(gca,'fontsize',14)
title('Log magnitude frequency response for filter 1')
xlabel('Normalized frequency')
ylabel('Log Magnitude (dB)')

% Zoom to passband ripple
subplot(3,2,5)
% axes('position',[0.596 0.25 0.231 0.128])
plot(-40:80/4096:40-80/4096,fft_h1,'linewidth',2)
grid on
axis([-6 6 -0.2 0.2])
title('Zoom to Passband Ripple','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('(dB)','fontsize',14)

% Pole Zero plot
set(gcf, 'WindowState', 'maximized');
subplot(3,2,6)
plot(0,0)
plot(exp(1i*2*pi*(0:0.01:1)),'b','linewidth',2)
hold on
plot(roots(h1),'ro','linewidth',2)
hold off
grid on
axis('equal')
axis([-1.2 1.2 -1.2 1.2])
title('Zero-pole: filter 1')

%% Filter 2 plots
figure(2)
subplot(3,1,1)

% Impulse Response
plot(0:length(h2)-1,h2,'linewidth',2)
grid on
% axis([0 50 -0.05 0.20])
set(gca,'fontsize',14)
set(gca,'xtick',[0:5:50])
title('Impulse response, filter 2','fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

% Log magnitude frequency response
subplot(3,1,2)
fft_h2 = 20*log10(abs(fftshift(fft(h2,4096))));
plot(-40:80/4096:40-80/4096,fft_h2,'linewidth',2)
grid on
axis([-20 20 -100 10])
set(gca,'fontsize',14)
title('Log magnitude frequency response for filter 2')
xlabel('Normalized frequency')
ylabel('Log Magnitude (dB)')

% Zoom to passband ripple
subplot(3,2,5)
% axes('position',[0.596 0.25 0.231 0.128])
plot(-40:80/4096:40-80/4096,fft_h2,'linewidth',2)
grid on
axis([-6 6 -0.2 0.2])
title('Zoom to Passband Ripple','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('(dB)','fontsize',14)

% Pole Zero plot
set(gcf, 'WindowState', 'maximized');
subplot(3,2,6)
plot(0,0)
plot(exp(1i*2*pi*(0:0.01:1)),'b','linewidth',2)
hold on
plot(roots(h2),'ro','linewidth',2)
hold off
grid on
axis('equal')
axis([-1.2 1.2 -1.2 1.2])
title('Zero-pole: filter 2')

%% Filter 3 plots
figure(3)
subplot(3,1,1)

% Impulse Response
plot(0:length(h3)-1,h3,'linewidth',2)
grid on
% axis([0 50 -0.05 0.20])
set(gca,'fontsize',14)
set(gca,'xtick',[0:5:50])
title('Impulse response, filter 3','fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

% Log magnitude frequency response
subplot(3,1,2)
fft_h3 = 20*log10(abs(fftshift(fft(h3,4096))));
plot(-40:80/4096:40-80/4096,fft_h3,'linewidth',2)
grid on
axis([-20 20 -100 10])
set(gca,'fontsize',14)
title('Log magnitude frequency response for filter 3')
xlabel('Normalized frequency')
ylabel('Log Magnitude (dB)')

% Zoom to passband ripple
subplot(3,2,5)
% axes('position',[0.596 0.25 0.231 0.128])
plot(-40:80/4096:40-80/4096,fft_h3,'linewidth',2)
grid on
axis([-6 6 -0.2 0.2])
title('Zoom to Passband Ripple','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('(dB)','fontsize',14)

% Pole Zero plot
set(gcf, 'WindowState', 'maximized');
subplot(3,2,6)
plot(0,0)
plot(exp(1i*2*pi*(0:0.01:1)),'b','linewidth',2)
hold on
plot(roots(h3),'ro','linewidth',2)
hold off
grid on
axis('equal')
axis([-1.2 1.2 -1.2 1.2])
title('Zero-pole: filter 3')

%% Filter 4 plots
figure(4)
subplot(3,1,1)

% Impulse Response
plot(0:length(h4)-1,h4,'linewidth',2)
grid on
% axis([0 50 -0.05 0.20])
set(gca,'fontsize',14)
set(gca,'xtick',[0:5:50])
title('Impulse response, filter 4','fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

% Log magnitude frequency response
subplot(3,1,2)
fft_h4 = 20*log10(abs(fftshift(fft(h4,4096))));
plot(-40:80/4096:40-80/4096,fft_h4,'linewidth',2)
grid on
axis([-20 20 -100 10])
set(gca,'fontsize',14)
title('Log magnitude frequency response for filter 4')
xlabel('Normalized frequency')
ylabel('Log Magnitude (dB)')

% Zoom to passband ripple
subplot(3,2,5)
% axes('position',[0.596 0.25 0.231 0.128])
plot(-40:80/4096:40-80/4096,fft_h4,'linewidth',2)
grid on
axis([-6 6 -0.2 0.2])
title('Zoom to Passband Ripple','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('(dB)','fontsize',14)

% Pole Zero plot
set(gcf, 'WindowState', 'maximized');
subplot(3,2,6)
plot(0,0)
plot(exp(1i*2*pi*(0:0.01:1)),'b','linewidth',2)
hold on
plot(roots(h4),'ro','linewidth',2)
hold off
grid on
axis('equal')
axis([-1.2 1.2 -1.2 1.2])
title('Zero-pole: filter 4')




