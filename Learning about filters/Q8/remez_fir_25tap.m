% Clearing data
clc;
close all;
clf reset;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up
% Constants:
f_s = 40e3;
f_snq = f_s/2;
Rs= 60;     % Stopband ripple (dB)
% Filter 1
f_pass = 7.5e3;
f_stop = 12e3;
% Filter Order
% N = (f_s/(f_stop-f_pass)) * attenuation(db) / 22
N = 25 - 1;
f_norm = [0, f_pass, f_stop, f_snq] / f_snq;
h1 = firpm(N, f_norm, [1 1 0 0], [1 10]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Impulse response
figure(1)
% Filter 1
subplot(4,1,1)
plot(0:100,[h1, zeros(1,101- length(h1),1)],'linewidth',2)
grid on
axis([0 50 -0.2 0.5])
set(gca,'fontsize',14)
set(gca,'xtick',[0:10:100])
title('Impulse response, filter 1','fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Log magnitude frequency response
% Filter 1
figure(2)
subplot(4,1,1)
fft_h1 = 20*log10(abs(fftshift(fft(h1,4096))));
plot(-20:40/4096:20-40/4096,fft_h1,'linewidth',2)
grid on
axis([-f_snq/1e3 f_snq/1e3 -100 10])
set(gca,'fontsize',14)
title('Log magnitude frequency response for filter 1')
xlabel('Normalized frequency')
ylabel('Log Magnitude (dB)')

axes('position',[0.596 0.25 0.231 0.128])
plot(-20:40/4096:20-40/4096,fft_h1,'linewidth',2)
hold on
plot([-f_pass/1e3 -f_pass/1e3 f_pass/1e3 f_pass/1e3],[-0.08 -0.05 -0.05 -0.08],'r','linewidth',2)
plot([-f_pass/1e3 -f_pass/1e3 f_pass/1e3 f_pass/1e3],[+0.08 +0.05 +0.05 +0.08],'r','linewidth',2)
hold off
grid on
axis([(-f_pass/1e3)-1 (f_pass/1e3)+1 -0.15 0.1])
title('Zoom to Passband Ripple','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('(dB)','fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Phase response and Group delay
figure(3)
[H, w] = freqz(h1, 1, 1024, f_s);

% Extend to negative frequencies
w_full = [-flip(w(2:end)); w];     % Frequency axis from -fs/2 to fs/2
H_full = [conj(flip(H(2:end))); H]; % Mirror and conjugate response

% Phase response (unwrap to make it continuous)
% Represents how much phase distortion the filter has at each frequency
phase = unwrap(angle(H_full));

% Group delay (manually compute for negative frequencies)
% Represents how much time delay (in samples) the signal components
% experience at each frequency
[gd, w_gd] = grpdelay(h1, 1, 1024, f_s);
gd_full = [flip(gd(2:end)); gd]; % Symmetric group delay for negative frequencies

% Plot Phase Response
subplot(2, 1, 1);
plot(w_full, phase);
title('Phase Response');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
grid on;

% Plot Group Delay
subplot(2, 1, 2);
plot([-flip(w_gd(2:end)); w_gd], gd_full);
title('Group Delay Response');
xlabel('Frequency (Hz)');
ylabel('Group Delay (samples)');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Pole Zero plot
figure(4)
plot(0,0)
plot(exp(1i*2*pi*(0:0.01:1)),'b','linewidth',2)
hold on
plot(roots(h1),'ro','linewidth',2)
hold off
grid on
axis('equal')
axis([-1.2 1.2 -1.2 1.2])
title('Zero-pole: 25-tap remez')


%% Match filter
figure(5)
vv=remez(158,[0 7.5 32.5 640]/640,[1 1 0 0]);
vv2=32*reshape([vv 0],32,5);
x3a = qpsk(1000);
gg = sqrt_nyq_y2(4,0.25,12,0);

% delayed version of matched filter output
x4a=filter(gg,1,x3a)/(gg*gg');
x4a=filter(vv2(31,:),1,x4a);
plot(0:4:400,x4a(1:4:401),'ro')
title('Real Time Series: Elliptic Filter, Match Filter, Small Delay')




