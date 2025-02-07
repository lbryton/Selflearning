% Clearing data
clc;
close all;
clf reset;
clear;

%%
% Set up

N=7;        % Order of filter
Rs= 60;     % Stopband ripple (dB)

% Math for stopband edge frequency
% W1 = (stopband frequency) / (sampling rate / 2)
sampling_rate = 80e3;
stop_freq = 10e3;
pass_freq = 5e3;

W1 = 2 * pass_freq / sampling_rate;
% W1 = 0.30;   % Normalized cutoff frequency (rounded up for margin)

% Designs a Chebyshev Type 2 digital filter
% Ripple in stopband, and ripple-free passband
% Outputs numerator (bb) and denominator (aa) coefficients
[bb,aa] =cheby1(N,0.1,W1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Second order sections (SOS): Represents higher-order filter
% as a cascade of second-order filters
% In our case, we represent a 7th order Chebyshev Type 2 digital filter
% as 4 2nd order filters
% Columns of the second order filter represents:
% c00, c10, c20, 1, c11, c22
% c00 - c20: represents coefficients for the numerator of the second-order
% filter
% 1 represents the scaling fators
% c11 - c22: represents coefficients for the demoninator of the second-order
% filter

[sos,g] = tf2sos(bb,aa);
% Plot magnitude and phase
% freqz(bb,aa)

% The produced second order section for cheby2 with 7-th degree
% 1.0000    1.0000         0    1.0000   -0.4406         0
% 1.0000    0.3187    1.0000    1.0000   -0.9782    0.2859
% 1.0000   -0.8076    1.0000    1.0000   -1.2103    0.5134
% 1.0000   -1.1418    1.0000    1.0000   -1.4942    0.8157

% Spliting by filters (also in order starting with c1 for cascade and 
% incrementing values for other 2nd order filters)
c1 = sos(1,:);
c2 = sos(2,:);
c3 = sos(3,:);
c4 = sos(4,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% Calculating stage factor (used to normalize gain of each sos)
filter_cnt = length(sos(:,1));
factor = zeros(1,filter_cnt);

% Gets factor for each second-order filter
for i=1:filter_cnt
    factor(i) = sum(sos(i,4:6))/sum(sos(i,1:3)); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% Calculating impulse response (100 samples)
x = [1 zeros(1,100)];
weights = zeros(1,filter_cnt*2-1);
disp(length(weights))
% goes through each sample
feed_forward = 0;
output = zeros(1,101);
for n=1:101
    % Cascades through each filter
    for m=1:filter_cnt
        % w = 0;
        if (m == 1)
            w = x(n)-weights(1) * sos(m,5);                         % feedback term
            feed_forward = (w + weights(1)*sos(1,2))*factor(1);     % scaled feed forward term
            weights(1) = w;                                         % shift feedback into register
        else
            w = feed_forward - weights(2*m-2) * sos(m,5) - ...      % feedback term
                weights(2*m-1) * sos(m,6);
            feed_forward = (w + weights(2*m-2) * sos(m,2) + ...     % scaled feed forward term
                weights(2*m-1) * sos(m,3)) * factor(m);
            weights(2*m-1) = weights(2*m-2);                        % shift feedback from first register to second
            weights(2*m-2) = w;                                     % shift feedback into first
        end
        if (m == filter_cnt)
            output(n) = feed_forward;
        end
    end
end

figure(1)
subplot(2,1,1)
plot(0:100,output,'linewidth',2)
grid on
axis([0 50 -0.05 0.20])
set(gca,'fontsize',14)
set(gca,'xtick',[0:5:100])
title('Filter Impulse response','fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

subplot(2,1,2)
fft_hh = fftshift(20*log10(abs(fft(output,1024))));
plot(-0.5:1/1024:0.5-1/1024,fft_hh,'linewidth',2)
grid on
axis([0 0.5 -80 10])
set(gca,'fontsize',14)
title('Frequency Response','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('Log Magnitude (dB)','fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% Zoom to passband ripple
% Notice how changes to Wn changes where it stops
figure(2)
pb_percent = pass_freq/sampling_rate;
sb_percent = stop_freq/sampling_rate;

plot(-0.5:1/1024:0.5-1/1024,fft_hh,'linewidth',2)
hold on
plot([sb_percent sb_percent 0.5],[-20 -60 -60],'r','linewidth',2)
plot([-sb_percent -sb_percent -0.5],[-20 -60 -60],'r','linewidth',2)
plot([-pb_percent -pb_percent pb_percent pb_percent],[-60 0 0 -60],'r','linewidth',2)
hold off
grid on
axis([-0.5 0.5 -80 10])
set(gca,'fontsize',14)
title('Frequency Response, Passband Ripple = 0.1 dB, Stopband Attenuation > 80 dB','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('Log Magnitude (dB)','fontsize',14)

axes('position',[0.596 0.25 0.231 0.128])
plot(-0.5:1/1024:0.5-1/1024,fft_hh,'linewidth',2)
hold on
plot([-0.0625 -0.0625 0.0625 0.0625],[-0.13 -0.1 -0.1 -0.13],'r','linewidth',2)
plot([-0.0625 -0.0625 0.0625 0.0625],[+0.13 +0.1 +0.1 +0.13],'r','linewidth',2)
hold off
grid on
axis([-0.1 0.1 -0.2 0.2])
title('Zoom to Passband Ripple','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('(dB)','fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Phase response and Group delay
figure(3)
[H, w] = freqz(bb, aa, 1024, sampling_rate);

% Extend to negative frequencies
w_full = [-flip(w(2:end)); w];     % Frequency axis from -fs/2 to fs/2
H_full = [conj(flip(H(2:end))); H]; % Mirror and conjugate response

% Phase response (unwrap to make it continuous)
% Represents how much phase distortion the filter has at each frequency
phase = unwrap(angle(H_full));

% Group delay (manually compute for negative frequencies)
% Represents how much time delay (in samples) the signal components
% experience at each frequency
[gd, w_gd] = grpdelay(bb, aa, 1024, sampling_rate);
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
% Pole-zero plot
% Used for understanding stability of system
% Zeros are the root of the numerators of the polynomial
% Poles are the root of the denominators of the polynomial

figure(4)
subplot(2,1,1)
plot(0,0)
plot(exp(1i*2*pi*(0:0.01:1)),'r','linewidth',2)
hold on

% Plotting zeroes
plot(-1,0,'ro','linewidth',2)
plot(roots(c2(1:3)),'ro','linewidth',2)
plot(roots(c3(1:3)),'ro','linewidth',2)
plot(roots(c4(1:3)),'ro','linewidth',2)

hold off
grid on
axis('equal')
axis([-1.2 1.2 -1.2 1.2])
title('Zero-pole: Zero section')

subplot(2,1,2)
plot(0,0)
plot(exp(1i*2*pi*(0:0.01:1)),'r','linewidth',2)
hold on

% Plotting poles
plot(-c1(5),0,'rx','linewidth',2)
plot(roots(c2(4:6)),'rx','linewidth',2)
plot(roots(c3(4:6)),'rx','linewidth',2)
plot(roots(c4(4:6)),'rx','linewidth',2)
hold off
grid on
axis('equal')
axis([-1.2 1.2 -1.2 1.2])
title('Zero-pole: Pole section')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency response of each SOS
% Ran for 100 samples
figure(5)
h1=filter(c1(1,2), c1(4:5), [1 zeros(1,100)])*sum(c1(4:5))/sum(c1(1:2));
h2=filter(c2(1,3), c2(4:6), [1 zeros(1,100)])*sum(c2(4:6))/sum(c2(1:3));
h3=filter(c3(1,3), c3(4:6), [1 zeros(1,100)])*sum(c3(4:6))/sum(c3(1:3));
h4=filter(c4(1,3), c4(4:6), [1 zeros(1,100)])*sum(c4(4:6))/sum(c4(1:3));

fh1=fftshift(20*log10(abs(fft(h1,2048))));
fh2=fftshift(20*log10(abs(fft(h2,2048))));
fh3=fftshift(20*log10(abs(fft(h3,2048))));
fh4=fftshift(20*log10(abs(fft(h4,2048))));

subplot(2,2,1)
plot(-0.5:1/2048:0.5-1/2048,fh1,'linewidth',2)
grid on
axis([-0.15 0.15 -90 20])
title('Frequency Response H_1')

subplot(2,2,2)
plot(-0.5:1/2048:0.5-1/2048,fh2,'linewidth',2)
grid on
axis([-0.15 0.15 -90 20])
title('Frequency Response H_2')

subplot(2,2,3)
plot(-0.5:1/2048:0.5-1/2048,fh3,'linewidth',2)
grid on
axis([-0.15 0.15 -90 20])
title('Frequency Response H_3')

subplot(2,2,4)
plot(-0.5:1/2048:0.5-1/2048,fh4,'linewidth',2)
grid on
axis([-0.15 0.15 -90 20])
title('Frequency Response H_4')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impulse response of each SOS
% Ran for 100 samples
figure(6)
h1=filter(c1(1,2), c1(4:5), [1 zeros(1,100)])*sum(c1(4:5))/sum(c1(1:2));
h2=filter(c2(1,3), c2(4:6), [1 zeros(1,100)])*sum(c2(4:6))/sum(c2(1:3));
h3=filter(c3(1,3), c3(4:6), [1 zeros(1,100)])*sum(c3(4:6))/sum(c3(1:3));
h4=filter(c4(1,3), c4(4:6), [1 zeros(1,100)])*sum(c4(4:6))/sum(c4(1:3));

subplot(2,2,1)
plot(0:100,h1,'linewidth',2)
grid on
% axis([-0.15 0.15 -90 20])
title('Impulse Response H_1')

subplot(2,2,2)
plot(0:100,h2,'linewidth',2)
grid on
% axis([-0.15 0.15 -90 20])
title('Impulse Response H_2')

subplot(2,2,3)
plot(0:100,h3,'linewidth',2)
grid on
% axis([-0.15 0.15 -90 20])
title('Impulse Response H_3')

subplot(2,2,4)
plot(0:100,h4,'linewidth',2)
grid on
% axis([-0.15 0.15 -90 20])
title('Impulse Response H_4')

