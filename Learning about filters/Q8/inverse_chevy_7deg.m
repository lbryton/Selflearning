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
f_s = 40e3;
f_stop = 12e3;
f_pass = 7.5e3;

W1 = 2 * f_stop / f_s;
% W1 = 0.30;   % Normalized cutoff frequency (rounded up for margin)

% Designs a Chebyshev Type 2 digital filter
% Ripple in stopband, and ripple-free passband
% Outputs numerator (bb) and denominator (aa) coefficients
[bb,aa] = cheby2(N,Rs,W1, 'low');

% Second order section
[sos,g] = tf2sos(bb,aa);

c1 = sos(1,:);
c2 = sos(2,:);
c3 = sos(3,:);
c4 = sos(4,:);

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
axis([0 50 -0.5 0.5])
set(gca,'fontsize',14)
set(gca,'xtick',[0:5:100])
title('Filter Impulse response','fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

subplot(2,1,2)
fft_hh = fftshift(20*log10(abs(fft(output,1024))));
plot(-0.5:1/1024:0.5-1/1024,fft_hh,'linewidth',2)
grid on
axis([-0.5 0.5 -80 10])
set(gca,'fontsize',14)
title('Frequency Response','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('Log Magnitude (dB)','fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% Zoom to passband ripple + Log mag frequency response
figure(2)
pb_percent = f_pass;
sb_percent = f_stop;

plot(-f_s/2:f_s/1024:f_s/2- f_s/1024,fft_hh,'linewidth',2)
hold on
plot([sb_percent sb_percent f_s/2],[-20 -60 -60],'r','linewidth',2)
plot([-sb_percent -sb_percent -f_s/2],[-20 -60 -60],'r','linewidth',2)
plot([-pb_percent -pb_percent pb_percent pb_percent],[-60 0 0 -60],'r','linewidth',2)
hold off
grid on
axis([-f_s/2 f_s/2 -80 10])
set(gca,'fontsize',14)
title('Frequency Response, Passband Ripple = 0.1 dB, Stopband Attenuation > 80 dB','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('Log Magnitude (dB)','fontsize',14)

axes('position',[0.596 0.25 0.231 0.128])
plot(-f_s/2:f_s/1024:f_s/2- f_s/1024,fft_hh,'linewidth',2)
hold on
plot([-pb_percent -pb_percent pb_percent pb_percent],[-0.13 -0.1 -0.1 -0.13],'r','linewidth',2)
plot([-pb_percent -pb_percent pb_percent pb_percent],[+0.13 +0.1 +0.1 +0.13],'r','linewidth',2)
hold off
grid on
axis([-7.5e3 7.5e3 -0.15 0.15])
title('Zoom to Passband Ripple','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('(dB)','fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Phase response and Group delay
figure(3)
[H, w] = freqz(bb, aa, 1024, f_s);

% Extend to negative frequencies
w_full = [-flip(w(2:end)); w];     % Frequency axis from -fs/2 to fs/2
H_full = [conj(flip(H(2:end))); H]; % Mirror and conjugate response

% Phase response (unwrap to make it continuous)
% Represents how much phase distortion the filter has at each frequency
phase = unwrap(angle(H_full));

% Group delay (manually compute for negative frequencies)
% Represents how much time delay (in samples) the signal components
% experience at each frequency
[gd, w_gd] = grpdelay(bb, aa, 1024, f_s);
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



%% 1000 QPSK
figure(5)
q = qpsk(1000);
qf1 = filter(bb,aa,q);
nyq_filt = sqrt_nyq_y2(4,0.25,12,0);
qf2 = filter(nyq_filt, 1 , qf1);
plot(0:4:400,qf2(1:4:401),'ro')
title('Real Time Series: Inverse cheby, Match Filter, Small Delay')
