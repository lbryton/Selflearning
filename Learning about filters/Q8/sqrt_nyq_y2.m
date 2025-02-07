function nyq_filt = sqrt_nyq_y2(sps, alpha, delay, plot_flag)
% sps: samples per symbol (4 input)
% alpha: excess BW
% delay: symbols from sample 1 to filter center
% plot_flag: 0 for no figures

% Need to get whole length from start of sample 1 to last sample of impulse
% Response

improved_res = 400;

%% Designing filter

% Getting filter length
arg = -delay : 1/sps : delay;
N = length(arg);
disp(N)
% Whole range
% NN/2 corresponds to 1/T
NN = improved_res * sps;                % High-resolution frequency grid
H1 = zeros(1, NN);                      % Frequency response initialized to zero
%% Passband
% Confused why we do (1 + NN/2)
% Fill
H1((NN/2 + 1) + (-NN/(2*sps) : NN/(2*sps))) = [0.5, ones(1, NN/sps - 1), 0.5];

% Getting alternate filter
M=alpha*400;

% Making M odd
if M-2*floor(M/2)==0
    M=M+1;
% Taper using kaiser
end
f1 = kaiser(M, 12);
f1 = f1 / sum(f1);                      % Scale
f2x = conv(H1, f1);
f2 = f2x(1+(M-1)/2:length(f2x)-(M-1)/2);
f3 = sqrt(f2);
f4 = fftshift(f3);
% Convert to time domain
f4 = 4 * fftshift((ifft(f4)));

f5 = f4((1+NN/2)+(-(N-1)/2:(N-1)/2));
f_real = real(f5);
if plot_flag == 1
    figure(1)
    plot(0:length(f_real)-1,f_real,'-o','linewidth',2)
    grid on
    title('SQRT Nyquist Filter Impulse Response')
    xlabel('Time Index')
    ylabel('Amplitude')
end
nyq_filt = f_real;
end