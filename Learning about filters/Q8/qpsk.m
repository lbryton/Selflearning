% clc;
% close all;
% clf reset;
% clear;
function q = qpsk(N)
numSymbols = N;
bits = randi([0 1], numSymbols * 2, 1);
symbols = (2*bits(1:2:end) - 1) + 1j * (2*bits(2:2:end) - 1);

scaling = 4;
q = upsample(symbols, scaling);
end
% sqrt_nyq_y2(4,0.25,12,1);
% example_sqrt_nyq_y2(4,0.25,12,1);
