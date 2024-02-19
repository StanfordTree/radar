clear all; close all; clc;

data = fopen('data.hw6');
signal = fread(data,[1750,20481],'int');
% 
% sig_even = mydat(2:2:end,:); % Get the even component
% sig_odd = mydat(1:2:end,:); % Get the odd component
% signal = sig_odd + 1i*sig_even;
signal = signal';
signal = signal(:,28:end);

imagesc(abs(signal))
xlabel('range'); ylabel('azimuth')


%%
clc; close all;

s = -5e11;
h = 696000;
prf = 2159.827;
fs = 16e6;
lam = 0.236057;

signal_fft = fft(signal,[],1);
% To get tau, want time bandwidth product around 100
tau = sqrt(100/abs(s));
fc = abs(s) + s*tau/2; %%%%%%%%%%%%%%%% MAYBE CHANGE
ref_fft = fft(makechirp(s,tau,fs,fc,1,1723));

range_compressed = zeros(size(signal_fft));
for i = 1:1723
    range_compressed(i,:) = signal_fft(i,:).*conj(ref_fft);
end

range_td = ifft(range_compressed,[],1);
imagesc(abs(range_td));