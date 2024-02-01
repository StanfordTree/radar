clear all
close all
clc

nhdr = 412;
nsamp = 10218;
nlines = 10100;

s = 4.189166e11;
tau = 37.12e-6;
fs = 18.96e6;
PRF = 1679.9;
r0 = 830000;
IQ_avg = 15.5;
p_vel = 7550;
lam = 0.0566;
fc = 0;
num = (nsamp-nhdr)/2;

chirp = makechirp(s,tau,fs,fc,1,num);
chp_fft = fft(chirp); % transform to frequency domain

datafile = fopen('ersdata.hw3');
data_hdr = fread(datafile,[nsamp,nlines], 'uint8');
data = data_hdr(nhdr+1:end,:);
%%
sig_even = data(2:2:end,:)-15.5; % Get the even component
sig_odd = data(1:2:end,:)-15.5; % Get the odd component
signal = sig_odd + 1i*sig_even; % Make the complex number
clear sig_even sig_odd

%%

signal_fft = fft(signal); % fft each azimuth line (column)
avg_signal = mean(abs(signal_fft),2); % Take the average of the spectra magnitudes
magnitude = 20*log10(avg_signal + 1e-30);
freq = linspace(-fs/2,fs/2,num); % Frequency axis
figure
plot(freq,fftshift(magnitude));
title('Average dB Magnitude Plot of Average Raw Data Spectrum');
xlabel('frequency (Hz)'); ylabel('Magnitude (dB)');
grid on

%% Range compress
for i = 1:10100
    new_signal(:,i) = signal_fft(:,i).*conj(chp_fft).';
end
%%
new_signal_td = (ifft(new_signal,[],1));
% new_signal_td = abs(new_signal_td);
img = imagesc(abs(new_signal_td))
%% Range compress again in a different way
% signal_spec = fft(curr_line,[], 1);
% signal_td = ifft(signal_spec.*conj(chirp_FFT));
az_spec = fft(new_signal_td,[],2);
avg_az = mean(abs(az_spec),1);
%%

%% Plotting for 1a
clc
figure(1)
az_avg = fftshift(abs(avg_az));
freq = linspace(-PRF/2, PRF/2, length(az_avg));
plot(freq, az_avg)

xlabel('PRF')
title('Azimuth spectrum averaged over range')
ylabel('Power')

%%
%f_centroid = 1370.9;


%% 1b

clc
sum = 0;
for k = 1:4903 % k = bin number
    for i = 2:length(chirp)
        curr_sum = signal(i,k)*conj(signal(i-1,k));
        sum = sum + curr_sum;
    end
    sum_per_bin(k) = sum;
end

phase = angle(sum_per_bin);
centroid = (PRF)*(mean(phase)/(2*pi));
%%