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
new_signal_td = (ifft(new_signal,[],1));
img = imagesc(abs(new_signal_td))
%% Range compress again in a different way
az_spec = fft(new_signal_td,[],2);
avg_az = mean(abs(az_spec),1);

%% Plotting for 1a
clc
figure(1)
az_avg = fftshift(abs(avg_az));
freq = linspace(-PRF/2, PRF/2, length(az_avg));
plot(freq, az_avg)

xlabel('PRF')
title('Azimuth spectrum averaged over range')
ylabel('Power')



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

%% Q2
f_centroid = -308.97;
clc
close all

for i = 1:size(new_signal_td,1)
    for j = 1:64
        adj_single_look(i,j) = new_signal_td(i,j)*exp(complex(0,2*pi*f_centroid-(j/PRF)));
    end
end
az_fft_single_look = fft(adj_single_look,[],2);

img = imagesc(abs(az_fft_single_look)');
xlabel('Azimuth')
ylabel('Range')
title('Single look azimuth processed image')
%%
clc

clear az_single n image_mat az_multi_look
tic
image_mat = zeros(size(new_signal_td));
counter = 2;
n = 2;
j = 1;
l = 0;
while (64*l + 64) <= 10100
    
    n = round((counter-1)*3.51);
    shift_no = n;
    for i = 1:size(new_signal_td,1)
        k = 1;
        for j = (64*l + 1):(64*l + 64)
            % we generate one range compressed image...not yet fft'd in
            % azimuth
            az_single(i,k) = new_signal_td(i,j)*exp(complex(0,2*pi*f_centroid-(j/PRF)));
            k = k+1;
        end
    end
     l = l+1;
     az_multi_look = fft(az_single,[],2);  %%fft'ing in azimuth dir. 
     image_mat(:,n:n+63) = image_mat(:,n:n+63) + az_multi_look;
     
     % image_mat = image_mat + comb_mat;
     counter = counter+1;
     az_single = zeros(4903,64);
     % comb_mat = zeros(4903,10100);
end
toc

%%
close all
image_mat = image_mat(:,1:613);

full_im = imagesc(2*abs(image_mat)')
colorbar
clim([1e4 0.8e5])
xlabel('Range')
ylabel('Elevation')
title('Unfocused range-dopplar processed image')

