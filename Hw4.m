clear all
clc

s = 4.189166e11;
tau = 37.12e-6;
fs = 18.96e6;

prf = 1679.9;
lam = 5.66e-2;
r0 = 830000;
v0 = 7550;
l = 10;
rad_E = 6378e3;
theta = deg2rad(23);
BW = s*tau;
fc = 0;


%Num of range bins
no_range_data = (10218-412)/2;

%n_valid calculation 
n = round(tau*fs);
n_valid = round(no_range_data - n);

%v_effective 
z = 770725; % calculated using wolfram from law of sines 
v_eff = v0*sqrt(rad_E/(rad_E + z)); %z 

%range_res = c/2(s*tau) = c/(2*BW)
slant_range_res = (3e8)/(2*BW)
slant_range_bin_spacing = 3e8/(2*fs)
rc = r0+(n_valid/2)*slant_range_bin_spacing;


beta = asin(rc*sin(theta)/rad_E);
i = theta + beta;
ground_range = slant_range_res/sin(i)

clc
nhdr = 412;
nsamp = 10218;
nlines = 10100;
num = (nsamp-nhdr)/2;
c = 3e8;

datafile = fopen('ersdata.hw3');
data_hdr = fread(datafile,[nsamp,nlines], 'uint8');
data = data_hdr(nhdr+1:end,:);

% Lets make chirp
chirp = makechirp(s,tau,fs,fc,1,4903);
chp_fft = fft(chirp);

%Convert data into real and complex signals
sig_even = data(2:2:end,:)-15.5; % Get the even component
sig_odd = data(1:2:end,:)-15.5; % Get the odd component
signal = sig_odd + 1i*sig_even; % Make the complex number
clear sig_even sig_odd

% signal = signal(1:4200,:);
signal_fft = fft(signal);
% range compress data
for i = 1:10100
    new_signal(:,i) = signal_fft(:,i).*conj(chp_fft).';
end
range_compressed = (ifft(new_signal,[],1));
other_dir = (fft(range_compressed,[],2));
% -308.97
f_centroid = -308.97;
img = zeros(4903,10100);
for j = 1:4903
    range = r0 + (j-1)*slant_range_bin_spacing;
    rd = sqrt(range^2 + (f_centroid*range*lam/(2*v_eff))^2);
    fr = -2*(v_eff^2)/(lam*rd);
    tau_az = rd*lam/(v_eff*l);
    ref = makechirp(fr, 0.8*tau_az, prf,f_centroid,1,10100);
    ref_fft = fft(ref);
    % plot(20*log10(abs(ref_fft)))

    img(j,:) = other_dir(j,:).*conj(ref_fft);
end

img = img(1:4200, :);
transformed = ifft(img, [], 2);

% for i = 1:10100
%     avg = mean(abs(transformed(:,i)));
%     transformed(:,i) = transformed(:,i)/avg;
% end
%%
close all
full_im = imagesc(abs(transformed)')
colormap gray
caxis([1000,60000])
xlabel('Range (m)')
ylabel('Azimuth (m)')
title('Processed SAR image')
%%
clc
k = 1;
clear sq_img
for j = 1:5:(10100-5)
    sum = sqrt(abs(transformed(:,j)).^2 + abs(transformed(:,j+1)).^2 + abs(transformed(:,j+2)).^2 + abs(transformed(:,j+3)).^2 + abs(transformed(:, j+4)).^2);
    sq_img(:,k) = sum;
    k = k+1;
end
%%
close all
full_im = imagesc(abs(sq_img)')
colormap gray
caxis([0,90000])

xlabel('Range (m)')

ylabel('Azimuth (m)')
title('Processed SAR image')



  