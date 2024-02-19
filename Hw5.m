clear all; close all; clc;

s = 10^12;
tau = 10e-6;
fs = 100e6;
n = 2048;
fc = (s*tau)/2;
signal = makechirp(s,tau,fs,fc,1,n);
signal_fft = fft(signal);
%%  s = 1e12
s1 = 1e12;
fc_prime1_1 = (s1*tau)/4;
fc_prime1_2 = 3*(s1*tau)/4;
ref1_part1 = fft(makechirp(s1,tau,fs,fc_prime1_1,1,n));
ref1_part1 = ref1_part1(1:n/2);
ref1_part2 = fft(makechirp(s1,tau,fs,fc_prime1_2,1,n));
ref1_part2 = ref1_part2(n/2 + 1:n);

for i = 1:length(ref1_part1)
    correlated_sig1(i) = signal_fft(i)*conj(ref1_part1(i));
end

for i = 1:length(ref1_part2)
    correlated_sig1(n/2+i) = signal_fft(i+500)*conj(ref1_part2(i));
end

correlated_sig1 = abs(ifft(correlated_sig1));
%% s= 1.01e12
s2 = 1.01e12;
fc_prime2_1 = (s2*tau)/4;
fc_prime2_2 = 3*(s2*tau)/4;
ref2_part1 = fft(makechirp(s2,tau,fs,fc_prime2_1,1,n));
ref2_part1 = ref2_part1(1:n/2);
ref2_part2 = fft(makechirp(s2,tau,fs,fc_prime2_2,1,n));
ref2_part2 = ref2_part2(n/2 + 1:n);

for i = 1:length(ref2_part1)
    correlated_sig2(i) = signal_fft(i)*conj(ref2_part1(i));
end

for i = 1:length(ref2_part2)
    correlated_sig2(n/2+i) = signal_fft(i+500)*conj(ref2_part2(i));
end

correlated_sig2 = abs(ifft(correlated_sig2));
%% s= 1.03e12
s3 = 1.03e12;
fc_prime3_1 = (s3*tau)/4;
fc_prime3_2 = 3*(s3*tau)/4;
ref3_part1 = fft(makechirp(s3,tau,fs,fc_prime3_1,1,n));
ref3_part1 = ref3_part1(1:n/2);
ref3_part2 = fft(makechirp(s3,tau,fs,fc_prime3_2,1,n));
ref3_part2 = ref3_part2(n/2 + 1:end);

for i = 1:length(ref3_part1)
    correlated_sig3(i) = signal_fft(i)*conj(ref3_part1(i));
end

for i = 1:length(ref3_part2)
    correlated_sig3(n/2+i) = signal_fft(i+500)*conj(ref3_part2(i));
end

correlated_sig3 = abs(ifft(correlated_sig3));
%% s = 0.98e12;
s4 = 0.98e12;
fc_prime4_1 = (s4*tau)/4;
fc_prime4_2 = 3*(s4*tau)/4;
ref4_part1 = fft(makechirp(s4,tau,fs,fc_prime4_1,1,n));
ref4_part1 = ref4_part1(1:n/2);
ref4_part2 = fft(makechirp(s4,tau,fs,fc_prime4_2,1,n));
ref4_part2 = ref4_part2(n/2 + 1:end);

for i = 1:length(ref4_part1)
    correlated_sig4(i) = signal_fft(i)*conj(ref4_part1(i));
end

for i = 1:length(ref4_part2)
    correlated_sig4(n/2+i) = signal_fft(i+500)*conj(ref4_part2(i));
end

correlated_sig4 = abs(ifft(correlated_sig4));
%%
figure(1); hold on;
t = linspace(0, tau, length(correlated_sig4));
plot(t,correlated_sig1);
plot(t,correlated_sig2);
plot(t,correlated_sig3);
plot(t,correlated_sig4);
grid on
legend('s = 1e12','s = 1.01e12','s = 1.03e12','s = 0.98e12')
xlabel('Time')
ylabel('Power')
title('Autofocus program for 4 different chirp slopes')
%%
clc
idx_sig = n - 249;
[max1,idx1] = max(correlated_sig1);
[max2,idx2] = max(correlated_sig2);
[max3,idx3] = max(correlated_sig3);
[max4,idx4] = max(correlated_sig4);

dt_1 = (idx_sig - idx4)/fs;
dS = 2*dt_1*(s4^2)/(s4*tau);
%%
clear all; close all; clc;
data = fopen('simlband.dat');
mydat = fread(data,[4096,2048],'float');

%range is columnsss

sig_even = mydat(2:2:end,:); % Get the even component
sig_odd = mydat(1:2:end,:); % Get the odd component
signal = sig_odd + 1i*sig_even;

signal = signal';
% imagesc(abs(signal))
% xlabel('range'); ylabel('azimuth')
s = 10^12;
tau = 10e-6;
fs = 24e6;


signal_fft = fft(signal,[],2);
ref_fft = fft(makechirp(s,tau,fs,0,1,2048));
%range compress
for i = 1:2048
    range_compressed(i,:) = signal_fft(i,:).*(ref_fft);
end
range_td = ifft(range_compressed,[],2);
imagesc(abs(range_td));
xlabel('Range'); ylabel('Azimuth')
prf = 250; v = 250; l = 2; lam = 0.25; r0 = 4653;
%%
fd = get_fd(2048,mydat,ref_fft,prf)
az_fft = fft(range_td, [], 1);
imagesc(abs(az_fft))
xlabel('Range'); ylabel('Azimuth')

%%
avg_az = fftshift(mean(abs(az_fft),2));
figure(1)
freq = linspace(-prf/2, prf/2, length(avg_az));
plot(freq, avg_az)
xlabel('Freq');ylabel('Power')
%%
clc; close all;
clear az_ref_fft
fd = -111.3; %hz
dR = 3e8/(2*fs);
img = zeros(2048,2048);
% 
% for i = 1:2048
%     range = r0 + (i-1)*dR;
%     rd = sqrt(range^2 + (fd*range*lam/(2*v))^2);
%     fr = -2*(v^2)/(lam*rd);
%     tau_az = rd*lam/(v*l);
%     az_ref = makechirp(fr, 0.8*tau_az, prf,fd,1,2048);
%     az_ref_fft = fft(az_ref);
%     img(:,i) = az_fft(:,i).*conj(az_ref_fft).';
% end
% 
% transformed = ifft(img, [], 1);
% full_im = imagesc(abs(transformed))
for i = 1:2048
    range = r0 + (i-1)*dR;
    rd = sqrt(range^2 + (fd*range*lam/(2*v))^2);
    fr = -2*(v^2)/(lam*rd);
    tau_az = rd*lam/(v*l);
    az_ref_fft = makechirp(fr,0.8*tau_az,prf,fd,1,2048);
    az_ref_fft = fft(az_ref_fft);
    img(:,i) = az_fft(:,i).*(az_ref_fft).';
end 

transformed = ifft(img,[],1);
imagesc(abs(transformed))
xlabel('Range');ylabel('Azimuth')
%% Range migration
clc
clear new_az_fft;
f = linspace(-prf/2 +fd,prf/2 + fd ,2048);
for i = 1:length(f)
    range = r0 + (i-1)*dR;
    deltaF(i) = (f(i)^2 - fd^2)*(lam^2*range/(8*v^2));
end
deltaF = round(deltaF/dR);
% plot(deltaF);
new_az_fft = az_fft;
new_az_fft(1:924,:) = circshift(az_fft(1:924,:),1,2);
new_az_fft(925:1105,:) = circshift(az_fft(925:1105,:),0,2);
new_az_fft(1106:1342,:) = circshift(az_fft(1106:1342,:),-1,2);
new_az_fft(1343:1512,:) = circshift(az_fft(1343:1512,:),-2,2);
new_az_fft(1513:1710,:) = circshift(az_fft(1513:1710,:),-3,2);
new_az_fft(1710:1820,:) = circshift(az_fft(1710:1820,:),-4,2);
new_az_fft(1821:end,:) = circshift(az_fft(1821:end,:),-5,2);

%%

clc; close all;
fd = -111.3; %hz
dR = 3e8/(2*fs);
new_img = zeros(2048,2048);

for i = 1:2048
    range = r0 + (i-1)*dR;
    rd = sqrt(range^2 + (fd*range*lam/(2*v))^2);
    fr = -2*(v^2)/(lam*rd);
    tau_az = rd*lam/(v*l);
    az_ref_fft = makechirp(fr,0.8*tau_az,prf,fd,1,2048);
    az_ref_fft = fft(az_ref_fft);
    new_img(:,i) = new_az_fft(:,i).*(az_ref_fft).';
end 

newtransformed = ifft(new_img,[],1);
imagesc(abs(newtransformed))
xlabel('Range');ylabel('Azimuth')