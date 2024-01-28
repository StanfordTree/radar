clear all
close all
clc

s = 4.189166e11;
tau = 37.12e-6;
fs = 18.96e6;
PRF = 1679.9e6;
r0 = 830000;
IQ_avg = 15.5;
p_vel = 7550;
lam = 0.0566;

n_samples = 10100/2;
data = fopen('ersdata.hw3');
data = fread(data,[10218,10100], 'uint8')';
data = data(:, 412:end);
imagesc(data)

chirp = make_chirp(s, tau, fs, n_samples);
chirp_FFT = get_FFT(chirp, n_samples);