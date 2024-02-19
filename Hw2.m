clear all
close all
clc

s = 10^12; 
tau = 10^-5; %length of chirp
fs = 10^8;  %sample rate

% step 1 : create a chirp signal
n_pts = fs*tau;
for i = 1:(n_pts+1)
    t = (i-500)/fs;
    phase = pi*s*(t^2);
    ref(i) = exp(complex(0,phase));
end
ref = [ref zeros(1,2048-1000)];
ref = circshift(ref, -500);

%step 2: calculate transform (going from time domain to frequency domain)
Y = fft(ref,2048);

Y = fftshift(Y);

%step 3: plot the power vs frequency 
% add 1e-30 for log purposes 
for i = 1:2048
    dB(i) = 20*log10(abs(Y(i))+1e-30);
end

%frequency s*tau is the bandwidth within the time the signal is sent 
freq = linspace(0, fs, 2048);

%plotting 
plot(freq,dB)
xlabel('Frequency')
ylabel('dB')
title('Chirp spectrum')





s = 1.01*10^12;
for i = 1:(n_pts+1)
    t = (i-500)/fs;
    phase = pi*s*(t^2);
    newref(i) = exp(complex(0,phase));
end
newref = [newref zeros(1,2048-1000)];
newref = circshift(newref, -500);

%step 2: calculate transform (going from time domain to frequency domain)
newref = fft(newref,2048);

newref = fftshift(newref);


%1a
for i = 1:2048
    signal(i) = Y(i)*conj(newref(i));
end

signal = fftshift(ifft(signal));
signal = 20*log10(abs(signal));
time = linspace(0, tau, 2048);
figure(2)
plot(time,signal)
xlabel('Time (s)')
ylabel('Power (dB)')
title('Compressed chirp for s = 1.03*10^{12} Hz')
%% Q2
clear all
close all
clc


s = 10^12; 
tau = 10^-5; %length of chirp
fs = 10^8;  %sample rate


%%%%%%%%%%%%%%%%%%%%%%%%% ref %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 1 : create a chirp signal
n_pts = fs*tau;
for i = 1:(n_pts+1)
    t = (i-500)/fs;
    phase = pi*s*(t^2);
    ref(i) = exp(complex(0,phase));
end
ref = [ref zeros(1,2048-1000)];
ref = circshift(ref, -500);

plot(1:2048, ref)

%step 2: calculate transform (going from time domain to frequency domain)
Y = fftshift(fft(ref,2048));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:1000
    t = (i-500)/fs;
    mixed_ref = pi*s*(t^2);
    chirp1(i) = exp(complex(0,mixed_ref));
    chirp2(i) = 5*exp(complex(0,mixed_ref));
    chirp3(i) = 2*exp(complex(0,mixed_ref));

end
chirp1 = [chirp1 zeros(1,2048-1000)];
chirp2 = [chirp2 zeros(1,2048-1000)];
chirp3 = [chirp3 zeros(1,2048-1000)];

newchrip1 = circshift(chirp1, -400);
newchirp2 = circshift(chirp2, -100);

mixed_chirp = newchrip1 + newchirp2 + chirp3;


figure(1)
plot(linspace(0,tau,2048), mixed_chirp)
xlabel('Time (s)')
ylabel('Amplitude')
title('Mixed chirp amplitude over time')


mixed_chirp_Y = fft(mixed_chirp,2048);
mixed_chirp_Y = fftshift(mixed_chirp_Y);

mixed_signal = mixed_chirp_Y.*conj(Y);

mixed_signal = (ifft(mixed_signal));

figure(4) 
plot(linspace(0,tau,2048),mixed_signal)
xlabel('Time (s)')
ylabel('Power')
title('Compresed signal')
%% Q3
clear all
close all
clc

data = fopen('ersdata');
data = fread(data,[10218,1024],'uint8').';
imagesc(data)
% xlabel('Samples')
% ylabel('Rows')
% title('ERS data image')
%3.1
s = 4.189166 * 10^11;
tau = 37.12e-6;
fs = 18.96*10^6;
n_pts = ceil(fs*tau);

for i = 1:n_pts
    t = (i-0.5*n_pts)/fs;
    phase = pi*s*(t^2);
    ref(i) = exp(complex(0,phase));
end
ref = [ref zeros(1,4903-n_pts)];
ref = circshift(ref, -n_pts*0.5);

Y = (fft(ref,4903));
dB = 20*log10(abs(Y)+1e-30);
freq = linspace(0, fs, length(dB));

figure(1)
plot(freq,dB)
xlabel('Frequency')
ylabel('dB')
title('Chirp spectrum for 3b')

%3c

data = data(:,413:end);
for i = 1:1024
    for k = 1:4903
        curr_line(k) = complex(data(i,2*k-1)-15.5, data(i,2*k)-15.5);
    end
    sig_td(i,:) = curr_line; 
    signal(i,:) = (fft(curr_line));
end

figure(2)
avg_val = mean(signal, 1);
plot(linspace(0, s*tau, 4903), abs(avg_val))
title('Averaged specta in the range direction')
xlabel('Frequency')
ylabel('Power (linear)')
%% %3d
close all
clc
for i = 1:1024
    new_signal(i,:) = signal(i,:).*conj(Y);
end

new_signal_td = fftshift(ifft(new_signal));
new_signal_td = abs(new_signal_td);
figure(1)
imagesc(new_signal_td)
xlabel('range')
ylabel('azimuth')
title('Range compressed image')



%% Q4
clear all
close all
clc

s = 10^11;
tau = 30*10^-6;
fs = 20e6;
fc = 10e6;

n_pts = fs*tau;
for i = 1:n_pts
    t = i/fs;
    phase = pi*s*(t^2) + 2*pi*fc*t; %account for new center freq
    ref(i) = exp(complex(0,phase));
end

freq = linspace(0, tau, n_pts);

Y = fft(ref);

for i = 1:length(Y)
    signal(i) = Y(i) * conj(Y(i));
end


impulse = fftshift(ifft(signal));
signal = 20*log10(abs(impulse)+1e-30);

for i = 1:length(signal)
    if signal(i) < -10
        signal(i) = 0;
    end
end
plot(freq, signal)
xlabel('Time (s)')
ylabel('Power (dB)')
title('Impulse response of I/Q system')


%double freq, take fft, take only positive side to get the real signal 
fs = 2*fs;
n_pts = fs*tau;

for i = 1:n_pts
    t = i/fs;
    phase = pi*s*(t^2) + 2*pi*fc*t; %account for new center freq
    ref(i) = exp(complex(0,phase));
end
real_ref = real(ref);
Y = fftshift(fft(real_ref));

newY = Y(600:end);
for i = 1:length(newY)
    signal(i) = newY(i) * conj(newY(i));
end

freq = linspace(0, tau, length(newY));

impulse = fftshift(ifft(signal));
signal = 20*log10(abs(impulse)+1e-30);
figure(2)
plot(freq, signal)
xlabel('Time (s)')
ylabel('Power (dB)')
title('Impulse response of Offset video system')

%% Q5
clear all
close all
clc

s = 10^11;
tau = 30*10^-6;
fs = 20e6;
fc = 10e6;

n_pts = fs*tau;
for i = 1:n_pts
    t = i/fs;
    phase = pi*s*(t^2) + 2*pi*fc*t; %account for new center freq
    ref(i) = exp(complex(0,phase));
end

freq = linspace(0, tau, n_pts)*fs;

w = (0.4:0.05:1);
for i = 1:length(w)
    weight = w(i) + (1-w(i))*cos(freq);
    new_ref = weight.*ref;

    weigthed_Y = fftshift(fft(new_ref));
    for j = 1:length(weigthed_Y)
        weighted_signal(j) = weigthed_Y(j) * conj(weigthed_Y(j));
    end

    impulse = fftshift(ifft(weighted_signal));
    weighted_signal = 20*log10(abs(impulse));

    peaks = findpeaks(weighted_signal); 
    [~,max_idx] = max(peaks);
    pslr_main = max(peaks);
    pslr_side = peaks(max_idx+1);

    pslr_val(i) = pslr_side-pslr_main;

    impulse = abs(impulse).^2;
    sidelobe = trapz(impulse(312:end));
    mainlobe = trapz(impulse(300:311));

    ISLR_val(i) = abs(sidelobe/mainlobe);

end
plot(w,pslr_val)  
hold on
plot(w,10*log10(abs(ISLR_val)))
xlabel('w')
ylabel('dB')
title('Peak and integrated sidelobe levels for w ranging from 0.4->1')

% 
% % Y = fft(ref);

% for i = 1:length(Y)
%     signal(i) = Y(i) * conj(Y(i));
% end
% impulse = fftshift(ifft(signal));
% signal = 20*log10(abs(impulse)+1e-30);
% 
% pslr_side = signal(311);
% pslr_main = max(signal);
% 
% pslr_val = pslr_side-pslr_main;
% 
% impulse = (abs(impulse)).^2;
% 
% sidelobe = trapz(abs(impulse(309:end)));
% mainlobe = trapz(abs(impulse(301:308)));
% 
% ISLR_val = abs(sidelobe/mainlobe);
% plot(freq, signal)
% xlabel('Time (s)')
% ylabel('Power (dB)')
% title('Impulse response of I/Q system')


%%
y = chirp(linspace(0,tau,2048),0,tau,s*tau);
z = fft(y);
z = fftshift(z);
plot(linspace(0,fs,2048), 20*log10(abs(z)+1e-30))