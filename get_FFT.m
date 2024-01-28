function FFT_chirp = get_FFT(ref, n_samples)
Y = fft(ref, n_samples);
FFT_chirp = fftshift(Y);
end