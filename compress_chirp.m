function signal = compress_chirp(Y, n_samples)
for i = 1:n_samples
    db(i) = 20*log10(abs(Y(i))+1e-30);
end
signal = dB;

end