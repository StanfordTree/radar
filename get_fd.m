function fd = get_fd(n,signal,chirp,PRF)

sum = 0;
for k = 1:n % k = bin number
    for i = 2:length(chirp)
        curr_sum = signal(i,k)*conj(signal(i-1,k));
        sum = sum + curr_sum;
    end
    sum_per_bin(k) = sum;
end

phase = angle(sum_per_bin);
fd = (PRF)*(mean(phase)/(2*pi));
end