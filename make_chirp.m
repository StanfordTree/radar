function chirp_spectrum = make_chirp(s, tau, fs, n_samples)

n_pts = fs*tau;
if mod(ceil(n_pts),2) ~= 0
    n_pts = floor(n_pts); 
else
    n_pts = ceil(n_pts);
end

for i = 1:n_pts
    t = (i-0.5*n_pts)/fs;
    phase = pi*s*(t^2);
    ref(i) = exp(complex(0, phase));
end

ref = [ref zeros(1, n_samples - n_pts)];
chirp_spectrum = circshift(ref, -0.5*n_pts);

end