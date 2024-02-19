function chirp = makechirp(s,tau,fs,fc,start,n)
    %Function to compute chirp - reused in all problems
    %s: slope
    %tau: pulse length
    %fs: sample rate
    %fc: center frequency
    %start: starting index of chirp
    %n: the length of the chirp including zero
    dt=1/fs;
    npts=tau*fs;
    t=[-npts/2:npts/2-1]*dt;
    phase=pi*s*t.^2+2*pi*fc*t;
    chirp=[zeros(1,start-1) exp(1i*phase) zeros(1,n-length(phase)-start+1)];
end