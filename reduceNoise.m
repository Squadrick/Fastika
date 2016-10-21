function [proc] = reduceNoise(f, fs)
N=size(f,1);
df = fs / N;
w = (-(N/2):(N/2)-1)*df;
n = 7;
beginFreq = 700 / (fs/2);
endFreq = 12000 / (fs/2);

[b,a] = butter(n, [beginFreq, endFreq]);
proc = filter(b,a,f);