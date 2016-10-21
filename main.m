[x1, Fs1, bits1] = wavread('mic1.wav');
[x2, Fs2, bits2] = wavread('mic2.wav');

x = [x1';x2'];

p = fastica(x);

p1 = p(1,:)';
p2 = p(2,:)';

wavwrite(p1,Fs1,bits1,'out1.wav');
wavwrite(p2,Fs2,bits2,'out2.wav');