nvar=1e-5;
CP=3;
x = [1, 3, -1, 7, 2];  %freq
xtemp = ifft(x);
xofdm = [xtemp(end-CP+1:end), xtemp];  %time 
h= [1, -2, 0.5];

y = conv(xofdm, h);
v = sqrt(nvar)*randn(size(y));
z=y+v;
ptrval=CP+1;
zfreq=fft(z(ptrval:ptrval-1+length(x)),length(x));

Dp=(CP+1)-ptrval;
zest=x.*fft(h,length(x)).*exp(-1i*(2*pi*(0:4)/length(x))*Dp);
[zfreq; zest]
chanout=ChanEst(xofdm,z,3);

