function [SNRout,r]=fcazT(RMAT, GF2, DelayEst, bbb1, LIA, UIA, NFFT2,ts2)

sincvec=fft(sinc((0:NFFT2-1)-DelayEst/ts2), NFFT2);
 bbb1p=conj(sincvec([NFFT2+LIA+1:NFFT2,2:UIA+1])).'.*bbb1;

g=bbb1p;
r=g.'*RMAT.';
SNRout=r(1)*conj(r(1))/sum(r(2:end).*conj(r(2:end)));
end
