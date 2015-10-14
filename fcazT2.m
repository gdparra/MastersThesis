function [SNRout,r]=fcazT(GF2, DelayEst, bbb1, LIA, UIA, NFFT2,ts2)

sincvec=fft(sinc((0:NFFT2-1)-DelayEst/ts2), NFFT2);
 bbb1p=conj(sincvec([NFFT2+LIA+1:NFFT2,2:UIA+1])).'.*bbb1;



ML=0;
MU=-LIA+UIA-1;
N=-LIA+UIA;
indval=zeros(1,N);
indval(1,:)=(0:N-1).';
%for rdex=1:cazacj.NREf
%g=[bbb2(-LIA+1:UIA-LIA);  bbb2(1:-LIA)];
g=bbb1p;
k=0;
for tau=ML:MU
    k=k+1;
    xax(k)=tau;
   % r(rdex,k)=(1/N)*sum(conj(gF.').*g([rem(tau+N,N)+1:N,1:rem(tau+N,N)]));
    r(k)=(1/N)*sum(conj(GF2.').*g([rem(tau+N,N)+1:N,1:rem(tau+N,N)]));

end
SNRout=r(1)*conj(r(1))/sum(r(2:end).*conj(r(2:end)));
end
