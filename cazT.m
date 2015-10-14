
                    sincvec=fft(sinc((0:NFFT2-1)-DelayToTarget(1)/ts2), NFFT2);
                    bbb1p=conj(sincvec([NFFT2+LIA+1:NFFT2,2:UIA+1])).'.*bbb1;
p=23;
N=-LIA+UIA;
gF=zeros(1,N);
switch rem(N,2)
    case 0
        for n=0:N-1
            gF(n+1)=exp(-1i*(2*pi/N)*p*n^2/2);
        end
        ML=-N/2;
        MU=N/2;
    case 1
        for n=0:N-1
            gF(n+1)=exp(-1i*(2*pi/N)*p*n*(n+1)/2);
        end
        ML=-(N+1)/2;
        MU=-ML-1;
end

GF2=[gF(end+LIA+1:end), gF(1:UIA)];


ML=0;
MU=-LIA+UIA-1;
N=-LIA+UIA;
indval=zeros(1,N);
indval(1,:)=(0:N-1).';
%for rdex=1:cazacj.NREf
for rdex=1:1
%g=[bbb2(-LIA+1:UIA-LIA);  bbb2(1:-LIA)];
g=bbb1p;
k=0;
for tau=ML:MU
    k=k+1;
    xax(k)=tau;
   % r(rdex,k)=(1/N)*sum(conj(gF.').*g([rem(tau+N,N)+1:N,1:rem(tau+N,N)]));
    r(rdex,k)=(1/N)*sum(conj(GF2.').*g([rem(tau+N,N)+1:N,1:rem(tau+N,N)]));

end
end
