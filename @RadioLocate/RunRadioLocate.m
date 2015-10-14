function [varargout] = RunCazacDespread(varargin)
dbstop if error

cazacj=varargin{1};
VecLongT=cazacj.VecLong;
%spn=37;
spn=cazacj.synchpt;
fbins=11;
veclongF=zeros(cazacj.NREf, cazacj.Nfft);
veclongFH=zeros(cazacj.NREf, cazacj.Nfft);

for k=0:cazacj.Nsymb-1
vectmp=VecLongT(spn+k*(cazacj.Nfft+cazacj.CP)+1:spn+k*(cazacj.Nfft+cazacj.CP)+cazacj.Nfft);
tmpF=fft(vectmp);
veclongF(1:cazacj.NREf,k+1)=tmpF(217:227);
veclongFH(1:cazacj.NREf,k+1)=tmpF(287:297);
end

ML=cazacj.ML;
MU=cazacj.MU;
N=cazacj.Nsymb;
indval=zeros(1,N);
indval(1,:)=[0:N-1].'
for rdex=1:cazacj.NREf
g=veclongF(rdex,:);
k=0;
for tau=ML:MU
    k=k+1;
    xax(k)=tau;
    r(rdex,k)=(1/N)*sum(conj(g).*g([rem(tau+N,N)+1:N,1:rem(tau+N,N)]));
end
end
varargout{1}=r;
varargout{2}=veclongF;
varargout{3}=veclongFH;
        
rdex=1;
switch varargin{2}
    case 'PlotAll'
figure(8)
subplot(2,1,1), plot(indval,real(r(rdex,:)));
xlabel('index');
ylabel('Real Amplitude');
titleSTR=strcat({'CAZAC Sequence Correlation,'},{', N='},{num2str(N)});
title(titleSTR);
subplot(2,1,2), plot(indval,imag(r(rdex,:)));
xlabel('index');
ylabel('Imag Amplitude');
    case 'NullPlot'
   
    otherwise
        disp('Not Valid Option');
end

    
    

       

        
     
        
        