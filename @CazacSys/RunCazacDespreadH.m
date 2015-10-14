function [varargout] = RunCazacDespreadH(varargin)
dbstop if error



HQFREQ=varargin{3};
cazacj=varargin{1};
gr=cazacj.vec;
ML=cazacj.ML;
MU=cazacj.MU;
N=cazacj.Nsymb;
spn=cazacj.CP;
%VecLongT=cazacj.VecLong;
fbins=11;
NFFT=cazacj.Nfft;
veclongF=zeros(cazacj.NREf, cazacj.Nsymb);
%veclongFH=zeros(cazacj.NREf, cazacj.Nfft);
HEST=zeros(cazacj.NREf,1);
GEST=conj(gr.');
coherent = varargin{4};
%cazacj.VecLong=freq2time(cazacj);

switch coherent
    case 'yes'
        %%%%%%%%%%%%%Channel Estimate%%%%%%%%%%%%%%%%%%%%%
 %       for k=0:cazacj.Nsymb-1
            vectmp=VecLongT(spn+1:spn+cazacj.Nfft);
            tmpF=fft(vectmp);
            veclongF(1:cazacj.NREf)=[tmpF(NFFT-4:NFFT,1); tmpF(1:6,1)];
            %veclongFH=[];
  %      end
        HEST=veclongF*GEST;
        HEST=HEST/cazacj.Nsymb;
        dbg77=1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k=0:cazacj.Nsymb-1
            vectmp=VecLongT(spn+k*(cazacj.Nfft+cazacj.CP)+1:spn+k*(cazacj.Nfft+cazacj.CP)+cazacj.Nfft);
            tmpF=fft(vectmp);
            veclongF(1:cazacj.NREf,k+1)=[tmpF(NFFT-4:NFFT,1); tmpF(1:6,1)]./HEST;
            veclongFH=[];
        end
    case 'no'
        for k=0:cazacj.Nsymb-1
            vectmp=VecLongT(spn+k*(cazacj.Nfft+cazacj.CP)+1:spn+k*(cazacj.Nfft+cazacj.CP)+cazacj.Nfft);
            tmpF=fft(vectmp);
            veclongF(1:cazacj.NREf,k+1)=[tmpF(NFFT-4:NFFT,1); tmpF(1:6,1)];
            veclongFH=[];
        end
    otherwise
        disp('invalid option');
        
end


indval=zeros(1,N);
indval(1,:)=[0:N-1].'
for rdex=1:cazacj.NREf
    g=veclongF(rdex,:);
    k=0;
    for tau=ML:MU
        k=k+1;
        xax(k)=tau;
        r(rdex,k)=(1/N)*sum(conj(gr).*g([rem(tau+N,N)+1:N,1:rem(tau+N,N)]));
        if tau==0
            dbg77=1;
        end
    end
end
varargout{1}=r;
varargout{2}=veclongF;
varargout{3}=veclongFH;

rdex=2;
switch varargin{2}
    case 'PlotAll'
        figure(8)
        subplot(2,1,1), plot(indval,real(r(rdex,:)));
        xlabel('index');
        ylabel('Real Amplitude');
        titleSTR=strcat({'CAZAC Sequence Correlation'},{', N='},{num2str(N)});
        title(titleSTR);
        subplot(2,1,2), plot(indval,imag(r(rdex,:)));
        xlabel('index');
        ylabel('Imag Amplitude');
    case 'NullPlot'
        
    otherwise
        disp('Not Valid Option');
end









