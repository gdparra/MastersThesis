%%%%%%%%%%%%%%%%%%%%%%%%%
%BasebandGen Script  EE4653
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%
dbstop if error  % stop if an error occurs 
clear all  %  clears work space between runs
rcrolloff=0.4;   %rolloff Factor for nyquist pulse
symbolrate=128e3;  %desired symbol rate
NumSamplesTxPulse=13; %Warning the program can modify this value!
SampleLaunchPeriod=8;  % period in samples between symbol launches
EbNodB=8;
%%%%%%%%%%%%%%%%%
NumBitsPerSymbol=2;  %2-Tupple bits
binary_sequence=[1,0,1,0,0,0,0,0,1,1,1,1,1,0,0,0,1,1,1,0,0,1,0,1];
NumSymbs=length(binary_sequence)/NumBitsPerSymbol;  %Number of symbols 
b=zeros(1, NumSymbs);
for k=1:NumSymbs
    %n-tupple: (e.g. 2-tupple b1,b0 --> 2*b1+b0 mapping to complex coeff.)
    k0k1=binary_sequence((k-1)*NumBitsPerSymbol+1:k*NumBitsPerSymbol)*[2,1]';
   switch k0k1;
% switch binary_sequence((k-1)*NumBitsPerSymbol+1:k*NumBitsPerSymbol)*[2,1]';

        % QPSK Mapping Defined
        case 0
            b(k)=exp(1i*((2*pi/4)*0+pi/4));
        case 1
            b(k)=exp(1i*((2*pi/4)*1+pi/4));
        case 3
            b(k)=exp(1i*((2*pi/4)*2+pi/4));
        case 2
            b(k)=exp(1i*((2*pi/4)*3+pi/4));
    end
end
%%%%%%%%%

%%%%%%%%%%%%%%%%%

Es=1;
p.rcrolloff=rcrolloff;
p.symbolrate=symbolrate;
p.NumSamplesTxPulse=NumSamplesTxPulse;
p.SampleLaunchPeriod=SampleLaunchPeriod;
a=BasebandGenNew(p, Es,b);

%%%%%%%%%%%%%%%%%%%%%%%%%%

%FOR BASEBAND MODE USE the AWGN Channel
EbNoLin=10^(EbNodB/10);
vn=(1/2)/EbNoLin;
yawgnbase=a.basebandsig+sqrt(vn/2)*randn(size(a.basebandsig))+1i*sqrt(vn/2)*randn(size(a.basebandsig));

figure(6);
xax1=a.Ts*(0:length(yawgnbase)-1);
subplot(2,1,1), plot(xax1, real(yawgnbase));
xlabel('sample time (sec)');
ylabel('Real');
titletext=strcat({'Baseband Signal in AWGN, EbNo='},{num2str(EbNodB)},...
    {'dB, SamplingFrequency= '}, {num2str(1/a.Ts)});
title(titletext);
subplot(2,1,2), plot(xax1, imag(yawgnbase));
xlabel('sample time (sec)');
ylabel('Imag');
%%%%%%%%%%%%%%%%%


figure(7);
fs=1/a.Ts;
NumPointsFFT=2048;
xax = ([0:NumPointsFFT-1]/NumPointsFFT)*fs;
yaxdB=20*log10(abs(fft(a.basebandsig,NumPointsFFT)));
plot([xax-xax(end)/2], [yaxdB(NumPointsFFT/2+1:NumPointsFFT),yaxdB(1:NumPointsFFT/2)], 'b');
titletext=strcat({'FFT of 8-PSK Signal, Rsymb='},{num2str(a.symbolrate)},{'Pulse Rolloff= '}, {num2str(a.rcrolloff)},...
    {' SamplingRate= '}, {num2str(fs)});
xtext='Frequency in Hz';
ytext='Baseband Signal Power in dB';
xlabel(xtext);
ylabel(ytext);
title(titletext);




