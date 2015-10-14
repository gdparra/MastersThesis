function [varargout] = InterfereGenScriptChanOnly(same_pdp,B,fldhomecombfig3,fldhomecombfig4,length_riceanh,loc4g,int_ratio,mm,TowInd)
%channel is defined at rate fs
dbstop if error  % stop if an error occurs
%clear all  %  clears work space between runs
fs=7.68e6;
Nfftsize=512;
Mchan=1;
%Line of site gain
tau0=(11*1/fs);
NLoop=1000;
yvec=zeros(NLoop,1);
channelgain=zeros(NLoop,1);
%for loopK=1:NLoop

%this is the average value over a suite of simulation trials (varies per trial)
interloc.riciankdB=8;

%scattering gain (power delay profile)
interloc.pdp    =    [  -8,   -6,  -4 ,   0,   0,   -4 ,     -8,     -9,   -10 ,  -12 ,  -14  ]; %dB scattering response!
interloc.PathDelays= [ 0.1,  0.3,  0.5, 0.7, 1.0,   1.3,   15.0,   15.2,   15.7,  17.2,   20.0]*0.2e-6;  %in seconds
interloc.ts=1/fs;

%not moving (zero velocity)
interloc.dopp=0;
%%%%%%%%%%%%%%%%%%%%%
%New Approach
% figure(2)
% hold on

for Tx=1:Mchan
    for Rx=1:Mchan
%         index = randperm(numel(interloc.pdp));
%         B = interloc.pdp(index);
%         same_pdp='no';
        
        switch same_pdp
            case 'yes'
                B=B;
            case 'no'
                index = randperm(numel(interloc.pdp));
                B = interloc.pdp(index);
        end
        
        
        [interloc.h00t] = GenRayleigh(B, interloc.PathDelays, interloc.ts,  interloc.dopp);
        %interloc.h00t=[0.102851428058823 - 0.254496286105682i,-0.0373485103931939 + 0.361391182791193i,0.0296117040699562 + 0.106418746549803i,-0.216334888506520 - 0.0308339681712780i,0.155918690292763 - 0.0353107544934772i,-0.618670913034324 - 0.0988670529293303i,0.0823037514768995 - 0.0548247778307555i,-0.121937988346907 + 0.00331107107926077i,-0.0382189044530078 + 0.00585350152982686i,0.0501112197893314 + 0.0748846982163497i,0.297874593322162 + 0.437736166752035i];
        %%%%%%%%%%%%%%%
        %forced_freq_H=fft(forced_time_h,length_riceanh);
        %forced_time_h=ifft(forced_freq_H,length_riceanh);
        %%%%%%%%%%%%%%%
        %convert the continuous rician channel to discretely sampled
        interloc.h00N0=MultiRes(interloc.h00t, interloc.PathDelays, interloc.ts, 'NonLinear');
        
        interloc.h00N = interloc.h00N0/sqrt(sum(interloc.h00N0.*conj(interloc.h00N0)));
        %interloc.h00N=forced_time_h;
        
        temph00N=fft(interloc.h00N,length_riceanh);
        interloc.h00N=ifft(temph00N,length_riceanh);
        
        
        radiochan(1:length(interloc.h00N),Rx+(Tx-1)*Mchan)=interloc.h00N;
%        plot(20*log10(abs(fft(radiochan(1:length(interloc.h00N),Rx+(Tx-1)*Mchan),Nfftsize))));
    end
end
% plot(20*log10(abs(fft(radiochan(1:3,1),64))));
%hold off

for k=1:Mchan^2
    [radiochan(1:3,k).'];
end
%Scattering
%print first five values

channelgain= sum(interloc.h00N.*conj(interloc.h00N));


%%%%%%%%%%%%%%%%%%%%%%%%%
%BasebandGen Script
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%

rcrolloff=0.5;   %rolloff Factor for nyquist pulse
symbolrate=960e3;  %desired symbol rate
NumSamplesTxPulse=13; %Warning the program can modify this value!
SampleLaunchPeriod=8;  % period in samples between symbol launches
EbNodB=8;
%%%%%%%%%%%%%%%%%
NumBitsPerSymbol=2;  %2-Tupple bits
%binary_sequence=[1,0,1,0,0,0,0,0,1,1,1,1,1,0,0,0,1,1,1,0,0,1,0,1];
binary_sequence=randi([0:1],1,24);
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
p.fldhomecombfig3=fldhomecombfig3;
p.fldhomecombfig4=fldhomecombfig4;
p.mm=mm;
p.TowInd=TowInd;
p.same_pdp=same_pdp;
a=BasebandGenNew(p, Es,b);

%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_freq_int=fft(a.txpulse,loc4g.Nfft);
temp_time_int=ifft(temp_freq_int,loc4g.Nfft);
% int_Tsymbol=zeros(1,loc4g.Nfft+loc4g.CP);
% int_Tsymbol=[temp_time_int(end-loc4g.CP+1:end),temp_time_int];
% a.txpulse=int_Tsymbol;
a.txpulse=temp_time_int;
%%%%%%%%%%%%%%%%%%%%%%%%%%

%FOR BASEBAND MODE USE the AWGN Channel
EbNoLin=10^(EbNodB/10);
vn=(1/2)/EbNoLin;

%Calculate nvar for interference
powSig=1;
nvar=powSig/10^(loc4g.SNR/10);
nvar_int=nvar*int_ratio;

pow_rcosine=sum(a.txpulse(1:end).*conj(a.txpulse(1:end)))/length(a.txpulse(1:end));

pow_adj_nvar=sqrt(nvar_int)/sqrt(pow_rcosine);

pow_adj_nvar2=sqrt(nvar)/sqrt(pow_rcosine); %Added after new WFG

adj_txpulse=a.txpulse*pow_adj_nvar;

adj_txpulse2=a.txpulse*pow_adj_nvar2; %Added after new WFG

pow_adj_txpulse=sum(adj_txpulse(1:end).*conj(adj_txpulse(1:end)))/length(adj_txpulse(1:end));

pow_adj_txpulse2=sum(adj_txpulse2(1:end).*conj(adj_txpulse2(1:end)))/length(adj_txpulse2(1:end)); %Added after new WFG

int_Tsymbol=zeros(1,loc4g.Nfft+loc4g.CP);

int_Tsymbol2=zeros(1,loc4g.Nfft+loc4g.CP); %Added after new WFG

int_Tsymbol=[adj_txpulse(end-loc4g.CP+1:end),adj_txpulse];

int_Tsymbol2=[adj_txpulse2(end-loc4g.CP+1:end),adj_txpulse2]; %Added after new WFG

pow_int_Tsymbol=sum(int_Tsymbol(1:end).*conj(int_Tsymbol(1:end)))/length(int_Tsymbol(1:end));

pow_int_Tsymbol2=sum(int_Tsymbol2(1:end).*conj(int_Tsymbol2(1:end)))/length(int_Tsymbol2(1:end)); %Added after new WFG

yawgnbase_pre=conv(int_Tsymbol,radiochan(:,1));
%yawgnbase=conv(a.txpulse,radiochan(:,1));

yawgnbase_pre2=conv(int_Tsymbol2,radiochan(:,1)); %Added after new WFG

pow_yawgnbase_pre_a=sum(yawgnbase_pre(1:end).*conj(yawgnbase_pre(1:end)))/length(yawgnbase_pre(1:end));

pow_yawgnbase_pre_a2=sum(yawgnbase_pre2(1:end).*conj(yawgnbase_pre2(1:end)))/length(yawgnbase_pre2(1:end)); %Added after new WFG

pow_yawgnbase_pre=sum(yawgnbase_pre(loc4g.CP+1:loc4g.Nfft+loc4g.CP).*conj(yawgnbase_pre(loc4g.CP+1:loc4g.Nfft+loc4g.CP)))/length(yawgnbase_pre(loc4g.CP+1:loc4g.Nfft+loc4g.CP));

pow_yawgnbase_pre2=sum(yawgnbase_pre2(loc4g.CP+1:loc4g.Nfft+loc4g.CP).*conj(yawgnbase_pre2(loc4g.CP+1:loc4g.Nfft+loc4g.CP)))/length(yawgnbase_pre2(loc4g.CP+1:loc4g.Nfft+loc4g.CP)); %Added after new WFG

%pow_adj_chan=sqrt(pow_int_Tsymbol)/sqrt(pow_yawgnbase_pre);
pow_adj_chan=sqrt(pow_adj_txpulse)/sqrt(pow_yawgnbase_pre);

pow_adj_chan2=sqrt(pow_adj_txpulse2)/sqrt(pow_yawgnbase_pre2); %Added after new WFG

yawgnbase=yawgnbase_pre*pow_adj_chan;

yawgnbase2=yawgnbase_pre2*pow_adj_chan2; %Added after new WFG

yawgnbase_freq=fft(yawgnbase(loc4g.CP+1:loc4g.Nfft+loc4g.CP),loc4g.Nfft);

yawgnbase_freq2=fft(yawgnbase2(loc4g.CP+1:loc4g.Nfft+loc4g.CP),loc4g.Nfft); %Added after new WFG

yawgnbase_freq_v=yawgnbase_freq.';

yawgnbase_freq_v2=yawgnbase_freq2.'; %Added after new WFG

pow_yawgnbase_a=sum(yawgnbase(1:end).*conj(yawgnbase(1:end)))/length(yawgnbase(1:end));

pow_yawgnbase_a2=sum(yawgnbase2(1:end).*conj(yawgnbase2(1:end)))/length(yawgnbase2(1:end)); %Added after new WFG

pow_yawgnbase=sum(yawgnbase(loc4g.CP+1:loc4g.Nfft+loc4g.CP).*conj(yawgnbase(loc4g.CP+1:loc4g.Nfft+loc4g.CP)))/length(yawgnbase(loc4g.CP+1:loc4g.Nfft+loc4g.CP));

pow_yawgnbase2=sum(yawgnbase2(loc4g.CP+1:loc4g.Nfft+loc4g.CP).*conj(yawgnbase2(loc4g.CP+1:loc4g.Nfft+loc4g.CP)))/length(yawgnbase2(loc4g.CP+1:loc4g.Nfft+loc4g.CP));

% figure(6);
% xax1=a.Ts*(0:length(yawgnbase)-1);
% subplot(2,1,1), plot(xax1, real(yawgnbase));
% xlabel('sample time (sec)');
% ylabel('Real');
% titletext=strcat({'Baseband Signal in AWGN, EbNo='},{num2str(EbNodB)},...
%     {'dB, SamplingFrequency= '}, {num2str(1/a.Ts)});
% title(titletext);
% subplot(2,1,2), plot(xax1, imag(yawgnbase));
% xlabel('sample time (sec)');
% ylabel('Imag');
%%%%%%%%%%%%%%%%%


% figure(7);
% fs=1/a.Ts;
% NumPointsFFT=64;
% xax = ([0:NumPointsFFT-1]/NumPointsFFT)*fs;
% yaxdB=20*log10(abs(fft(yawgnbase,NumPointsFFT)));
% plot([xax-xax(end)/2], [yaxdB(NumPointsFFT/2+1:NumPointsFFT),yaxdB(1:NumPointsFFT/2)], 'b');
% titletext=strcat({'FFT of 8-PSK Signal, Rsymb='},{num2str(a.symbolrate)},{'Pulse Rolloff= '}, {num2str(a.rcrolloff)},...
%     {' SamplingRate= '}, {num2str(fs)});
% xtext='Frequency in Hz';
% ytext='Baseband Signal Power in dB';
% xlabel(xtext);
% ylabel(ytext);
% title(titletext);
varargout{1}=yawgnbase;
varargout{2}=radiochan;
varargout{3}=B;
varargout{4}=yawgnbase_freq_v;
varargout{5}=yawgnbase2;
varargout{6}=yawgnbase_freq_v2;


end




