function [varargout]=wfg(adhocobj,loc4g,bins_av,Intbaseband,Interchantime,Intbaseband_freq, Intbaseband2,Intbaseband_freq2,nvar_ratio,TowInd,mm,simode,smode,fldhomecombfig2)
Nfft=loc4g.Nfft;%64;
h=loc4g.h00N;%[0.8423i,  0.5391];  %Channel Time
Hk=fft(h,Nfft);        %Channel Freq
binsinmat=bins_av;%[2:32,34:64];%Non-zero bins contemplated
zerobin=[1,(loc4g.Nfft/2)+1];%[1,33];

%Calculate nvar for interference
powSig=1;
Nvar=powSig/10^(loc4g.SNR/10);%0.001;            %variance of the noise

interference_power=Intbaseband_freq.*conj(Intbaseband_freq);
if strcmp(smode,'AWGN+F+I')==true
    Ipow=interference_power.';%3*interference_power.';%3*Nvar*ones(size(Hk)); %Interference power
else
    Ipow(1,1:loc4g.Nfft)=0;
end
Txp=1;                 %Tx Power
Srx=0.1*ones(size(Hk));               %Rx power (typically coupled via path loss)
Srx(zerobin)=zeros(size(zerobin));
Nvar=((powSig/10^(loc4g.SNR/10)))*ones(size(Hk))*nvar_ratio;%0.001*ones(size(Hk));
ifft_Nvar=(ifft(Nvar,Nfft));
ifft_power_Nvar=sum(ifft_Nvar.*conj(ifft_Nvar))/length(ifft_Nvar);

gval=(1/length(binsinmat))*sum(Srx(binsinmat)./...
    (Nvar(binsinmat)+Ipow(binsinmat))); %Average SNIR at RX
gamk=zeros(size(Hk));
gamk(binsinmat)=(Hk(binsinmat).*conj(Hk(binsinmat)))*...
    sum(Srx(binsinmat)./(Nvar(binsinmat)+Ipow(binsinmat)));  %Variation in signal power gain at Rx
oldlen=length(binsinmat);
gam0=oldlen/(1+sum((gamk(binsinmat)).^-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamk(1)=0;
gamk((loc4g.Nfft2/2)+1)=0;

Freq_IndexN=[((-loc4g.fs/2)/loc4g.DF):((loc4g.fs/2)/loc4g.DF)-1];
Freq_IndexHz=(Freq_IndexN*loc4g.DF)/1e6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%New WFG 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [atemp, bind0a]=find(1./gamk(binsinmat) <= 1/gam0);   %waterfilling bins
% bind=binsinmat(bind0a);
% 
% sort_bins=sort(gamk(bind),'descend');
% gamma_0p=sort_bins(loc4g.SegSize+1)+1e-10;
% bind1=bind;
% 
% [atemp,bind0a]=find(1./gamk(binsinmat) <= 1/gamma_0p);
% bind1=binsinmat(bind0a);
% gamma_1p=gamma_0p; %added
% Pt=length(bind1)/gamma_1p-sum(1./gamk(bind1));
% ind2 = length(find(((1/gamma_1p - 1./gamk(bind1))*Pt) > 0));
% 
% bind=bind1;
% gam0=gamma_0p;

%%%%%%%%%%%
%Sigma Improvement %Mod W45 V3
%%%%%%%%%%%
gamklog=10*log10(gamk);

Std_1=std(gamklog(binsinmat));
Mean_gamklog=mean(gamklog(binsinmat));

errvec2=0;
errvec3=0;
kit=0;
kit2=0;
for k=1:length(binsinmat)
    if abs(gamklog(binsinmat(k))-Mean_gamklog) < Std_1;
        kit=kit+1;
        errvec2(kit)=binsinmat(k);
    end
    
    if gamklog(binsinmat(k))> Std_1+Mean_gamklog;
        kit2=kit2+1;
        errvec3(kit2)=binsinmat(k);
    end
end
%Revise minimum value
segment_values=[4,8,16,32,48,64,80,96,112,128]; %Mod W45 V3
[tmpval,mval]=min(abs(segment_values-length(errvec3))); %Mod W45 V3

adhocobj.SegSize=segment_values(mval); %Mod W45 V3
loc4g.SegSize=segment_values(mval); %Mod W45 V3
loc4g.cazacobj.SegSize=segment_values(mval); %Mod W45 V3

sort_bins=sort(gamklog(binsinmat),'descend'); %Mod W45 V3
gamklog_value=sort_bins(loc4g.SegSize+1); %Mod W45 V3
[tmpval2,mval2]=min(abs(gamklog-gamklog_value)); %Mod W45 V3
gamma_0p=gamk(mval2)+1e-10; %Mod W45 V3

%gamma_0p=length(errvec3)/(1+sum((gamk(errvec3)).^-1));
%[atemp,bind0a]=find(1./gamk(errvec3) <= 1/gamma_0p);
[atemp,bind0a]=find(1./gamk(binsinmat) <= 1/gamma_0p); %Mod W45 V3
bind1=binsinmat(bind0a);
%gamma_1p=length(bind1)/(1+sum((gamk(bind1)).^-1));%gamma_0p; %added
gamma_1p=gamma_0p; %Mod W45 V3
Pt=length(bind1)/gamma_1p-sum(1./gamk(bind1));

bind=bind1;
gam0=gamma_0p;

%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End - New WFG 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gainv=gamk(bind);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Proportion Keepers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gam0a=length(bind)/(1+sum((gamk(bind)).^-1));
% xgam0=gam0/gam0a;
% xgain1=((length(bind)/gam0)-1)/sum(1./(gamk(bind)));
% xgain2=1/xgain1;
% gainv2=gamk(bind)*xgain2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End - Proportion Keepers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Pt=length(bind)/gam0-sum(1./gamk(bind));
%Pt=length(bind)/gam0-sum(1./(gamk(bind)*xgain2)); %xgain2 adeed to keep Pt=1
PowerGain=zeros(1,Nfft);
PowerGain_norm=zeros(1,Nfft);
PowerGain(bind)=(1/gamma_1p - 1./gainv)*Pt;        %Linear power, not in dB!
%PowerGain(bind)=(1/gam0 - 1./gainv)*Pt; %Previous gainv2 results in neg WFG
Pt_one=1/sum(PowerGain); %Mod W45 V3
PowerGain(bind)=PowerGain(bind)*Pt_one; %Mod W45 V3


if strcmp(simode,'AWGN')==true;
    adhocobj.SegSize=128; %Mod W45 V4
    loc4g.SegSize=128; %Mod W45 V4
    loc4g.cazacobj.SegSize=128; %Mod W45 V4
    bind=bins_av(1:loc4g.SegSize);
    PowerGain_norm=zeros(1,Nfft);
    PowerGain(bind)=1;
    gam0=1;
    gainv=gamk(bind);
else
% This Test Eliminates previous WFG
%    PowerGain(bind)=1;
end
% adhocobj.SegSize=length(bind); %Mod W45 V3
% loc4g.SegSize=length(bind); %Mod W45 V3
% loc4g.cazacobj.SegSize=length(bind); %Mod W45 V3
adhocobj.NumbSeg=ceil(adhocobj.NumberSubCarr/adhocobj.SegSize); %Mod W45 V3
loc4g.NumbSeg=ceil(loc4g.NREf/loc4g.SegSize); %Mod W45 V3
loc4g.cazacobj.NumbSeg=ceil(loc4g.cazacobj.NREf/loc4g.cazacobj.SegSize); %Mod W45 V3

test_FSymbolO=zeros(1,Nfft);
if loc4g.SegSize>loc4g.NREf
    zadoffchu=loc4g.cazacobj.REarr(1:loc4g.NREf).';
    test_FSymbolO(1,bind(1:loc4g.NREf))=zadoffchu(1:loc4g.NREf).*PowerGain(bind(1:loc4g.NREf));
else
    zadoffchu=loc4g.cazacobj.REarr(1:loc4g.SegSize).';
    test_FSymbolO(1,bind(1:loc4g.SegSize))=zadoffchu(1:length(bind)).*PowerGain(bind); %Mod W45 V3
end
Frequency_Power_G=sum(test_FSymbolO(1,:).*conj(test_FSymbolO(1,:)));
%PowerGain2=sqrt((Nfft^2)/Frequency_Power_G);
PowerGain2=sqrt((1)/Frequency_Power_G);

PowerGain_norm(bind)=PowerGain(bind)*PowerGain2;
PowerGain_sum=sum(PowerGain_norm(bind));

comp_matrix(1:length(bind),1)=1/gam0;
comp_matrix(:,2)=(1./gamk(bind));
comp_matrix(:,3)=(1/gam0 - 1./gainv);
ind = find(comp_matrix(:,3) > 0);

gamrecip=zeros(1,Nfft);
gamrecip(binsinmat)=1./gamk(binsinmat);
Hk_binsinmat=zeros(1,Nfft);
Hk_binsinmat(binsinmat)=real(Hk(binsinmat));

% figure(30)
% plot(0:Nfft-1, PowerGain_norm, 'b', 0:Nfft-1, gamrecip, 'r')
% xlabel('Frequency Bin Index')
% ylabel('Linear Signal Power')
% title('Blue: WaterFill Gain, Red: (N+I)/S')

%figure(3)
figure('visible','off');
switch simode
    
    case 'AWGN_Fadding_Interference'
        
        %plot(Freq_IndexHz, PowerGain_norm, 'b', Freq_IndexHz, gamrecip, 'r')
        stem(Freq_IndexHz, PowerGain_norm, 'b','Marker','none')
        hold on
        plot(Freq_IndexHz, gamrecip, 'r')
        hold off;        
        titlefig=strcat({['Frequency Response of Blue: Water Filling Gain and Red: (N+I)/S'];['Sim=',num2str(mm),'   Tower=',num2str(TowInd)]});
        title(titlefig)
        xlabel('Frequency Index MHz')
        ylabel('Linear Signal Power')
        
    case'AWGN_Fadding'
        
        stem(Freq_IndexHz, PowerGain_norm, 'b','Marker','none')
        hold on
        plot(Freq_IndexHz, gamrecip, 'r')
        hold off;
        titlefig=strcat({['Frequency Response of Blue: Water Filling Gain and Red: N/S'];['Sim=',num2str(mm),'   Tower=',num2str(TowInd)]});
        title(titlefig)
        xlabel('Frequency Index MHz')
        ylabel('Linear Signal Power')
        
    case 'AWGN'
        
        stem(Freq_IndexHz, PowerGain_norm, 'b','Marker','none')
        hold on
        plot(Freq_IndexHz, gamrecip, 'r')
        hold off;
        xlabel('Frequency Index MHz')
        ylabel('Linear Signal Power')
        titlefig=strcat({['Frequency Response of Blue: AWGN Gain and Red: N/S'];['Sim=',num2str(mm),'   Tower=',num2str(TowInd)]});
        title(titlefig)
        %axis([-inf,inf,0,1.2])
        
end

saveas(gcf, fldhomecombfig2);

varargout{1}=PowerGain_norm;
varargout{2}=bind;
end



