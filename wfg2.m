function [varargout]=wfg2(adhocobj,loc4g,bins_av,Intbaseband,Interchantime,Intbaseband_freq, Intbaseband2,Intbaseband_freq2,nvar_ratio)
Nfft=loc4g.Nfft;%64;
h=loc4g.h00N;%[0.8423i,  0.5391];  %Channel Time
Hk=fft(h,Nfft);        %Channel Freq
binsinmat=bins_av;%[2:32,34:64];%Non-zero bins contemplated
zerobin=[1,(loc4g.Nfft/2)+1];%[1,33];

%Calculate nvar for interference
powSig=1;
Nvar=powSig/10^(loc4g.SNR/10);%0.001;            %variance of the noise

interference_power=Intbaseband_freq.*conj(Intbaseband_freq);
Ipow=interference_power.';%3*interference_power.';%3*Nvar*ones(size(Hk)); %Interference power
Txp=1;                 %Tx Power
Srx=0.1*ones(size(Hk));               %Rx power (typically coupled via path loss)
Srx(zerobin)=zeros(size(zerobin));
Nvar=(powSig/10^(loc4g.SNR/10))*ones(size(Hk))*nvar_ratio;%0.001*ones(size(Hk));
gval=(1/length(binsinmat))*sum(Srx(binsinmat)./...
    (Nvar(binsinmat)+Ipow(binsinmat))); %Average SNIR at RX
gamk=zeros(size(Hk));
gamk(binsinmat)=(Hk(binsinmat).*conj(Hk(binsinmat)))*...
    sum(Srx(binsinmat)./(Nvar(binsinmat)+Ipow(binsinmat)));  %Variation in signal power gain at Rx
oldlen=length(binsinmat);
gam0=oldlen/(1+sum((gamk(binsinmat)).^-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Old WFG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indv=0;
SL=1;
bindold=length(binsinmat);
while SL==1
    indv=indv+1;
    gam(indv)=gam0;
    [atemp, bind0a]=find(1./gamk(binsinmat) <= 1/gam(indv));   %waterfilling bins
    bind=binsinmat(bind0a);
    gam0=length(bind)/(1+sum((gamk(bind).^-1)));  %update threshold
    Pt=length(bind)/gam0-sum(1./gamk(bind))
    if indv >= 2
        if length(bind)==bindold                      %stop if the iterations converge
            SL=0;
        end
    end
    bindold=length(bind);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End - Old WFG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%New WFG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [atemp, bind0a]=find(1./gamk(binsinmat) <= 1/gam0);   %waterfilling bins
% bind=binsinmat(bind0a);
% 
% bind1=bind;
% inc_factor=min(gamk(bind1))/1000;
% gamma_0p=gam0;
% while length(bind1)>(loc4g.SegSize)
%     gamma_0p_pre=gamma_0p;
%     gamma_0p=gamma_0p+inc_factor;
%     [atemp,bind0a]=find(1./gamk(binsinmat) <= 1/gamma_0p);
%     %xbind1=find(bind1>(min(bins_av)-1)&bind1<(max(bins_av)+1));%Added
%     %bind1=bind1(xbind1);%Added
%     bind1=binsinmat(bind0a);
%     gamma_1p=length(bind1)/(1+sum((gamk(bind1).^-1)));
%     Pt=length(bind1)/gamma_1p-sum(1./gamk(bind1));
% end
% 
% if length(bind1)<(loc4g.SegSize)
%     gamma_0p=gamma_0p_pre;
%     [atemp,bind0a]=find(1./gamk(binsinmat) <= gamma_0p);
%     %xbind1=find(bind1>(min(bins_av)-1)&bind1<(max(bins_av)+1));%Added
%     %bind1=bind1(xbind1);%Added
%     bind1=binsinmat(bind0a);
%     gamma_1p=length(bind1)/(1+sum((gamk(bind1).^-1)));
%     Pt=length(bind1)/gamma_1p-sum(1./gamk(bind1));
% end
% bind=bind1;
% gam0=length(bind)/(1+sum((gamk(bind).^-1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End - New WFG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gainv=gamk(bind);
Pt=length(bind)/gam0-sum(1./gamk(bind));
PowerGain=zeros(1,Nfft);%zeros(1,length(Nfft));
PowerGain(bind)=(1/gam0 - 1./gainv)*Pt;        %Linear power, not in dB!

comp_matrix(1:length(bind),1)=1/gam0;
comp_matrix(:,2)=(1./gamk(bind));
comp_matrix(:,3)=(1/gam0 - 1./gainv);
ind = find(comp_matrix(:,3) > 0);

gamrecip=zeros(1,Nfft);
gamrecip(binsinmat)=1./gamk(binsinmat);
Hk_binsinmat=zeros(1,Nfft);
Hk_binsinmat(binsinmat)=real(Hk(binsinmat));
figure(31)
plot(0:Nfft-1, PowerGain, 'k', 0:Nfft-1, gamrecip, 'r', 0:Nfft-1, Hk_binsinmat, 'b')
xlabel('Frequency Bin Index')
ylabel('Linear Signal Power')
title('Blue: Ricean Channel Freq. Response, Red: (N+I)/S,  Black: WaterFill Gain')

varargout{1}=PowerGain
end



