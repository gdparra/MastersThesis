%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main Simulation Script - Simulation W5 V1A Complete & Simple
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varargout] = delayest(timedataout,loc4g,cazacobj,zcobj,vn,SignalPower,PowerGain,defdel,nvar_ratio,Total_pow_vn,fldhomecombfig5,fldhomecombfig8,mm,TowInd,simode,comp_table_h,RiceanEst,fldhomecombfig9)
cp_inc_samples=(loc4g.CP/0.1)+1;
max_samples=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Water Filling Equalizer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for tkt=1:loc4g.Nfft2
    W_gain(tkt,1)=(conj(PowerGain(tkt))/((conj(PowerGain(tkt))*(PowerGain(tkt)))));
end
W_gain(isnan(W_gain))=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference Symbols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ocazac(1:loc4g.NREf,1)=cazacobj.REarr(1:loc4g.NREf,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OFDM Symbol Rx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for CntNumbSeg=1:zcobj.NumbSeg
    %Received Signal Power Estimation
    rx_tsignal=timedataout(1:zcobj.NumbSeg,loc4g.CP2+1:loc4g.CP2+loc4g.Nfft2);
    pow_rx_tsignal(CntNumbSeg,1)=sum(rx_tsignal(CntNumbSeg,1:end).*conj(rx_tsignal(CntNumbSeg,1:end)))/length(rx_tsignal(CntNumbSeg,1:end));
    %Convert received segments to frequency
    Rxp(:,CntNumbSeg)=fft(timedataout(CntNumbSeg,loc4g.CP2+1:loc4g.CP2+loc4g.Nfft2).',loc4g.Nfft2);
    Frequency_Power_AVG_Rxp(CntNumbSeg,:)=sum(Rxp(:,CntNumbSeg).*conj(Rxp(:,CntNumbSeg)));
    %Computed normalized noise according to %25 ratio
    vn_r(CntNumbSeg,1)=vn(CntNumbSeg,1)*nvar_ratio;
    %Computed noise with no gain according to %25 ratio
    nvar_ng(CntNumbSeg,1)=cazacobj.nvar_ng(1,CntNumbSeg)*nvar_ratio;
    %Equlizer for normalized signal
    PowerGain2_eq(CntNumbSeg,1)=(conj(cazacobj.PowerGain2(CntNumbSeg,1))/((conj(cazacobj.PowerGain2(CntNumbSeg,1))*cazacobj.PowerGain2(CntNumbSeg,1))));
end
saved_Rxp(1:loc4g.Nfft2,1:zcobj.NumbSeg)=Rxp(1:loc4g.Nfft2,1:zcobj.NumbSeg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimated Channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h00NDel=RiceanEst.h00NDel(TowInd+1,1:end);
HDel=fft(h00NDel,length(h00NDel));
Hlosorg=fft(RiceanEst.h00Nlos,loc4g.Nfft2).';
Hscatorg=fft(RiceanEst.h00Nscat,loc4g.Nfft2).';
Horg=fft(RiceanEst.h00N,loc4g.Nfft2).';
HDel2=fft(h00NDel,loc4g.Nfft2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars cnt0 current timedelay delayf1 delayi1 Hest HGain h00Nest hGain 
clearvars h00Nscatest Hscatest HscatGain HscatPG Hscat HscatGain2
clearvars h00Nlosest Hlosest HlosGain HlosPG Hlos HlosGain2
clearvars Rxp_temp Rxp_temp2 adj_Rxp adj_Rx Fg0p Fg0 tempc
clearvars adj_RxGain adj_RxPG adj_Rxpn adj_RxGain2 adj_RxPG2 saved_adj_Rxp saved_Fg0
clearvars zc2 tau xax cr sumvec saved_cr cnta AxValue BxIndex

cnt0=0;
for timedelay=[0:floor(defdel),defdel,ceil(defdel):length(h00NDel)-1+20]%W6 V2A
    cnt0=cnt0+1;
    for loopN=0:loc4g.Nfft2-1
        delayf1(loopN+1,cnt0)=exp(-1i*2*pi*loopN*(timedelay/loc4g.Nfft2));
        delayi1(loopN+1,cnt0)=exp(1i*2*pi*loopN*(timedelay/loc4g.Nfft2));
    end
%%%%%W6 V2A     
%     for loopN=0:length(h00NDel)-1
%         Hest(loopN+1,cnt0)=HDel(loopN+1)*exp(1i*2*pi*loopN*(timedelay/length(h00NDel)));
%     end
    Hest(:,cnt0)=HDel(1,:);%W6 V2A
%%%%%W6 V2A     
    HGain(cnt0,1)=sum(Hest(:,cnt0).*conj(Hest(:,cnt0)));
    
    h00Nest=ifft(Hest(:,cnt0),length(h00NDel));
%%%%%W6 V2A      
    if strcmp(simode,'AWGN')
        h00Nest=h00NDel.';
    end
%%%%%W6 V2A      
    hGain(:,cnt0)=h00Nest.*conj(h00Nest);
    h00Nscatest=[0;h00Nest(2:end,1)];
    
    Hscatest=fft(h00Nscatest,loc4g.Nfft2);
    HscatGain(cnt0,1)=sum(Hscatest.*conj(Hscatest));
    HscatPG=sqrt((70.045)/HscatGain(cnt0,1));
%%%%%W6 V2A      
    if strcmp(simode,'AWGN')
        HscatPG=1;
    end
%%%%%W6 V2A 
    Hscat(:,cnt0)=Hscatest*HscatPG;
    HscatGain2(cnt0,1)=sum(Hscat(:,cnt0).*conj(Hscat(:,cnt0)));

    h00Nlosest=[h00Nest(1,1),zeros(1,length(h00Nest)-1)];
    
    Hlosest=fft(h00Nlosest,loc4g.Nfft2);
    HlosGain(cnt0,1)=sum(Hlosest.*conj(Hlosest));
    HlosPG=sqrt((441.95)/HlosGain(cnt0,1));
%%%%%W6 V2A      
    if strcmp(simode,'AWGN')
        HlosPG=1;
    end
%%%%%W6 V2A 
    Hlos(:,cnt0)=(Hlosest.')*HlosPG;
    HlosGain2(cnt0,1)=sum(Hlos(:,cnt0).*conj(Hlos(:,cnt0)));
end

% Graphs
% figure(100)
% plot(0:cnt0-1,HGain)
% xlabel('H ^T_o');
% ylabel('Gain Amplitude');
% title('Ricean H Gain per ^T_o');
% axis([0,160,0,200])
% 
% figure(101)
% plot(0:cnt0-1,HlosGain)
% xlabel('LOS ^T_o');
% ylabel('Gain Amplitude');
% title('LOS Gain per ^T_o');
% axis([0,160,0,550])
% 
% figure(102)
% plot(0:cnt0-1,HlosGain2)
% xlabel('LOS ^T_o');
% ylabel('Norm. Gain Amplitude');
% title('LOS Normalized Gain per T_o');
% axis([0,160,0,550])
% 
% figure(103)
% plot(0:cnt0-1,HscatGain)
% xlabel('ICI ^T_o');
% ylabel('Gain Amplitude');
% title('ICI Gain per ^T_o');
% axis([0,160,0,550])
% 
% figure(104)
% plot(0:cnt0-1,HscatGain2)
% xlabel('ICI ^T_o');
% ylabel('Norm. Gain Amplitude');
% title('ICI Normalized Gain per ^T_o');
% axis([0,160,0,100])

%%%%%W6 V2A 
for CntNumbSeg=1:zcobj.NumbSeg
    Rxp(:,CntNumbSeg)=Rxp(:,CntNumbSeg).*W_gain(:,1)*PowerGain2_eq(CntNumbSeg,1);
    cazacobj.FSymbol_g0(CntNumbSeg,:)=(cazacobj.FSymbol_g0(CntNumbSeg,:).').*W_gain(:,1)*PowerGain2_eq(CntNumbSeg,1);
end
%%%%%W6 V2A

for cnta=1:cnt0
    for CntNumbSeg=1:zcobj.NumbSeg
        Rxp_temp(:,CntNumbSeg)=(Rxp(:,CntNumbSeg)-((cazacobj.FSymbol_g0(CntNumbSeg,:).').*Hscat(:,cnta).*delayf1(:,cnta)));      
        %Rxp_temp2(:,CntNumbSeg)=Rxp_temp(:,CntNumbSeg).*delayi1(:,cnta).*W_gain(:,1)*(1/Hlos(1,cnta))*PowerGain2_eq(CntNumbSeg,1);%W6 V2A
        Rxp_temp2(:,CntNumbSeg)=Rxp_temp(:,CntNumbSeg).*delayi1(:,cnta)*(1/Hlos(1,cnta));%W6 V2A
        adj_Rxp(1:zcobj.SegSize,CntNumbSeg)=Rxp_temp2(cazacobj.bins_gain(1:zcobj.SegSize),CntNumbSeg);
        adj_Rx(((CntNumbSeg-1)*zcobj.SegSize)+1:CntNumbSeg*zcobj.SegSize,1)=adj_Rxp(1:zcobj.SegSize,CntNumbSeg);
        %Fg0p(1:zcobj.SegSize,CntNumbSeg)=zcobj.FSymbol_g0(CntNumbSeg,cazacobj.bins_gain(1:zcobj.SegSize))*Hlosorg(1,1).*(W_gain(cazacobj.bins_gain(1:zcobj.SegSize),1).')*(1/Hlosorg(1,1))*PowerGain2_eq(CntNumbSeg,1);%W6 V2A
        Fg0p(1:zcobj.SegSize,CntNumbSeg)=zcobj.FSymbol_g0(CntNumbSeg,cazacobj.bins_gain(1:zcobj.SegSize))*Hlosorg(1,1)*(1/Hlosorg(1,1));%W6 V2A
        Fg0(((CntNumbSeg-1)*zcobj.SegSize)+1:CntNumbSeg*zcobj.SegSize,1)=Fg0p(1:zcobj.SegSize,CntNumbSeg);
    end
    tempc=adj_Rx(1:loc4g.NREf,1);
    adj_RxGain(1,cnta)=sum(tempc.*conj(tempc));
    adj_RxPG=sqrt((loc4g.NREf)/adj_RxGain(1,cnta));
    
    if adj_RxGain(1,cnta)>loc4g.NREf
        adj_Rxpn(1:loc4g.NREf,1)=adj_Rx(1:loc4g.NREf,1)*adj_RxPG;
    else
        adj_Rxpn(1:loc4g.NREf,1)=adj_Rx(1:loc4g.NREf,1);
    end
    adj_RxGain2(:,cnta)=adj_Rxpn(1:loc4g.NREf,1).*conj(adj_Rxpn(1:loc4g.NREf,1));
    adj_RxPG2=sqrt((1)./adj_RxGain2(:,cnta));
    
    if adj_RxGain(1,cnta)>loc4g.NREf
        adj_Rxpn(1:loc4g.NREf,1)=adj_Rx(1:loc4g.NREf,1)*adj_RxPG.*adj_RxPG2;
    else
        adj_Rxpn(1:loc4g.NREf,1)=adj_Rx(1:loc4g.NREf,1);
    end
    
    Temp_mat=[Fg0(1:loc4g.NREf),adj_Rxpn(1:loc4g.NREf)];
    saved_adj_Rxp(:,cnta)=adj_Rxpn;
    saved_Fg0(:,cnta)=Fg0(1:loc4g.NREf,1);
    zc2=adj_Rxpn;
    tau=0;
    xax(cnta)=tau;
    cr=sum(conj(Fg0(1:loc4g.NREf)).*zc2([rem(tau+loc4g.NREf,loc4g.NREf)+1:loc4g.NREf,1:rem(tau+loc4g.NREf,loc4g.NREf)]));
    sumvec(cnta)=real(cr);
    saved_cr(cnta,1:length(cr))=cr(1:length(cr));
    cnta=cnta+1;
end

[AxValue BxIndex]=sort(sumvec,'descend');
index_max_phase=find(sumvec(1,:)==max(sumvec(1,:)),1,'last');
comp_mat=[saved_Fg0(:,round(defdel)+1),Ocazac,saved_adj_Rxp(:,round(defdel)+1),saved_adj_Rxp(:,index_max_phase)];
firstphase=index_max_phase-2;
% figure(105)
% stem(0,saved_cr(round(defdel)+1,:))
% xlabel('-N/2:N/2+1 Correlation Resource Element');
% ylabel('Real Component Amplitude');
% title('Actual Delay Correlation');
% axis([-loc4g.NREf/2 loc4g.NREf/2 0 max(sumvec)])
% 
% figure(106)
% stem(0,saved_cr(index_max_phase,:))
% xlabel('-N/2:N/2+1 Correlation Resource Element');
% ylabel('Gain Amplitude');
% title('Estimated Delay Correlation');
% axis([-loc4g.NREf/2 loc4g.NREf/2 0 max(sumvec)])
% 
% figure(107)
% plot(0:cnta-2,sumvec)
% xlabel('Correlation Zero Flag for every ^To');
% ylabel('Zero Flag Amplitude');
% title('Zero Flag Real Value Amplitude');

if strcmp(simode,'AWGN')
    %figure(9)
    figure('visible','off');
    plot(1:(cnta-1),sumvec(1:(cnta-1)));
    xlabel('Correlation Zero Flag Value per Phase Change');
    ylabel('Amplitude - Zero Flag');
    title({['1st Est. Phase - Phase=',num2str(firstphase),',   Real CR Samp=',num2str(floor(defdel))];['Def. Phase=',num2str(defdel),'   Sim=',num2str(mm),'   Tower=',num2str(TowInd)]})
    saveas(gcf, fldhomecombfig5);
    %print(gcf, '-djpeg', fldhomecombfig11,'-r500');
else
    %figure(9)
    figure('visible','off');
    plot(1:(cnta-1),sumvec(1:(cnta-1)));
    xlabel('Correlation Zero Flag Value per Phase Change');
    ylabel('Amplitude - Zero Flag');
    title({['1st Est. Phase - Phase=',num2str(firstphase),',   Real CR Samp=',num2str(floor(defdel))];['Def. Phase=',num2str(defdel),'   Sim=',num2str(mm),'   Tower=',num2str(TowInd)]})
    saveas(gcf, fldhomecombfig5);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars cnt0 current timedelay delayf1 delayi1 Hest HGain h00Nest hGain 
clearvars h00Nscatest Hscatest HscatGain HscatPG Hscat HscatGain2
clearvars h00Nlosest Hlosest HlosGain HlosPG Hlos HlosGain2
clearvars Rxp_temp Rxp_temp2 adj_Rxp adj_Rx Fg0p Fg0 tempc
clearvars adj_RxGain adj_RxPG adj_Rxpn adj_RxGain2 adj_RxPG2 saved_adj_Rxp saved_Fg0
clearvars zc2 tau xax cr sumvec saved_cr cnta AxValue BxIndex

cnt0=0;
current=(index_max_phase-1);
for timedelay=current-1:0.1:current+1
    cnt0=cnt0+1;
    for loopN=0:loc4g.Nfft2-1
        delayf1(loopN+1,cnt0)=exp(-1i*2*pi*loopN*(timedelay/loc4g.Nfft2));
        delayi1(loopN+1,cnt0)=exp(1i*2*pi*loopN*(timedelay/loc4g.Nfft2));
    end
%%%%%W6 V2A     
%     for loopN=0:length(h00NDel)-1
%         Hest(loopN+1,cnt0)=HDel(loopN+1)*exp(1i*2*pi*loopN*(timedelay/length(h00NDel)));
%     end
    Hest(:,cnt0)=HDel(1,:);%W6 V2A
%%%%%W6 V2A     
    HGain(cnt0,1)=sum(Hest(:,cnt0).*conj(Hest(:,cnt0)));
    
    h00Nest=ifft(Hest(:,cnt0),length(h00NDel));
%%%%%W6 V2A      
    if strcmp(simode,'AWGN')
        h00Nest=h00NDel.';
    end
%%%%%W6 V2A
    hGain(:,cnt0)=h00Nest.*conj(h00Nest);
    h00Nscatest=[0;h00Nest(2:end,1)];
    
    Hscatest=fft(h00Nscatest,loc4g.Nfft2);
    HscatGain(cnt0,1)=sum(Hscatest.*conj(Hscatest));
    HscatPG=sqrt((70.045)/HscatGain(cnt0,1));
%%%%%W6 V2A      
    if strcmp(simode,'AWGN')
        HscatPG=1;
    end
%%%%%W6 V2A 
    Hscat(:,cnt0)=Hscatest*HscatPG;
    HscatGain2(cnt0,1)=sum(Hscat(:,cnt0).*conj(Hscat(:,cnt0)));

    h00Nlosest=[h00Nest(1,1),zeros(1,length(h00Nest)-1)];
    
    Hlosest=fft(h00Nlosest,loc4g.Nfft2);
    HlosGain(cnt0,1)=sum(Hlosest.*conj(Hlosest));
    HlosPG=sqrt((441.95)/HlosGain(cnt0,1));
%%%%%W6 V2A      
    if strcmp(simode,'AWGN')
        HlosPG=1;
    end
%%%%%W6 V2A 
    Hlos(:,cnt0)=(Hlosest.')*HlosPG;
    HlosGain2(cnt0,1)=sum(Hlos(:,cnt0).*conj(Hlos(:,cnt0)));
end

% Graphs
% figure(100)
% plot(current-1:0.1:current+1,HGain)
% xlabel('H ^T_o');
% ylabel('Gain Amplitude');
% title('Ricean H Gain per ^T_o');
% axis([current-1,current+1,0,200])
% 
% figure(101)
% plot(current-1:0.1:current+1,HlosGain)
% xlabel('LOS ^T_o');
% ylabel('Gain Amplitude');
% title('LOS Gain per ^T_o');
% axis([current-1,current+1,0,550])
% 
% figure(102)
% plot(current-1:0.1:current+1,HlosGain2)
% xlabel('LOS ^T_o');
% ylabel('Norm. Gain Amplitude');
% title('LOS Normalized Gain per T_o');
% axis([current-1,current+1,0,550])
% 
% figure(103)
% plot(current-1:0.1:current+1,HscatGain)
% xlabel('ICI ^T_o');
% ylabel('Gain Amplitude');
% title('ICI Gain per ^T_o');
% axis([current-1,current+1,0,550])
% 
% figure(104)
% plot(current-1:0.1:current+1,HscatGain2)
% xlabel('ICI ^T_o');
% ylabel('Norm. Gain Amplitude');
% title('ICI Normalized Gain per ^T_o');
% axis([current-1,current+1,0,100])

for cnta=1:cnt0
    for CntNumbSeg=1:zcobj.NumbSeg
        Rxp_temp(:,CntNumbSeg)=(Rxp(:,CntNumbSeg)-((cazacobj.FSymbol_g0(CntNumbSeg,:).').*Hscat(:,cnta).*delayf1(:,cnta)));      
        %Rxp_temp2(:,CntNumbSeg)=Rxp_temp(:,CntNumbSeg).*delayi1(:,cnta).*W_gain(:,1)*(1/Hlos(1,cnta))*PowerGain2_eq(CntNumbSeg,1);%W6 V2A  
        Rxp_temp2(:,CntNumbSeg)=Rxp_temp(:,CntNumbSeg).*delayi1(:,cnta)*(1/Hlos(1,cnta));%W6 V2A
        adj_Rxp(1:zcobj.SegSize,CntNumbSeg)=Rxp_temp2(cazacobj.bins_gain(1:zcobj.SegSize),CntNumbSeg);
        adj_Rx(((CntNumbSeg-1)*zcobj.SegSize)+1:CntNumbSeg*zcobj.SegSize,1)=adj_Rxp(1:zcobj.SegSize,CntNumbSeg);
        %Fg0p(1:zcobj.SegSize,CntNumbSeg)=zcobj.FSymbol_g0(CntNumbSeg,cazacobj.bins_gain(1:zcobj.SegSize))*Hlosorg(1,1).*(W_gain(cazacobj.bins_gain(1:zcobj.SegSize),1).')*(1/Hlosorg(1,1))*PowerGain2_eq(CntNumbSeg,1);%W6 V2A
        Fg0p(1:zcobj.SegSize,CntNumbSeg)=zcobj.FSymbol_g0(CntNumbSeg,cazacobj.bins_gain(1:zcobj.SegSize))*Hlosorg(1,1)*(1/Hlosorg(1,1));%W6 V2A
        Fg0(((CntNumbSeg-1)*zcobj.SegSize)+1:CntNumbSeg*zcobj.SegSize,1)=Fg0p(1:zcobj.SegSize,CntNumbSeg);
    end
    tempc=adj_Rx(1:loc4g.NREf,1);
    adj_RxGain(1,cnta)=sum(tempc.*conj(tempc));
    adj_RxPG=sqrt((loc4g.NREf)/adj_RxGain(1,cnta));
    
    if adj_RxGain(1,cnta)>loc4g.NREf
        adj_Rxpn(1:loc4g.NREf,1)=adj_Rx(1:loc4g.NREf,1)*adj_RxPG;
    else
        adj_Rxpn(1:loc4g.NREf,1)=adj_Rx(1:loc4g.NREf,1);
    end
    adj_RxGain2(:,cnta)=adj_Rxpn(1:loc4g.NREf,1).*conj(adj_Rxpn(1:loc4g.NREf,1));
    adj_RxPG2=sqrt((1)./adj_RxGain2(:,cnta));
    
    if adj_RxGain(1,cnta)>loc4g.NREf
        adj_Rxpn(1:loc4g.NREf,1)=adj_Rx(1:loc4g.NREf,1)*adj_RxPG.*adj_RxPG2;
    else
        adj_Rxpn(1:loc4g.NREf,1)=adj_Rx(1:loc4g.NREf,1);
    end
    
    Temp_mat=[Fg0(1:loc4g.NREf),adj_Rxpn(1:loc4g.NREf)];
    saved_adj_Rxp(:,cnta)=adj_Rxpn;
    saved_Fg0(:,cnta)=Fg0(1:loc4g.NREf,1);
    zc2=adj_Rxpn;
    tau=0;
    xax(cnta)=tau;
    cr=sum(conj(Fg0(1:loc4g.NREf)).*zc2([rem(tau+loc4g.NREf,loc4g.NREf)+1:loc4g.NREf,1:rem(tau+loc4g.NREf,loc4g.NREf)]));
    sumvec(cnta)=real(cr);
    saved_cr(cnta,1:length(cr))=cr(1:length(cr));
    cnta=cnta+1;
end

[AxValue BxIndex]=sort(sumvec,'descend');
index_max_phase2=find(sumvec(1,:)==max(sumvec(1,:)),1,'last');
comp_mat=[saved_Fg0(:,10+1),Ocazac,saved_adj_Rxp(:,10+1),saved_adj_Rxp(:,index_max_phase2)];
firstphase=(index_max_phase-2)+(((index_max_phase2-1)*0.1)-0.1);
% figure(105)
% stem(0,saved_cr(10+1,:))
% xlabel('-N/2:N/2+1 Correlation Resource Element');
% ylabel('Real Component Amplitude');
% title('Actual Delay Correlation');
% axis([-loc4g.NREf/2 loc4g.NREf/2 0 max(sumvec)])
% 
% figure(106)
% stem(0,saved_cr(index_max_phase2,:))
% xlabel('-N/2:N/2+1 Correlation Resource Element');
% ylabel('Gain Amplitude');
% title('Estimated Delay Correlation');
% axis([-loc4g.NREf/2 loc4g.NREf/2 0 max(sumvec)])
% 
% figure(107)
% plot(current-1:0.1:current+1,sumvec)
% xlabel('Correlation Zero Flag for every ^To');
% ylabel('Zero Flag Amplitude');
% title('Zero Flag Real Value Amplitude');
% axis([current-1 current+1 0 max(sumvec)])

if strcmp(simode,'AWGN')
    %figure(10)
    figure('visible','off');
    plot(1:(cnta-1),sumvec(1:(cnta-1)));
    xlabel('Correlation Zero Flag Value per Phase Change');
    ylabel('Amplitude - Zero Flag');
    title({['2nd Est. Phase - Phase=',num2str(firstphase),',   Real CR Samp=',num2str(floor(defdel*10)/10)];['Def. Phase=',num2str(defdel),'   Sim=',num2str(mm),'   Tower=',num2str(TowInd)]})
    saveas(gcf, fldhomecombfig8);
    %print(gcf, '-djpeg', fldhomecombfig11,'-r500');
else
    %figure(10)
    figure('visible','off');
    plot(1:(cnta-1),sumvec(1:(cnta-1)));
    xlabel('Correlation Zero Flag Value per Phase Change');
    ylabel('Amplitude - Zero Flag');
    title({['2nd Est. Phase - Phase=',num2str(firstphase),',   Real CR Samp=',num2str(floor(defdel*10)/10)];['Def. Phase=',num2str(defdel),'   Sim=',num2str(mm),'   Tower=',num2str(TowInd)]})
    saveas(gcf, fldhomecombfig8);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars cnt0 current timedelay delayf1 delayi1 Hest HGain h00Nest hGain 
clearvars h00Nscatest Hscatest HscatGain HscatPG Hscat HscatGain2
clearvars h00Nlosest Hlosest HlosGain HlosPG Hlos HlosGain2
clearvars Rxp_temp Rxp_temp2 adj_Rxp adj_Rx Fg0p Fg0 tempc
clearvars adj_RxGain adj_RxPG adj_Rxpn adj_RxGain2 adj_RxPG2 saved_adj_Rxp saved_Fg0
clearvars zc2 tau xax cr sumvec saved_cr cnta AxValue BxIndex

cnt0=0;
current=(index_max_phase-2)+((index_max_phase2-1)*0.1);
for timedelay=current-0.1:0.001:current+0.1
    cnt0=cnt0+1;
    for loopN=0:loc4g.Nfft2-1
        delayf1(loopN+1,cnt0)=exp(-1i*2*pi*loopN*(timedelay/loc4g.Nfft2));
        delayi1(loopN+1,cnt0)=exp(1i*2*pi*loopN*(timedelay/loc4g.Nfft2));
    end
%%%%%W6 V2A     
%     for loopN=0:length(h00NDel)-1
%         Hest(loopN+1,cnt0)=HDel(loopN+1)*exp(1i*2*pi*loopN*(timedelay/length(h00NDel)));
%     end
    Hest(:,cnt0)=HDel(1,:);%W6 V2A
%%%%%W6 V2A     
    HGain(cnt0,1)=sum(Hest(:,cnt0).*conj(Hest(:,cnt0)));
    
    h00Nest=ifft(Hest(:,cnt0),length(h00NDel));
%%%%%W6 V2A      
    if strcmp(simode,'AWGN')
        h00Nest=h00NDel.';
    end
%%%%%W6 V2A
    hGain(:,cnt0)=h00Nest.*conj(h00Nest);
    h00Nscatest=[0;h00Nest(2:end,1)];
    
    Hscatest=fft(h00Nscatest,loc4g.Nfft2);
    HscatGain(cnt0,1)=sum(Hscatest.*conj(Hscatest));
    HscatPG=sqrt((70.045)/HscatGain(cnt0,1));
%%%%%W6 V2A      
    if strcmp(simode,'AWGN')
        HscatPG=1;
    end
%%%%%W6 V2A 
    Hscat(:,cnt0)=Hscatest*HscatPG;
    HscatGain2(cnt0,1)=sum(Hscat(:,cnt0).*conj(Hscat(:,cnt0)));

    h00Nlosest=[h00Nest(1,1),zeros(1,length(h00Nest)-1)];
    
    Hlosest=fft(h00Nlosest,loc4g.Nfft2);
    HlosGain(cnt0,1)=sum(Hlosest.*conj(Hlosest));
    HlosPG=sqrt((441.95)/HlosGain(cnt0,1));
%%%%%W6 V2A      
    if strcmp(simode,'AWGN')
        HlosPG=1;
    end
%%%%%W6 V2A 
    Hlos(:,cnt0)=(Hlosest.')*HlosPG;
    HlosGain2(cnt0,1)=sum(Hlos(:,cnt0).*conj(Hlos(:,cnt0)));
end

% Graphs
% figure(100)
% plot(current-0.1:0.001:current+0.1,HGain)
% xlabel('H ^T_o');
% ylabel('Gain Amplitude');
% title('Ricean H Gain per ^T_o');
% %axis([current-0.1,current+0.1,0,200])
% 
% figure(101)
% plot(current-0.1:0.001:current+0.1,HlosGain)
% xlabel('LOS ^T_o');
% ylabel('Gain Amplitude');
% title('LOS Gain per ^T_o');
% xlim([current-0.1,current+0.1])
% axis 'auto y'
% 
% figure(102)
% plot(current-0.1:0.001:current+0.1,HlosGain2)
% xlabel('LOS ^T_o');
% ylabel('Norm. Gain Amplitude');
% title('LOS Normalized Gain per T_o');
% %axis([current-0.1,current+0.1,0,550])
% 
% figure(103)
% plot(current-0.1:0.001:current+0.1,HscatGain)
% xlabel('ICI ^T_o');
% ylabel('Gain Amplitude');
% title('ICI Gain per ^T_o');
% xlim([current-0.1,current+0.1])
% axis 'auto y'
% 
% figure(104)
% plot(current-0.1:0.001:current+0.1,HscatGain2)
% xlabel('ICI ^T_o');
% ylabel('Norm. Gain Amplitude');
% title('ICI Normalized Gain per ^T_o');
% axis([current-0.1,current+0.1,0,100])

for cnta=1:cnt0
    for CntNumbSeg=1:zcobj.NumbSeg
        Rxp_temp(:,CntNumbSeg)=(Rxp(:,CntNumbSeg)-((cazacobj.FSymbol_g0(CntNumbSeg,:).').*Hscat(:,cnta).*delayf1(:,cnta)));      
        %Rxp_temp2(:,CntNumbSeg)=Rxp_temp(:,CntNumbSeg).*delayi1(:,cnta).*W_gain(:,1)*(1/Hlos(1,cnta))*PowerGain2_eq(CntNumbSeg,1);%W6 V2A
        Rxp_temp2(:,CntNumbSeg)=Rxp_temp(:,CntNumbSeg).*delayi1(:,cnta)*(1/Hlos(1,cnta));%W6 V2A
        adj_Rxp(1:zcobj.SegSize,CntNumbSeg)=Rxp_temp2(cazacobj.bins_gain(1:zcobj.SegSize),CntNumbSeg);
        adj_Rx(((CntNumbSeg-1)*zcobj.SegSize)+1:CntNumbSeg*zcobj.SegSize,1)=adj_Rxp(1:zcobj.SegSize,CntNumbSeg);
        %Fg0p(1:zcobj.SegSize,CntNumbSeg)=zcobj.FSymbol_g0(CntNumbSeg,cazacobj.bins_gain(1:zcobj.SegSize))*Hlosorg(1,1).*(W_gain(cazacobj.bins_gain(1:zcobj.SegSize),1).')*(1/Hlosorg(1,1))*PowerGain2_eq(CntNumbSeg,1);%W6 V2A
        Fg0p(1:zcobj.SegSize,CntNumbSeg)=zcobj.FSymbol_g0(CntNumbSeg,cazacobj.bins_gain(1:zcobj.SegSize))*Hlosorg(1,1)*(1/Hlosorg(1,1));%W6 V2A
        Fg0(((CntNumbSeg-1)*zcobj.SegSize)+1:CntNumbSeg*zcobj.SegSize,1)=Fg0p(1:zcobj.SegSize,CntNumbSeg);
    end
    tempc=adj_Rx(1:loc4g.NREf,1);
    adj_RxGain(1,cnta)=sum(tempc.*conj(tempc));
    adj_RxPG=sqrt((loc4g.NREf)/adj_RxGain(1,cnta));
    
    if adj_RxGain(1,cnta)>loc4g.NREf
        adj_Rxpn(1:loc4g.NREf,1)=adj_Rx(1:loc4g.NREf,1)*adj_RxPG;
    else
        adj_Rxpn(1:loc4g.NREf,1)=adj_Rx(1:loc4g.NREf,1);
    end
    adj_RxGain2(:,cnta)=adj_Rxpn(1:loc4g.NREf,1).*conj(adj_Rxpn(1:loc4g.NREf,1));
    adj_RxPG2=sqrt((1)./adj_RxGain2(:,cnta));
    
    if adj_RxGain(1,cnta)>loc4g.NREf
        adj_Rxpn(1:loc4g.NREf,1)=adj_Rx(1:loc4g.NREf,1)*adj_RxPG.*adj_RxPG2;
    else
        adj_Rxpn(1:loc4g.NREf,1)=adj_Rx(1:loc4g.NREf,1);
    end
    
    Temp_mat=[Fg0(1:loc4g.NREf),adj_Rxpn(1:loc4g.NREf)];
    saved_adj_Rxp(:,cnta)=adj_Rxpn;
    saved_Fg0(:,cnta)=Fg0(1:loc4g.NREf,1);
    zc2=adj_Rxpn;
    tau=0;
    xax(cnta)=tau;
    cr=sum(conj(Fg0(1:loc4g.NREf)).*zc2([rem(tau+loc4g.NREf,loc4g.NREf)+1:loc4g.NREf,1:rem(tau+loc4g.NREf,loc4g.NREf)]));
    sumvec(cnta)=real(cr);
    saved_cr(cnta,1:length(cr))=cr(1:length(cr));
    cnta=cnta+1;
end

[AxValue BxIndex]=sort(sumvec,'descend');
index_max_phase3=find(sumvec(1,:)==max(sumvec(1,:)),1,'last');
comp_mat=[saved_Fg0(:,100+1),Ocazac,saved_adj_Rxp(:,100+1),saved_adj_Rxp(:,index_max_phase3)];
firstphase=(index_max_phase-2)+(((index_max_phase2-1)*0.1)-0.1)+((index_max_phase3-1)*0.001);
% figure(105)
% stem(0,saved_cr(100+1,:))
% xlabel('-N/2:N/2+1 Correlation Resource Element');
% ylabel('Real Component Amplitude');
% title('Actual Delay Correlation');
% axis([-loc4g.NREf/2 loc4g.NREf/2 0 max(sumvec)])
% 
% figure(106)
% stem(0,saved_cr(index_max_phase3,:))
% xlabel('-N/2:N/2+1 Correlation Resource Element');
% ylabel('Gain Amplitude');
% title('Estimated Delay Correlation');
% axis([-loc4g.NREf/2 loc4g.NREf/2 0 max(sumvec)])
% 
% figure(107)
% plot(current-0.1:0.001:current+0.1,sumvec)
% xlabel('Correlation Zero Flag for every ^To');
% ylabel('Zero Flag Amplitude');
% title('Zero Flag Real Value Amplitude');
% xlim([current-0.1,current+0.1])
% axis 'auto y'

if strcmp(simode,'AWGN')
    %figure(12)
    figure('visible','off');
    plot(1:(cnta-1),sumvec(1:(cnta-1)));
    xlabel('Correlation Zero Flag Value per Phase Change');
    ylabel('Amplitude - Zero Flag');
    title({['2nd Est. Phase - Phase=',num2str(firstphase),',   Real CR Samp=',num2str(defdel)];['Def. Phase=',num2str(defdel),'   Sim=',num2str(mm),'   Tower=',num2str(TowInd)]})
    saveas(gcf, fldhomecombfig9);
    %print(gcf, '-djpeg', fldhomecombfig11,'-r500');
else
    %figure(12)
    figure('visible','off');
    plot(1:(cnta-1),sumvec(1:(cnta-1)));
    xlabel('Correlation Zero Flag Value per Phase Change');
    ylabel('Amplitude - Zero Flag');
    title({['2nd Est. Phase - Phase=',num2str(firstphase),',   Real CR Samp=',num2str(defdel)];['Def. Phase=',num2str(defdel),'   Sim=',num2str(mm),'   Tower=',num2str(TowInd)]})
    saveas(gcf, fldhomecombfig9);
end

calculated_delay=(index_max_phase-2)+(((index_max_phase2-1)*0.1)-0.1)+((index_max_phase3-1)*0.001);
calculated_time=calculated_delay*loc4g.ts2;

varargout{1}=calculated_time;
varargout{2}=calculated_delay;