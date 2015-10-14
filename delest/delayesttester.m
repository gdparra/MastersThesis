%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main Simulation Script (Complete Scrpit with testers Simulation W5 V1A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varargout] = delayest(timedataout,loc4g,cazacobj,zcobj,vn,SignalPower,PowerGain,defdel,nvar_ratio,Total_pow_vn,fldhomecombfig11,fldhomecombfig12,mm,TowInd,simode,comp_table_h,RiceanEst)
cp_inc_samples=(loc4g.CP/0.1)+1;
max_samples=1;

% Water Filling Equalizer
for tkt=1:loc4g.Nfft2
    W_gain(tkt,1)=(conj(PowerGain(tkt))/((conj(PowerGain(tkt))*(PowerGain(tkt)))));
end
W_gain(isnan(W_gain))=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference Symbols
Ocazac(1:loc4g.NREf,1)=cazacobj.REarr(1:loc4g.NREf,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OFDM Symbol Rx
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Estimated Channel
h00NDel=RiceanEst.h00NDel(TowInd+1,1:end);
HDel=fft(h00NDel,length(h00NDel));
Hlosorg=fft(RiceanEst.h00Nlos,loc4g.Nfft2).';
Hscatorg=fft(RiceanEst.h00Nscat,loc4g.Nfft2).';
Horg=fft(RiceanEst.h00N,loc4g.Nfft2).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt0=0;
for timedelay=0:length(h00NDel)-1
    cnt0=cnt0+1;
    for loopN=0:loc4g.Nfft2-1
        delayf1(loopN+1,cnt0)=exp(-1i*2*pi*loopN*(timedelay/loc4g.Nfft2));
        delayi1(loopN+1,cnt0)=exp(1i*2*pi*loopN*(timedelay/loc4g.Nfft2));
    end
    for loopN=0:length(h00NDel)-1
        Hest(loopN+1,cnt0)=HDel(loopN+1)*exp(1i*2*pi*loopN*(timedelay/length(h00NDel)));
    end
    
    % Section added to Normalize Channel and set Kdb on LOS on est. channel
%     h00Ntemp=ifft(Hest(:,cnt0),length(h00NDel));
%     h00N01(1)=h00Ntemp(1);
%     h00Nscattemp = [0;h00Ntemp(2:end)];
%     GLos=sqrt(10^(loc4g.riciankdB/10)*(sum(h00Nscattemp.*conj(h00Nscattemp)))/...
%         (h00Ntemp(1)*conj(h00Ntemp(1))));
%     h00Nlostemp =  [GLos*h00Ntemp(1);zeros(size(h00Ntemp(2:end)))];
%     h00Ntemp(1)=h00Nlostemp(1);
%     h00Ntemp2 = h00Ntemp/sqrt(sum(h00Ntemp.*conj(h00Ntemp)));
%     h00Ntemp2(1)=h00N01(1);
%     Hest(:,cnt0)=fft(h00Ntemp2,length(h00NDel));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    HGain(cnt0,1)=sum(Hest(:,cnt0).*conj(Hest(:,cnt0)));
    test=mean(abs(Hest(:,cnt0)-RiceanEst.H.'));
    
    h00Nest=ifft(Hest(:,cnt0),length(h00NDel));
    hGain(:,cnt0)=h00Nest.*conj(h00Nest);
    h00Nscatest=[0;h00Nest(2:end,1)];
    test=mean(abs(h00Nscatest-RiceanEst.h00Nscat.'));
    
    Hscatest=fft(h00Nscatest,loc4g.Nfft2);
    HscatGain(cnt0,1)=sum(Hscatest.*conj(Hscatest));
    HscatPG=sqrt((70.045)/HscatGain(cnt0,1));
    Hscat(:,cnt0)=Hscatest*HscatPG;
    HscatGain2(cnt0,1)=sum(Hscat(:,cnt0).*conj(Hscat(:,cnt0)));

    h00Nlosest=[h00Nest(1,1),zeros(1,length(h00Nest)-1)];
    test=mean(abs(h00Nlosest-RiceanEst.h00Nlos));
    
    Hlosest=fft(h00Nlosest,loc4g.Nfft2);
    HlosGain(cnt0,1)=sum(Hlosest.*conj(Hlosest));
    HlosPG=sqrt((441.95)/HlosGain(cnt0,1));
    Hlos(:,cnt0)=(Hlosest.')*HlosPG;
    HlosGain2(cnt0,1)=sum(Hlos(:,cnt0).*conj(Hlos(:,cnt0)));
end

% Graphs
figure(100)
plot(0:cnt0-1,HGain)
xlabel('H ^T_o');
ylabel('Gain Amplitude');
title('Ricean H Gain per ^T_o');
axis([0,160,0,200])

figure(101)
plot(0:cnt0-1,HlosGain)
xlabel('LOS ^T_o');
ylabel('Gain Amplitude');
title('LOS Gain per ^T_o');
axis([0,160,0,550])

figure(102)
plot(0:cnt0-1,HlosGain2)
xlabel('LOS ^T_o');
ylabel('Norm. Gain Amplitude');
title('LOS Normalized Gain per T_o');
axis([0,160,0,550])

figure(103)
plot(0:cnt0-1,HscatGain)
xlabel('ICI ^T_o');
ylabel('Gain Amplitude');
title('ICI Gain per ^T_o');
axis([0,160,0,550])

figure(104)
plot(0:cnt0-1,HscatGain2)
xlabel('ICI ^T_o');
ylabel('Norm. Gain Amplitude');
title('ICI Normalized Gain per ^T_o');
axis([0,160,0,100])

for cnta=1:cnt0
    
    for CntNumbSeg=1:zcobj.NumbSeg
        Rxp_temp(:,CntNumbSeg)=(Rxp(:,CntNumbSeg)-((cazacobj.FSymbol_g0(CntNumbSeg,:).').*Hscat(:,cnta).*delayf1(:,cnta)));%.*delayi1(:,cnta);%.*W_gain(:,1);
        ICI(:,CntNumbSeg)=((cazacobj.FSymbol_g0(CntNumbSeg,:).').*Hscat(:,cnta).*delayf1(:,cnta));
        %Testers
        Rxporg(:,CntNumbSeg)=((cazacobj.FSymbol_g0(CntNumbSeg,:).').*Hlos(:,cnta).*delayf1(:,cnta))+...
               ((cazacobj.FSymbol_g0(CntNumbSeg,:).').*Hscat(:,cnta).*delayf1(:,cnta));
        matrix0=[Rxp(zcobj.bins_gain(1:zcobj.SegSize),CntNumbSeg),Rxporg(zcobj.bins_gain(1:zcobj.SegSize),CntNumbSeg)];
        lostest=((cazacobj.FSymbol_g0(CntNumbSeg,:).').*Hlos(:,cnta).*delayf1(:,cnta));
        matrix1=[lostest(zcobj.bins_gain(1:zcobj.SegSize),1),Rxp_temp(zcobj.bins_gain(1:zcobj.SegSize),1)];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Rxp_temp2(:,CntNumbSeg)=Rxp_temp(:,CntNumbSeg).*delayi1(:,cnta).*W_gain(:,1)*(1/Hlos(1,cnta))*PowerGain2_eq(CntNumbSeg,1);
        
        % Testers
        Rxpsim(:,CntNumbSeg)=(cazacobj.FSymbol_g0(CntNumbSeg,:).').*Hlos(:,defdel+1).*delayf1(:,defdel+1)+...
                             (cazacobj.FSymbol_g0(CntNumbSeg,:).').*Hscat(:,defdel+1).*delayf1(:,defdel+1);
        mtvlim=zcobj.SegSize;%zcobj.SegSize%zcobj.NREf
        mtv(1:mtvlim,1)=(cazacobj.FSymbol_g0(CntNumbSeg,zcobj.bins_gain(1:mtvlim)).');%CAZAC
        mtv(1:mtvlim,2)=Hlos(zcobj.bins_gain(1:mtvlim),defdel+1);%LOS
        mtv(1:mtvlim,3)=delayf1(zcobj.bins_gain(1:mtvlim),defdel+1);%Defdel
        mtv(1:mtvlim,4)=mtv(1:mtvlim,1);%CAZAC
        mtv(1:mtvlim,5)=Hscat(zcobj.bins_gain(1:mtvlim),defdel+1);%SCAT
        mtv(1:mtvlim,6)=mtv(1:mtvlim,3);%Defdel
        mtv(1:mtvlim,7)=mtv(1:mtvlim,1).*mtv(1:mtvlim,2).*mtv(1:mtvlim,3);%LOS Component
        mtv(1:mtvlim,8)=mtv(1:mtvlim,4).*mtv(1:mtvlim,5).*mtv(1:mtvlim,6);%ICI
        mtv(1:mtvlim,9)=mtv(1:mtvlim,7)+mtv(1:mtvlim,8);%LOS Component+ICI = Rxp
        mtv(1:mtvlim,10)=(cazacobj.FSymbol_g0(CntNumbSeg,zcobj.bins_gain(1:mtvlim)).');%CAZAC
        mtv(1:mtvlim,11)=Hscat(zcobj.bins_gain(1:mtvlim),cnta);%SCAT
        mtv(1:mtvlim,12)=delayf1(zcobj.bins_gain(1:mtvlim),cnta);%^Tf
        mtv(1:mtvlim,13)=mtv(1:mtvlim,10).*mtv(1:mtvlim,11).*mtv(1:mtvlim,12);
        mtv(1:mtvlim,14)=mtv(1:mtvlim,9)-mtv(1:mtvlim,13);%Rxp-RxpSCAT Component ^Tf
        mtv(1:mtvlim,15)=mtv(1:mtvlim,7)-mtv(1:mtvlim,14);%LOS Component - RxpLOS Component ^Tf
        mtv(1:mtvlim,16)=mtv(1:mtvlim,14).*delayi1(zcobj.bins_gain(1:mtvlim),cnta);%(RxpLOS Component ^Tf)*^Ti
        mtv(1:mtvlim,17)=(cazacobj.FSymbol_g0(CntNumbSeg,zcobj.bins_gain(1:mtvlim)).').*Hlos(zcobj.bins_gain(1:mtvlim),defdel+1);%CAZAC*LOS
        mtv(1:mtvlim,18)=mtv(1:mtvlim,16)-mtv(1:mtvlim,17);%(RxpLOS Component ^Tf)*^Ti-CAZAC*LOS
        mtv2(1,cnta)=sum(mtv(1:mtvlim,1).*conj(mtv(1:mtvlim,1)));
        mtv2(2,cnta)=sum(mtv(1:mtvlim,2).*conj(mtv(1:mtvlim,2)));
        mtv2(3,cnta)=sum(mtv(1:mtvlim,3).*conj(mtv(1:mtvlim,3)));
        mtv2(4,cnta)=sum(mtv(1:mtvlim,4).*conj(mtv(1:mtvlim,4)));
        mtv2(5,cnta)=sum(mtv(1:mtvlim,5).*conj(mtv(1:mtvlim,5)));
        mtv2(6,cnta)=sum(mtv(1:mtvlim,6).*conj(mtv(1:mtvlim,6)));
        mtv2(7,cnta)=sum(mtv(1:mtvlim,7).*conj(mtv(1:mtvlim,7)));
        mtv2(8,cnta)=sum(mtv(1:mtvlim,8).*conj(mtv(1:mtvlim,8)));
        mtv2(9,cnta)=sum(mtv(1:mtvlim,9).*conj(mtv(1:mtvlim,9)));
        mtv2(10,cnta)=sum(mtv(1:mtvlim,10).*conj(mtv(1:mtvlim,10)));
        mtv2(11,cnta)=sum(mtv(1:mtvlim,11).*conj(mtv(1:mtvlim,11)));
        mtv2(12,cnta)=sum(mtv(1:mtvlim,12).*conj(mtv(1:mtvlim,12)));
        mtv2(13,cnta)=sum(mtv(1:mtvlim,13).*conj(mtv(1:mtvlim,13)));
        mtv2(14,cnta)=sum(mtv(1:mtvlim,14).*conj(mtv(1:mtvlim,14)));
        mtv2(15,cnta)=sum(mtv(1:mtvlim,15).*conj(mtv(1:mtvlim,15)));
        mtv2(16,cnta)=sum(mtv(1:mtvlim,16).*conj(mtv(1:mtvlim,16)));
        mtv2(17,cnta)=sum(mtv(1:mtvlim,17).*conj(mtv(1:mtvlim,17)));
        mtv2(18,cnta)=sum(mtv(1:mtvlim,18).*conj(mtv(1:mtvlim,18)));
        mtv3(1,1)=find(mtv2(1,:)==min(mtv2(1,:)),1,'first');
        mtv3(1,2)=find(mtv2(1,:)==max(mtv2(1,:)),1,'last');
        mtv3(2,1)=find(mtv2(2,:)==min(mtv2(2,:)),1,'first');
        mtv3(2,2)=find(mtv2(2,:)==max(mtv2(2,:)),1,'last');
        mtv3(3,1)=find(mtv2(3,:)==min(mtv2(3,:)),1,'first');
        mtv3(3,2)=find(mtv2(3,:)==max(mtv2(3,:)),1,'last');
        mtv3(4,1)=find(mtv2(4,:)==min(mtv2(4,:)),1,'first');
        mtv3(4,2)=find(mtv2(4,:)==max(mtv2(4,:)),1,'last');
        mtv3(5,1)=find(mtv2(5,:)==min(mtv2(5,:)),1,'first');
        mtv3(5,2)=find(mtv2(5,:)==max(mtv2(5,:)),1,'last');
        mtv3(6,1)=find(mtv2(6,:)==min(mtv2(6,:)),1,'first');
        mtv3(6,2)=find(mtv2(6,:)==max(mtv2(6,:)),1,'last');
        mtv3(7,1)=find(mtv2(7,:)==min(mtv2(7,:)),1,'first');
        mtv3(7,2)=find(mtv2(7,:)==max(mtv2(7,:)),1,'last');
        mtv3(8,1)=find(mtv2(8,:)==min(mtv2(8,:)),1,'first');
        mtv3(8,2)=find(mtv2(8,:)==max(mtv2(8,:)),1,'last');
        mtv3(9,1)=find(mtv2(9,:)==min(mtv2(9,:)),1,'first');
        mtv3(9,2)=find(mtv2(9,:)==max(mtv2(9,:)),1,'last');
        mtv3(10,1)=find(mtv2(10,:)==min(mtv2(10,:)),1,'first');
        mtv3(10,2)=find(mtv2(10,:)==max(mtv2(10,:)),1,'last');
        mtv3(11,1)=find(mtv2(11,:)==min(mtv2(11,:)),1,'first');
        mtv3(11,2)=find(mtv2(11,:)==max(mtv2(11,:)),1,'last');
        mtv3(12,1)=find(mtv2(12,:)==min(mtv2(12,:)),1,'first');
        mtv3(12,2)=find(mtv2(12,:)==max(mtv2(12,:)),1,'last');
        mtv3(13,1)=find(mtv2(13,:)==min(mtv2(13,:)),1,'first');
        mtv3(13,2)=find(mtv2(13,:)==max(mtv2(13,:)),1,'last');
        mtv3(14,1)=find(mtv2(14,:)==min(mtv2(14,:)),1,'first');
        mtv3(14,2)=find(mtv2(14,:)==max(mtv2(14,:)),1,'last');
        mtv3(15,1)=find(mtv2(15,:)==min(mtv2(15,:)),1,'first');
        mtv3(15,2)=find(mtv2(15,:)==max(mtv2(15,:)),1,'last');
        mtv3(16,1)=find(mtv2(16,:)==min(mtv2(16,:)),1,'first');
        mtv3(16,2)=find(mtv2(16,:)==max(mtv2(16,:)),1,'last');
        mtv3(17,1)=find(mtv2(17,:)==min(mtv2(17,:)),1,'first');
        mtv3(17,2)=find(mtv2(17,:)==max(mtv2(17,:)),1,'last');
        mtv3(18,1)=find(mtv2(18,:)==min(mtv2(18,:)),1,'first');
        mtv3(18,2)=find(mtv2(18,:)==max(mtv2(18,:)),1,'last');
        SCATcomp(:,1)=mtv(:,11);
        SCATcomp(:,2)=mtv(:,11).*mtv(:,12);
        SCATgaincomp(1,1)=sum(SCATcomp(:,1).*conj(SCATcomp(:,1)));
        SCATgaincomp(1,2)=sum(SCATcomp(:,2).*conj(SCATcomp(:,2)));
        
        matrix2=[Rxp(zcobj.bins_gain(1:zcobj.SegSize),CntNumbSeg),Rxpsim(zcobj.bins_gain(1:zcobj.SegSize),CntNumbSeg),Rxporg(zcobj.bins_gain(1:zcobj.SegSize),CntNumbSeg)];
        
        indest=Rxp_temp2(zcobj.bins_gain(1:zcobj.SegSize),CntNumbSeg);
        indorg=(cazacobj.FSymbol_g0(CntNumbSeg,zcobj.bins_gain(1:zcobj.SegSize)).')*Hlos(1,cnta).*W_gain(zcobj.bins_gain(1:zcobj.SegSize),1)*(1/Hlos(1,cnta));%*PowerGain2_eq(CntNumbSeg,1);
        matrix3=[indest,indorg];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        adj_Rxp(1:zcobj.SegSize,CntNumbSeg)=Rxp_temp2(cazacobj.bins_gain(1:zcobj.SegSize),CntNumbSeg);
        adj_Rx(((CntNumbSeg-1)*zcobj.SegSize)+1:CntNumbSeg*zcobj.SegSize,1)=adj_Rxp(1:zcobj.SegSize,CntNumbSeg);
        
        Fg0p(1:zcobj.SegSize,CntNumbSeg)=zcobj.FSymbol_g0(CntNumbSeg,cazacobj.bins_gain(1:zcobj.SegSize))*Hlosorg(1,1).*(W_gain(cazacobj.bins_gain(1:zcobj.SegSize),1).')*(1/Hlosorg(1,1))*PowerGain2_eq(CntNumbSeg,1);
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
    %adj_Rxpn(1:loc4g.NREf,1)=adj_Rx(1:loc4g.NREf,1);
    adj_RxPG2=sqrt((1)./adj_RxGain2(:,cnta));
    
    if adj_RxGain(1,cnta)>loc4g.NREf
        adj_Rxpn(1:loc4g.NREf,1)=adj_Rx(1:loc4g.NREf,1)*adj_RxPG.*adj_RxPG2;
    else
        adj_Rxpn(1:loc4g.NREf,1)=adj_Rx(1:loc4g.NREf,1);
    end
    Temp_mat=[Fg0(1:loc4g.NREf),adj_Rxpn(1:loc4g.NREf)];
    Abs_mat(cnta,1)=mean(abs(Temp_mat(:,1)-Temp_mat(:,2)));
    saved_adj_Rxp(:,cnta)=adj_Rxpn;
    saved_Fg0(:,cnta)=Fg0(1:loc4g.NREf,1);
    zc2=adj_Rxpn;
    
%    k=0;
%     for tau=-loc4g.NREf/2:loc4g.NREf/2
%         k=k+1;
%         xax(k)=tau;
%         cr(k)=sum(conj(Fg0(1:loc4g.NREf)).*zc2([rem(tau+loc4g.NREf,loc4g.NREf)+1:loc4g.NREf,1:rem(tau+loc4g.NREf,loc4g.NREf)]));
%     end %Fg0, zc2

    tau=0;
    xax(cnta)=tau;
    cr=sum(conj(Fg0(1:loc4g.NREf)).*zc2([rem(tau+loc4g.NREf,loc4g.NREf)+1:loc4g.NREf,1:rem(tau+loc4g.NREf,loc4g.NREf)]));
    
    %sumvec(cnta)=real(cr(loc4g.NREf/2+1));
    sumvec(cnta)=real(cr);
    saved_cr(cnta,1:length(cr))=cr(1:length(cr));
    cnta=cnta+1;
end

[AxValue BxIndex]=sort(sumvec,'descend');
index_max_phase=find(sumvec(1,:)==max(sumvec(1,:)),1,'last');
comp_mat=[saved_Fg0(:,round(defdel)+1),Ocazac,saved_adj_Rxp(:,round(defdel)+1),saved_adj_Rxp(:,index_max_phase)];

figure(105)
stem(0,saved_cr(round(defdel)+1,:))
xlabel('-N/2:N/2+1 Correlation Resource Element');
ylabel('Real Component Amplitude');
title('Actual Delay Correlation');
axis([-loc4g.NREf/2 loc4g.NREf/2 0 max(sumvec)])

figure(106)
stem(0,saved_cr(index_max_phase,:))
xlabel('-N/2:N/2+1 Correlation Resource Element');
ylabel('Gain Amplitude');
title('Estimated Delay Correlation');
axis([-loc4g.NREf/2 loc4g.NREf/2 0 max(sumvec)])

figure(107)
plot(0:cnta-2,sumvec)
xlabel('Correlation Zero Flag for every ^To');
ylabel('Zero Flag Amplitude');
title('Zero Flag Real Value Amplitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars cnt0 timedelay delayf1 delayi1 Hest HGain h00Nest hGain 
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
    for loopN=0:length(h00NDel)-1
        Hest(loopN+1,cnt0)=HDel(loopN+1)*exp(1i*2*pi*loopN*(timedelay/length(h00NDel)));
    end
       
    HGain(cnt0,1)=sum(Hest(:,cnt0).*conj(Hest(:,cnt0)));
    
    h00Nest=ifft(Hest(:,cnt0),length(h00NDel));
    hGain(:,cnt0)=h00Nest.*conj(h00Nest);
    h00Nscatest=[0;h00Nest(2:end,1)];
    
    Hscatest=fft(h00Nscatest,loc4g.Nfft2);
    HscatGain(cnt0,1)=sum(Hscatest.*conj(Hscatest));
    HscatPG=sqrt((70.045)/HscatGain(cnt0,1));
    Hscat(:,cnt0)=Hscatest*HscatPG;
    HscatGain2(cnt0,1)=sum(Hscat(:,cnt0).*conj(Hscat(:,cnt0)));

    h00Nlosest=[h00Nest(1,1),zeros(1,length(h00Nest)-1)];
    
    Hlosest=fft(h00Nlosest,loc4g.Nfft2);
    HlosGain(cnt0,1)=sum(Hlosest.*conj(Hlosest));
    HlosPG=sqrt((441.95)/HlosGain(cnt0,1));
    Hlos(:,cnt0)=(Hlosest.')*HlosPG;
    HlosGain2(cnt0,1)=sum(Hlos(:,cnt0).*conj(Hlos(:,cnt0)));
end

% Graphs
figure(100)
plot(current-1:0.1:current+1,HGain)
xlabel('H ^T_o');
ylabel('Gain Amplitude');
title('Ricean H Gain per ^T_o');
axis([current-1,current+1,0,200])

figure(101)
plot(current-1:0.1:current+1,HlosGain)
xlabel('LOS ^T_o');
ylabel('Gain Amplitude');
title('LOS Gain per ^T_o');
axis([current-1,current+1,0,550])

figure(102)
plot(current-1:0.1:current+1,HlosGain2)
xlabel('LOS ^T_o');
ylabel('Norm. Gain Amplitude');
title('LOS Normalized Gain per T_o');
axis([current-1,current+1,0,550])

figure(103)
plot(current-1:0.1:current+1,HscatGain)
xlabel('ICI ^T_o');
ylabel('Gain Amplitude');
title('ICI Gain per ^T_o');
axis([current-1,current+1,0,550])

figure(104)
plot(current-1:0.1:current+1,HscatGain2)
xlabel('ICI ^T_o');
ylabel('Norm. Gain Amplitude');
title('ICI Normalized Gain per ^T_o');
axis([current-1,current+1,0,100])

for cnta=1:cnt0
    for CntNumbSeg=1:zcobj.NumbSeg
        Rxp_temp(:,CntNumbSeg)=(Rxp(:,CntNumbSeg)-((cazacobj.FSymbol_g0(CntNumbSeg,:).').*Hscat(:,cnta).*delayf1(:,cnta)));      
        Rxp_temp2(:,CntNumbSeg)=Rxp_temp(:,CntNumbSeg).*delayi1(:,cnta).*W_gain(:,1)*(1/Hlos(1,cnta))*PowerGain2_eq(CntNumbSeg,1);     
        adj_Rxp(1:zcobj.SegSize,CntNumbSeg)=Rxp_temp2(cazacobj.bins_gain(1:zcobj.SegSize),CntNumbSeg);
        adj_Rx(((CntNumbSeg-1)*zcobj.SegSize)+1:CntNumbSeg*zcobj.SegSize,1)=adj_Rxp(1:zcobj.SegSize,CntNumbSeg);
        Fg0p(1:zcobj.SegSize,CntNumbSeg)=zcobj.FSymbol_g0(CntNumbSeg,cazacobj.bins_gain(1:zcobj.SegSize))*Hlosorg(1,1).*(W_gain(cazacobj.bins_gain(1:zcobj.SegSize),1).')*(1/Hlosorg(1,1))*PowerGain2_eq(CntNumbSeg,1);
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

figure(105)
stem(0,saved_cr(10+1,:))
xlabel('-N/2:N/2+1 Correlation Resource Element');
ylabel('Real Component Amplitude');
title('Actual Delay Correlation');
axis([-loc4g.NREf/2 loc4g.NREf/2 0 max(sumvec)])

figure(106)
stem(0,saved_cr(index_max_phase2,:))
xlabel('-N/2:N/2+1 Correlation Resource Element');
ylabel('Gain Amplitude');
title('Estimated Delay Correlation');
axis([-loc4g.NREf/2 loc4g.NREf/2 0 max(sumvec)])

figure(107)
plot(current-1:0.1:current+1,sumvec)
xlabel('Correlation Zero Flag for every ^To');
ylabel('Zero Flag Amplitude');
title('Zero Flag Real Value Amplitude');
axis([current-1 current+1 0 max(sumvec)])

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
    for loopN=0:length(h00NDel)-1
        Hest(loopN+1,cnt0)=HDel(loopN+1)*exp(1i*2*pi*loopN*(timedelay/length(h00NDel)));
    end
       
    HGain(cnt0,1)=sum(Hest(:,cnt0).*conj(Hest(:,cnt0)));
    
    h00Nest=ifft(Hest(:,cnt0),length(h00NDel));
    hGain(:,cnt0)=h00Nest.*conj(h00Nest);
    h00Nscatest=[0;h00Nest(2:end,1)];
    
    Hscatest=fft(h00Nscatest,loc4g.Nfft2);
    HscatGain(cnt0,1)=sum(Hscatest.*conj(Hscatest));
    HscatPG=sqrt((70.045)/HscatGain(cnt0,1));
    Hscat(:,cnt0)=Hscatest*HscatPG;
    HscatGain2(cnt0,1)=sum(Hscat(:,cnt0).*conj(Hscat(:,cnt0)));

    h00Nlosest=[h00Nest(1,1),zeros(1,length(h00Nest)-1)];
    
    Hlosest=fft(h00Nlosest,loc4g.Nfft2);
    HlosGain(cnt0,1)=sum(Hlosest.*conj(Hlosest));
    HlosPG=sqrt((441.95)/HlosGain(cnt0,1));
    Hlos(:,cnt0)=(Hlosest.')*HlosPG;
    HlosGain2(cnt0,1)=sum(Hlos(:,cnt0).*conj(Hlos(:,cnt0)));
end

% Graphs
figure(100)
plot(current-0.1:0.001:current+0.1,HGain)
xlabel('H ^T_o');
ylabel('Gain Amplitude');
title('Ricean H Gain per ^T_o');
axis([current-0.1,current+0.1,0,200])

figure(101)
plot(current-0.1:0.001:current+0.1,HlosGain)
xlabel('LOS ^T_o');
ylabel('Gain Amplitude');
title('LOS Gain per ^T_o');
xlim([current-0.1,current+0.1])
axis 'auto y'

figure(102)
plot(current-0.1:0.001:current+0.1,HlosGain2)
xlabel('LOS ^T_o');
ylabel('Norm. Gain Amplitude');
title('LOS Normalized Gain per T_o');
axis([current-0.1,current+0.1,0,550])

figure(103)
plot(current-0.1:0.001:current+0.1,HscatGain)
xlabel('ICI ^T_o');
ylabel('Gain Amplitude');
title('ICI Gain per ^T_o');
xlim([current-0.1,current+0.1])
axis 'auto y'

figure(104)
plot(current-0.1:0.001:current+0.1,HscatGain2)
xlabel('ICI ^T_o');
ylabel('Norm. Gain Amplitude');
title('ICI Normalized Gain per ^T_o');
axis([current-0.1,current+0.1,0,100])

for cnta=1:cnt0
    for CntNumbSeg=1:zcobj.NumbSeg
        Rxp_temp(:,CntNumbSeg)=(Rxp(:,CntNumbSeg)-((cazacobj.FSymbol_g0(CntNumbSeg,:).').*Hscat(:,cnta).*delayf1(:,cnta)));      
        Rxp_temp2(:,CntNumbSeg)=Rxp_temp(:,CntNumbSeg).*delayi1(:,cnta).*W_gain(:,1)*(1/Hlos(1,cnta))*PowerGain2_eq(CntNumbSeg,1);     
        adj_Rxp(1:zcobj.SegSize,CntNumbSeg)=Rxp_temp2(cazacobj.bins_gain(1:zcobj.SegSize),CntNumbSeg);
        adj_Rx(((CntNumbSeg-1)*zcobj.SegSize)+1:CntNumbSeg*zcobj.SegSize,1)=adj_Rxp(1:zcobj.SegSize,CntNumbSeg);
        Fg0p(1:zcobj.SegSize,CntNumbSeg)=zcobj.FSymbol_g0(CntNumbSeg,cazacobj.bins_gain(1:zcobj.SegSize))*Hlosorg(1,1).*(W_gain(cazacobj.bins_gain(1:zcobj.SegSize),1).')*(1/Hlosorg(1,1))*PowerGain2_eq(CntNumbSeg,1);
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

figure(105)
stem(0,saved_cr(100+1,:))
xlabel('-N/2:N/2+1 Correlation Resource Element');
ylabel('Real Component Amplitude');
title('Actual Delay Correlation');
axis([-loc4g.NREf/2 loc4g.NREf/2 0 max(sumvec)])

figure(106)
stem(0,saved_cr(index_max_phase3,:))
xlabel('-N/2:N/2+1 Correlation Resource Element');
ylabel('Gain Amplitude');
title('Estimated Delay Correlation');
axis([-loc4g.NREf/2 loc4g.NREf/2 0 max(sumvec)])

figure(107)
plot(current-0.1:0.001:current+0.1,sumvec)
xlabel('Correlation Zero Flag for every ^To');
ylabel('Zero Flag Amplitude');
title('Zero Flag Real Value Amplitude');
xlim([current-0.1,current+0.1])
axis 'auto y'

calculated_delay=(index_max_phase-2)+(((index_max_phase2-1)*0.1)-0.1)+((index_max_phase3-1)*0.001);
calculated_time=calculated_delay*loc4g.ts2;

varargout{1}=calculated_time;
varargout{2}=calculated_delay;
        
