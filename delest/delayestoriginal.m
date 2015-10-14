%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main Simulation Script (Not Working Properly, Return to last week of 2014
%delest for great results or Simulation W5 V1A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varargout] = delayest(timedataout,loc4g,cazacobj,zcobj,vn,SignalPower,PowerGain,defdel,nvar_ratio,Total_pow_vn,fldhomecombfig11,fldhomecombfig12,mm,TowInd,simode,comp_table_h,RiceanEst)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of Window Frame Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cp_inc_samples=(loc4g.CP/0.1)+1;
max_samples=1; %Perfect Detector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of Window Frame Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Registry Allocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saved_sumvec=zeros(max_samples,cp_inc_samples);
saved_abs_sumvec=zeros(max_samples,cp_inc_samples);
max_sumvec=zeros(1,max_samples);
saved_cr=zeros(cp_inc_samples,loc4g.NREf+1);
sumvec=zeros(1,cp_inc_samples);
abs_sumvec=zeros(1,cp_inc_samples);
kind=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of Registry Allocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tkt=1:loc4g.Nfft2
    W_gain2(tkt,1)=(conj(PowerGain(tkt)*cazacobj.PowerGain2(1,1))/((conj(PowerGain(tkt)*cazacobj.PowerGain2(1,1))*(PowerGain(tkt)*cazacobj.PowerGain2(1,1)))));
end

for tkt=1:loc4g.Nfft2
    W_gain(tkt,1)=(conj(PowerGain(tkt))/((conj(PowerGain(tkt))*(PowerGain(tkt)))));
end
W_gain(isnan(W_gain))=0;
W_gain2(isnan(W_gain))=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reference Symbols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ocazac(1:loc4g.NREf,1)=cazacobj.REarr(1:loc4g.NREf,1);% % Change NREf for SegSize
OIcazac(1:loc4g.NREf,1)=[cazacobj.REarr(loc4g.NREf+loc4g.LI+1:loc4g.NREf,1);cazacobj.REarr(1:loc4g.UI,1)];
%g0=zeros(1,loc4g.Nfft);
%g0(cazacobj.bins_gain(1:cazacobj.NREf))=cazacobj.REarr(1:cazacobj.NREf);% Change NREf for SegSize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of Reference Symbols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

saved_Rxp(1:loc4g.Nfft2,1:zcobj.NumbSeg)=Rxp(1:loc4g.Nfft2,1:zcobj.NumbSeg);% CntNumbSeg replaces i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Original CAZAC Sequence Retrieval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cnt=0;
for timedelay=0:0.1:233.1
    cnt=cnt+1;
    for loopN=0:loc4g.Nfft2-1
        delaytf(loopN+1,cnt)=exp(-1i*2*pi*loopN*(timedelay/loc4g.Nfft2));
        delayti(loopN+1,cnt)=exp(1i*2*pi*loopN*(timedelay/loc4g.Nfft2));
    end
    
end

for timedelay=0:0.1:233.1
    kind=kind+1;
    k=0;
    
    for CntNumbSeg=1:zcobj.NumbSeg
        
        Hscat=fft(loc4g.h00Nscat(1:end),loc4g.Nfft2);
        h00NDel=RiceanEst.h00NDel(TowInd+1,1:end);
        HDel=fft(h00NDel,length(h00NDel));
        
        for loopN=0:length(h00NDel)-1
            Hest(loopN+1)=HDel(loopN+1)*exp(1i*2*pi*loopN*(timedelay/length(h00NDel)));
        end
        
        h00Nest=ifft(Hest,length(h00NDel));
        h00Nscatest=[0,h00Nest(1,2:end)];
        Hscatest=fft(h00Nscatest,loc4g.Nfft2);
        Hscat=Hscatest;
        
        for loopN=0:loc4g.Nfft2-1
            Hscat_temp(loopN+1)=Hscat(loopN+1)*exp(-1i*2*pi*loopN*(timedelay/loc4g.Nfft2));
        end
        
        Hscat=Hscat_temp;
        
        for Rk=1:loc4g.Nfft2
            Rxp_temp(Rk,CntNumbSeg)=Rxp(Rk,CntNumbSeg)-(cazacobj.FSymbol_g0(CntNumbSeg,Rk)*Hscat(Rk));
             test_fad_ref(Rk,CntNumbSeg)=(cazacobj.FSymbol_g0(CntNumbSeg,Rk)*Hscat(Rk));
%             ratio_symb(Rk,CntNumbSeg)=((sqrt(zcobj.NREf^2)/zcobj.SegSize)/(PowerGain(1,Rk)*zcobj.PowerGain2(CntNumbSeg,1)));
%             Rxp_temp_ref(Rk,CntNumbSeg)=Rxp_temp(Rk,CntNumbSeg)*ratio_symb(Rk,CntNumbSeg);
%             Rxp_temp(Rk,CntNumbSeg)=Rxp_temp_ref(Rk,CntNumbSeg);
        end
        
%         Rxp_temp_ref(isinf(Rxp_temp_ref))=0;
%         Rxp_temp(isinf(Rxp_temp))=0;
%         
        comp_table_fad(:,1)=test_fad_ref(:,CntNumbSeg);
        comp_table_fad(:,2)=comp_table_h(1:512,4);
        comp_table_fad(:,3)=[comp_table_fad(:,2)-comp_table_fad(:,1)];
        comp_table_fad(:,5)=[comp_table_h(:,5)+comp_table_h(:,7)];
        comp_table_fad(:,6)=Rxp_temp(:,CntNumbSeg);
        comp_table_fad(:,7)=comp_table_fad(:,5)-comp_table_fad(:,6);
%         comp_table_fad(:,8)=Rxp_temp_ref(:,CntNumbSeg);
        comp_table_fad(zcobj.bins_gain(1:zcobj.SegSize),9)=1;
        
        tmpvec_Rxp_temp(:,CntNumbSeg)=ifft(Rxp_temp(:,CntNumbSeg),loc4g.Nfft);
        powSig_Rxp_temp(CntNumbSeg,:)=sum(tmpvec_Rxp_temp(1:end,CntNumbSeg).*conj(tmpvec_Rxp_temp(1:end,CntNumbSeg)))/length(tmpvec_Rxp_temp(1:end,CntNumbSeg));
        Frequency_Power_AVG_Rxp_temp(CntNumbSeg,:)=sum((Rxp_temp(:,CntNumbSeg).*conj(Rxp_temp(:,CntNumbSeg))));
        
%         for Rk=1:loc4g.Nfft2
%             ratio_symb(Rk,CntNumbSeg)=((sqrt(zcobj.NREf^2)/zcobj.SegSize)/(PowerGain(1,Rk)*zcobj.PowerGain2(CntNumbSeg,1)));
%         end
%         
        
        Hlos=fft(loc4g.h00Nlos(1:end),loc4g.Nfft2);
%        Hlos=ones(1,loc4g.Nfft2);
        h00Nlosest=[h00Nest(1,1),zeros(1,length(h00Nest)-1)];
        Hlosest=fft(h00Nlosest,loc4g.Nfft2);
        Hlos=Hlosest;
        
        Hscat_indicator_estimated(:,CntNumbSeg)=cazacobj.FSymbolG2(CntNumbSeg,:)*Hlos(1);
        Hscat_indicator_calculated(:,CntNumbSeg)=Rxp_temp(zcobj.bins_gain(1:zcobj.SegSize),CntNumbSeg);
        Hscat_mat=[Hscat_indicator_estimated(:,CntNumbSeg),Hscat_indicator_calculated(:,CntNumbSeg)];
        for loopN=0:loc4g.Nfft2-1
            Hlos_temp(loopN+1)=Hlos(loopN+1)*exp(-1i*2*pi*loopN*(timedelay/loc4g.Nfft2));
        end
        
%        Hlos=Hlos_temp;
% revision Needed         
        vn2(CntNumbSeg,1)=vn(CntNumbSeg,1)*loc4g.Nfft;%sqrt(262144/2)+1i*sqrt(262144/2);
        if strcmp(simode,'AWGN') %Mod W45 V3
            powSig_Rxp_tempv2(CntNumbSeg,1)=zcobj.powSig(CntNumbSeg,1)*loc4g.Nfft;
        elseif strcmp(simode,'AWGN_Fadding') %Mod W45 V3
            powSig_Rxp_tempv2(CntNumbSeg,1)=powSig_Rxp_temp(CntNumbSeg,:)*loc4g.Nfft;
        elseif strcmp(simode,'AWGN_Fadding_Interference') %Mod W45 V3
            powSig_Rxp_tempv2(CntNumbSeg,1)=powSig_Rxp_temp(CntNumbSeg,:)*loc4g.Nfft;
        end
        
        for Rk2=1:loc4g.Nfft2
            Hlos_eq(Rk2,CntNumbSeg)=(conj(Hlos(Rk2))/((vn2(CntNumbSeg,1)/powSig_Rxp_tempv2(CntNumbSeg,1))+(conj(Hlos(Rk2))*Hlos(Rk2)))); %Added N Gain Equalizer %*((conj(loc4g.Nfft))/((conj(loc4g.Nfft)*loc4g.Nfft)))
        end
% revision Needed        
%         for Rk2=1:loc4g.Nfft2
%             Rxp_temp2(Rk2,CntNumbSeg)=Rxp_temp(Rk2,CntNumbSeg)*Hlos_eq(Rk2,CntNumbSeg);
%         end
        
        for Rk2=0:loc4g.Nfft2-1
            Rxp_temp2(Rk2+1,CntNumbSeg)=Rxp_temp(Rk2+1,CntNumbSeg)*exp(1i*2*pi*Rk2*(timedelay/loc4g.Nfft2));
        end
        
        Hlos_indicator_estimated(CntNumbSeg,:)=cazacobj.FSymbol_g0(CntNumbSeg,zcobj.bins_gain(:))*Hlos(1);
        Hlos_indicator_calculated(:,CntNumbSeg)=Rxp_temp2(zcobj.bins_gain(1:zcobj.SegSize),CntNumbSeg);
        
%        Hlos_eq_estimated(CntNumbSeg,:)=Rxp_temp(:,CntNumbSeg)*(conj(Hlos(Rk2))/((Total_pow_vn(CntNumbSeg,1)/powSig_Rxp_temp(CntNumbSeg,1))+(conj(Hlos(Rk2))*Hlos(Rk2))));
        Hlos_eq_mat=[Rxp_temp2((zcobj.bins_gain(1:3)),CntNumbSeg),Hlos_indicator_estimated(CntNumbSeg,1:3).'];
        
        tmpvec_Rxp_temp2(:,CntNumbSeg)=ifft(Rxp_temp2(:,CntNumbSeg),loc4g.Nfft);
        powSig_Rxp_temp2(CntNumbSeg,1)=sum(tmpvec_Rxp_temp2(1:end,CntNumbSeg).*conj(tmpvec_Rxp_temp2(1:end,CntNumbSeg)))/length(tmpvec_Rxp_temp2(1:end,CntNumbSeg));
        Frequency_Power_AVG_Rxp_temp2(CntNumbSeg,:)=sum((Rxp_temp2(:,CntNumbSeg).*conj(Rxp_temp2(:,CntNumbSeg))));
        
        for Rk3=1:loc4g.Nfft2
            Rxp_temp3(Rk3,CntNumbSeg)=Rxp_temp2(Rk3,CntNumbSeg)*W_gain(Rk3,1);%*(1/Hlos(1))*(1/cazacobj.PowerGain2(CntNumbSeg,1));
        end
        
        PowerGain2_indicator_estimated(:,CntNumbSeg)=zcobj.FSymbolG(CntNumbSeg,:).';
        PowerGain2_indicator_calculated(:,CntNumbSeg)=Rxp_temp3(zcobj.bins_gain(1:zcobj.SegSize),CntNumbSeg);
        
        PowerGain2_indicator_mat=[PowerGain2_indicator_calculated(1:3,CntNumbSeg),PowerGain2_indicator_estimated(1:3,CntNumbSeg)];
        
        tmpvec_Rxp_temp3(:,CntNumbSeg)=ifft(Rxp_temp3(:,CntNumbSeg),loc4g.Nfft);
        powSig_Rxp_temp3(CntNumbSeg,1)=sum(tmpvec_Rxp_temp3(1:end,CntNumbSeg).*conj(tmpvec_Rxp_temp3(1:end,CntNumbSeg)))/length(tmpvec_Rxp_temp3(1:end,CntNumbSeg));
        Frequency_Power_AVG_Rxp_temp3(CntNumbSeg,:)=sum((Rxp_temp3(:,CntNumbSeg).*conj(Rxp_temp3(:,CntNumbSeg))));
        
        for Rk4=1:loc4g.Nfft2
            Rxp_temp4(Rk4,CntNumbSeg)=Rxp_temp3(Rk4,CntNumbSeg);%*W_gain(Rk4);%/PowerGain(Rk4)
        end
        
        tmpvec_Rxp_temp4(:,CntNumbSeg)=ifft(Rxp_temp4(:,CntNumbSeg),loc4g.Nfft);
        powSig_Rxp_temp4(CntNumbSeg,1)=sum(tmpvec_Rxp_temp4(1:end,CntNumbSeg).*conj(tmpvec_Rxp_temp4(1:end,CntNumbSeg)))/length(tmpvec_Rxp_temp4(1:end,CntNumbSeg));
        Frequency_Power_AVG_Rxp_temp4(CntNumbSeg,:)=sum((Rxp_temp4(:,CntNumbSeg).*conj(Rxp_temp4(:,CntNumbSeg))));
        
        adj_Rxp(1:zcobj.SegSize,CntNumbSeg)=Rxp_temp4(cazacobj.bins_gain(1:zcobj.SegSize),CntNumbSeg);
        adj_Rx(((CntNumbSeg-1)*zcobj.SegSize)+1:CntNumbSeg*zcobj.SegSize,1)=adj_Rxp(1:zcobj.SegSize,CntNumbSeg);
        
        Fg0p(1:zcobj.SegSize,CntNumbSeg)=zcobj.FSymbol_g0(CntNumbSeg,cazacobj.bins_gain(1:zcobj.SegSize))*Hlos(1).*W_gain(cazacobj.bins_gain(1:zcobj.SegSize),1).';
        Fg0(((CntNumbSeg-1)*zcobj.SegSize)+1:CntNumbSeg*zcobj.SegSize,1)=Fg0p(1:zcobj.SegSize,CntNumbSeg);
        
    end
    
    adj_Rxpn(1:loc4g.NREf,1)=adj_Rx(1:loc4g.NREf,1);
    Temp_mat=[Fg0(1:loc4g.NREf),adj_Rxpn(1:loc4g.NREf)];
    Abs_mat(kind,1)=mean(abs(Temp_mat(:,1)-Temp_mat(:,2)));
    
    saved_adj_Rxp{:,kind}=adj_Rxpn;
    zc2=adj_Rxpn;
    
    %     Alternative way to create a cyclic autocorrelation
    %     taux(1:loc4g.NREf+1,1)=-loc4g.NREf/2:loc4g.NREf/2;
    %     rang_rx(:,1)=rem(taux(1:loc4g.NREf+1)+loc4g.NREf,loc4g.NREf)+1;
    %     rang_rx(:,2)=rem(taux(1:loc4g.NREf+1)+loc4g.NREf,loc4g.NREf);
    %     for sampdel=1:loc4g.NREf+1
    %         final_table(1:loc4g.NREf,sampdel)=zc2([rang_rx(sampdel,1):loc4g.NREf,1:rang_rx(sampdel,2)]);
    %     end
    %     for correlation=1:loc4g.NREf+1
    %         cr2(correlation)=(1/loc4g.NREf)*sum(conj(Ocazac).*final_table(:,correlation));
    %     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Autocorrelation between signals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for tau=-loc4g.NREf/2:loc4g.NREf/2
        k=k+1;
        xax(k)=tau;
        cr(k)=sum(conj(Fg0(1:loc4g.NREf)).*zc2([rem(tau+loc4g.NREf,loc4g.NREf)+1:loc4g.NREf,1:rem(tau+loc4g.NREf,loc4g.NREf)]));
    end %Ocazac, %Fg0(1:loc4g.NREf)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %End of Autocorrelation between signals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Registry Allocation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sumvec(kind)=real(cr(loc4g.NREf/2+1));
    abs_sumvec(kind)=abs(cr(loc4g.NREf/2+1));
    saved_sumvec(1,1:length(sumvec))=sumvec(1,1:length(sumvec));
    saved_abs_sumvec(1,1:length(abs_sumvec))=abs_sumvec(1,1:length(abs_sumvec));
    saved_cr(kind,1:length(cr))=cr(1:length(cr));
    sumcr_0(kind)=sum(abs(cr(1:loc4g.NREf+1)));
    sumcr_real(kind)=sum(real(cr(1:loc4g.NREf+1)));
    max_sumvec(1,1)=max(sumvec(1,1:length(sumvec)));
    
    proportion_denominator_sum(kind)=sum(real(cr(1:(loc4g.NREf/2))))+sum(real(cr((loc4g.NREf/2)+2:loc4g.NREf+1)));
    proportion_denominator_avg(kind)=(sum(real(cr(1:(loc4g.NREf/2))))+sum(real(cr((loc4g.NREf/2)+2:loc4g.NREf+1))))/loc4g.NREf;

    proportion_denominator_sum_abs(kind)=sum(abs(cr(1:(loc4g.NREf/2))))+sum(abs(cr((loc4g.NREf/2)+2:loc4g.NREf+1)));
    proportion_denominator_avg_abs(kind)=(sum(abs(cr(1:(loc4g.NREf/2))))+sum(abs(cr((loc4g.NREf/2)+2:loc4g.NREf+1))))/loc4g.NREf;  
    
    proportion_sum(kind)=sumvec(kind)/proportion_denominator_sum(kind);
    proportion_avg(kind)=sumvec(kind)/proportion_denominator_avg(kind);
    
    proportion_sum_abs(kind)=abs_sumvec(kind)/proportion_denominator_sum_abs(kind);
    proportion_avg_abs(kind)=abs_sumvec(kind)/proportion_denominator_avg_abs(kind);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %End of Registry Allocation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Comparators
    comparator1(kind,5)=sumvec(kind);
    comparator1(kind,6)=abs_sumvec(kind);
    comparator1(kind,7)=sumcr_0(kind);
    comparator1(kind,8)=sumcr_real(kind);
    comparator1(kind,9)=proportion_sum(kind);
    comparator1(kind,10)=proportion_avg(kind);
    comparator1(kind,11)=proportion_sum_abs(kind);
    comparator1(kind,12)=proportion_avg_abs(kind);
    if kind>1
        comparator1(kind,1)=sumvec(kind)-sumvec(kind-1);
        comparator1(kind,2)=abs_sumvec(kind)-abs_sumvec(kind-1);
        comparator1(kind,3)=sumcr_0(kind)-sumcr_0(kind-1);
        comparator1(kind,4)=sumcr_real(kind)-sumcr_real(kind-1);
    end
    
end

comparatora(1,1)=find(comparator1(:,5)==max(comparator1(:,5)));
comparatora(1,2)=find(comparator1(:,6)==max(comparator1(:,6)));
comparatora(1,3)=find(comparator1(:,7)==min(comparator1(:,7)));
comparatora(1,4)=find(comparator1(:,8)==min(comparator1(:,8)));
comparatora(1,5)=find(comparator1(:,9)==max(comparator1(:,9)));
comparatora(1,6)=find(comparator1(:,10)==max(comparator1(:,10)));
comparatora(1,7)=find(comparator1(:,11)==max(comparator1(:,11)));
comparatora(1,8)=find(comparator1(:,12)==max(comparator1(:,12)));
comparatora(1,9)=floor(defdel*10)+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End - Original CAZAC Sequence Retrieval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of Window Frame Samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Registry Allocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saved_max_sumvec(1,1:length(max_sumvec))=max_sumvec(1,1:length(max_sumvec));
max_sumvec;
index_value=max(max_sumvec);
index=find(max_sumvec==index_value,1,'last');
index=max(index);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End - Registry Allocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cyclic Autocorrelation 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(simode,'AWGN')
    index_max_phase=find(sumvec(1,:)==max(sumvec(1,:)),1,'last');
    firtphase=(index_max_phase-2)*0.1;
    %Figure 9
    target_cr=sumvec;
    figure(9)
    plot(1:kind,sumvec(1:kind));
    xlabel('Correlation Bin 0 Value per Phase Change');
    ylabel('Amplitude - Bin 0');
    %titleSTR=strcat({'1st Est. Phase - Phase='},{num2str(firtphase)},{', Smallest |CR| Samp='},{num2str(index_max_phase)});
    %title(titleSTR);
    title({['1st Est. Phase - Phase=',num2str(firtphase),',   Real CR Samp=',num2str(index_max_phase)];['Def. Phase=',num2str(defdel),'   Sim=',num2str(mm),'   Tower=',num2str(TowInd)]})
    saveas(gcf, fldhomecombfig11);
    %print(gcf, '-djpeg', fldhomecombfig11,'-r500');
else
    index_max_phase=find(sumcr_0(1,:)==min(sumcr_0(1,:)),1,'last');
    firtphase=(index_max_phase-2)*0.1;
    
     index_max_phase_sumvec=find(sumvec(1,:)==max(sumvec(1,:)),1,'last');
     firtphase_sumvec=(index_max_phase_sumvec-2)*0.1;
     
     %index_max_phase=index_max_phase_sumvec; %Mod W45 V4
     %firtphase=firtphase_sumvec; %Mod W45 V4
    %%%%%%%%%%
    %Hypothesis
    %%%%%%%%%%
    
    %END - Hypothesis
    
    %Figure 9
    target_cr=sumcr_0;
    figure(9)
    subplot(2,1,1)
    plot(1:kind,sumvec(1:kind));
    xlabel('Correlation Bin 0 Value per Phase Change');
    ylabel('Amplitude - Bin 0');
    title({['1st Est. Phase - Phase=',num2str(firtphase_sumvec),',   Real CR Samp=',num2str(index_max_phase_sumvec)];['Def. Phase=',num2str(defdel),'   Sim=',num2str(mm),'   Tower=',num2str(TowInd)]})
    subplot(2,1,2)
    plot(1:kind,sumcr_0(1:kind));
    xlabel('Sum(|CR|) - Phase Range');
    ylabel('Amplitude - Sum (|CR|)');
    %titleSTR=strcat({'1st Est. Phase - Phase='},{num2str(firtphase)},{', Smallest |CR| Samp='},{num2str(index_max_phase)});
    %title(titleSTR);
    title({['1st Est. Phase - Phase=',num2str(firtphase),',   Smallest Sum(|CR|) Samp=',num2str(index_max_phase)];['Def. Phase=',num2str(defdel),'   Sim=',num2str(mm),'   Tower=',num2str(TowInd)]})
    saveas(gcf, fldhomecombfig11);
    %print(gcf, '-djpeg', fldhomecombfig11,'-r500');
    
    index_max_phase=index_max_phase_sumvec; %Mod W45 V4
    firtphase=firtphase_sumvec; %Mod W45 V4
%     figure(20)
%     subplot(4,1,1)
%     plot(1:kind,sumvec(1:kind));
%     title({['Real Phase=',num2str(comparatora(1,5))]})
%     ylabel('sumvec')
%     subplot(4,1,2)
%     plot(1:kind,abs_sumvec(1:kind));
%     ylabel('abs sumvec')
%     subplot(4,1,3)
%     plot(1:kind,sumcr_0(1:kind));
%     ylabel('sumcr 0')
%     subplot(4,1,4)
%     plot(1:kind,sumcr_real(1:kind));
%     ylabel('sumcr real')
    
    %figure(21)
    %plot(-loc4g.NREf/2:loc4g.NREf/2,saved_cr(floor(defdel*10)+1,:));
    %ylabel('Correlation Amplitude');
    %xlabel('Tau Index');
    %title({['Correlation Plot of Estimated Phase=',num2str((index_max_phase_sumvec-1)/10),'   Real Phase=',num2str(defdel)];['SNIRdB=',num2str(cazacobj.SNRdB),'   Mode=AWGN+F+I'];['Zadoff-Chu=',num2str(cazacobj.NREf),'   Number of Symbols=',num2str(cazacobj.NumbSeg)]})
    %title({['Correlation Plot of Def. Phase=',num2str(defdel)];['SNRdB=',num2str(cazacobj.SNRdB),'   Mode=',num2str(simode)];['Zadoff-Chu=',num2str(cazacobj.NREf),'   Number of Symbols=',num2str(cazacobj.NumbSeg)]})
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Registry Allocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saved_cr2=zeros(cp_inc_samples,loc4g.NREf+1);
kind2=0;
timedelay2=0;
kvec=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End - Registry Allocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Hypothesis
[m mi] = sort(sumcr_0);
lowest7index = mi(1:7);
sumvec(lowest7index);
[row,col] = find(sumvec>max(sumvec)-25);;

[m2 mi2]=sort(sumvec);
highest7index = mi2(end-15:end);
sumvec(highest7index);

Rxp(1:loc4g.Nfft2,1:zcobj.NumbSeg)=saved_Rxp(1:loc4g.Nfft2,1:zcobj.NumbSeg);
for timedelay2=(index_max_phase-2-3)*0.1:0.0001:(index_max_phase+1)*0.1
    kind2=kind2+1;
    loopN=0;
    adj_Rx=zeros(loc4g.NREf,1);
    adj_Rxp=zeros(zcobj.SegSize,zcobj.NumbSeg);
    adj_Rxpn=zeros(loc4g.NREf,1);
    zc2=zeros(loc4g.NREf,1);
    tau=0;
    xax=0;
    cr=zeros(1,loc4g.NREf+1);
    k=0;
    
    for CntNumbSeg=1:zcobj.NumbSeg
        
        Hscat=fft(loc4g.h00Nscat(1:end),loc4g.Nfft2);
        h00NDel=RiceanEst.h00NDel(TowInd+1,1:end);
        HDel=fft(h00NDel,length(h00NDel));
        
        for loopN=0:length(h00NDel)-1
            Hest(loopN+1)=HDel(loopN+1)*exp(1i*2*pi*loopN*(timedelay2/length(h00NDel)));
        end
        
        h00Nest=ifft(Hest,length(h00NDel));
        h00Nscatest=[0,h00Nest(1,2:end)];
        Hscatest=fft(h00Nscatest,loc4g.Nfft2);
        Hscat=Hscatest;
        
        for loopN=0:loc4g.Nfft2-1
            Hscat_temp(loopN+1)=Hscat(loopN+1)*exp(-1i*2*pi*loopN*(timedelay2/loc4g.Nfft2));
        end
        
        Hscat=Hscat_temp;
        
        for Rk=1:loc4g.Nfft2
            Rxp_temp(Rk,CntNumbSeg)=Rxp(Rk,CntNumbSeg)-(cazacobj.FSymbol_g0(CntNumbSeg,Rk)*Hscat(Rk));
%             test_fad_ref(Rk,CntNumbSeg)=(cazacobj.FSymbol_g0(CntNumbSeg,Rk)*Hscat(Rk));
%             ratio_symb(Rk,CntNumbSeg)=((sqrt(zcobj.NREf^2)/zcobj.SegSize)/(PowerGain(1,Rk)*zcobj.PowerGain2(CntNumbSeg,1)));
%             Rxp_temp_ref(Rk,CntNumbSeg)=Rxp_temp(Rk,CntNumbSeg)*ratio_symb(Rk,CntNumbSeg);
%             Rxp_temp(Rk,CntNumbSeg)=Rxp_temp_ref(Rk,CntNumbSeg);
        end
        
        tmpvec_Rxp_temp(:,CntNumbSeg)=ifft(Rxp_temp(:,CntNumbSeg),loc4g.Nfft);
        powSig_Rxp_temp(CntNumbSeg,:)=sum(tmpvec_Rxp_temp(1:end,CntNumbSeg).*conj(tmpvec_Rxp_temp(1:end,CntNumbSeg)))/length(tmpvec_Rxp_temp(1:end,CntNumbSeg));
        
        Hlos=fft(loc4g.h00Nlos(1:end),loc4g.Nfft2);
%        Hlos=ones(1,loc4g.Nfft2);
        h00Nlosest=[h00Nest(1,1),zeros(1,length(h00Nest)-1)];
        Hlosest=fft(h00Nlosest,loc4g.Nfft2);
        Hlos=Hlosest;
        
        Hscat_indicator_estimated(:,CntNumbSeg)=cazacobj.FSymbolG2(CntNumbSeg,:)*Hlos(1);
        Hscat_indicator_calculated(:,CntNumbSeg)=Rxp_temp(zcobj.bins_gain(1:zcobj.SegSize),CntNumbSeg);
        Hscat_mat=[Hscat_indicator_estimated(:,CntNumbSeg),Hscat_indicator_calculated(:,CntNumbSeg)];
        for loopN=0:loc4g.Nfft2-1
            Hlos_temp(loopN+1)=Hlos(loopN+1)*exp(-1i*2*pi*loopN*(timedelay2/loc4g.Nfft2));
        end
        
%        Hlos=Hlos_temp;
        
        for Rk2=1:loc4g.Nfft2
            Hlos_eq(Rk2,CntNumbSeg)=(conj(Hlos(Rk2))/((vn2(CntNumbSeg,1)/powSig_Rxp_tempv2(CntNumbSeg,1))+(conj(Hlos(Rk2))*Hlos(Rk2)))); %Added N Gain Equalizer %*((conj(loc4g.Nfft))/((conj(loc4g.Nfft)*loc4g.Nfft)))
        end
        
%         for Rk2=1:loc4g.Nfft2
%             Rxp_temp2(Rk2,CntNumbSeg)=Rxp_temp(Rk2,CntNumbSeg)*Hlos_eq(Rk2,CntNumbSeg);
%         end

        for Rk2=0:loc4g.Nfft2-1
            Rxp_temp2(Rk2+1,CntNumbSeg)=Rxp_temp(Rk2+1,CntNumbSeg)*exp(1i*2*pi*Rk2*(timedelay2/loc4g.Nfft2));
        end
        
        Hlos_indicator_estimated(CntNumbSeg,:)=cazacobj.FSymbol_g0(CntNumbSeg,zcobj.bins_gain(:))*Hlos(1);
        Hlos_indicator_calculated(:,CntNumbSeg)=Rxp_temp2(zcobj.bins_gain(1:zcobj.SegSize),CntNumbSeg);
        
        Hlos_eq_estimated(CntNumbSeg,:)=Rxp_temp(:,CntNumbSeg)*(conj(Hlos(Rk2))/((Total_pow_vn(CntNumbSeg,1)/powSig_Rxp_temp(CntNumbSeg,1))+(conj(Hlos(Rk2))*Hlos(Rk2))));
        Hlos_eq_mat=[Rxp_temp2((zcobj.bins_gain(1:3)),CntNumbSeg),Hlos_indicator_estimated(CntNumbSeg,1:3).',Hlos_eq_estimated(CntNumbSeg,zcobj.bins_gain(1:3)).'];
        
        tmpvec_Rxp_temp2(:,CntNumbSeg)=ifft(Rxp_temp2(:,CntNumbSeg),loc4g.Nfft);
        powSig_Rxp_temp2(CntNumbSeg,1)=sum(tmpvec_Rxp_temp2(1:end,CntNumbSeg).*conj(tmpvec_Rxp_temp2(1:end,CntNumbSeg)))/length(tmpvec_Rxp_temp2(1:end,CntNumbSeg));
        
        for Rk3=1:loc4g.Nfft2
            Rxp_temp3(Rk3,CntNumbSeg)=Rxp_temp2(Rk3,CntNumbSeg)*W_gain(Rk3,1);%*(1/Hlos(1))*(1/cazacobj.PowerGain2(CntNumbSeg,1));
        end
        
%        PowerGain2_indicator_estimated(:,CntNumbSeg)=zcobj.FSymbolG(CntNumbSeg,:).';
        PowerGain2_indicator_estimated(:,CntNumbSeg)=(cazacobj.FSymbol_g0(CntNumbSeg,zcobj.bins_gain(1:zcobj.SegSize)).')*Hlos(1).*W_gain(zcobj.bins_gain((1:zcobj.SegSize)),1);
        PowerGain2_indicator_calculated(:,CntNumbSeg)=Rxp_temp3(zcobj.bins_gain(1:zcobj.SegSize),CntNumbSeg);
        
        PowerGain2_indicator_mat=[PowerGain2_indicator_calculated(:,CntNumbSeg),PowerGain2_indicator_estimated(:,CntNumbSeg)];
        PowerGain2_comp(kind2,CntNumbSeg)=mean(abs(PowerGain2_indicator_mat(:,1)-PowerGain2_indicator_mat(:,2)));
        tmpvec_Rxp_temp3(:,CntNumbSeg)=ifft(Rxp_temp3(:,CntNumbSeg),loc4g.Nfft);
        powSig_Rxp_temp3(CntNumbSeg,1)=sum(tmpvec_Rxp_temp3(1:end,CntNumbSeg).*conj(tmpvec_Rxp_temp3(1:end,CntNumbSeg)))/length(tmpvec_Rxp_temp3(1:end,CntNumbSeg));
        
        for Rk4=1:loc4g.Nfft2
            Rxp_temp4(Rk4,CntNumbSeg)=Rxp_temp3(Rk4,CntNumbSeg);%*W_gain(Rk4);%/PowerGain;
        end
        
        tmpvec_Rxp_temp4(:,CntNumbSeg)=ifft(Rxp_temp4(:,CntNumbSeg),loc4g.Nfft);
        powSig_Rxp_temp4(CntNumbSeg,1)=sum(tmpvec_Rxp_temp4(1:end,CntNumbSeg).*conj(tmpvec_Rxp_temp4(1:end,CntNumbSeg)))/length(tmpvec_Rxp_temp4(1:end,CntNumbSeg));
        
        adj_Rxp(1:zcobj.SegSize,CntNumbSeg)=Rxp_temp4(cazacobj.bins_gain(1:zcobj.SegSize),CntNumbSeg);
        adj_Rx(((CntNumbSeg-1)*zcobj.SegSize)+1:CntNumbSeg*zcobj.SegSize,1)=adj_Rxp(1:zcobj.SegSize,CntNumbSeg);
    
        Fg0p(1:zcobj.SegSize,CntNumbSeg)=zcobj.FSymbol_g0(CntNumbSeg,cazacobj.bins_gain(1:zcobj.SegSize))*Hlos(1).*W_gain(cazacobj.bins_gain(1:zcobj.SegSize),1).';
        Fg0(((CntNumbSeg-1)*zcobj.SegSize)+1:CntNumbSeg*zcobj.SegSize,1)=Fg0p(1:zcobj.SegSize,CntNumbSeg);
    
    end
    
    adj_Rxpn(1:loc4g.NREf,1)=adj_Rx(1:loc4g.NREf,1);
    Temp_mat=[Fg0(1:loc4g.NREf),adj_Rxpn(1:loc4g.NREf)];
    Abs_mat2(kind2,1)=mean(abs(Temp_mat(:,1)-Temp_mat(:,2)));
    
    saved_adj_Rxpk{:,kind2}=adj_Rxpn;
    zc2=adj_Rxpn;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Autocorrelation between signals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for tau=-loc4g.NREf/2:loc4g.NREf/2
        k=k+1;
        xax(k)=tau;
        cr(k)=sum(conj(Fg0(1:loc4g.NREf)).*zc2([rem(tau+loc4g.NREf,loc4g.NREf)+1:loc4g.NREf,1:rem(tau+loc4g.NREf,loc4g.NREf)]));
    end %Ocazac,%Fg0(1:loc4g.NREf)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %End of Autocorrelation between signals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Registry Allocation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sumvec2(kind2)=real(cr(loc4g.NREf/2+1));
    abs_sumvec2(kind2)=abs(cr(loc4g.NREf/2+1));
    saved_cr2(kind2,1:length(cr))=cr(1:length(cr));
    sumcr(kind2)=sum(abs(cr(1:loc4g.NREf-1)));
    
    comparator2(kind2,4)=sumvec2(kind2);
    comparator2(kind2,5)=abs_sumvec2(kind2);
    comparator2(kind2,6)=sumcr(kind2);
    
    %Comparators
    if kind2>1
        comparator2(kind2,1)=sumvec2(kind2-1)-sumvec2(kind2);
        comparator2(kind2,2)=abs_sumvec2(kind2-1)-abs_sumvec2(kind2);
        comparator2(kind2,3)=sumcr(kind2-1)-sumcr(kind2);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if strcmp(simode,'AWGN')
    %[diffval index_max_phase3] = min(sumcr-1);%MMSE
    %secondphase=(index_max_phase3-1)*0.0001;
    index_max_phase3=find(sumvec2(1,:)==max(sumvec2(1,:)),1,'last'); %Testing
    secondphase=(index_max_phase3-1)*0.0001; %Testing
    
    calculated_delay=(index_max_phase-2)*0.1+(index_max_phase3-1)*0.0001-(3)*0.1;
    calculated_time=calculated_delay*loc4g.ts2;
    
    %Figure 10
    target_cr=sumvec2;
    figure(10)
    plot(1:kind2,sumvec2(1:kind2));
    xlabel('Correlation Bin 0 Value per Phase Change');
    ylabel('Amplitude - Bin 0');
    %titleSTR=strcat({'2nd Est. Phase - Phase='},{num2str(secondphase)},{' Total Est. Phase='},{num2str(calculated_delay)},{' MMSE |CR| to 1='},{num2str(index_max_phase3)});
    %title(titleSTR);
    title({['2nd Est. Phase - Phase=',num2str(secondphase),', Maximum value found at Phase=',num2str(index_max_phase3),', Sim=',num2str(mm),', Tower=',num2str(TowInd)];['Total Est. Phase=',num2str(calculated_delay),'  -  Def. Phase=',num2str(defdel)]})
    saveas(gcf, fldhomecombfig12);
    %print(gcf, '-djpeg', fldhomecombfig12,'-r500');
else
    %index_max_phase3=find(sumcr(1,:)==min(sumcr(1,:)),1,'last');%Look for value closer to 1
    index_max_phase3=find(sumvec2(1,:)==max(sumvec2(1,:)),1,'last'); % Replacement
    secondphase=(index_max_phase3-1)*0.0001;
    
    calculated_delay=(index_max_phase-2)*0.1+(index_max_phase3-1)*0.0001-(3)*0.1;
    calculated_time=calculated_delay*loc4g.ts2;
    
    %Figure 10
    target_cr=sumcr;
    figure(10)
    plot(1:kind2,sumvec2(1:kind2)); %Replaced sumcr with sumvec2
    xlabel('Sum(|CR|) Sample Number');
    ylabel('Sum of |CR|');
    title({['2nd Est. Phase - Phase=',num2str(secondphase),', Minimum value found at Phase=',num2str(index_max_phase3),', Sim=',num2str(mm),', Tower=',num2str(TowInd)];['Total Est. Phase=',num2str(calculated_delay),'  -  Def. Phase=',num2str(defdel)]})
   saveas(gcf, fldhomecombfig12);
    %print(gcf, '-djpeg', fldhomecombfig12,'-r500');
end
varargout{1}=calculated_time;
varargout{2}=calculated_delay;
end
