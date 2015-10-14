clear all
close all
clear classes
dbstop if error
areaval        = (40*1609.344)^2;  %miles diameter to meters square: rectangular

%%%%%%%%%%%%%%% Tower Locations %%%%%%%%%%%%%
InitVecSize    = 64;  %used in hyperbolic solver initial guess
DelaySpread    = 20e-6;
InitMethod     = 'locus';  %choices:  'random' or 'locus'
NumIter        = 8; %Number of times we iterate in hyperbolic solution solverNumSymbols;
frametype      = 1;     %location<--,  information, synchronization
numtargets     = 1;   %currently only one target supported
noiseval       = 0;   % rather than use timing estimate to create range error, we can use an explic noise variance
Numsim         = 35;  % reruns with new noise statistics  and a new fading channel, should choose new target location
Fsval          = 7.68e6;   %sampling frequency of multicarrier modulation  (subcarrier space: 15e3)
Fsval2         = 4.995e6;   %Resource elements used for data
CPsize         = 2*161+1;  %size of cyclic prefix
simode         = 'AWGN_Fadding';% Modes: AWGN_Fadding_Interference, AWGN_Fadding, AWGN
MonteCarloSims = 16;
RangeV=5e-6;
SimResP1=  10000;  %defines the division or resolution for the phase estimate associated with the (range delay)
fldhome        = 'C:\Users\ejl334\Desktop\Sample2B\Results';
mkdir(fldhome)

if strcmp(simode,'AWGN_Fadding_Interference')
    nvar_ratio=(0.5);
    int_ratio=(0.5);
    smode='AWGN+F+I';
    titleSNR='SNIR=';
elseif strcmp(simode,'AWGN_Fadding')
    nvar_ratio=(1);
    int_ratio=(1e-50);
    smode='AWGN+F';
    titleSNR='SNR=';
elseif strcmp(simode,'AWGN')
    nvar_ratio=(1);
    int_ratio=(1e-50);
    smode='AWGN';
    titleSNR='SNR=';
end

SynchWindow=(CPsize-round(DelaySpread/(1/Fsval)))*(1/Fsval);
errvec=zeros(Numsim,1);
format shortg;
avgNumbSeg=zeros(Numsim,MonteCarloSims); %Mod W45 V4
for SimScen=1%1:MonteCarloSims
    switch SimScen
        
        case 1
            SNRvaldB=15; RicianKdB=8; NumSymbols=1; Numref=7; NumberSubCarr  =128; SegSize=[]; areaval=(40*1609.344)^2;
        case 2
            SNRvaldB=0; RicianKdB=8; NumSymbols=1; Numref=7; NumberSubCarr  =128; SegSize=[]; areaval=(40*1609.344)^2;
        case 3
            SNRvaldB=-5; RicianKdB=8; NumSymbols=1; Numref=7; NumberSubCarr  =128; SegSize=[]; areaval=(40*1609.344)^2;
        case 4
            SNRvaldB=-15; RicianKdB=8; NumSymbols=1; Numref=7; NumberSubCarr  =128; SegSize=[]; areaval=(40*1609.344)^2;
        case 5
            SNRvaldB=15; RicianKdB=8; NumSymbols=1; Numref=7; NumberSubCarr  =512; SegSize=[]; areaval=(40*1609.344)^2;
        case 6
            SNRvaldB=0; RicianKdB=8; NumSymbols=1; Numref=7; NumberSubCarr  =512; SegSize=[]; areaval=(40*1609.344)^2;
        case 7
            SNRvaldB=-5; RicianKdB=8; NumSymbols=1; Numref=7; NumberSubCarr  =512; SegSize=[]; areaval=(40*1609.344)^2;
        case 8
            SNRvaldB=-15; RicianKdB=8; NumSymbols=1; Numref=7; NumberSubCarr  =512; SegSize=[]; areaval=(40*1609.344)^2;
        case 9
            SNRvaldB=15; RicianKdB=8; NumSymbols=1; Numref=7; NumberSubCarr  =2048; SegSize=[]; areaval=(40*1609.344)^2;
        case 10
            SNRvaldB=0; RicianKdB=8; NumSymbols=1; Numref=7; NumberSubCarr  =2048; SegSize=[]; areaval=(40*1609.344)^2;
        case 11
            SNRvaldB=-5; RicianKdB=8; NumSymbols=1; Numref=7; NumberSubCarr  =2048; SegSize=[]; areaval=(40*1609.344)^2;
        case 12
            SNRvaldB=-15; RicianKdB=8; NumSymbols=1; Numref=7; NumberSubCarr  =2048; SegSize=[]; areaval=(40*1609.344)^2; 
        case 13
            SNRvaldB=15; RicianKdB=8; NumSymbols=1; Numref=7; NumberSubCarr  =8192; SegSize=[]; areaval=(40*1609.344)^2;
        case 14
            SNRvaldB=0; RicianKdB=8; NumSymbols=1; Numref=7; NumberSubCarr  =8192; SegSize=[]; areaval=(40*1609.344)^2;
        case 15
            SNRvaldB=-5; RicianKdB=8; NumSymbols=1; Numref=7; NumberSubCarr  =8192; SegSize=[]; areaval=(40*1609.344)^2;
        case 16
            SNRvaldB=-15; RicianKdB=8; NumSymbols=1; Numref=7; NumberSubCarr  =8192; SegSize=[]; areaval=(40*1609.344)^2;    
    end
    
    warning('off','MATLAB:uitreenode:DeprecatedFunction');
    warning('off','MATLAB:uitab:DeprecatedFunction');
    warning('off','MATLAB:uitabgroup:DeprecatedFunction');
    warning('off','MATLAB:RandStream:ReadingInactiveLegacyGeneratorState');
    warning('off','MATLAB:uitree:DeprecatedFunction');
    warning('off','MATLAB:JavaComponentThreading');
    warning('off','MATLAB:uiflowcontainer:DeprecatedFunction');
    warning('off','MATLAB:JavaEDTAutoDelegation');
    warning('off','MATLAB:RandStream:ActivatingLegacyGenerators');
    warning('off','MATLAB:class:DynPropDuplicatesMethod');
    warning('off','MATLAB:Debugger:BreakpointSuppressed');
    warning('off','MATLAB:uitable:DeprecatedFunction');
    warning('off','MATLAB:uigridcontainer:DeprecatedFunction');
    warning('off','MATLAB:mir_warning_unrecognized_pragma');
    warning('off','MATLAB:class:InvalidDynamicPropertyName');
    warning('off','MATLAB:MKDIR:DirectoryExists');
    warning('off','comm:channel_channelfilter_initialize:channeltaps');
    warning('off','MATLAB:plot:IgnoreImaginaryXYPart');
    
    %NumbSeg=ceil(NumberSubCarr/SegSize); %Mod W45 V3
    NumbSeg=[]; %Mod W45 V3
    BW  = 15e3*NumberSubCarr;
    DelEstLongFreq2=zeros(NumSymbols,Numref);
    DelSamples_Saved=zeros(Numsim,Numref);
    Target_Locations_MC=zeros(SimScen,Numsim);
    
    ExpandedEq=1;
    mm=1;
    fldhomecombfig1=char(strcat(fldhome, {'\'},{'fig1.jpg'}));
    
    adhocobj=AdHocNodes(Numref, areaval, noiseval, numtargets, NumIter,...
        ExpandedEq, InitVecSize, InitMethod, SynchWindow,fldhomecombfig1,...
        SNRvaldB,RicianKdB,NumberSubCarr,SegSize,NumbSeg,...
        nvar_ratio,int_ratio,smode,titleSNR);  %constructor;
    
    Target_Locations_MC(SimScen,mm)=adhocobj.target;
    DelayToTarget=get(adhocobj, 'TDelay');
    TotalNumbSeg=zeros(Numsim,Numref); %Mod W45 V4
    for mm = 1:Numsim
        %Coordinate Locations of Numref Towers and one random Target are generated
        ExpandedEq=1;  % Redundant hyperbolic equation parameter facter: multiplier by number of reference toweres
        
        fldstr3=char(strcat({'SimResult'},{num2str(SimScen)}));
        fldstr4=char(strcat({'NumSim'},{num2str(mm)}));
        fldhome2=char(strcat(fldhome, {'\'},fldstr3));
        fldhomecomb1=char(strcat(fldhome, {'\'},fldstr3));
        fldhomecombfig1=char(strcat(fldhome2, {'\'},{'fig1.jpg'}));
        mkdir(fldhome,fldstr3)
        mkdir(fldhome2,fldstr4)
        
        if mm==1
            saveas(gcf,fldhomecombfig1);
        end
        
        for TowInd=0:Numref-1
            same_pdp='no';
            pdp_value=0;
            fldstr5=char(strcat({'TowInd'},{num2str(TowInd)}));
            fldhome3=char(strcat(fldhome2, {'\'},fldstr4));
            fldhomecombfig0=char(strcat(fldhome3, {'\'},fldstr5, {'\'}, {'fig1.jpg'}));
            fldhomecombfig2=char(strcat(fldhome3, {'\'},fldstr5, {'\'}, {'fig2.jpg'})); %AV
            fldhomecombfig3=char(strcat(fldhome3, {'\'},fldstr5, {'\'}, {'fig3.jpg'}));
            fldhomecombfig4=char(strcat(fldhome3, {'\'},fldstr5, {'\'}, {'fig4.jpg'}));
            fldhomecombfig5=char(strcat(fldhome3, {'\'},fldstr5, {'\'}, {'fig5.jpg'})); %AV
            fldhomecombfig8=char(strcat(fldhome3, {'\'},fldstr5, {'\'}, {'fig6.jpg'}));
            fldhomecombfig9=char(strcat(fldhome3, {'\'},fldstr5, {'\'}, {'fig7.jpg'}));
            fldhomecombfig10=char(strcat(fldhome3, {'\'},fldstr5, {'\'}, {'fig8.jpg'})); %AV
            fldhomecombfig11=char(strcat(fldhome3, {'\'},fldstr5, {'\'}, {'fig9.jpg'})); %AV
            fldhomecombfig12=char(strcat(fldhome3, {'\'},fldstr5, {'\'}, {'fig10.jpg'})); %AV
            mkdir(fldhome3,fldstr5)
            TowerIndex=TowInd;
            loc4g =RadioLocate(SNRvaldB, RicianKdB, Numref, NumSymbols, frametype,...
                DelayToTarget, TowerIndex+1, CPsize, Fsval, BW, NumberSubCarr,fldhomecombfig1,SegSize,NumbSeg,mm,simode);
            DelayToTarget_Samples=DelayToTarget/loc4g.ts;
            length_riceanh=length(loc4g.h00N);
            [Intbaseband, Interchantime,pdp_value,Intbaseband_freq, Intbaseband2,Intbaseband_freq2]=InterfereGenScriptChanOnly(same_pdp,pdp_value,fldhomecombfig3,fldhomecombfig4,length_riceanh,loc4g,int_ratio,mm,TowInd);
            same_pdp='yes';
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WFG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Determine Useful Bins
            min_re=floor(((-loc4g.fs/2)-(-Fsval2/2))/loc4g.DF);%Added
            max_re=floor((((loc4g.fs/2)-(Fsval2/2))/loc4g.DF))+1;%Added
            %New available bins range allos data only in the mid 5 MHz
            bins_av=[max_re:loc4g.Nfft2/2,(loc4g.Nfft2/2)+2:loc4g.Nfft2+min_re];%Added
            [PowerGain,bind1]=wfg(adhocobj,loc4g,bins_av,Intbaseband,Interchantime,Intbaseband_freq,...
                                  Intbaseband2,Intbaseband_freq2,nvar_ratio,TowInd,mm,simode,smode,fldhomecombfig2);
            CSymbols_PT=length(bind1);                   
            Freq_IndexN=[((-loc4g.fs/2)/loc4g.DF):((loc4g.fs/2)/loc4g.DF)-1];
            Freq_IndexHz=(Freq_IndexN*loc4g.DF)/1e6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End - WFG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            zcobj=loc4g.cazacobj;
            
            for indexval=0:get(zcobj, 'Nsymb')-1
                if rem(indexval,5)==0
                    disp(strcat({'Symbol='},{num2str(indexval)}, {', Tower='},{num2str(TowInd)}, {', Simulation='},{num2str(mm)}));
                end
                h00N=get(loc4g, 'DiscreteChan');
                zcobj=RunCazacSys(zcobj, indexval+1,TowInd+1,DelayToTarget,CSymbols_PT,PowerGain,bind1,Intbaseband,h00N);
                RiceanEst=ChannEst(h00N,DelayToTarget_Samples,loc4g.Nfft,loc4g.SNR,zcobj.powSig,nvar_ratio,simode);%W6 V2A simode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimation Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                timedatain=get(zcobj, 'veclong');
                timedatainv=[timedatain].';
                TxSig00=[];%zeros(zcobj.NumbSeg,length(h00N)+length(timedatainv)-1);
                h00N=get(loc4g, 'DiscreteChan');
                TxSig01=[];%zeros(zcobj.NumbSeg,zcobj.Nfft+zcobj.CP+length(h00N)-1);
                timedataout=[];%zeros(zcobj.NumbSeg,length(TxSig01));
                
                for CntNumbSeg=1:zcobj.NumbSeg
                    powSig=zcobj.powSig;
                    nvar_i(CntNumbSeg,1)=powSig(CntNumbSeg,1)/10^(get(loc4g,'SNR')/10);
                    nvar(CntNumbSeg,1)=nvar_i(CntNumbSeg,1)*nvar_ratio;
                    TxSig00(CntNumbSeg,:)=conv(h00N, timedatainv(CntNumbSeg,1:end));
                    len_Tx_Default=length(TxSig00(CntNumbSeg,:));
                    
                    cmp_sigs=zeros(loc4g.Nfft,3);
                    fft_TxSig00=fft(TxSig00(CntNumbSeg,loc4g.CP+1:loc4g.CP+loc4g.Nfft));
                    cmp_sigs(zcobj.bins_gain(1:zcobj.SegSize),1)=fft_TxSig00(zcobj.bins_gain(1:zcobj.SegSize));
                    fft_interference=fft(Intbaseband(loc4g.CP+1:loc4g.CP+loc4g.Nfft));
                    cmp_sigs(zcobj.bins_gain(1:zcobj.SegSize),2)=fft_interference(zcobj.bins_gain(1:zcobj.SegSize));
                    cmp_sigs(:,3)=Intbaseband_freq;
                    cmp_sigs(zcobj.bins_gain(1:zcobj.SegSize),4)=1;
                    fft_ricean=fft(h00N,loc4g.Nfft);
                    
                    if strcmp(simode,'AWGN_Fadding_Interference')
                        TxSig01(CntNumbSeg,:)=TxSig00(CntNumbSeg,:)+Intbaseband;
                    else
                        TxSig01(CntNumbSeg,:)=TxSig00(CntNumbSeg,:);
                    end
                    
                    noise_var(CntNumbSeg,:)=(sqrt(nvar(CntNumbSeg,1)/2)*randn(size(TxSig01(CntNumbSeg,1:end)))+1i*sqrt(nvar(CntNumbSeg,1)/2)*randn(size(TxSig01(CntNumbSeg,1:end))));
                    %noise_var(CntNumbSeg,1:length(TxSig01))=sqrt(nvar(CntNumbSeg,1)/2)+1i*sqrt(nvar(CntNumbSeg,1)/2);
                    noise_var_fft(CntNumbSeg,:)=fft(noise_var(loc4g.CP2+1:loc4g.CP2+loc4g.Nfft2),loc4g.Nfft2);
                    noise_var_fft_Power_AVG(CntNumbSeg,1)=sum(noise_var_fft(CntNumbSeg,:).*conj(noise_var_fft(CntNumbSeg,:)));
                    
                    fft_allnoise=noise_var_fft(1,:)+Intbaseband_freq(:,1).';
                    N_o_S=fft_allnoise./fft_ricean;
                    fft_timedatainv=fft(timedatainv(1,loc4g.CP2+1:loc4g.CP2+loc4g.Nfft2),loc4g.Nfft2);
                    fft_timedatainv2=zeros(1,loc4g.Nfft2);
                    fft_timedatainv2(1,zcobj.bins_gain(1:zcobj.SegSize))=fft_timedatainv(1,zcobj.bins_gain(1:zcobj.SegSize));
                    
                    pow_TxSig00(CntNumbSeg,:)=sum(TxSig00(CntNumbSeg,loc4g.CP2+1:loc4g.CP2+loc4g.Nfft2).*conj(TxSig00(CntNumbSeg,loc4g.CP2+1:loc4g.CP2+loc4g.Nfft2)))/length(TxSig00(CntNumbSeg,loc4g.CP2+1:loc4g.CP2+loc4g.Nfft2)); %I
                    pow_TxSig01(CntNumbSeg,:)=sum(TxSig01(CntNumbSeg,loc4g.CP2+1:loc4g.CP2+loc4g.Nfft2).*conj(TxSig01(CntNumbSeg,loc4g.CP2+1:loc4g.CP2+loc4g.Nfft2)))/length(TxSig01(CntNumbSeg,loc4g.CP2+1:loc4g.CP2+loc4g.Nfft2));
                    pow_noise_var(CntNumbSeg,:)=sum(noise_var(CntNumbSeg,1:end).*conj(noise_var(CntNumbSeg,1:end)))/length(noise_var(CntNumbSeg,1:end));
                    pow_noise_var_fft(CntNumbSeg,:)=sum(noise_var(CntNumbSeg,loc4g.CP2+1:loc4g.CP2+loc4g.Nfft2).*conj(noise_var(CntNumbSeg,loc4g.CP2+1:loc4g.CP2+loc4g.Nfft2)))/length(noise_var(CntNumbSeg,loc4g.CP2+1:loc4g.CP2+loc4g.Nfft2)); %I
                    pow_Intbaseband=sum(Intbaseband(loc4g.CP2+1:loc4g.CP2+loc4g.Nfft2).*conj(Intbaseband(loc4g.CP2+1:loc4g.CP2+loc4g.Nfft2)))/length(Intbaseband(loc4g.CP2+1:loc4g.CP2+loc4g.Nfft2));
                    Total_pow_sig(CntNumbSeg,:)=pow_TxSig00(CntNumbSeg,1)+pow_noise_var_fft(CntNumbSeg,1)+pow_Intbaseband;
                    Total_pow_vn(CntNumbSeg,:)=pow_noise_var_fft(CntNumbSeg,1)+pow_Intbaseband;
                    
                    if CntNumbSeg==1
                        if strcmp(simode,'AWGN_Fadding_Interference')
                            %figure(29)
                            figure('visible','off');
                            subplot(2,1,1)
                            stem(Freq_IndexHz,real(fft_timedatainv2),'marker','none');
                            hold on;
                            plot(Freq_IndexHz,real(N_o_S),'r');
                            %plot(Freq_IndexHz,real(cmp_sigs(:,2)).','g');
                            hold off;
                            titlefig=strcat({['Amplitude Components per Sample Comparison'];['Blue: TxSig(n), Red: N+I/S'];['Symbol=',num2str(CntNumbSeg),'   Tower=',num2str(TowInd)]});
                            title(titlefig)
                            xlabel('Frequency Index MHz')
                            ylabel('Real Component')
                            
                            subplot(2,1,2)
                            stem(Freq_IndexHz,imag(fft_timedatainv2),'marker','none');
                            hold on;
                            plot(Freq_IndexHz,imag(N_o_S),'r');
                            %plot(Freq_IndexHz,imag(cmp_sigs(:,2).'),'g');
                            hold off;
                            xlabel('Frequency Index MHz')
                            ylabel('Imaginary Component')
                        else
                            subplot(2,1,1)
                            stem(Freq_IndexHz,real(fft_timedatainv2),'marker','none');
                            hold on;
                            plot(Freq_IndexHz,real(N_o_S),'r');
                            hold off;
                            titlefig=strcat({['Amplitude Components per Sample Comparison'];['Blue: TxSig(n), Red: N/S'];['Symbol=',num2str(CntNumbSeg),'   Tower=',num2str(TowInd)]});
                            title(titlefig)
                            xlabel('Frequency Index MHz')
                            ylabel('Real Component')
                            
                            subplot(2,1,2)
                            stem(Freq_IndexHz,imag(fft_timedatainv2),'marker','none');
                            hold on;
                            plot(Freq_IndexHz,imag(N_o_S),'r');
                            hold off;
                            xlabel('Frequency Index MHz')
                            ylabel('Imaginary Component')
                        end
                        saveas(gcf, fldhomecombfig0);
                    end
                    
%                     h00N_fad_test=[0,h00N(2:end)];
%                     h00N_los_test=h00N(1);
%                     fft_timedatainv(CntNumbSeg,:)=fft(timedatainv(CntNumbSeg,loc4g.CP+1:loc4g.CP+loc4g.Nfft),loc4g.Nfft);
%                     fft_fading(1,:)=fft(h00N_fad_test,loc4g.Nfft);
%                     fft_los(1,:)=fft(h00N_los_test,loc4g.Nfft);
%                     test_fading(CntNumbSeg,:)=fft_timedatainv(CntNumbSeg,:).*fft_fading(1,:);
%                     test_hlos(CntNumbSeg,:)=fft_timedatainv(CntNumbSeg,:).*fft_los(1,:);
%                     freq_noise_var_fft2(CntNumbSeg,:)=fft(noise_var(CntNumbSeg,loc4g.CP+1:loc4g.CP+loc4g.Nfft),loc4g.Nfft2);
                    timedataout(CntNumbSeg,:)=TxSig01(CntNumbSeg,1:end)+noise_var(CntNumbSeg,1:end);
%                     test_TxSig00(CntNumbSeg,:)=fft(TxSig00(CntNumbSeg,loc4g.CP+1:loc4g.CP+loc4g.Nfft));
%                     fft_timedata_out(CntNumbSeg,:)=fft(timedataout(CntNumbSeg,loc4g.CP+1:loc4g.CP+loc4g.Nfft),loc4g.Nfft);
%                     
%                     comp_table_h(:,1:8)=[fft_timedatainv(1,:).',fft_fading(1,:).',...
%                         fft_los(1,:).',test_fading(1,:).',test_hlos(1,:).',...
%                         test_TxSig00(1,:).',freq_noise_var_fft2(1,:).',...
%                         fft_timedata_out(1,:).'];
%                     comp_table_h(zcobj.bins_gain(1:zcobj.SegSize),9)=1;
                    comp_table_h=[];
                end
                TotalNumbSeg(mm,TowInd+1)=loc4g.NumbSeg; %Mod W45 V4
                [DelEstLongFreqC,DelSamples]=delayest(timedataout,loc4g,loc4g.cazacobj,zcobj,nvar_i,powSig,PowerGain,DelayToTarget_Samples(TowerIndex+1),nvar_ratio,Total_pow_vn,fldhomecombfig5,fldhomecombfig8,mm,TowInd,simode,comp_table_h,RiceanEst,fldhomecombfig9);
                close all
                DelEstLongFreq2(indexval+1,TowInd+1)=DelEstLongFreqC;
                DelSamples_Saved(mm,TowInd+1)=DelSamples;
                
            end
            
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End - Estimation Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        avgNumbSeg(mm,SimScen)=ceil(sum(TotalNumbSeg(mm,:))/loc4g.NumRefNodes); %Mod W45 V4
        
        Numbas=4;   %Basis nodes
        noiseval=0.1; %variance from true location in random; variance of the 2-Dim Gausian PDF in Cluster
        numtargets=1; %right now we only solve for  1 target
        
        InitLoc.nodes  = 16;  %Controles the number of different initial guesses
        InitLoc.method = 'locus';  %how the initial guess is selected: 'random', 'locus', 'offset' (offset: will force only 1 initial guess)
        InitLoc.offset = 1;  %offset value in offset initial guess mode
        
        initvec=InitVecSize;
        spatialvec=Numref;
        
        zerrmat=zeros(length(spatialvec), length(initvec));
        for xvind=1:length(initvec)
            for yvind=1:length(spatialvec)
                disp(yvind);
                
                DelEstForXYCoord=DelEstLongFreq2;
                
                [adhocobj, outputdata,longoutdata, outval]=RunAdHocNode(adhocobj,DelEstForXYCoord);  %run methods
                
                EZ = outval(end,1);
                OZ = outval(end,3);
                ErrVal= EZ-OZ;
                
            end
            
        end
        errvec(mm)=ErrVal;
        
    end
    avgNumbCarrSymb(SimScen,1)=ceil((NumberSubCarr*Numsim*Numref)/sum(sum(TotalNumbSeg(:,:)))); %Mod W45 V4
    avgNumbSegT(SimScen,1)=ceil(sum(sum(TotalNumbSeg(:,:)))/(loc4g.NumRefNodes*Numsim)); %Mod W45 V4
    
    ErrorStd=std(errvec);
    ErrorMean=mean(abs(errvec));
    errvec2=0;
    kit=0;
    for k=1:length(errvec)
        if abs(errvec(k)-ErrorMean) < ErrorStd;
            kit=kit+1;
            errvec2(kit)=errvec(k);
        end
    end
    if sum(abs(errvec2))==0
        errvec2=errvec;
    end
    
    ErrorStd=std(errvec2);
    ErrorMean=mean(abs(errvec2));
    errvec3=0;
    kit=0;
    for k=1:length(errvec2)
        if abs(errvec2(k)-ErrorMean) < ErrorStd;
            kit=kit+1;
            errvec3(kit)=errvec2(k);
        end
    end
    
    if sum(abs(errvec3))==0
        errvec3=errvec2;
    end
    
    ErrorStd=std(errvec3);
    ErrorMean=mean(abs(errvec3));
    errvec4=0;
    kit=0;
    for k=1:length(errvec3)
        if abs(errvec3(k)-ErrorMean) < ErrorStd;
            kit=kit+1;
            errvec4(kit)=errvec3(k);
        end
    end
    
    if sum(abs(errvec4))==0
        errvec4=errvec3;
    end
    
    ErrorStd=std(errvec4);
    ErrorMean=mean(abs(errvec4));
    errvec5=0;
    kit=0;
    for k=1:length(errvec4)
        if abs(errvec4(k)-ErrorMean) < ErrorStd;
            kit=kit+1;
            errvec5(kit)=errvec4(k);
        end
    end
    
    if sum(abs(errvec5))==0
        errvec5=errvec4;
    end
    
    %errvec2=errvec4; %Added
    
    fldstr3=char(strcat({'SimResult'},{num2str(SimScen)}));
    fldhomecomb1=char(strcat(fldhome, {'\'},fldstr3));
    fldhomecombfig6=char(strcat(fldhome, {'\'},fldstr3, {'\'}, {'fig3.jpg'}));
    fldhomecombfig7=char(strcat(fldhome, {'\'},fldstr3, {'\'}, {'fig2.jpg'}));
    fldhomecombdat3=char(strcat(fldhome, {'\'},fldstr3, {'\'}, {'data.mat'}));
    
    mkdir(fldhome,fldstr3)
    %%%%%%%%%%%%%%%%
    TargXY=get(adhocobj,'TargetVec');
    TowerXY0=get(adhocobj, 'TowerLoc');
    TowerXY=[TowerXY0, TowerXY0(1)];
    %%%%%%%%%%%%%%%%
    
    %figure(11)
    figure('visible','off');
    ErrorStd2=std(errvec5);
    ErrorMean2=mean(abs(errvec5));
    if Numsim==1 %Mod W45 V4
        ErrorAbs=abs(errvec5); %Mod W45 V4
    else %Mod W45 V4
        %ErrorAbs=sum(abs(errvec2)/Numsim); %Mod W45 V4
        ErrorAbs=ErrorMean2; %Mod W46 V1
    end %Mod W45 V4
    errvec6=errvec5+TargXY;
    plot(real(errvec6), imag(errvec6),'k.',real(TargXY), imag(TargXY), 'ro')
    legstr1=strcat({'SNR= '},{num2str(SNRvaldB)}, {'dB, RicK='},{num2str(RicianKdB)},....
        {'dB, NSymb='},{num2str(NumSymbols)}, {', AbsErr='},{num2str(ErrorAbs)},{', STD'},{num2str(ErrorStd2)});
    ylabel('Error y');
    xlabel('Error x');
    
    if strcmp(simode,'AWGN_Fadding_Interference')
        title({['Mode=',smode,'   ',titleSNR,num2str(SNRvaldB),'   Ref.Nodes=',num2str(Numref),'   Ricean KdB=',num2str(RicianKdB)];['Sub.Carriers=',num2str(NumberSubCarr),'   Avg.Sub.Carriers/Symb=',num2str(avgNumbCarrSymb(SimScen,1)),'   Avg.Numb.Symbs=',num2str(avgNumbSegT(SimScen,1)),];['Interference Ratio=',num2str(int_ratio),'  Noise Ratio=',num2str(nvar_ratio)];['Mean=',num2str(ErrorAbs),'m','  STD=',num2str(ErrorStd2)]})
    else                                                                                                                                                                                                   %loc4g.SegSize                                        %loc4g.NumbSeg %Mod W45 V4
        title({['Mode=',smode,'   ',titleSNR,num2str(SNRvaldB),'   Ref.Nodes=',num2str(Numref),'   Ricean KdB=',num2str(RicianKdB)];['Sub.Carriers=',num2str(NumberSubCarr),'   Avg.Sub.Carriers/Symb=',num2str(avgNumbCarrSymb(SimScen,1)),'   Avg.Numb.Symbs=',num2str(avgNumbSegT(SimScen,1)),];['Noise Ratio=',num2str(nvar_ratio),'   Mean=',num2str(ErrorAbs),'m','  STD=',num2str(ErrorStd2)]})
    end                                                                                                                                                                                                    %loc4g.SegSize                                        %loc4g.NumbSeg %Mod W45 V4
    
    saveas(gcf, fldhomecombfig6);
    
    TimeWindow=NumbSeg*(CPsize+512)*(1/Fsval);
    legstr1=strcat({'SNR= '},{num2str(SNRvaldB)}, {'dB, RicK='},{num2str(RicianKdB)},....
        {', uErr='},{num2str(ErrorMean2)},{', STD'},{num2str(ErrorStd2)});
    legstr2=strcat({'X Err(m), total Area='},{num2str(areaval)}, {'m^2'});
    legstr3=strcat({'Time Window='},{num2str(TimeWindow)},{' s'});
    %figure(9)
    figure('visible','off');
    plot(real(errvec5+TargXY), imag(errvec5+TargXY),'k.', real(TowerXY), imag(TowerXY), 'k', real(TargXY), imag(TargXY), 'ro');
    axis([-4e04 4e04 -4e04 4e04])
    ylabel(legstr3);
    xlabel(legstr2);

    if strcmp(simode,'AWGN_Fadding_Interference')
        title({['Mode=',smode,'   ',titleSNR,num2str(SNRvaldB),'   Ref.Nodes=',num2str(Numref),'   Ricean KdB=',num2str(RicianKdB)];['Sub.Carriers=',num2str(NumberSubCarr),'   Avg.Sub.Carriers/Symb=',num2str(avgNumbCarrSymb(SimScen,1)),'   Avg.Numb.Symbs=',num2str(avgNumbSegT(SimScen,1)),];['Interference Ratio=',num2str(int_ratio),'  Noise Ratio=',num2str(nvar_ratio)];['Mean=',num2str(ErrorAbs),'m','  STD=',num2str(ErrorStd2)]})
    else                                                                                                                                                                                                   %loc4g.SegSize                                        %loc4g.NumbSeg %Mod W45 V4
        title({['Mode=',smode,'   ',titleSNR,num2str(SNRvaldB),'   Ref.Nodes=',num2str(Numref),'   Ricean KdB=',num2str(RicianKdB)];['Sub.Carriers=',num2str(NumberSubCarr),'   Avg.Sub.Carriers/Symb=',num2str(avgNumbCarrSymb(SimScen,1)),'   Avg.Numb.Symbs=',num2str(avgNumbSegT(SimScen,1)),];['Noise Ratio=',num2str(nvar_ratio),'   Mean=',num2str(ErrorAbs),'m','  STD=',num2str(ErrorStd2)]})
    end                                                                                                                                                                                                    %loc4g.SegSize                                        %loc4g.NumbSeg %Mod W45 V4
    
    saveas(gcf, fldhomecombfig7);
    %print(gcf, '-djpeg', fldhomecombfig7,'-r500');
    temp_errvec2=abs(errvec5);
    %saved_abs_error(SimScen,1:length(temp_errvec2))=temp_errvec2(1,1:end); %Mod W45 V4
    saved_abs_error(SimScen,1:length(temp_errvec2))=temp_errvec2; %Mod W45 V4
    save(fldhomecombdat3, 'TargXY', 'TowerXY', 'errvec', 'errvec2', 'errvec3', 'errvec4', 'errvec5','errvec6','ErrorMean', 'ErrorMean2', 'ErrorStd', 'ErrorStd2','saved_abs_error','Target_Locations_MC','DelayToTarget_Samples','DelSamples_Saved','DelayToTarget','TotalNumbSeg','avgNumbSeg','avgNumbCarrSymb','avgNumbSegT');
    
    clearvars -except areaval InitVecSize DelaySpread InitMethod NumIter frametype numtargets noiseval Numsim Fsval Fsval2 CPsize simode MonteCarloSims RangeV SimResP1 fldhome...
              nvar_ratio int_ratio smode titleSNR SynchWindow format avgNumbSeg DelSamples_Saved DelayToTarget DelayToTarget_Samples TargXY Target_Locations_MC TotalNumbSeg...
              TowerXY avgNumbCarrSymb avgNumbSeg SimScen avgNumbSegT
    close all      
    noiseval=0;
    format shortg;
    dbg77=1;
end
