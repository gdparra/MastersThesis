function [cazac] = freq2time(varargin)

cazac=varargin{1};
TowInd=varargin{2};
DelayToTarget=varargin{3};
CSymbols_PT=varargin{4};
PowerGain=varargin{5};
bind1=varargin{6};
Intbaseband=varargin{7};
h00N=varargin{8};
indexval=cazac.symbind;
Nfft=cazac.Nfft*cazac.USR;

FSymbol=zeros(Nfft,cazac.NumbSeg);
FSymbolO=zeros(Nfft,cazac.NumbSeg);
FSymbolG=zeros(cazac.NumbSeg,Nfft);
FSymbol_R=zeros(cazac.NumbSeg,Nfft);
FSymbolG2=zeros(cazac.NumbSeg,Nfft);
FSymbol_g0=zeros(cazac.NumbSeg,Nfft);
FSymbol_ga=zeros(cazac.NumbSeg,Nfft);
FSymbol_gb=zeros(cazac.NumbSeg,Nfft);

CP=cazac.CP*cazac.USR;
%TsymbolLong=zeros((Nfft+CP)*cazac.Nsymb,1);
lind=0;
nsymbind=cazac.symbind;


indval= [cazac.Nfft+cazac.ML+1:cazac.Nfft,2:cazac.MU+1];

%Remove bins 1 and 257 availability for use
A = 1;
B=(Nfft/2)+1;
Z1 = bind1(find(bind1~=A));
Z2 = Z1(find(Z1~=B));

%=========================================================
%Limiting Parameters
%=========================================================

for CntNumbSeg=1:cazac.NumbSeg;
    if CntNumbSeg<cazac.NumbSeg
        minrange=((CntNumbSeg-1)*cazac.SegSize)+1;
        maxrange=CntNumbSeg*cazac.SegSize;
    else
        minrange=((CntNumbSeg-1)*cazac.SegSize)+1;
        maxrange=cazac.NREf;
    end
    
    
    NumRE=maxrange-(minrange-1);
    %=========================================================
    %End of Limiting Parameters
    %=========================================================
    %cazac.NREf replaced by cazac.SegSize
    %FSymbolO(Z2(1:cazac.NREf))=cazac.REarr(1:cazac.NREf,indexval);
    FSymbolO(Z2(1:NumRE),CntNumbSeg)=cazac.REarr(minrange:maxrange,indexval); %cazac.NREf replaced by cazac.SegSize
    %FSymbolO(Z2(1:cazac.NREf))=1;
    %FSymbolO(Z2(1:NumRE),CntNumbSeg)=0;
    DelaySamples=DelayToTarget(1,TowInd)/cazac.ts2;
    %DelaySamples=0;
    %=========================================================
    %Delay Modifier placed due to Ricean cancelation on delest
    %=========================================================
    for k=0:Nfft-1
        FSymbol(k+1,CntNumbSeg)=FSymbolO(k+1,CntNumbSeg)*exp(-1i*2*pi*k*((DelaySamples/Nfft)));
    end
    
    %=========================================================
    %Section added to calculate power without gain
    %=========================================================
    Frequency_Power_RE(:,CntNumbSeg)=abs(FSymbol(:,CntNumbSeg));
%     Frequency_Power_AVG(:,CntNumbSeg)=sum((FSymbol(:,CntNumbSeg).*conj(FSymbol(:,CntNumbSeg)))/length(FSymbol(:,CntNumbSeg)));
    Frequency_Power_AVG(CntNumbSeg,1)=sum(FSymbol(:,CntNumbSeg).*conj(FSymbol(:,CntNumbSeg)));
    tmpvec_ng(:,CntNumbSeg)=ifft(FSymbol(1:end,CntNumbSeg),Nfft);
    powSig_ng(1,CntNumbSeg)=sum(tmpvec_ng(1:end,CntNumbSeg).*conj(tmpvec_ng(1:end,CntNumbSeg)))/length(tmpvec_ng(1:end,CntNumbSeg)); %I
    nvar_ng(1,CntNumbSeg)=powSig_ng(1,CntNumbSeg)/10^(cazac.SNRdB/10); %I
    
    %=========================================================
    %FSymbolG=zeros(cazac.NumbSeg,Nfft);
    for k=0:NumRE-1
        FSymbolG(CntNumbSeg,Z2(k+1))=FSymbol(Z2(k+1),CntNumbSeg)*PowerGain(Z2(k+1));
    end
    %=========================================================
    
    %=========================================================
    %Section added to calculate power with gain
    %=========================================================
    Frequency_Power_RE_G(CntNumbSeg,:)=abs(FSymbolG(CntNumbSeg,:));
    Frequency_Power_AVG_G(CntNumbSeg,:)=sum(FSymbolG(CntNumbSeg,:).*conj(FSymbolG(CntNumbSeg,:)));
    
    %=========================================================
    %Change to Time (Feedback to Freq)
    %=========================================================
    tmpvec(CntNumbSeg,:)=ifft(FSymbolG(CntNumbSeg,:),Nfft);
    %Tsymbol(CntNumbSeg,:)=[tmpvec(CntNumbSeg,end-CP+1:end);tmpvec(CntNumbSeg,:)].';
    powSig_g(CntNumbSeg,1)=sum(tmpvec(CntNumbSeg,1:end).*conj(tmpvec(CntNumbSeg,1:end)))/length(tmpvec(CntNumbSeg,1:end)); %I
    nvar_g(CntNumbSeg,1)=powSig_g(CntNumbSeg,1)/10^(cazac.SNRdB/10); %I
%     Previous Power Gain 2    
%     measure1_powSig(CntNumbSeg,1)=(1/powSig_g(CntNumbSeg,1));
%     PowerGain2(CntNumbSeg,1)=sqrt(measure1_powSig(CntNumbSeg,1));
    PowerGain2(CntNumbSeg,1)=sqrt((Nfft^2)/Frequency_Power_AVG_G(CntNumbSeg,:));
    W_gain=(conj(PowerGain)./((conj(PowerGain).*PowerGain)));
    %FSymbol_R=zeros(cazac.NumbSeg,Nfft);
    for k=0:NumRE-1 %replacement of cazac.NREf by NumRE
        FSymbol_R(CntNumbSeg,Z2(k+1))=FSymbolG(CntNumbSeg,Z2(k+1)).*(W_gain(Z2(k+1)));
    end
    
    FSymbol_mat=[FSymbol(Z2(1:NumRE),CntNumbSeg),FSymbol_R(CntNumbSeg,Z2(1:NumRE)).'];
    %=========================================================
    %Section added to calculate power with gain 2
    %=========================================================
    %FSymbolG2=zeros(cazac.NumbSeg,Nfft);
    FSymbolG2(CntNumbSeg,Z2)=FSymbolG(CntNumbSeg,Z2)*PowerGain2(CntNumbSeg,1);
    Frequency_Power_RE_G2(CntNumbSeg,:)=abs(FSymbolG2(CntNumbSeg,:));
    Frequency_Power_AVG_G2(CntNumbSeg,1)=sum(FSymbolG2(CntNumbSeg,:).*conj(FSymbolG2(CntNumbSeg,:)));
    
    %=========================================================
    %Change to Time (After Feedback)
    %=========================================================
    %Reference Only
    tmpvec2_ref(CntNumbSeg,:)=ifft(FSymbolG2(CntNumbSeg,:),Nfft);
    Tsymbol2_ref(:,CntNumbSeg)=[tmpvec2_ref(CntNumbSeg,end-CP+1:end),tmpvec2_ref(CntNumbSeg,:)].';
    powSig_g2_ref(CntNumbSeg,1)=sum(tmpvec2_ref(CntNumbSeg,1:end).*conj(tmpvec2_ref(CntNumbSeg,1:end)))/length(tmpvec2_ref(CntNumbSeg,1:end));
    %end of Reference
    tmpvec2(CntNumbSeg,:)=ifft(FSymbolG2(CntNumbSeg,:),Nfft); % New Addition *Nfft
    Tsymbol2(:,CntNumbSeg)=[tmpvec2(CntNumbSeg,end-CP+1:end),tmpvec2(CntNumbSeg,:)].';
    powSig_g2(CntNumbSeg,1)=sum(tmpvec2(CntNumbSeg,1:end).*conj(tmpvec2(CntNumbSeg,1:end)))/length(tmpvec2(CntNumbSeg,1:end)); %I
    
    %FSymbol_g0=zeros(cazac.NumbSeg,Nfft);
    for k=0:NumRE-1 %replacement of cazac.NREf by NumRE
        FSymbol_g0(CntNumbSeg,Z2(k+1))=FSymbolO(Z2(k+1),CntNumbSeg)*PowerGain(Z2(k+1))*PowerGain2(CntNumbSeg,1); % New Addition *Nfft
    end
                                      
    temp_ifft_g0=ifft(FSymbol_g0(CntNumbSeg,:),Nfft);
    FSymbol_g0(CntNumbSeg,:)=fft(temp_ifft_g0,Nfft);
    
    test(:,1)=FSymbol(:,1);
    test(:,2)=FSymbolO(:,1);
    test(:,3)=test(:,1)-test(:,2);
    
    for k=0:NumRE-1 %replacement of cazac.NREf by NumRE
        FSymbol_ga(CntNumbSeg,Z2(k+1))=FSymbolO(Z2(k+1),CntNumbSeg)*PowerGain(Z2(k+1)); % New Addition *Nfft
    end
    
    test(:,5)=FSymbolG(1,:).';
    test(:,6)=FSymbol_ga(1,:).';
    test(:,7)=test(:,5)-test(:,6);
    
    for k=0:NumRE-1 %replacement of cazac.NREf by NumRE
        FSymbol_gb(CntNumbSeg,Z2(k+1))=FSymbolO(Z2(k+1),CntNumbSeg)*PowerGain(Z2(k+1))*PowerGain2(CntNumbSeg,1); % New Addition *Nfft
    end
    
    test(:,9)=FSymbolG2(1,:).';
    test(:,10)=FSymbol_gb(1,:).';
    test(:,11)=test(:,9)-test(:,10);
    
    FSymbolG2_Norm(CntNumbSeg,:)=fft(Tsymbol2(CP+1:end,CntNumbSeg),Nfft);
    Frequency_Power_AVG_G2_Norm(CntNumbSeg,1)=sum(FSymbolG2_Norm(CntNumbSeg,:).*conj(FSymbolG2_Norm(CntNumbSeg,:)));
    FSymbolG2_mat=[FSymbolG2_Norm(CntNumbSeg,Z2(1:NumRE)).',FSymbol_g0(CntNumbSeg,Z2(1:NumRE)).'];
    
%     for k=0:Nfft-1
%         FSymbol_gx(k+1,CntNumbSeg)=FSymbol_g0(CntNumbSeg,k+1)*exp(-1i*2*pi*k*((DelaySamples/Nfft)));
%     end
%     
%     FSymbolG2_mat_gx=zeros(Nfft,3);
%     FSymbolG2_mat_gx(:,1)=FSymbolG2_Norm(CntNumbSeg,:).';
%     FSymbolG2_mat_gx(:,2)=FSymbol_gx(:,1);
%     FSymbolG2_mat_gx(:,3)=FSymbolG2_mat_gx(:,1)-FSymbolG2_mat_gx(:,2);
%     
%     tmpvec2(CntNumbSeg,:)=ifft(FSymbol_gx(:,CntNumbSeg),Nfft); % New Addition *Nfft
%     Tsymbol2(:,CntNumbSeg)=[tmpvec2(CntNumbSeg,end-CP+1:end),tmpvec2(CntNumbSeg,:)].';
%     powSig_g2(CntNumbSeg,1)=sum(tmpvec2(CntNumbSeg,1:end).*conj(tmpvec2(CntNumbSeg,1:end)))/length(tmpvec2(CntNumbSeg,1:end)); %I

    
end
cazac.VecLong=Tsymbol2;
cazac.bins_gain=Z2;
cazac.powSig=powSig_g2;
cazac.PowerGain2=PowerGain2;
cazac.nvar_ng=nvar_ng;
cazac.nvar_g=nvar_g;
cazac.powSig_ng=powSig_ng;
cazac.powSig_g=powSig_g;
cazac.FSymbol_g0=FSymbol_g0;
cazac.FSymbolG(1:cazac.NumbSeg,1:cazac.SegSize)=FSymbolG(1:cazac.NumbSeg,Z2(1:cazac.SegSize));
cazac.FSymbolG2(1:cazac.NumbSeg,1:cazac.SegSize)=FSymbolG2(1:cazac.NumbSeg,Z2(1:cazac.SegSize));

dbg77=1;

