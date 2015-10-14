classdef ChannEst < handle
    %ChannEst Estimates Ricean Channel at the Receiver based on the h00N
    %previously defined
    properties
        h00N=1;
        DelayToTarget_Samples=[];
        Nfft=512;
        h00Nscat=[];
        h00Nlos=[];
        H=[];
        Hscat=[];
        Hlos=[];
        HDel=[];
        HscatDel=[];
        HlosDel=[];
        h00NDel=[];
        fading=[];
        los=[];
        SNR=[];
        powSig=[];
        nvar_ratio=[];
        simode=[];
    end
    
    methods
        function obj= ChannEst(varargin)
            for k=1:nargin
                switch k
                    case 1
                        obj.h00N=varargin{1};
                    case 2
                        obj.DelayToTarget_Samples=varargin{2};
                    case 3
                        obj.Nfft=varargin{3};
                    case 4
                        obj.SNR=varargin{4};
                    case 5
                        obj.powSig=varargin{5};
                    case 6
                        obj.nvar_ratio=varargin{6};
                    case 7
                        obj.simode=varargin{7};
                end
            end
            %Functions are Set here
            
            %Time Domain
            obj.h00Nscat = [0,obj.h00N(2:end)];
            obj.h00Nlos =  [obj.h00N(1),zeros(size(obj.h00N(2:end)))];
            
            obj.Nfft=length(obj.h00N);
            
            nvar_i=obj.powSig(1,1)/10^(obj.SNR/10);
            nvar(1,1)=nvar_i(1,1)*obj.nvar_ratio;
            noise_var(1,:)=(sqrt(nvar/2)*randn(size(obj.h00N))+1i*sqrt(nvar(1,1)/2)*randn(size(obj.h00N)));
            temph00N=obj.h00N;%+noise_var;
            
            %Frequency DOmain
            H=fft(temph00N,obj.Nfft);
            Hscat=fft(obj.h00Nscat(1:end),obj.Nfft);
            Hlos=fft(obj.h00Nlos(1:end),obj.Nfft);
            Hsum=Hlos+Hscat;
            res_mat1=[H.',Hsum.'];
            
            for m=1:length(obj.DelayToTarget_Samples)
                
                phaseshift=obj.DelayToTarget_Samples(m);
%%%%%W6 V2A                 
%                 for k=0:length(H)-1
%                     HDel(m,k+1)=H(k+1)*exp(-1i*2*pi*k*(phaseshift/obj.Nfft));
%                 end
                  HDel(m,:)=H(:);%W6 V2A 
%                 for k=0:length(H)-1
%                     ChannUnDel(m,k+1)=HDel(m,k+1)*exp(1i*2*pi*k*(phaseshift/obj.Nfft));
%                 end
%%%%%W6 V2A                 
                h00NDel(m,:)=ifft(HDel(m,:), length(obj.h00N))+noise_var;
                TrainVec = [1, 3, -1, 7, 2];
                zvec = conv(TrainVec, h00NDel(m,:));
                xsize=length(TrainVec);
                ChanSize=length(h00NDel(m,:));
                TrainMat=zeros(xsize+ChanSize-1,ChanSize);
                for k=0:size(TrainMat,2)-1
                    TrainMat(k+1:k+xsize,k+1)=TrainVec.';
                end
                XMMSE=TrainMat;
                difval=length(zvec)-size(XMMSE,1);
                if difval <= 0
                    newh00NDel=inv(XMMSE'*XMMSE)*XMMSE'*[zvec, zeros(1,-difval)].';
                else
                    newh00NDel=inv(XMMSE'*XMMSE)*XMMSE'*zvec(1:end-difval).';
                end
                h00NDel(m,:)=newh00NDel;
%%%%%W6 V2A                 
                if strcmp(obj.simode,'AWGN')
                    h00NDel(m,:)=obj.h00N;
                end
%%%%%W6 V2A                 
                %Theory 1: Not Valid
                for k=0:length(H)-1
                    HscatDel(m,k+1)=Hscat(k+1)*exp(-1i*2*pi*k*(phaseshift/obj.Nfft));
                end
                
                for k=0:length(H)-1
                    HlosDel(m,k+1)=Hlos(k+1)*exp(-1i*2*pi*k*(phaseshift/obj.Nfft));
                end
                
                HsumDel(m,:)=HlosDel(m,:)+HscatDel(m,:);
                res_mat2=[HDel(1,:).',HsumDel(1,:).'];
                %%%%%%%%%%%%%%%%%%%%
                
                %Theory 2: Valid.
                h00NscatDel(m,:)=[0,h00NDel(m,2:end)];
                h00NlosDel(m,:)=[h00NDel(m,1),zeros(1,length(h00NDel)-1)];
                
                fading(m,:)=fft(h00NscatDel(m,:),obj.Nfft);
                los(m,:)=fft(h00NlosDel(m,:),obj.Nfft);
                %%%%%%%%%%%%%%%%%%%%
                res_mat3=[HscatDel(m,:).',fading(m,:).'];
                res_mat4=[HlosDel(m,:).',los(m,:).'];
                
            end
            
            obj.H=H;
            obj.Hscat=Hscat;
            obj.Hlos=Hlos;
            obj.HDel=HDel;
            obj.HscatDel=fading;
            obj.HlosDel=los;
            obj.h00NDel=h00NDel;
            obj.fading=fading;
            obj.los=los;
        end
        
    end
    
end

