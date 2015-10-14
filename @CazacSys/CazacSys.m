classdef CazacSys < handle
    
    %dbstop if error
    
    properties
        Band   = 170e3;
        DF     = 15e3;
        NREf   =[];
        P      = 23;  %Relative prime
        Nsymb  = 1023; %Number of time domain OFDM systems
        CP     = 37;
        REarr  =[];
        ML     = [];
        MU     = [];
        fs     = 7680000;
        Nfft   = 512;
        vec    = [];
        sdr    = 'yes';
        fc     = 3.315e6;
        delay  = 0;
        USR    = 1;
        fs2    = [];
        ts2    = [];
        VecLong= [];
        VecTone= [];
        vecLT  = [];
        symbind= 53535353;
        Nfft2  = [];
        CP2    = [];
        synchpt= [];
        frametype=[];
        bins_gain=[];
        powSig=[];
        PowerGain2=[];
        SNRdB=[];
        nvar_ng=[];
        nvar_g=[];
        powSig_ng=[];
        powSig_g=[];
        FSymbol_g0=[];
        FSymbolG=[];
        FSymbolG2=[];
        SegSize=[];
        NumbSeg=[];
    end
    
    methods   %constructor method
        function obj=CazacSys(varargin)
            valcaz=varargin{1};
            obj.NREf   = floor(obj.Band/obj.DF); %Number of cazac Subcarriers
            obj.REarr  = zeros(obj.NREf,obj.Nsymb);
            obj.fs2    = obj.fs*obj.USR;
            obj.ts2    = 1/obj.fs2;
            obj.VecLong= zeros((obj.Nfft+obj.CP)*obj.USR,1);
            obj.VecTone= zeros((obj.Nfft+obj.CP)*obj.USR,1) ;
            obj.vecLT  = zeros((obj.Nfft+obj.CP)*obj.USR,1);
            obj.Nfft2  = obj.Nfft*obj.USR;
            obj.CP2    = obj.CP  *obj.USR;
            obj.synchpt= obj.CP;
            obj.vec    = zeros(obj.Nfft,1);
            obj.Band  =valcaz.Band;
            obj.DF    =valcaz.DF;
            obj.NREf  =valcaz.NREf;
            obj.P     =valcaz.P;
            obj.Nsymb =valcaz.Nsymb;
            obj.CP    =valcaz.CP;
            obj.fs    =valcaz.fs;
            obj.Nfft  =valcaz.Nfft;
            obj.fc    =valcaz.fc;
            obj.delay =valcaz.delay;
            obj.USR   =valcaz.USR;
            obj.fs2   =valcaz.fs2;
            obj.ts2   =valcaz.ts2;
            obj.Nfft2  = obj.Nfft*obj.USR;
            obj.CP2    = obj.CP  *obj.USR;
            obj.frametype=valcaz.frametype;
            obj.SNRdB=valcaz.SNRdB;
            obj.SegSize=valcaz.SegSize;
            obj.NumbSeg=valcaz.NumbSeg;
            
            
            g=zeros(1,obj.Nfft);
            gL=zeros(1,obj.Nfft2+obj.CP2);
            N=obj.Nsymb;
            p=obj.P;
            g=ones(1,N);
%             switch rem(N,2)
%                 case 0
%                     for n=0:N-1
%                         g(n+1)=exp(-1i*(2*pi/N)*p*n^2/2);
%                     end
%                 case 1
%                     for n=0:N-1
%                         g(n+1)=exp(-1i*(2*pi/N)*p*n*(n+1)/2);
%                     end
%                     
%             end
            
            N=obj.Nfft2+obj.CP2;
            switch rem(obj.Nfft2,2)
                case 0
                    for n=0:N-1
                        gL(n+1)=exp(-1i*(2*pi/N)*p*n^2/2);
                    end
                case 1
                    for n=0:N-1
                        gL(n+1)=exp(-1i*(2*pi/N)*p*n*(n+1)/2);
                    end
            end
            
            
            %N=11;
            N=obj.NREf;
            gF=zeros(1,N);
            if obj.frametype==2
                gF(n+1)=exp(-1i*(2*pi/4)*randi([0 3]));
            elseif obj.frametype==0 || obj.frametype==1
                switch rem(N,2)
                    case 0
                        for n=0:N-1
                            gF(n+1)=exp(-1i*(2*pi/N)*p*n^2/2);
                        end
                        obj.ML=-N/2;
                        obj.MU=N/2;
                    case 1
                        for n=0:N-1
                            gF(n+1)=exp(-1i*(2*pi/N)*p*n*(n+1)/2);
                        end
                        obj.ML=-(N+1)/2;
                        obj.MU=-obj.ML-1;
                end
            end

            obj.REarr=zeros(obj.NREf, obj.Nsymb);
            for k=1:obj.Nsymb
                for p=1:obj.NREf
                    obj.REarr(p,k)=gF(p)*g(k);
                end
            end
            
            %cazacsys.REarr=repmat(g, cazacsys.NREf,1);
            obj.vec=g;
            obj.vecLT=gL;
            
            
        end
    end
end
