classdef RadioLocate < handle
%function [radioloc, cazacobj] = RadioLocate(varargin)
%dbstop if error
    
properties 
clusternodes =[];
NumRefNodes =[];
cazacobj = [];
Band   = 11*15e3;
DF     = 15e3;
NREf = [];
freqvec = [];
Nsymb  = 1023; %Number of time domain OFDM systems
CP     = 2*161+1;
GRA=[];
REarr=[];
fs     = 7.68e6;
ts = [];
Nfft = [];
fc     = 3.315e6;
USR    = 1;
rprime = 23;
staticdelay = 1e-5;
Nfft2 = [];
CP2 = [];
%radioloc.cazacobj = [];
ts2  = [];
%%%%%%%%%%%%%%%%%%%%%%%%
LI=[];
UI=[];
SNR    = 0;
dopp   = 0; 
frame   = [1];  %0: Chan. Est., 1: Radio Location, 2: Information
riciankdB= 8;  %Rician K factor in dB
pdp    = [   -8  , -6  , -4  ,   0  ,   0  ,  -4  ,  -8   ,   -9  ,  -10  ,  -12  ,  -14  ]; %dB scattering response!
ricipathdel= [0,  0.1,  0.3,  0.5,   0.7,   1.0,   1.3,   15.0,   15.2,   15.7,   17.2,   20.0]*1e-6;  %in seconds
% radioloc.pdp             =  [      0   ,  0    ,  0   ,   0     ,   0   ,   0   ,      0    ,    0    ,     0    ,      0   ,      0    ]; %dB scattering response!
% radioloc.ricipathdel= [0,  0.1,  0.3,  0.5,   0.7,   1.0,   1.3,   15.0, 15.2,   15.7,   17.2,   20.0]*1e-6;  %in seconds
chanest= [];
h00t      =  [];
h00N      =  [];
h00Nscat  = [];
h00Nlos =[];
spatialvec= [5];
Numref =[];
pathdel     =  [];
HQFREQ2  = [];
HQFREQscat = [];
HQFREQLos     =[];
HQFREQLosB    =[];
HQFREQLosC    =[];
FREQISI   = [];   
ChanDbg=[];
ChanDbg2=[];
h00Nsinc=[];    
end

methods   %constructor method

function obj=RadioLocate(varargin)
obj.freqvec= [-(obj.NREf-1)/2:(obj.NREf-1)/2];
obj.NREf     = floor(obj.Band/obj.DF); 
obj.REarr   = zeros(obj.NREf,obj.Nsymb);
obj.ts           =1/obj.fs;
obj.Nfft       = obj.fs/obj.DF;
obj.Nfft2    = obj.Nfft*obj.USR;
obj.CP2      = obj.CP  *obj.USR;
obj.ts2        = obj.ts/obj.USR;
obj.GRA      =zeros(obj.NREf, obj.Nsymb);
obj.Numref=obj.spatialvec(1);  %Number of nodes (excluding basis) that have a distinct location geography;             

for k=1:nargin
    switch k
        case 1
           obj.SNR           =varargin{1};
        case 2
           obj.riciankdB     =varargin{2};
        case 3
           obj.NumRefNodes   =varargin{3};
        case 4
           obj.Nsymb         =varargin{4};
        case 5
           obj.frame         =varargin{5};
        case 6
           obj.pathdel       =varargin{6};
        case 7
            TowerIndex             =varargin{7};
        case 8
           obj.CP            =varargin{8};
        case 9
           obj.fs            =varargin{9};
           obj.ts     =1/obj.fs;
           obj.Nfft   =obj.fs/obj.DF;
          obj.ts2      =obj.ts/obj.USR;
        case 10
           obj.Band     = varargin{10};
        case 11
           obj.NREf=varargin{11};
    end
end

%radioloc.NREf   = floor(radioloc.Band/radioloc.DF); 
switch rem(obj.NREf,2)
    case 0
       obj.LI=-obj.NREf/2;
       obj.UI= obj.NREf/2;
    case 1
        obj.LI=  -(obj.NREf+1)/2;
        obj.UI=   (obj.NREf-1)/2;
end
obj.clusternodes  =obj.NREf-1;

Valcaz.Band    =obj.Band;
Valcaz.DF      =obj.DF;
Valcaz.NREf    =obj.NREf;
Valcaz.P       =obj.rprime;
Valcaz.Nsymb   =obj.Nsymb;
Valcaz.CP      =obj.CP;
Valcaz.fs      =obj.fs;
Valcaz.Nfft    =obj.fs/obj.DF;
Valcaz.fc      =obj.fc;
Valcaz.delay   =obj.staticdelay;
Valcaz.USR     =obj.USR;
Valcaz.fs2     =1/obj.ts2;
Valcaz.ts2     =obj.ts2;
Valcaz.frametype=obj.frame;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj.cazacobj=CazacSys(Valcaz);
obj.GRA=get(obj.cazacobj,'REArr');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%
%Prior Approach
%%%%%%%%%%%
PathDelays=obj.pathdel(TowerIndex)+obj.ricipathdel;

%%%%%%%%%%%%%%%%%%%%%
%New Approach
[obj.h00t] = GenRice(obj.riciankdB,obj.pdp, PathDelays, obj.ts,  obj.dopp);
  %%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%
  obj.h00N=MultiRes(obj.h00t, PathDelays, obj.ts, 'NonLinear');
  obj.h00Nscat=MultiRes([0, obj.h00t(2:end)], PathDelays, obj.ts, 'NonLinear');
  [obj.h00Nlos, obj.h00Nsinc, ToDel]=MultiRes([obj.h00t(1), zeros(size(obj.h00t(2:end)))], PathDelays, obj.ts, 'NonLinear');


  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
NFFT2=obj.Nfft2;
HQFREQp=fft(obj.h00N, NFFT2);  %
obj.HQFREQ2=[HQFREQp(1,NFFT2+obj.LI+1:NFFT2), HQFREQp(1,2:obj.UI+1)].';
                         
HQFREQtmp=fft(obj.h00Nscat,NFFT2);
obj.HQFREQscat=[HQFREQtmp(1,NFFT2+obj.LI+1:NFFT2), HQFREQtmp(1,2:obj.UI+1)].';
obj.FREQISI=repmat(obj.HQFREQscat,1,size(obj.GRA,2)).*...
                                  [obj.GRA(end+obj.LI+1:end,:); obj.GRA(1:obj.UI, :)];
                                        
HQFREQtmp=fft(obj.h00Nlos,NFFT2);
obj.HQFREQLos=[HQFREQtmp(1,NFFT2+obj.LI+1:NFFT2), ...
                                             HQFREQtmp(1,2:obj.UI+1)].';
 obj.ChanDbg=repmat(obj.HQFREQLos,1,size(obj.GRA,2)).*...
                                  [obj.GRA(end+obj.LI+1:end,:); obj.GRA(1:obj.UI, :)];

HQFREQtmp=fft(obj.h00Nsinc,NFFT2);                              
obj.HQFREQLosC=obj.h00t(1)*[HQFREQtmp(1,NFFT2+obj.LI+1:NFFT2), ...
                                             HQFREQtmp(1,2:obj.UI+1)].';
obj.HQFREQLosB=[HQFREQtmp(1,NFFT2+obj.LI+1:NFFT2), ...
                                             HQFREQtmp(1,2:obj.UI+1)].';
obj.ChanDbg2=repmat(obj.HQFREQLosC,1,size(obj.GRA,2)).*...
                                  [obj.GRA(end+obj.LI+1:end,:); obj.GRA(1:obj.UI, :)];
end
end
end

