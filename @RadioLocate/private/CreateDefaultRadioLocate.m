function radioloc = CreateDefaultRadioLocate

%radioloc.Band   = 170e3;
Band   = 11*15e3;
DF     = 15e3;
NREf   = floor(obj.Band/obj.DF); 
freqvec= [-(obj.NREf-1)/2:(obj.NREf-1)/2];
Nsymb  = 1023; %Number of time domain OFDM systems
CP     = 2*161+1;
REarr  = zeros(obj.NREf,obj.Nsymb);
fs     = 7.68e6;
ts     =1/obj.fs;
Nfft   = obj.fs/obj.DF;
fc     = 3.315e6;
USR    = 1;
rprime = 23;
staticdelay = 1e-5;
Nfft2  = obj.Nfft*obj.USR;
CP2    = obj.CP  *obj.USR;
%radioloc.cazacobj = [];
ts2      = obj.ts/obj.USR;
GRA      =zeros(obj.NREf, obj.Nsymb);
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
Numref     =obj.spatialvec(1);  %Number of nodes (excluding basis) that have a distinct location geography; 
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
%%%%%%%%%%%%%%%%%%%%Coordinate System Solver%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
