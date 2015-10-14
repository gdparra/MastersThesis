function cazacsys = CreateDefaultCazacSys

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

Band   = 170e3;
DF     = 15e3;
NREf=[];
P      = 23;  %Relative prime
Nsymb  = 1023; %Number of time domain OFDM systems
CP     = 37;
REarr =[];
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
symbind= 1;
Nfft2  = [];
CP2    = [];
synchpt= [];
frametype=[];



