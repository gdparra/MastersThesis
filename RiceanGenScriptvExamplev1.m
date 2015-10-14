%channel is defined at rate fs
fs=960e3;
%Line of site gain
tau0=(11*1/fs);
NLoop=1000;
yvec=zeros(NLoop,1);
channelgain=zeros(NLoop,1);
%for loopK=1:NLoop
    
%this is the average value over a suite of simulation trials (varies per trial)    
radioloc.riciankdB=8;

%scattering gain (power delay profile)
radioloc.pdp    = [       -8  , -6  , -4  ,   0  ,   0  ,  -4  ,  -8   ,   -9  ,  -10  ,  -12  ,  -14  ]; %dB scattering response!
radioloc.PathDelays= [0,  0.1,  0.3,  0.5,   0.7,   1.0,   1.3,   15.0,   15.2,   15.7,   17.2,   20.0]*1e-6;  %in seconds
radioloc.ts=1/fs;

%not moving (zero velocity)
radioloc.dopp=0;
%%%%%%%%%%%%%%%%%%%%%
%New Approach
[radioloc.h00t] = GenRice(radioloc.riciankdB,radioloc.pdp, radioloc.PathDelays, radioloc.ts,  radioloc.dopp);

%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%
  %convert the continuous rician channel to discretely sampled
  radioloc.h00N0=MultiRes(radioloc.h00t, radioloc.PathDelays, radioloc.ts, 'NonLinear');
  
  radioloc.h00N = radioloc.h00N0/sqrt(sum(radioloc.h00N0.*conj(radioloc.h00N0)));

  
  %Scattering 
  radioloc.h00Nscat = [0,radioloc.h00N(2:end)];
  GLos=sqrt(10^(radioloc.riciankdB/10)*(sum(radioloc.h00Nscat.*conj(radioloc.h00Nscat)))/...
      (radioloc.h00N(1)*conj(radioloc.h00N(1))));
    
  %Line of site
  radioloc.h00Nlos =  [GLos*radioloc.h00N(1),zeros(size(radioloc.h00N(2:end)))];
  
  %print first five values
  [radioloc.h00N(1:5).', radioloc.h00Nlos(1:5).', radioloc.h00Nscat(1:5).'];

  
  [radioloc.h00N(1:5).', radioloc.h00Nlos(1:5).', radioloc.h00Nscat(1:5).'];
  Kdbval=10*log10(sum(radioloc.h00Nlos.*conj(radioloc.h00Nlos))/...
      sum(radioloc.h00Nscat.*conj(radioloc.h00Nscat)));
  
  yvec=Kdbval;
  channelgain= sum(radioloc.h00N.*conj(radioloc.h00N));
  
  
