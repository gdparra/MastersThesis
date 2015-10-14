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
radioloc.pdp    =    [  -8,   -6,  -4 ,   0,   0,   -4 ,     -8,     -9,   -10 ,  -12 ,  -14  ]; %dB scattering response!
radioloc.PathDelays= [ 0.1,  0.3,  0.5, 0.7, 1.0,   1.3,   15.0,   15.2,   15.7,  17.2,   20.0]*0.2e-6;  %in seconds
radioloc.ts=1/fs;

%not moving (zero velocity)
radioloc.dopp=0;
%%%%%%%%%%%%%%%%%%%%%
%New Approach
figure(1)
hold on
Mchan=5;
for Tx=1:Mchan
    for Rx=1:Mchan
index = randperm(numel(radioloc.pdp));
B = radioloc.pdp(index);
[radioloc.h00t] = GenRayleigh(B, radioloc.PathDelays, radioloc.ts,  radioloc.dopp);

%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%
  %convert the continuous rician channel to discretely sampled
  radioloc.h00N0=MultiRes(radioloc.h00t, radioloc.PathDelays, radioloc.ts, 'NonLinear');
  
  radioloc.h00N = radioloc.h00N0/sqrt(sum(radioloc.h00N0.*conj(radioloc.h00N0)));
  radiochan(1:3,Rx+(Tx-1)*Mchan)=radioloc.h00N;
  plot(20*log10(abs(fft(radiochan(1:3,Rx+(Tx-1)*Mchan),64))));
    end
end

for k=1:Mchan^2
        [radiochan(1:3,k).']
end
  %Scattering 
  %print first five values
 
  channelgain= sum(radioloc.h00N.*conj(radioloc.h00N));
  
  
