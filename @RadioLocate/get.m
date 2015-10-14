function value = get(radioloc, property)

switch property
% case 'cazac'
%     value=radioloc.cazacobj;
case 'DeltaF'
    value = radioloc.DF;
case 'NREf'
    value = radioloc.NREf;   
case 'Nsymb'
    value = radioloc.Nsymb;
case 'GRA'
    value = radioloc.GRA;
case 'Nfft'
     value=radioloc.Nfft;         
case 'gr'
     value=radioloc.vec;
case 'Nfft2'
     value=radioloc.Nfft2;
case 'CP2'
     value=radioloc.CP2;
case 'FreqISI'
     value=radioloc.FREQISI;
case 'SNR'
     value=radioloc.SNR;
case 'DiscreteChan'
     value =radioloc.h00N;
case 'TS2'
     value=radioloc.ts2; 
case 'H00T'
     value=radioloc.h00t;
case 'HQFREQ'
     value=radioloc.HQFREQ2;
case 'HQFREQscat'
     value=radioloc.HQFREQscat;
case 'HQFREQLos'
     value=radioloc.HQFREQLos;     
case 'GeoDelays'
     value=radioloc.DelayToTarget;
case 'LOSDelay'
     value=radioloc.pathdel;
case 'synchoff'
    value=radioloc.deoff;
case 'Nscat'
     value=radioloc.h00Nscat;        
case 'LI'
     value=radioloc.LI;
case 'UI'
    value=radioloc.UI;
case 'ChannelDbg'
    value=radioloc.ChanDbg;
case 'ChannelDbg2'
    value=radioloc.HQFREQLosB;
case 'ChannelHC'
    value=radioloc.HQFREQLosC;
case 'SincF'
     value=radioloc.h00Nsinc;
case 'HOONLos'
     value=radioloc.h00Nlos;
end
