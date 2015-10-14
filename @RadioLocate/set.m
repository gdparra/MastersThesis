function radioloc = set(radioloc,property,value, value2)

GRA=value2;
switch property 
    case 'AcquireChannel'
        TowerIndex=value;
PathDelays=radioloc.pathdel(TowerIndex)+radioloc.ricipathdel;
NFFT2=radioloc.Nfft2;
[radioloc.h00t] = GenRice(radioloc.riciankdB,radioloc.pdp, PathDelays, radioloc.ts,  radioloc.dopp);
  
  radioloc.h00N=MultiRes(radioloc.h00t, PathDelays, radioloc.ts, 'NonLinear');
  radioloc.h00Nscat=MultiRes([0, radioloc.h00t(2:end)], PathDelays, radioloc.ts, 'NonLinear');
  [radioloc.h00Nlos, radioloc.h00Nsinc, ToDel]=MultiRes([radioloc.h00t(1), zeros(size(radioloc.h00t(2:end)))], PathDelays, radioloc.ts, 'NonLinear');




%====================
HQFREQp=fft(radioloc.h00N, NFFT2);  %relevance ??
radioloc.HQFREQ2=[HQFREQp(1,NFFT2+radioloc.LI+1:NFFT2), ...
                                             HQFREQp(1,2:radioloc.UI+1)].';
%radioloc.ChanDbg=repmat(radioloc.HQFREQ2,1,size(GRA,2)).*...
 %                                 [GRA(end+radioloc.LI+1:end,:); GRA(1:radioloc.UI, :)];
                                         
HQFREQtmp=fft(radioloc.h00Nscat,NFFT2);
radioloc.HQFREQscat=[HQFREQtmp(1,NFFT2+radioloc.LI+1:NFFT2),... 
                                              HQFREQtmp(1,2:radioloc.UI+1)].';
radioloc.FREQISI=repmat(radioloc.HQFREQscat,1,size(GRA,2)).*...
                                  [GRA(end+radioloc.LI+1:end,:); GRA(1:radioloc.UI, :)];
                                        
HQFREQtmp=fft(radioloc.h00Nlos,NFFT2);
radioloc.HQFREQLos=[HQFREQtmp(1,NFFT2+radioloc.LI+1:NFFT2), ...
                                             HQFREQtmp(1,2:radioloc.UI+1)].';
 radioloc.ChanDbg=repmat(radioloc.HQFREQLos,1,size(GRA,2)).*...
                                  [GRA(end+radioloc.LI+1:end,:); GRA(1:radioloc.UI, :)];

HQFREQtmp=fft(radioloc.h00Nsinc,NFFT2);                              
radioloc.HQFREQLosC=radioloc.h00t(1)*[HQFREQtmp(1,NFFT2+radioloc.LI+1:NFFT2), ...
                                             HQFREQtmp(1,2:radioloc.UI+1)].';
radioloc.HQFREQLosB=[HQFREQtmp(1,NFFT2+radioloc.LI+1:NFFT2), ...
                                             HQFREQtmp(1,2:radioloc.UI+1)].';
radioloc.ChanDbg2=repmat(radioloc.HQFREQLosC,1,size(GRA,2)).*...
                                  [GRA(end+radioloc.LI+1:end,:); GRA(1:radioloc.UI, :)];
%====================
%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
end

