function [nvec] = TDelayToXYerr(TRangeEstimates,zval0,clusternodes, targetval,Cspeed, CoarseDelay)
%dbstop in TDelayToXYerr.m at 25

nvec=zeros(length(zval0)*clusternodes,1);
DRange=Cspeed*TRangeEstimates;
Rabs=(targetval-zval0);

for k=0:length(zval0)-1
    if size(DRange,3) > 1
    tmpDRange=reshape(DRange(:,:,k+1),clusternodes,1)+Cspeed*CoarseDelay(k+1);
    else
    tmpDRange=reshape(DRange(:,k+1), clusternodes, 1)+Cspeed*CoarseDelay(k+1);
    end
    
    %nvecprime=(Rabs(k+1)/abs(Rabs(k+1)))*sign(abs(tmpDRange)-abs(Rabs(k+1))).*abs(abs(tmpDRange)-abs(Rabs(k+1)));
    nvecprime=(Rabs(k+1)/abs(Rabs(k+1)))*abs(tmpDRange-abs(Rabs(k+1)));
%     maxval=0;
%     for p=1:length(nvecprime)
%      if abs(nvecprime(p)) > abs(maxval)
%          maxind=p;
%          maxval=nvecprime(p);
%      end
%      zval=ones(1:length(nvecprime), 1);
%      figure(k+10)
%      plot(1:length(nvecprime), abs(nvecprime), 'k',  1:length(nvecprime), zval*abs(Rabs(k+1)), 'r');
%     dbg77=1;
%     end
    nvec(k*clusternodes+1:(k+1)*clusternodes)=nvecprime;
end
dbg77=1;

