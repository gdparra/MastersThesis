function [varargout] = GenRice(riceankdB,scatterdB,PathDelays,ts,fd)

% riceanklin = 10^(riceankdB/10);
% scatterlin=sum(10.^(scatterdB/10));
% kdb=10*log10(riceanklin*scatterlin);
AvgPathGaindB=scatterdB;
chan = rayleighchan(ts,fd, PathDelays, AvgPathGaindB);
chan.NormalizePathGains = 1;
ChannelOut=chan.PathGains/norm(chan.PathGains);
kdB2=10*log10(abs(ChannelOut(1)^2)/sum(abs(ChannelOut(2:end).^2)));
varargout{1}=ChannelOut;
varargout{2}=kdB2;
varargout{3}=ChannelOut(2:end);
