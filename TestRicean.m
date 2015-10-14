clear all
fs=7.68e6;
ts=1/fs;
%ts=2.5e-6;
fd=0;
delayval=50e-6;
riceankdB = 8; %dB
riceanklin = 10^(riceankdB/10);
Nit=500;
sval=zeros(Nit,1);

for n=1:Nit
PathDelays=delayval+[0, 0.1, 0.3, 0.5, 0.7, 1.0, 1.3, 15.0, 15.2, 15.7, 17.2, 20.0]*1e-6;
scatterdB=[-8,  -6,  -4,   0,   0,  -4,   -8,   -9,  -10,  -12, -14];
scatterlin=mean(10.^(scatterdB/10));
kdb=10*log10(riceanklin*scatterlin);
AvgPathGaindB=[kdb, scatterdB];
chan = ricianchan(ts,fd, 1, PathDelays, AvgPathGaindB);
chan.NormalizePathGains = 1;
h00N=chan.PathGains/norm(chan.PathGains);
sval(n)=20*log10(abs(h00N(1))/mean(abs(h00N(2:end))));
end

dbg77=1;