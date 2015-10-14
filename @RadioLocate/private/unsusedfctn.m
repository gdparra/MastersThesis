function [cazac] = freq2time(varargin)

cazac=varargin{1};
indexval=cazac.symbind;
Nfft=cazac.Nfft*cazac.USR;
FSymbol=zeros(Nfft,1);

CP=cazac.CP*cazac.USR;
%TsymbolLong=zeros((Nfft+CP)*cazac.Nsymb,1);
lind=0;
nsymbind=cazac.symbind;


indval= [cazac.Nfft2-4:cazac.Nfft2,1:6];




FSymbol(1:6)=cazac.REarr(6:11, indexval);
FSymbol(end-4:end)=cazac.REarr(1:5,indexval);

tmpvec=ifft(FSymbol(:),Nfft);
Tsymbol=[tmpvec(end-CP+1:end);tmpvec].';

cazac.VecLong=Tsymbol;

dbg77=1;

