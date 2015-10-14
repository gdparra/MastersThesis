function [ ChanEstOut ] = ChanEst(varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%NumEq=10;
%NumEq=varargin{1};
TrainVec = varargin{1};
zvec = varargin{2};
xsize=length(TrainVec);
ChanSize=varargin{3};
 
TrainMat=zeros(xsize+ChanSize-1,ChanSize);
for k=0:size(TrainMat,2)-1
TrainMat(k+1:k+xsize,k+1)=TrainVec.';
end
 
XMMSE=TrainMat;
difval=length(zvec)-size(XMMSE,1);
if difval <= 0
ChanEstOut=inv(XMMSE'*XMMSE)*XMMSE'*[zvec, zeros(1,-difval)].';
else
ChanEstOut=inv(XMMSE'*XMMSE)*XMMSE'*zvec(1:end-difval).';
end
 
figure(1)
plot(conv(TrainVec.',ChanEstOut))
title('Convolution of ChanEstOut and h(n)');
ylabel('Amplitude');
xlabel('Index');

end

