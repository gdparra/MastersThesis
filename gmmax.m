function [ output_args ] = gmmax( DataLong)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% MU1 = [4];
% SIGMA1 = [0.2 ];
% MU2 = [-2 ];
% SIGMA2 = [1];
% X = [mvnrnd(MU1,SIGMA1,1000);mvnrnd(MU2,SIGMA2,100)];
X=reshape(DataLong, size(DataLong,1)*size(DataLong,2),1);
%plot(1:length(X), X, '.')
%hold on
options=statset('Display', 'final');
obj = gmdistribution.fit(X,2, 'Options', options);
obj.mu;
obj.PComponents;
obj.NlogL;
properties('gmdistribution');
h = ezplot(@(x)pdf(obj,[x]),[-4e-6 20e-6]);
P = posterior(obj,X);
[a,b]=max(get(h,'Ydata'));
Dval0=get(h,'Xdata');
output_args=Dval0(b);
end

