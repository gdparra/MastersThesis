function [ output_args ] = histmax( DataLong, BinSize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[a,b]=hist(DataLong);
[e,f]=max(a);
colsel=median(f);

end

