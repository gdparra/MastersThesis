function [mean, stdev ] = ExampleFunction( x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
           %STAT Interesting statistics.
           n = length(x);
           mean = avg(x,n);
           stdev = sqrt(sum((x-avg(x,n)).^2)/n);
 
 
end

          %-------------------------
           function mean = avg(x,n)
           %AVG subfunction
           mean = sum(x)/n;
           end