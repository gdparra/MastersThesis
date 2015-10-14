%abc=FreqDelEst-mean(FreqDelEst);
sigval=std(FreqDelEst)
lll=1;
for kkk=1:length(FreqDelEst)
    if  abs(FreqDelEst(kkk)-mean(FreqDelEst)) < sigval
          abc1(lll)=FreqDelEst(kkk);
          lll=lll+1;
    end
end