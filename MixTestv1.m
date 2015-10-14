MU1 = [4];
SIGMA1 = [0.2 ];
MU2 = [-2 ];
SIGMA2 = [1];
X = [mvnrnd(MU1,SIGMA1,1000);mvnrnd(MU2,SIGMA2,100)];

plot(1:length(X), X, '.')
hold on
options=statset('Display', 'final');
obj = gmdistribution.fit(X,2, 'Options', options);
obj.mu
obj.PComponents
obj.NlogL
properties('gmdistribution')
h = ezplot(@(x)pdf(obj,[x]),[-12 10]);
P = posterior(obj,X);
[a,b]=max(get(h,'Ydata'));
Dval0=get(h,'Xdata');
Dval=Dval0(b);
% delete(h)
% scatter(X(:,1),X(:,2),10,P(:,1),'.')
% hb = colorbar;
% ylabel(hb,'Component 1 Probability')