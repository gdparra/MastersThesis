nvar=1e-5;
x = [1, 3, -1, 7, 2];
h= [1, -2, 0.5];
y = conv(x, h);
v = sqrt(nvar)*randn(size(y));
z=y+v;
chanout=ChanEst(x,z,3);
