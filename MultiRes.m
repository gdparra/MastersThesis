function [ varargout ] = MultiRes( xdatain, Tin, To,  TypeMeth)
dbstop if error
ap=size(xdatain);
if ap(1)==1
    xdatain=xdatain.';
end
ta=size(Tin);
if ta(1)==1
    Tin=Tin.';
end

switch TypeMeth
    case 'Linear'
        N=floor(length(xdatain)*(To/Tin));
        xdataout=zeros(N,1);
        
        for n=0:N-1
            xdataout(n+1)=xdatain.'*sinc(n-(1/To)*Tin*[0:M-1].');
        end
    case 'NonLinear'
        N=max(floor(Tin(end)/To),1);
        xdataout=zeros(N,1);
        sdatain=zeros(size(xdatain));
        sdatain(1)=1;
        sdataout=zeros(N,1);
        for n=0:N-1
            xdataout(n+1)=xdatain.'*sinc(n-(Tin/To));
            %sdataout(n+1)=sdatain.'*sinc(n-(1/To)*Tin);
        end
            sdataout=(sinc((0:N-1)-Tin(1)/To)).';
            dbg77=1;
end
varargout{1}=xdataout.';%/norm(xdataout);
varargout{2}=sdataout.';
varargout{3}=Tin(1);
end