function [varargout] = RunCazacSys(varargin)
dbstop if error

cazacj=varargin{1};
indexval=varargin{3};
cazacj.symbind= indexval;
cazacj=freq2time(cazacj);
varargout{1}=cazacj;
switch varargin{2}
    case 'PlotAll'
        g2=cazacj.REarr;
        g=g2;
        ML=cazacj.ML;
        MU=cazacj.MU;
        N=cazacj.Nsymb;
        indval=[0:N-1];
        r=zeros(cazacj.NREf,cazacj.Nsymb);
        for rdex=1:cazacj.NREf
            k=0;
            for tau=ML:MU
                k=k+1;
                xax(k)=tau;
                r(rdex, k)=(1/N)*sum(conj(g(rdex,:)).*g2(rdex,[rem(tau+N,N)+1:N,1:rem(tau+N,N)]));
            end
        end
        
        rdex=1;
        
        figure(1)
        subplot(2,1,1), plot(indval,real(g(rdex,:)));
        xlabel('index')
        ylabel('Real Amplitude')
        titleSTR=strcat({'CAZAC Sequence Enabling Spreading,'},{', N='},{num2str(N)});
        title(titleSTR);
        subplot(2,1,2), plot(indval,imag(g(rdex,:)));
        xlabel('index')
        ylabel('Imag Amplitude')
        
        figure(2)
        subplot(2,1,1), plot(indval,real(r(rdex,:)));
        xlabel('index')
        ylabel('Real Amplitude')
        titleSTR=strcat({'CAZAC Sequence Correlation'},{', N='},{num2str(N)});
        title(titleSTR);
        subplot(2,1,2), plot(indval,imag(r(rdex,:)));
        xlabel('index')
        ylabel('Imag Amplitude')
        
        figure(3)
        [xax, yax]=meshgrid(0:cazacj.NREf-1,0:cazacj.Nsymb-1);
        mesh(xax', yax', real(cazacj.REarr));
        
        figure(4)
        mesh(xax', yax', imag(cazacj.REarr));
    case 'NullPlot'
        
        
        
        
    otherwise
        disp('Unknown Option');
end




