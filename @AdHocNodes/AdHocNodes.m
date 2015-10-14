classdef AdHocNodes < handle
    
    properties
        
        NumRefNodes   =10; %Number of Reference Nodes
        NumBasisNodes =4;   %Number of Basis Nodes
        AreaVal       =10; %Area of the system
        NoiseVarOTDOA =0.01; %Variance add in the range space
        NumTargetval  =3;  %Number of nodes that we wish to location
        InitNumNodes  =16;   %number of initial guesses
        InitMethod    ='locus'; %initial guess method: locus, random, offset
        NumIterations =1;    %number of iterations for converence
        RefnodeModel  ='gpsdenied';
        clusternodes  =1;
        bas           =[];  %  acutal basis vector locations
        target        =[]; % x actual target location values
        mnodes0       =[]; % x node locations without noise
        mnodes        =[]; % x node locations with noise
        Nvec          =[]; % x acutal noise statistics drawn from the distribution
        InitNodes     =[];  % x
        Initoffset    =1.5;  % x
        TDelay        =[];
        TDelayCoarse  =[];
        SynchWind     =0.95;
        SynchWindow   =[]%
        offset        =[];
        zval0         =[];
        Cspeed        =3e8;
        fldhomecombfig1=[];
        SNRvaldB      =[];
        RicianKdB     =[];
        NumberSubCarr =[];
        SegSize       =[];
        NumbSeg       =[];
        nvar_ratio    =[];
        int_ratio     =[];
        smode         =[];
        titleSNR      =[];
    end
    
    methods   %constructor method
        function obj= AdHocNodes(varargin)
            for k=1:nargin
                switch k
                    case 1
                        obj.NumRefNodes   =varargin{1};
                    case 2
                        obj.AreaVal       =varargin{2};
                    case 3
                        obj.NoiseVarOTDOA =varargin{3};
                    case 4
                        obj.NumTargetval  =varargin{4};
                    case 5
                        obj.NumIterations =varargin{5};
                    case 6
                        obj.clusternodes  =varargin{6};
                    case 7
                        obj.InitNumNodes  =varargin{7};
                    case 8
                        obj.InitMethod    =varargin{8};
                    case 9
                        obj.SynchWindow   =varargin{9};
                        % case 10
                        % obj.zval0         =varargin{10};
                    case 10
                        obj.fldhomecombfig1=varargin{10};
                    case 11
                        obj.SNRvaldB      =varargin{11};
                    case 12
                        obj.RicianKdB     =varargin{12};
                    case 13
                        obj.NumberSubCarr =varargin{13};
                    case 14
                        obj.SegSize       =varargin{14};
                    case 15
                        obj.NumbSeg       =varargin{15};
                    case 16
                        obj.nvar_ratio    =varargin{16};
                    case 17
                        obj.int_ratio     =varargin{17};
                    case 18
                        obj.smode         =varargin{18};
                    case 19
                        obj.titleSNR      =varargin{19};
                        
                end
            end
            
            a= sqrt(obj.AreaVal)/2;
            b= sqrt(obj.AreaVal)/2;
            dummyval      ='phold';
            
            obj.bas   =(b/2)*exp(1i*2*pi*[0:obj.NumBasisNodes-1]./obj.NumBasisNodes); %From (b/2) to (b)
            %obj.target=(a + (b-a).*rand(obj.NumTargetval,1))+1i*(a + (b-a).*rand(obj.NumTargetval,1));
            obj.target=(-a + (b-(-a)).*rand(obj.NumTargetval,1))+1i*(-a + (b-(-a)).*rand(obj.NumTargetval,1));
            %obj.target=[20068.1914336508 + 25164.7161064206i;];%[-32187-32187i];
            %obj.target=(a+i*b);
            %obj.target=(a)*rand(obj.NumTargetval,1)+1i*(b)*rand(obj.NumTargetval,1);
            %Plot Location of Basis Nodes
            %figure(1)
            figure('visible','off');
            plot(obj.bas,'d');
            hold on;
            plot(obj.target,'p');
            %             titlefiga=strcat({'SNRdB='},{num2str(obj.SNRvaldB)}, {' Ref.Nodes='},{num2str(obj.NumRefNodes)},...
            %                 {' RicKdB='},{num2str(obj.RicianKdB)}, {' #Sub.Carr='},{num2str(obj.NumberSubCarr)});
            %             titlefigb=strcat({'#Carr/Symb='},{num2str(obj.SegSize)}, {' #Symbs='},{num2str(obj.NumbSeg)},...
            %                 {' Int.Ratio='},{num2str(obj.int_ratio)}, {' VN.Ratio='},{num2str(obj.nvar_ratio)});
            axis([-64374 64374 -64374 64374]) %Mod W46 V2
            ylabel(' Y Meters ');
            xlabel(' X Meters ');
            %title({['SNIRdB=',num2str(obj.SNRvaldB),'   Ref.Nodes=',num2str(obj.NumRefNodes),'   Ricean KdB=',num2str(obj.RicianKdB)];['Numb.Sub.Carriers=',num2str(obj.NumberSubCarr),'   Numb.Carriers/Symb=',num2str(obj.SegSize)];['Numb.Symbs=',num2str(obj.NumbSeg),'  Interference Ratio=',num2str(obj.int_ratio),'  Noise Ratio=',num2str(obj.nvar_ratio)]})
            if strcmp(obj.smode,'AWGN_Fadding_Interference')
                title({['Mode=',obj.smode,'   ',obj.titleSNR,num2str(obj.SNRvaldB),'   Ref.Nodes=',num2str(obj.NumRefNodes),'   Ricean KdB=',num2str(obj.RicianKdB)];['Sub.Carriers=',num2str(obj.NumberSubCarr),'   Sub.Carriers/Symb=',num2str(obj.SegSize),'   Numb.Symbs=',num2str(obj.NumbSeg),];['Interference Ratio=',num2str(obj.int_ratio),'  Noise Ratio=',num2str(obj.nvar_ratio)]})
            else
                title({['Mode=',obj.smode,'   ',obj.titleSNR,num2str(obj.SNRvaldB),'   Ref.Nodes=',num2str(obj.NumRefNodes),'   Ricean KdB=',num2str(obj.RicianKdB)];['Sub.Carriers=',num2str(obj.NumberSubCarr),'   Sub.Carriers/Symb=',num2str(obj.SegSize),'   Numb.Symbs=',num2str(obj.NumbSeg),];['Noise Ratio=',num2str(obj.nvar_ratio)]})
            end
            
            %title('Location of Bas. Nodes(D), Target(Star), Ref.
            %Nodes(T)');% Move Back to SNR
            
            
            %obj.target=(a + (b-a).*0.2+1i*(a + (b-a).*(-0.7)));
            % xval= (a + (b-a).*rand(obj.NumRefNodes,1));
            % yval= (a + (b-a).*rand(obj.NumRefNodes,1));
            % zval=  xval + 1i*yval;
            
            %zval=(0.67*b)*exp(1i*(2*pi*(([0:obj.NumRefNodes-1]./obj.NumRefNodes)+rand(1))));
            %zval=(0.67*b)*exp(1i*(2*pi*(([0:obj.NumRefNodes-1]./obj.NumRefNodes)+ 1/7)));
            zval1=(b)*exp(1i*(2*pi*(([0:obj.NumRefNodes-1]./obj.NumRefNodes)+ 1/7)));
            
            %%%%%%%%%%%%%%%%dithering code
            varZ=(b/4)^2;
            ditherv=sqrt(varZ/2)*randn(size(0:obj.NumRefNodes-1))+...
                    1i*sqrt(varZ/2)*randn(size(0:obj.NumRefNodes-1));
            bas0   =ditherv+zval1; %From (b/2) to (b)
            gamb=2*b;
            Var = real(bas0) <= gamb;
            Dar = real(bas0) > gamb;
            zvalR=real(bas0).*Var+ b*ones(1,obj.NumRefNodes).*Dar;
            Vbr = real(bas0) >= -gamb;
            Dbr = real(bas0) < -gamb;
            zvalR=zvalR.*Vbr- gamb*ones(1,obj.NumRefNodes).*Dbr;
            
            Vai = imag(bas0) <= gamb;
            Dai = imag(bas0) > gamb;
            zvalI=imag(bas0).*Vai+ gamb*ones(1,obj.NumRefNodes).*Dai;
            Vbi = imag(bas0) >= -gamb;
            Dbi = imag(bas0) < -gamb;
            zvalI=zvalI.*Vbi- gamb*ones(1,obj.NumRefNodes).*Dbi;
            zval= zvalR+1i*zvalI;
            %%%%%%%%%%%%%%%%
            
            %zval(1)=[32187+32187i];
            obj.zval0 =zval;
            plot(obj.zval0,'^');
%Test       obj.target=obj.zval0(1)-(19.53125)*17.555; %Test W1 %17.555
%           plot(obj.target,'p'); %Test W1
            hold off;
            saveas(gcf, obj.fldhomecombfig1);
            
            ZMAT=ones(obj.clusternodes, obj.NumRefNodes)*diag(obj.zval0);
            obj.mnodes0=reshape(ZMAT, size(ZMAT,1)*size(ZMAT,2),1);
            obj.TDelayCoarse=floor(((abs(obj.zval0-obj.target))/obj.Cspeed)/(obj.SynchWindow));
            obj.TDelay=(abs(obj.zval0-obj.target))/obj.Cspeed-obj.TDelayCoarse*obj.SynchWindow;
            %obj.TDelay(1)=(abs(obj.zval0(1)-obj.target))/obj.Cspeed
            
            switch obj.InitMethod
                case 'locus'
                    %                      obj.InitNodes=(b/3)*exp(1i*2*pi*[0:obj.InitNumNodes-1]./obj.InitNumNodes);
                    obj.InitNodes=(b/2)*exp(1i*2*pi*[0:obj.InitNumNodes-1]./obj.InitNumNodes);
                case 'random'
                    obj.InitNodes= (a + (b-a).*rand(obj.InitNumNodes,1))+1i*(a + (b-a).*rand(obj.InitNumNodes,1));
                case 'offset'
                    obj.Initoffset=varargin{6}.offset;
                otherwise
                    error('init guess method not supported');
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
    end
end



