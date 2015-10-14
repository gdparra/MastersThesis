function [adhocobj, outdata, outdatalong, outval] = RunAdHocNode(varargin)
dbstop if error
%dbstop in RunAdHocNode.m at 42;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
for k=1:nargin
  switch k         
      case 1
adhocobj=varargin{1};
      case 2
RangeEstimates=mean(varargin{2},1);    
  end 
end


rangemethod='RefNodes';

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch adhocobj.RefnodeModel
    case 'random'
    adhocobj.Nvec    =  (1/sqrt(2))*sqrt(adhocobj.NoiseVarOTDOA).*(randn(adhocobj.NumRefNodes,1)+1i*randn(adhocobj.NumRefNodes,1));
    adhocobj.mnodes  = zval +adhocobj.Nvec;
    case 'cluster'
         adhocobj.mnodes=zeros(adhocobj.NumRefNodes*adhocobj.clusternodes,1);
         ZMAT=ones(adhocobj.clusternodes, adhocobj.NumRefNodes)*diag(zval);
         adhocobj.Nvec=sqrt(adhocobj.NoiseVarOTDOA/2).*randn(size(ZMAT,1)*size(ZMAT,2),1) +...
              1i*sqrt(adhocobj.NoiseVarOTDOA/2).*randn(size(ZMAT,1)*size(ZMAT,2),1);
         adhocobj.mnodes0=reshape(ZMAT, size(ZMAT,1)*size(ZMAT,2),1);
         adhocobj.mnodes =adhocobj.mnodes0+adhocobj.Nvec;
         adhocobj.mnodes =adhocobj.mnodes0;
     case 'gpsdenied'  %RadioLocation%
         TDelayTime=adhocobj.SynchWindow*adhocobj.TDelayCoarse;
         adhocobj.Nvec=TDelayToXYerr(RangeEstimates,adhocobj.zval0,adhocobj.clusternodes,  adhocobj.target, adhocobj.Cspeed, TDelayTime); %XY lOCATION error converted to distance
         adhocobj.mnodes =adhocobj.mnodes0+adhocobj.Nvec.*exp(1i*2*pi*[rand(size(adhocobj.Nvec))]);
    otherwise
        error('not supported');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   switch rangemethod
       case 'RefNodes'
        
           Klong=get(adhocobj, 'RefVec0');      %Nodes in the system (noiseless)
           Kvec=real(Klong(1:end).*conj(Klong(1:end)));  %Euclidean distance squared: Ki
           Kvecp=Kvec(2:end);                            %[K2,K3  ...KM]'
           K1vec=Kvec(1)*ones(size(Kvecp));              %[K1,K1, ...K1]'
           
           dvec=get(adhocobj, 'RefVec');             %noisey referenc location:  Ri
           rvec=get(adhocobj, 'TargetVec');          %target value locations:  (x,y)
           Ri1=abs(dvec(2:end)-rvec(1))-abs(dvec(1)-rvec(1));            %noise range estimates: Ri1  i=2,3...N
           R1 =abs(rvec(1)-dvec(1));

          

%Chan Methods          

           %rvec(1)  %1st target location 
           %Taylor Method
           %%%%%%%%%%presume a guess close to the real solution
           Mloops= get(adhocobj, 'NumInitNodes');
               switch get(adhocobj, 'InitMethod')                
                   case 'offset'
                   u=get(adhocobj, 'InitOffset');  % offset from true value
                   zg=real(rvec(1))-u +1i*(imag(rvec(1))+u);
                   Mloops=1;
                   case 'locus'
                   zg =get(adhocobj, 'InitNodes');
                   case 'random'
                   zg =get(adhocobj, 'InitNodes');
                   otherwise
                   error('init guess method not supported');
               end
               
          pmin=1e6;
          pminind=1;
          for p=1:Mloops
           targN=zg(p);
           xg=real(targN);
           yg=imag(targN);   
           for k=1:get(adhocobj,'NumIterations')
           targN=xg+1i*yg;
           A(:,1)=-2*xg+2*real(Klong(2:end));
           A(:,2)=-2*yg+2*imag(Klong(2:end));
%           CHAN calibration used to verify A matrix calculation
           %delt0a=-(abs(rvec(1)-Klong)).^2+ (abs(real(Klong)-real(rvec(1)))).^2+(abs(imag(Klong)-imag(rvec(1)))).^2;
           %delt0b=                                    -(abs(targN(1)-Klong)).^2+ abs(Klong).^2 -2*real(Klong)       *real(targN)        -2*imag(Klong)       *imag(targN)        +(abs(targN))^2;
           %delt0d=                                    -(abs(rvec(1)-dvec)).^2+ abs(dvec).^2   -2*real(dvec)        *real(targN)        -2*imag(dvec)        *imag(targN)        +(abs(targN))^2;
            delt0d=-((Ri1).^2+2*R1*Ri1+R1^2*ones(size(Ri1)))                     + Kvecp       +(-2*real(Klong(2:end))*real(targN)        -2*imag(Klong(2:end))*imag(targN)        +(abs(targN))^2);
           %delt0d=delt0d(1:end-1);             
           di=abs(Klong(2:end));
           DEL2= (A'*A)\A'*delt0d;
           
           nvalx=xg+DEL2(1);  %estimated target location x
           nvaly=yg+DEL2(2);  %estimated target location y
           %rvec(1): actual location   DEL2: increment   nval: new guess
           outdata(k,:)=[rvec(1), DEL2(1)+1i*DEL2(2), nvalx+1i*nvaly];
           xg=nvalx;
           yg=nvaly;
          
           end

           outdatalong{p}=outdata;
           if abs(outdata(end,2)) < abs(pmin)
           pmin = outdata(end,2);
           pminind=p;
           end

          end
           dbg88=1;
           
       case 'BasisNodes'
           
       otherwise
           error('Not a valid option: RunAdHocNode');
   end
        outval=outdatalong{pminind};
        
     
        
        