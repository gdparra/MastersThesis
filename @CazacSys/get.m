function value = get(cazacsys, property)

switch property
case 'Bandwidth'
    value = cazacsys.Band;
case 'DeltaF'
    value = cazacsys.DF;
case 'NREf'
    value = cazacsys.NREf;
case 'NRETot'        
    value = cazacsys.Nsymb*cazacsys.NREf;    
case 'Nsymb'
    value = cazacsys.Nsymb;
case 'CPsize'
    value = cazacsys.CP;
case 'ML'
    value = cazacsys.ML;
case 'MU'
    value = cazacsys.MU;
case 'REArr'
    value = cazacsys.REarr;
case 'veclong'
     value=cazacsys.VecLong;
case 'Nfft'
     value=cazacsys.Nfft;         
case 'gr'
     value=cazacsys.vec;
case 'Nfft2'
     value=cazacsys.Nfft2;
case 'CP2'
     value=cazacsys.CP2;
case 'Index'
     value=cazacsys.symbind;
case 'T2SamplePt'
     value=cazacsys.synchpt;
case  'vec2'
     value=cazacsys.VecTone;
case 'vecLT';
     value=cazacsys.vecLT;
end
