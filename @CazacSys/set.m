function cazacsys = set(cazacsys,property,value)

switch property
case 'Bandwidth'
    cazacsys.Band =value;
case 'DeltaF'
    cazacsys.DF =value;
case 'NREf'
    cazacsys.NREF=value;   
case 'Nsymb'
    cazacsys.Nsymb=value;
case 'CPsize'
    cazacsys.CP=value;
case 'ML'
    cazacsys.ML=value;
case 'MU'
    cazacsys.MU=value;
case 'REArr'
    cazacsys.REarr=value;
case 'veclong'
    cazacsys.VecLong=value;
case 'T2SamplePt'
    cazacsys.synchpt=value;
end
