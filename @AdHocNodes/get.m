function value = get(ahnobj, property)

switch property
case 'NumRefNodes'
    value = ahnobj.NumRefNodes;
case 'RefVec0'
    value = ahnobj.mnodes0;
case 'RefVec'
    value = ahnobj.mnodes;
case 'NoiseVec'        
    value = ahnobj.Nvec;    
case 'NoiseVar'
    value=ahnobj.NoiseVarOTDOA;
case 'NumBasisNodes'
    value = ahnobj.NumBasisNodes;
case 'BaseVec'
    value = ahnobj.Bas;
case 'TargetVec'
    value = ahnobj.target;
case 'NumTargetNodes'
    value = ahnobj.NumTargetval;
case 'NumInitNodes'
    value = ahnobj.InitNumNodes;
case 'InitMethod'
    value = ahnobj.InitMethod;
case 'InitOffset'
    value = ahnobj.Initoffset;
case 'InitNodes'
    value = ahnobj.InitNodes;
case 'NumIterations'
    value = ahnobj.NumIterations;
case 'TDelay'
    value = ahnobj.TDelay;
case 'cluster'
    value=ahnobj.clusternodes;
case 'CoarseDelay'
    value=ahnobj.TDelayCoarse;
case 'TowerLoc'
    value=ahnobj.zval0;  
end
