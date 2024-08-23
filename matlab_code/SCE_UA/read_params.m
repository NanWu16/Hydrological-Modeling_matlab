function [PARAM, param_bound] = read_params(CON)
    AllFilpath.LumpedParamFilepath = strcat("../../Data/Lumpara_", CON.catid, ".xlsx");
    if ~exist(AllFilpath.LumpedParamFilepath, 'file'), disp('The Lumped Parameters file is missing'); return; end
    PARAM       = xlsread(AllFilpath.LumpedParamFilepath, "Sheet1");
    param_bound = xlsread(AllFilpath.LumpedParamFilepath, "boundary_m");
end

