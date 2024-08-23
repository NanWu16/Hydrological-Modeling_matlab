function UL = read_underlying(CON)

    AllFilpath.Soiltype_Filepath = "../../Data/Soiltype.dat";
    if ~exist(AllFilpath.Soiltype_Filepath, 'file'), disp('The soil type file is missing'); return; end
    ctemp       = readtable(AllFilpath.Soiltype_Filepath);
    UL.STSWC    = ctemp.SWC;        
    UL.STFC     = ctemp.FC;       
    UL.STWP     = ctemp.WP;        
    UL.STSG     = ctemp.SG;        
    UL.STSHC    = ctemp.SHC_cm_h_;  
    UL.STSF     = ctemp.SF_cm_;     


    AllFilpath.CNtype_Filepath = "../../Data/CN.dat";
    if ~exist(AllFilpath.CNtype_Filepath, 'file'), disp('The CN type file is missing'); return; end
    UL.CN = readtable(AllFilpath.CNtype_Filepath);
    UL.CN = table2array(UL.CN(:,3:6));


    AllFilpath.Landcover_Filepath = "../../Data/Landcover.dat";
    if ~exist(AllFilpath.Landcover_Filepath, 'file'), disp('The landcover file is missing'); return; end
    ctemp       = readtable(AllFilpath.Landcover_Filepath);
    UL.LAI      = table2array(ctemp(:, 3 : 14));    
    UL.MaxLAI   = ctemp.MaxLAI;                     
    UL.CTopH    = ctemp.Canopy_TopHeight;           

    AllFilpath.DEM_Filepath = strcat("../../Data/ASC_file/", CON.catid, '_DEM.asc');
    if ~exist(AllFilpath.DEM_Filepath, 'file'), disp('The dem file is missing'); return; end
    [UL.Ny, UL.Nx, UL.XllCorner, UL.YllCorner, UL.DDem, UL.Nodata] = read_ascinformation(AllFilpath.DEM_Filepath);
    UL.RasterData(:, :, 1) = textread(AllFilpath.DEM_Filepath,'', 'headerlines', 6);

    UL.Dy = 6370.997 * UL.DDem * pi / 180 * 1000;
    UL.Dx = 6370.997 * cos(UL.YllCorner *pi / 180) * UL.DDem * pi / 180 * 1000;   
    UL.Cellsize = int32(UL.DDem / (30 / 3600) * 1000);

    UL.Garea = UL.Dx * UL.Dy/1000000;

    AllFilpath.RiverChannel_Filepath = strcat("../../Data/ASC_file/", CON.catid, '_RiverChannel.asc');
    if ~exist(AllFilpath.RiverChannel_Filepath, 'file'), disp('The river system file is missing'); return; end
    UL.RasterData(:, :, 2) = textread(AllFilpath.RiverChannel_Filepath,'', 'headerlines', 6);

    AllFilpath.TensionWater_Filepath = strcat("../../Data/ASC_file/", CON.catid, '_TensionWaterCapacity.asc');
    if ~exist(AllFilpath.TensionWater_Filepath, 'file'), disp('The tension water file is missing'); return; end
    UL.RasterData(:, :, 3) = textread(AllFilpath.TensionWater_Filepath,'', 'headerlines', 6);


    AllFilpath.FreeWater_Filepath = strcat("../../Data/ASC_file/", CON.catid, '_FreedomWaterCapacity.asc');
    if ~exist(AllFilpath.FreeWater_Filepath, 'file'), disp('The Free water file is missing'); return; end
    UL.RasterData(:, :, 4) = textread(AllFilpath.FreeWater_Filepath,'', 'headerlines', 6);

    AllFilpath.Soil_0_30_Filepath = strcat("../../Data/ASC_file/", CON.catid, '_0-30cmSoilType.asc');
    if ~exist(AllFilpath.Soil_0_30_Filepath, 'file'), disp('The soil type from 0 to 30cm file is missing'); return; end
    UL.RasterData(:, :, 5) = textread(AllFilpath.Soil_0_30_Filepath,'', 'headerlines', 6);

    AllFilpath.Soil_30_100_Filepath = strcat("../../Data/ASC_file/", CON.catid, '_30-100cmSoilType.asc');
    if ~exist(AllFilpath.Soil_30_100_Filepath, 'file'), disp('The soil type from 30 to 100cm file is missing'); return; end
    UL.RasterData(:, :, 6) = textread(AllFilpath.Soil_30_100_Filepath,'', 'headerlines', 6);


    AllFilpath.Humus_Filepath = strcat("../../Data/ASC_file/", CON.catid, '_HumusSoilDepth.asc');
    if ~exist(AllFilpath.Humus_Filepath, 'file'), disp('The Humus thickness file is missing'); return; end
    UL.RasterData(:, :, 7) = textread(AllFilpath.Humus_Filepath,'', 'headerlines', 6);


    AllFilpath.VegType_Filepath = strcat("../../Data/ASC_file/", CON.catid, '_VegetationType.asc');
    if ~exist(AllFilpath.VegType_Filepath, 'file'), disp('The vegetation type file is missing'); return; end
    UL.RasterData(:, :, 8) = textread(AllFilpath.VegType_Filepath,'', 'headerlines', 6);


    AllFilpath.FlowAcc_Filepath = strcat("../../Data/ASC_file/", CON.catid, '_FlowAccumulationArea.asc');
    if ~exist(AllFilpath.FlowAcc_Filepath, 'file'), disp('The flow accumulation file is missing'); return; end
    UL.RasterData(:, :, 9) = textread(AllFilpath.FlowAcc_Filepath,'', 'headerlines', 6);


    AllFilpath.RunoffProportion_Filepath = strcat("../../Data/ASC_file/", CON.catid, '_RunoffDistributionRatio.asc');
    if ~exist(AllFilpath.RunoffProportion_Filepath, 'file'), disp('The runoff propotion file is missing'); return; end
    UL.RasterData(:, :, 10) = textread(AllFilpath.RunoffProportion_Filepath,'', 'headerlines', 6);

    AllFilpath.FlowFirection_Filepath = strcat("../../Data/ASC_file/", CON.catid, '_GridFlowDirection.asc');
    if ~exist(AllFilpath.FlowFirection_Filepath, 'file'), disp('The Flow Direction file is missing'); return; end
    UL.RasterData(:, :, 11) = textread(AllFilpath.FlowFirection_Filepath,'', 'headerlines', 6);


    AllFilpath.VadoseZone_Filepath = strcat("../../Data/ASC_file/", CON.catid, '_VadoseZoneDepth.asc');
    if ~exist(AllFilpath.VadoseZone_Filepath, 'file'), disp('The vadose zone file is missing'); return; end
    UL.RasterData(:, :, 12) = textread(AllFilpath.VadoseZone_Filepath,'', 'headerlines', 6);

    AllFilpath.landuse_Filepath = strcat("../../Data/ASC_file/", CON.catid, '_landuse.asc');
    if ~exist(AllFilpath.landuse_Filepath, 'file'), disp('The landuse file is missing'); return; end
    UL.landuse = textread(AllFilpath.landuse_Filepath,'', 'headerlines', 6);


    AllFilpath.TopographicIndex_Filepath = strcat("../../Data/ASC_file/", CON.catid, '_TopographicIndex.asc');
    if ~exist(AllFilpath.TopographicIndex_Filepath, 'file'), disp('The TopographicIndex file is missing'); return; end
    UL.TopographicIndex = textread(AllFilpath.TopographicIndex_Filepath,'', 'headerlines', 6);

    UL.RasterData(:, :, 4)  = UL.RasterData(:, :, 4)  * CON.FreeWaterCoef;
    UL.RasterData(:, :, 3)  = UL.RasterData(:, :, 3)  * CON.TensionWaterCoef;
    UL.RasterData(:, :, 7)  = UL.RasterData(:, :, 7)  * CON.FreeWaterCoef;
    UL.RasterData(:, :, 12) = UL.RasterData(:, :, 12) * CON.TensionWaterCoef;
    UL.RasterData(:, 178, :) = [];
    UL.TopographicIndex(:, 178) = [];

    AllFilpath.SnowDepth_Filepath = strcat("../../Data/ASC_file/", CON.catid, '_Snowdepth.000');
    if ~exist(AllFilpath.SnowDepth_Filepath, 'file'), disp('The snow depth file is missing'); return; end
    ctemp = readtable(AllFilpath.SnowDepth_Filepath, 'headerlines', 6, 'FileType', 'text');
    UL.GridDs = table2array(ctemp);

    AllFilpath.OutletInfoFilepath = "../../Data/Qobs_Station_of_Subbasin.txt";
    if ~exist(AllFilpath.OutletInfoFilepath, 'file'), disp('The outlet info file is missing'); return; end
    UL.lonlat = textread(AllFilpath.OutletInfoFilepath, '', 'headerlines', 1);


    AllFilpath.CalOrderFilepath = strcat("../../Data/CalSort_01.txt");
    if ~exist(AllFilpath.CalOrderFilepath, 'file'), disp('The calculating order file is missing'); return; end
    ctemp           = textread(AllFilpath.CalOrderFilepath, '', 'headerlines', 1);
    UL.Ncalsort     = size(ctemp, 1);
    UL.SortingOrder = ctemp(:, 2);
    UL.SortingRow   = ctemp(:, 3);
    UL.SortingCol   = ctemp(:, 4);
end