%%% Program: GXAJ

clc;
clearvars -except Prec PET Temp

addpath(genpath('./read'))
addpath(genpath('./write'))
addpath(genpath('./initialize'))

AllFilpath.CatchID_Filepath = '../Data/CatID.dat';
if ~exist(AllFilpath.CatchID_Filepath, 'file'), disp('The basin information file is missing'); return; end

[catid, STCD, NoBasin, NoStation] = read_catchment(AllFilpath.CatchID_Filepath);
disp(strcat("Catchment name:", STCD));


AllFilpath.Configure_Filepath = "../Data/Configure.dat";
if ~exist(AllFilpath.Configure_Filepath, 'file'), disp('The global congifure file is missing'); return; end
[IYB, IMB, IDB, IYE, IME, IDE, Nwarmday, Nsm, SCF, FrozenCoef, FreeWaterCoef, TensionWaterCoef, Div, save_state, load_state, Inrunmethod, DT, TimeStepCmp] = read_configure(AllFilpath.Configure_Filepath);

sim_startday = datetime([num2str(IYB),'-',num2str(IMB),'-',num2str(IDB)], 'Inputformat', 'yyyy-MM-dd');
sim_endday = datetime([num2str(IYE),'-',num2str(IME),'-',num2str(IDE)], 'Inputformat', 'yyyy-MM-dd');

data_startday = datetime('2000-01-01', 'InputFormat', 'yyyy-MM-dd');   
loc_start = days(sim_startday - data_startday) + 1;
loc_end = days(sim_endday - data_startday) + 1;
load ../Data/Prec/Precipitation.mat
Prec = Precipitation(:, :, loc_start : loc_end);
clear Precipitation;
load ../Data/PET/PET0.mat
PET = PET0(:, :, loc_start : loc_end);
clear PET0;
load ../Data/Temp/Temperature.mat
Temp = Temperature(:, :, loc_start : loc_end);
clear Temperature;
disp(strcat("Calibration period: ", num2str(IYB), "年", num2str(IMB), "月", num2str(IDB), "日 --- ", num2str(IYE), "年", num2str(IME), "月", num2str(IDE), "日")) % 打印模型运行的起止年月日

[IY1, IM1, ID1, IH1, IY2, IM2, ID2, IH2] = data_length(IYB, IMB, IDB, 0, IYE, IME, IDE, 23);
ntime = int32(length(IY1));
ntime1 = int32(length(IY2));

AllFilpath.Discharge_Filepath = strcat("../Data/QH_", catid, ".xlsx");
if ~exist(AllFilpath.Discharge_Filepath, 'file'), disp('The observation runoff file is missing'); return; end
TotalQobs = zeros(ntime1, NoStation);
TotalQ = xlsread(AllFilpath.Discharge_Filepath);
r = size(TotalQ, 1);
for ii = 1 : ntime1
    for jj = 1 : r
        if TotalQ(jj, 1) == IY2(ii) && TotalQ(jj, 2) == IM2(ii) && TotalQ(jj, 3) == ID2(ii)
            if TimeStepCmp == 1
                for kk = 1 : NoStation
                    TotalQobs(ii, kk) = TotalQ(jj, kk + 3);
                end
            elseif TimeStepCmp == 0
                for kk = 1:NoStation
                    TotalQtmp = 0;
                    for i = jj:(DT + jj - 1)
                        TotalQtmp = TotalQtmp + TotalQ(i, kk + 3);
                    end
                    TotalQobs(ii, kk) = TotalQtmp / DT;
                end
            end
            break;
        end
    end
end
disp(strcat("Total length of observed discharge:", num2str(ntime1)))

AllFilpath.Soiltype_Filepath = "../Data/Soiltype.dat";
if ~exist(AllFilpath.Soiltype_Filepath, 'file'), disp('The soil type file is missing'); return; end
ctemp = readtable(AllFilpath.Soiltype_Filepath);
STSWC = ctemp.SWC;
STFC = ctemp.FC;
STWP = ctemp.WP;
STSG = ctemp.SG;
STSHC = ctemp.SHC_cm_h_;
STSF = ctemp.SF_cm_;

AllFilpath.CNtype_Filepath = "../Data/CN.dat";
if ~exist(AllFilpath.CNtype_Filepath, 'file'), disp('The CN type file is missing'); return; end
CN = readtable(AllFilpath.CNtype_Filepath);
CN = table2array(CN(:,3:6));

AllFilpath.Landcover_Filepath = "../Data/Landcover.dat";
if ~exist(AllFilpath.Landcover_Filepath, 'file'), disp('The landcover file is missing'); return; end
ctemp = readtable(AllFilpath.Landcover_Filepath);
LAI = table2array(ctemp(:, 3 : 14));
MaxLAI = ctemp.MaxLAI;
CTopH = ctemp.Canopy_TopHeight;

AllFilpath.DEM_Filepath = strcat("../Data/ASC_file/", catid, '_DEM.asc');
if ~exist(AllFilpath.DEM_Filepath, 'file'), disp('The dem file is missing'); return; end
[Ny, Nx, XllCorner, YllCorner, DDem, Nodata] = read_ascinformation(AllFilpath.DEM_Filepath);
RasterData(:, :, 1) = textread(AllFilpath.DEM_Filepath,'', 'headerlines', 6);

Dy = 6370.997 * DDem * pi / 180 * 1000;
Dx = 6370.997 * cos(YllCorner *pi / 180) * DDem * pi / 180 * 1000;
Cellsize = int32(DDem / (30 / 3600) * 1000);

Garea = Dx * Dy/1000000;


AllFilpath.RiverChannel_Filepath = strcat("../Data/ASC_file/", catid, '_RiverChannel.asc');
if ~exist(AllFilpath.RiverChannel_Filepath, 'file'), disp('The river system file is missing'); return; end
RasterData(:, :, 2) = textread(AllFilpath.RiverChannel_Filepath,'', 'headerlines', 6);

AllFilpath.TensionWater_Filepath = strcat("../Data/ASC_file/", catid, '_TensionWaterCapacity.asc');
if ~exist(AllFilpath.TensionWater_Filepath, 'file'), disp('The tension water file is missing'); return; end
RasterData(:, :, 3) = textread(AllFilpath.TensionWater_Filepath,'', 'headerlines', 6);

AllFilpath.FreeWater_Filepath = strcat("../Data/ASC_file/", catid, '_FreedomWaterCapacity.asc');
if ~exist(AllFilpath.FreeWater_Filepath, 'file'), disp('The Free water file is missing'); return; end
RasterData(:, :, 4) = textread(AllFilpath.FreeWater_Filepath,'', 'headerlines', 6);

AllFilpath.Soil_0_30_Filepath = strcat("../Data/ASC_file/", catid, '_0-30cmSoilType.asc');
if ~exist(AllFilpath.Soil_0_30_Filepath, 'file'), disp('The soil type from 0 to 30cm file is missing'); return; end
RasterData(:, :, 5) = textread(AllFilpath.Soil_0_30_Filepath,'', 'headerlines', 6);

AllFilpath.Soil_30_100_Filepath = strcat("../Data/ASC_file/", catid, '_30-100cmSoilType.asc');
if ~exist(AllFilpath.Soil_30_100_Filepath, 'file'), disp('The soil type from 30 to 100cm file is missing'); return; end
RasterData(:, :, 6) = textread(AllFilpath.Soil_30_100_Filepath,'', 'headerlines', 6);

AllFilpath.Humus_Filepath = strcat("../Data/ASC_file/", catid, '_HumusSoilDepth.asc');
if ~exist(AllFilpath.Humus_Filepath, 'file'), disp('The Humus thickness file is missing'); return; end
RasterData(:, :, 7) = textread(AllFilpath.Humus_Filepath,'', 'headerlines', 6);

AllFilpath.VegType_Filepath = strcat("../Data/ASC_file/", catid, '_VegetationType.asc');
if ~exist(AllFilpath.VegType_Filepath, 'file'), disp('The vegetation type file is missing'); return; end
RasterData(:, :, 8) = textread(AllFilpath.VegType_Filepath,'', 'headerlines', 6);

AllFilpath.FlowAcc_Filepath = strcat("../Data/ASC_file/", catid, '_FlowAccumulationArea.asc');
if ~exist(AllFilpath.FlowAcc_Filepath, 'file'), disp('The flow accumulation file is missing'); return; end
RasterData(:, :, 9) = textread(AllFilpath.FlowAcc_Filepath,'', 'headerlines', 6);

AllFilpath.RunoffProportion_Filepath = strcat("../Data/ASC_file/", catid, '_RunoffDistributionRatio.asc');
if ~exist(AllFilpath.RunoffProportion_Filepath, 'file'), disp('The runoff propotion file is missing'); return; end
RasterData(:, :, 10) = textread(AllFilpath.RunoffProportion_Filepath,'', 'headerlines', 6);

AllFilpath.FlowFirection_Filepath = strcat("../Data/ASC_file/", catid, '_GridFlowDirection.asc');
if ~exist(AllFilpath.FlowFirection_Filepath, 'file'), disp('The Flow Direction file is missing'); return; end
RasterData(:, :, 11) = textread(AllFilpath.FlowFirection_Filepath,'', 'headerlines', 6);

AllFilpath.VadoseZone_Filepath = strcat("../Data/ASC_file/", catid, '_VadoseZoneDepth.asc');
if ~exist(AllFilpath.VadoseZone_Filepath, 'file'), disp('The vadose zone file is missing'); return; end
RasterData(:, :, 12) = textread(AllFilpath.VadoseZone_Filepath,'', 'headerlines', 6);

AllFilpath.landuse_Filepath = strcat("../Data/ASC_file/", catid, '_landuse.asc');
if ~exist(AllFilpath.landuse_Filepath, 'file'), disp('The landuse file is missing'); return; end
landuse = textread(AllFilpath.landuse_Filepath,'', 'headerlines', 6);

AllFilpath.TopographicIndex_Filepath = strcat("../Data/ASC_file/", catid, '_TopographicIndex.asc');
if ~exist(AllFilpath.TopographicIndex_Filepath, 'file'), disp('The TopographicIndex file is missing'); return; end
TopographicIndex = textread(AllFilpath.TopographicIndex_Filepath,'', 'headerlines', 6);
TopographicIndex = TopographicIndex(:,1:177);

RasterData(:, :, 4) = RasterData(:, :, 4) * FreeWaterCoef;
RasterData(:, :, 3) = RasterData(:, :, 3) * TensionWaterCoef;
RasterData(:, :, 7) = RasterData(:, :, 7) * FreeWaterCoef;
RasterData(:, :, 12) = RasterData(:,:,12) * TensionWaterCoef;

if load_state == 0
    AllFilpath.SnowDepth_Filepath = strcat("../Data/ASC_file/", catid, '_Snowdepth.000');
else
    AllFilpath.SnowDepth_Filepath = strcat("../Data/Initial_state/", catid, '_Snowdepth.asc');
end
if ~exist(AllFilpath.SnowDepth_Filepath, 'file'), disp('The snow depth file is missing'); return; end
ctemp = readtable(AllFilpath.SnowDepth_Filepath, 'headerlines', 6, 'FileType', 'text');
GridDs = table2array(ctemp);

AllFilpath.OutletInfoFilepath = "../Data/Qobs_Station_of_Subbasin.txt";
if ~exist(AllFilpath.OutletInfoFilepath, 'file'), disp('The outlet info file is missing'); return; end
lonlat = textread(AllFilpath.OutletInfoFilepath, '', 'headerlines', 1);

disp('Read lumped parameters')
AllFilpath.LumpedParamFilepath = strcat("../Data/Lumpara_", catid, ".xlsx");
if ~exist(AllFilpath.LumpedParamFilepath, 'file'), disp('The Lumped Parameters file is missing'); return; end
TotalLumPara = xlsread(AllFilpath.LumpedParamFilepath);

GridQSim = zeros(Nx, Ny, ntime1);
initial_GridQSim = ini_zero_2dim(Nx, Ny, 1);
[GridQi, GridQg, GridQs] = ini_zero_2dim(Nx, Ny, 3);
[GridPObs, GridEObs, GridTc, GridTf] = ini_zero_2dim(Nx, Ny, 4);
[GridW, GridS, GridFLC, GridWU, GridWL, GridWD] = ini_zero_2dim(Nx, Ny, 6);
GridRP = zeros(Nx, Ny);
GridCN = zeros(Nx, Ny);
ThitaS= zeros(Nx, Ny);
ThitaF = zeros(Nx, Ny);
ThitaW = zeros(Nx, Ny);
GridHIni = zeros(Nx, Ny, 4);
GWWP = zeros(Nx, Ny);
GWSWC = zeros(Nx, Ny);
GridF = ones(Nx, Ny) * 0.1;
[He, WaterArea, RiverPoint, FlowDirection, DRMOrderNo, DRMWM] = ini_zero_2dim(Nx, Ny, 6);
[DRMSM, DRMMNfC, DRMAlpha, DRMMNfS, NextNoIJ, DRMBmax, DRMSSlope, DRMCSlope, DRMfc] = ini_zero_2dim(Nx, Ny, 9);
[SType030, SType30100, TopIndex, LCover, DRMKg, DRMKi, HumousT, ThickoVZ, GridWM, GridSM, BasBulk] = ini_zero_2dim(Nx, Ny, 11);
TotalPre_R = zeros(ntime1, 5);

for kk = 1 : NoBasin
    tic
    InputBasin = '00';
    if kk < 10
        InputBasin(2:2) = num2str(kk, '%1d');
    else
        InputBasin(1:2) = num2str(kk, '%2d');
    end
    for i = 1 : ntime1
        Qobs(i,1) = TotalQobs(i, 1);
    end

    AllFilpath.CalOrderFilepath = strcat("../Data/CalSort_", InputBasin, ".txt");
    if ~exist(AllFilpath.CalOrderFilepath, 'file'), disp('The calculating order file is missing'); return; end
    ctemp = textread(AllFilpath.CalOrderFilepath, '', 'headerlines', 1);
    Ncalsort = size(ctemp, 1);
    SortingOrder = ctemp(:, 2);
    SortingRow = ctemp(:, 3);
    SortingCol = ctemp(:, 4);

    OC = TotalLumPara(1, 2);
    ROC = TotalLumPara(1, 3);
    KEpC = TotalLumPara(1, 4);
    DeeperC = TotalLumPara(1, 5);
    AlUpper = TotalLumPara(1, 6);
    AlLower = TotalLumPara(1, 7);
    CCg = TotalLumPara(1, 8);
    CCi = TotalLumPara(1, 9);
    CCS = TotalLumPara(1, 10);
    LagTime = TotalLumPara(1, 11);
    CCS_HM = TotalLumPara(1, 12);
    LagTime_HM = TotalLumPara(1, 13);
    MKch = TotalLumPara(1, 14);
    MKs = TotalLumPara(1, 15);
    MKi = TotalLumPara(1, 16);
    MKg = TotalLumPara(1, 17);
    MXch = TotalLumPara(1, 18);
    MXs = TotalLumPara(1, 19);
    MXi = TotalLumPara(1, 20);
    MXg = TotalLumPara(1, 21);
    AlDeeper = 1 - AlUpper - AlLower;

    DRMGCNo = 0;
    SumKgKi = 0;
    for j = 1 : Ncalsort
        ii = SortingRow(j);
        jj = SortingCol(j);

        He(ii, jj) = RasterData(ii, jj, 1);
        RiverPoint(ii, jj) = RasterData(ii, jj, 2);
        DRMWM(ii, jj) = RasterData(ii, jj, 3);
        DRMSM(ii, jj) = RasterData(ii, jj, 4);
        SType030(ii, jj) = RasterData(ii, jj, 5) + 1;
        SType30100(ii, jj) = RasterData(ii, jj, 6) + 1;
        HumousT(ii, jj) = RasterData(ii, jj, 7);
        LCover(ii, jj) = RasterData(ii, jj, 8) + 1;
        DRMfc(ii, jj) = RasterData(ii, jj, 10);
        FlowDirection(ii, jj) = RasterData(ii, jj, 11);
        ThickoVZ(ii, jj) = RasterData(ii, jj, 12);
        if (He(ii, jj) ~= Nodata)
            DRMGCNo = DRMGCNo + 1;
            NextNoIJ(ii, jj) = DRMGCNo;
            if (HumousT(ii, jj) <= 300)
                ThitaS(ii, jj) = STSWC(SType030(ii, jj));
                ThitaF(ii, jj) = STFC(SType030(ii, jj));
                ThitaW(ii, jj) = STWP(SType030(ii, jj));
            else
                ThitaS(ii, jj) = STSWC(SType030(ii, jj)) * (300 / HumousT(ii, jj)) + STSWC(SType30100(ii, jj)) * (1 - 300 / HumousT(ii, jj));
                ThitaF(ii, jj) = STFC(SType030(ii, jj)) * (300 / HumousT(ii, jj)) + STFC(SType30100(ii, jj)) * (1 - 300 / HumousT(ii, jj));
                ThitaW(ii, jj) = STWP(SType030(ii, jj)) * (300 / HumousT(ii, jj)) + STWP(SType30100(ii, jj)) * (1 - 300 / HumousT(ii, jj));
            end
            DRMKi(ii, jj) = ((ThitaF(ii, jj) / ThitaS(ii, jj)) ^ OC) / (1 + ROC / (1 + 2 * (1 - ThitaW(ii, jj))));
            DRMKg(ii, jj) = (ThitaF(ii, jj) / ThitaS(ii, jj)) ^ OC - DRMKi(ii, jj);

            if isnan(DRMKi(ii, jj))
                DRMKi(ii, jj) =0.5;
            end
            if isnan(DRMKg(ii, jj))
                DRMKg(ii, jj) = 0.2;
            end
            KgKi_tmp = (ThitaF(ii, jj) / ThitaS(ii, jj)) ^ OC;
            if isnan(KgKi_tmp)
                KgKi_tmp = 0.7;
            end
            SumKgKi = SumKgKi + KgKi_tmp;
        end
    end
    KgKi = SumKgKi / DRMGCNo;
    disp(strcat("Ki + Kg=:", num2str(KgKi)));

    JMonth = IM2(1);
    ETKcbmin = 0.175;

    if load_state == 0
        for j = 1 : DRMGCNo
            ii = SortingRow(j);
            jj = SortingCol(j);
            ETKcb = 1.07 * (1 - exp(-0.84 * LAI(LCover(ii, jj), JMonth)));
            ETKcbmax = 1.07 * (1 - exp(-0.84 * MaxLAI(LCover(ii, jj)))) + 0.05;
            if (ETKcb - ETKcbmin < 0)
                Flc = 0;
            else
                Flc = ((ETKcb - ETKcbmin) / (ETKcbmax - ETKcbmin)) ^ (1 + 0.5 * CTopH(LCover(ii, jj)));
            end

            GridFLC(ii, jj) = Flc;

            Dp = DRMfc(ii, jj);
            GWM = DRMWM(ii, jj);
            GSM = DRMSM(ii, jj);
            Kg = DRMKg(ii, jj);
            Ki = DRMKi(ii, jj);
            ZUpper = AlUpper * ThickoVZ(ii, jj);
            ZLower = AlLower * ThickoVZ(ii, jj);
            ZDeeper = AlDeeper * ThickoVZ(ii, jj);

            if (ZUpper > 300)
                GWUM = (STFC(SType030(ii, jj)) - STWP(SType030(ii, jj))) * 300 + (STFC(SType30100(ii, jj)) - STWP(SType30100(ii, jj))) * (ZUpper - 300);
                GWLM = (STFC(SType30100(ii, jj)) - STWP(SType30100(ii, jj))) * ZLower;
            else
                GWUM = (STFC(SType030(ii, jj)) - STWP(SType030(ii, jj))) * ZUpper;

                if (ZUpper + ZLower > 300.0)
                    GWLM = (STFC(SType030(ii, jj)) - STWP(SType030(ii, jj))) * (300 - ZUpper) + (STFC(SType30100(ii, jj)) - STWP(SType30100(ii, jj))) * (ZLower - 300 + ZUpper);
                else
                    GWLM = (STFC(SType030(ii, jj)) - STWP(SType030(ii, jj))) * ZLower;
                end
            end

            if (ThickoVZ(ii, jj) > 300)
                GWWP(ii, jj) = STWP(SType030(ii, jj)) * 300 + STWP(SType30100(ii, jj)) * (ThickoVZ(ii, jj) - 300);
            else
                GWWP(ii, jj) = STWP(SType030(ii, jj)) * ThickoVZ(ii, jj);
            end

            if (ThickoVZ(ii, jj) > 300)
                GWSWC(ii, jj) = (STSWC(SType030(ii, jj)) - STWP(SType030(ii, jj))) * 300 + (STSWC(SType30100(ii, jj)) - STWP(SType30100(ii, jj))) * (ThickoVZ(ii, jj) - 300);
            else
                GWSWC(ii, jj) = (STSWC(SType030(ii, jj)) - STWP(SType030(ii, jj))) * ThickoVZ(ii, jj);
            end

            GWDM = GWM - GWUM - GWLM;
            GridWU(ii, jj) = GWUM / 10 * Nsm;
            GridWL(ii, jj) = GWLM / 10 * Nsm;
            GridWD(ii, jj) = GWDM / 10 * Nsm;
            GridW(ii, jj) = GridWU(ii, jj) + GridWL(ii, jj) + GridWD(ii, jj);
            GridS(ii, jj) = DRMSM(ii, jj) / 10 * Nsm;

            if (Qobs(1) > 0)
                GridQi(ii, jj) = Qobs(1) / double(DRMGCNo) / 2;
                GridQg(ii, jj) = Qobs(1) / double(DRMGCNo) / 2;
                GridQs(ii, jj) = Qobs(1) / double(DRMGCNo) / 2;
            else
                GridQi(ii, jj) = 0;
                GridQg(ii, jj) = 0;
                GridQs(ii, jj) = 0;
            end
            GridHIni(ii, jj, 1) = GridW(ii, jj);
            GridHIni(ii, jj, 2) = GridS(ii, jj);
            GridHIni(ii, jj, 3) = GridWU(ii, jj);
            GridHIni(ii, jj, 4) = GridWL(ii, jj);
        end
    else
        AllFilpath.GridW_Filepath = "../Data/Initial_state/GridW.asc";
        AllFilpath.GridS_Filepath = "../Data/Initial_state/GridS.asc";
        AllFilpath.GridWU_Filepath = "../Data/Initial_state/GridWU.asc";
        AllFilpath.GridWL_Filepath = "../Data/Initial_state/GridWL.asc";
        AllFilpath.GridQi_Filepath = "../Data/Initial_state/GridQi.asc";
        AllFilpath.GridQg_Filepath = "../Data/Initial_state/GridQg.asc";
        AllFilpath.GridQs_Filepath = "../Data/Initial_state/GridQs.asc";

        if ~exist(AllFilpath.GridW_Filepath, 'file'), disp('The GridW file is missing'); return; end
        if ~exist(AllFilpath.GridS_Filepath, 'file'), disp('The GridS file is missing'); return; end
        if ~exist(AllFilpath.GridWU_Filepath, 'file'), disp('The GridWU file is missing'); return; end
        if ~exist(AllFilpath.GridWL_Filepath, 'file'), disp('The GridWL file is missing'); return; end
        if ~exist(AllFilpath.GridQi_Filepath, 'file'), disp('The GridQi file is missing'); return; end
        if ~exist(AllFilpath.GridQg_Filepath, 'file'), disp('The GridQg file is missing'); return; end
        if ~exist(AllFilpath.GridQs_Filepath, 'file'), disp('The GridQs file is missing'); return; end
        ctemp = readtable(AllFilpath.GridW_Filepath, 'headerlines', 6, 'FileType', 'text');
        GridHIni(:, :, 1) = table2array(ctemp);
        ctemp = readtable(AllFilpath.GridS_Filepath, 'headerlines', 6, 'FileType', 'text');
        GridHIni(:, :, 2) = table2array(ctemp);
        ctemp = readtable(AllFilpath.GridWU_Filepath, 'headerlines', 6, 'FileType', 'text');
        GridHIni(:, :, 3) = table2array(ctemp);
        ctemp = readtable(AllFilpath.GridWL_Filepath, 'headerlines', 6, 'FileType', 'text');
        GridHIni(:, :, 4) = table2array(ctemp);
        ctemp = readtable(AllFilpath.GridQi_Filepath, 'headerlines', 6, 'FileType', 'text');
        GridQi = table2array(ctemp);
        ctemp = readtable(AllFilpath.GridQg_Filepath, 'headerlines', 6, 'FileType', 'text');
        GridQg = table2array(ctemp);
        ctemp = readtable(AllFilpath.GridQs_Filepath, 'headerlines', 6, 'FileType', 'text');
        GridQs = table2array(ctemp);
    end

    TSteps = ntime1;
    GridW0 = GridW;

    if load_state == 0
        if Qobs(1) > 0
            QSim = Qobs(1) * ones(DRMGCNo, LagTime);
        else
            QSim = zeros(DRMGCNo, LagTime);
        end
    else
        AllFilpath.SimGridQFilepath = "../Initial_state/SimGridQ.asc";
        if ~exist(AllFilpath.GridW_Filepath, 'file'), disp('The initial simulated runoff file is missing'); return; end
        ctemp = readtable(AllFilpath.GridW_Filepath, 'headerlines', 6, 'FileType', 'text');
        initial_GridQSim = table2array(ctemp);
        for ii = 1:DRMGCNo
            for i = 1 : LagTime
                QSim(ii, i) = initial_GridQSim(SortingRow(ii), SortingCol(jj));
            end
        end
    end

    for j = 1 : DRMGCNo
        ii = SortingRow(j);
        jj = SortingCol(j);
        if landuse(ii, jj) == 0
            landuse(ii, jj) = 4;
        end
        GridCN(ii, jj) = CN(landuse(ii, jj),STSG(SType030(ii, jj)));
        if GridCN(ii, jj) < 70
            GridRP(ii, jj) = 1;
        else
            GridRP(ii, jj) = 2;
        end
        if TopographicIndex(ii, jj)>15
            GridRP(ii, jj) = 1;
        elseif TopographicIndex(ii, jj) < 5
            GridRP(ii, jj) = 2;
        end
    end

    DRMPO = zeros(DRMGCNo, 1); DRMET = zeros(DRMGCNo, 1); Dis = 0;
    Pcum = zeros(Nx, Ny); Icum = zeros(2, DRMGCNo);
    DRMIca = zeros(DRMGCNo, 1); DRMIch = zeros(DRMGCNo, 1);
    SimQ = zeros(TSteps+1, 1);
    RFC = zeros(Nx, Ny);
    Pnet = zeros(DRMGCNo, 1);
    InflowQs = zeros(2, DRMGCNo); OutflowQs = zeros(2, DRMGCNo);
    InflowQi = zeros(2, DRMGCNo); OutflowQi = zeros(2, DRMGCNo);
    InflowQg = zeros(2, DRMGCNo); OutflowQg = zeros(2, DRMGCNo);
    InflowQch = zeros(2, DRMGCNo); OutflowQch = zeros(2, DRMGCNo);
    NoWM = zeros(TSteps, 1); NoSM = zeros(TSteps, 1); Qout = zeros(TSteps, 1);
    GridQ = zeros(DRMGCNo, 1);
    SumQs = zeros(DRMGCNo, 1); SumQi = zeros(DRMGCNo, 1); SumQg = zeros(DRMGCNo, 1); SumQch = zeros(DRMGCNo, 1);

    SumDRMIca  =  0;
    Ct = Garea / DT / 3.6;
    Cg = CCg ^ (DT / 24.);
    Ci = CCi ^ (DT / 24.);
    GridWM = Nodata * ones(Nx, Ny);
    GridSM = Nodata * ones(Nx, Ny);
    GridW = GridHIni(:, :, 1);
    GridS = GridHIni(:, :, 2);
    GridWU = GridHIni(:, :, 3);
    GridWL = GridHIni(:, :, 4);
    GridD1 = 0 * ones(Nx, Ny);
    GridWi = GridDs;
    GridWq = 0 * ones(Nx, Ny);
    GridS1 = 0 * ones(Nx, Ny);
    GridATI1 = 0 * ones(Nx, Ny);
    GridDensity_x1 = 0.15 * ones(Nx, Ny);
    GridTs = 0 * ones(Nx, Ny);
    GridTx = 0 * ones(Nx, Ny);
    FIa = 0 * ones(Nx, Ny);
    SFD = 0 * ones(Nx, Ny);
    SFDi = 0 * ones(Nx, Ny);
    SFDg = 0 * ones(Nx, Ny);
    GS0 = 0 * ones(Nx, Ny);
    GWU0 = 0 * ones(Nx, Ny);
    GWL0 = 0 * ones(Nx, Ny);
    GWD0 = 0 * ones(Nx, Ny);
    GridTn = 0 * ones(Nx, Ny);
    E_act = 0 * ones(Nx, Ny);

    for i = 1 : TSteps
        NoWM(i) = 0;
        NoSM(i) = 0;

        LumPara = TotalLumPara(IM2(i), 2 : 30);
        OC = LumPara(1);
        ROC = LumPara(2);
        KEpC = LumPara(3);
        DeeperC = LumPara(4);
        AlUpper = LumPara(5);
        AlLower = LumPara(6);
        CCg = LumPara(7);
        CCi = LumPara(8);
        CCS = LumPara(9);
        LagTime = LumPara(10);
        CCS_HM = LumPara(11);
        LagTime_HM = LumPara(12);
        MKch = LumPara(13);
        MKs = LumPara(14);
        MKi = LumPara(15);
        MKg = LumPara(16);
        MXch = LumPara(17);
        MXs = LumPara(18);
        MXi = LumPara(19);
        MXg = LumPara(20);
        UADJ = LumPara(21);
        MBASE = LumPara(22);
        MFMAX = LumPara(23);
        MFMIN = LumPara(24);
        TIPM = LumPara(25);
        NMF = LumPara(26);
        PLWHC = LumPara(27);
        DAYGM = LumPara(28);
        R1 = LumPara(29);
        AlDeeper = 1 - AlUpper - AlLower;

        Jday = IM2(i);
        if Jday ~= JMonth
            JMonth = Jday;
            for jNow = 1:DRMGCNo
                ii = SortingRow(jNow);
                jj = SortingCol(jNow);
                ETKcb = 1.07 * (1 - exp(-0.84 * LAI(LCover(ii, jj), JMonth)));
                ETKcbmax = 1.07 * (1 - exp(-0.84 * MaxLAI(LCover(ii, jj)))) + 0.05;
                if ETKcb - ETKcbmin < 0
                    Flc = 0;
                else
                    Flc = ((ETKcb - ETKcbmin) / (ETKcbmax - ETKcbmin)) ^ (1 + 0.5 * CTopH(LCover(ii, jj)));
                end
                GridFLC(ii, jj) = Flc;
            end
        end

        Inputdate1 = '00000000';
        Inputdate1(1:4) = num2str(IY2(i));
        if IM2(i) >= 10
            Inputdate1(5:6) = num2str(IM2(i));
        else
            Inputdate1(6:6) = num2str(IM2(i));
        end

        if ID2(i) >= 10
            Inputdate1(7:8) = num2str(ID2(i));
        else
            Inputdate1(8:8) = num2str(ID2(i));
        end

        GridPObs = Prec(:, :, i);
        GridEObs = PET(:, :, i);
        GridTObs = Temp(:, :, i);
        GridTf = GridTObs * 1.8 + 32;

        if load_state == 0
            SumQs = zeros(DRMGCNo, 1);
            SumQi = zeros(DRMGCNo, 1);
            SumQg = zeros(DRMGCNo, 1);
            SumQch = zeros(DRMGCNo, 1);
        else
            if load_state == 1
                AllFilpath.SumQFilepath = "../Data/Initial_state/SumQ.txt";
                ctemp = textread(AllFilpath.SumQFilepath, '', 'headerlines', 0);
                SumQch = ctemp(:, 1); SumQs = ctemp(:, 2);
                SumQi = ctemp(:, 3); SumQg = ctemp(:, 4);
            end
        end
        TotalPre_rain = 0;
        TotalPre_snow = 0;
        Total_GRs = 0;
        Total_GRi = 0;
        Total_GRg = 0;
        for k = 1 : DRMGCNo

            Lf = 80;
            ci = 0.05;
            ii = SortingRow(k);
            jj = SortingCol(k);
            if GridTf(ii, jj) <= 28.0
                fs = 1;
            else
                fs = exp(-(GridTf(ii,jj)/8-3.5)^5);
            end
            if GridTf(ii, jj) >= 40.0
                fs = 0;
            end
            Pn = GridPObs(ii,jj) * fs * SCF ;
            Wi = GridWi(ii,jj);
            Wi = Wi + Pn;
            fr = 1 - fs;

            if GridTObs(ii,jj) < -15
                Density = 0.05;
            elseif GridTObs(ii,jj) > 0
                Density = 0.05 + 0.0017 * (GridTObs(ii,jj))^1.5;
            else
                Density = 0.15;
            end
            Hn = 0.1 * Pn / Density;

            if GridTObs(ii,jj) < 0
                Tn = GridTObs(ii,jj);
            else
                Tn = 0;
            end
            if Tn < 0
                Dp = -(Tn * Pn)/(Lf / ci);
            else
                Dp = 0;
            end


            if Wi > 0

                start_date_0321 = datenum(IY2(i), 3, 21);
                end_date_0321 = datenum(IY2(i), IM2(i), ID2(i));
                N = end_date_0321 - start_date_0321;
                Sv = 0.5 * sin((N * 2 *  pi)/ 366) + 0.5;
                Mf = (DT/6) * (Sv * (MFMAX - MFMIN) + MFMIN);
                esat = 2.7489 * 10^8 * exp(-4278.63/(GridTObs(ii,jj) + 242.792));
                Pa = 33.86 * (29.9 - 0.335 * He(ii, jj) + 0.00022 * He(ii, jj)^2.4);
                if GridTObs(ii,jj) > 0
                    Tr = GridTObs(ii,jj);
                else
                    Tr = 0;
                end

                if GridPObs(ii,jj) * fr > 6
                    Mr = 6.12 * 10^(-10) * DT * ((GridTObs(ii,jj) + 273)^4-273^4) + 0.0125 * GridPObs(ii,jj) * fr * Tr + 8.5 * UADJ * (DT/6) * ((0.9 * esat - 6.11) + 0.00057 * Pa * GridTObs(ii,jj));
                else
                    Mr = 0;
                end
                if Mr > Wi
                    Mr = Wi;
                end

                if GridTObs(ii,jj) >= MBASE && GridPObs(ii,jj) * fr <= 6
                    Mnr = Mf * (GridTObs(ii,jj) - MBASE) + 0.0125 * GridPObs(ii,jj) * fr * Tr ;
                else
                    Mnr = 0;
                end
                if Mnr > Wi
                    Mnr = Wi;
                end

                if GridTObs(ii,jj) < 0
                    Tsur = GridTObs(ii,jj);
                else
                    Tsur = 0;
                end

                TIPMt = 1 - (1 - TIPM)^(DT/6);
                ATI1 = GridATI1(ii,jj);
                ATI2 = ATI1 + TIPMt * (GridTObs(ii,jj) - ATI1);
                if  ATI2 > 0
                    ATI2 = 0;
                end
                if Pn > 1.5 * DT
                    ATI2 = Tn;
                end
                if ATI2 > Tsur
                    Dt = NMF * (DT/6) * (Mf/MFMAX) * (ATI2 - Tsur);
                else
                    Dt = 0;
                end
                GridATI1(ii,jj) = ATI2;

                D1 = GridD1(ii,jj);
                Wq = GridWq(ii,jj);
                Wix = Wi;
                Wqt = Wq;
                Qw = Mr + Mnr + GridPObs(ii,jj) * fr;
                D2 = D1 + Dp + Dt;

                if D2 > Qw
                    Wi = Wi - Mr - Mnr + Qw ;
                    Wqx = PLWHC * Wi;
                    Wq = Wq;
                    D2 = D2 - Qw;
                    Qf = Qw;
                    E = 0;
                else
                    Wi = Wi - Mr - Mnr + D2;
                    Wqx = PLWHC * Wi;
                    Qf = D2;
                    Wq =  Wq + Qw -D2;
                    D2 = 0;
                    if Wq > Wqx
                        E = Wq + Qw -D2 - Wqx;
                        Wq = Wqx;
                        Ts = 0;
                    else
                        E = 0;
                    end
                end
                GridD1(ii,jj) = D2;
                GridWq(ii,jj) = Wq;
                if Wi < 0
                    Wi = 0;
                end

                S1 = GridS1(ii,jj);
                Omr = (S1 + E) * R1;
                S2 = S1 + E - Omr;
                Omr = E;
                GridS1(ii,jj) = S2;

                Mg = DAYGM * (DT/24);
                if Wi == 0
                    Og = 0;
                else
                    if Mg > Wi
                        Mg = Wi;
                    end
                    Og = Mg + (Mg / Wi) * Wq;
                end
                Wi = Wi - Og;

                if Wi < 0
                    Wi = 0;
                end
                GridWi(ii,jj) = Wi;
                Os = Omr + Og;
                if Wi < Os
                    Os = Wi;
                end
                Pr = 0;
                TotalPre_rain = TotalPre_rain + Pr;
                TotalPre_snow = TotalPre_snow + Os;
            else
                Os = 0;
                Pr = GridPObs(ii,jj) * fr;
                TotalPre_rain = TotalPre_rain + Pr;
                TotalPre_snow = TotalPre_snow + Os;
            end

            DRMPO(NextNoIJ(ii, jj)) = Pr;

            DRMET(NextNoIJ(ii, jj)) = KEpC * GridEObs(ii, jj);
            GLAI = LAI(LCover(ii, jj), JMonth);
            Scmax = 0.935 + 0.498 * GLAI - 0.00575 * GLAI^2;
            Cp = GridFLC(ii, jj);
            Cvd = 0.046 * GLAI;
            Pcum(ii, jj) = Pcum(ii, jj) + DRMPO(NextNoIJ(ii, jj));
            Icum(2, NextNoIJ(ii, jj)) = Cp * Scmax * (1 - exp((-Cvd) * Pcum(ii, jj) / Scmax));
            DRMIca(NextNoIJ(ii, jj)) = Icum(2, NextNoIJ(ii, jj)) - Icum(1, NextNoIJ(ii, jj));

            if DRMET(NextNoIJ(ii, jj)) >= DRMIca(NextNoIJ(ii, jj))
                DRMET(NextNoIJ(ii, jj)) = DRMET(NextNoIJ(ii, jj)) - DRMIca(NextNoIJ(ii, jj));
                DRMIca(NextNoIJ(ii, jj)) = 0.0;
                if DRMET(NextNoIJ(ii, jj)) >= SumDRMIca
                    DRMET(NextNoIJ(ii, jj)) = DRMET(NextNoIJ(ii, jj)) - SumDRMIca;
                    SumDRMIca = 0.0;
                else
                    SumDRMIca = SumDRMIca - DRMET(NextNoIJ(ii, jj));
                    DRMET(NextNoIJ(ii, jj)) = 0.0;
                end
            else
                DRMIca(NextNoIJ(ii, jj)) = DRMIca(NextNoIJ(ii, jj)) - DRMET(NextNoIJ(ii, jj));
                SumDRMIca = SumDRMIca + DRMIca(NextNoIJ(ii, jj));
                DRMET(NextNoIJ(ii, jj)) = 0.0;
            end


            Pnet(NextNoIJ(ii, jj)) = DRMPO(NextNoIJ(ii, jj)) - (Icum(2, NextNoIJ(ii, jj)) - Icum(1, NextNoIJ(ii, jj))) - DRMET(NextNoIJ(ii, jj));
            Icum(1, NextNoIJ(ii, jj)) = Icum(2, NextNoIJ(ii, jj));

            Pe = Pnet(NextNoIJ(ii, jj)) + Os;
            Ek = DRMET(NextNoIJ(ii, jj));
            GW = GridW(ii, jj);
            GS = GridS(ii, jj);
            GWU = GridWU(ii, jj);
            GWL = GridWL(ii, jj);
            GWD = GW - GWU - GWL;
            GWM = DRMWM(ii, jj);
            GSM = DRMSM(ii, jj);

            F = GridF(ii, jj);
            if isnan(F)==1
                F = 0.1;
            end
            Kg = DRMKg(ii, jj);
            Ki = DRMKi(ii, jj);
            KKi = (1. - ((1 - (Kg + Ki)) ^ (DT / 24))) / (1 + Kg / Ki);
            KKg = KKi * Kg / Ki;
            Kg = KKg;
            Ki = KKi;
            Ks = STSHC(SType030(ii, jj)) * 10*24;
            Sf = STSF(SType030(ii, jj)) * 10;
            f = Ks * (1 + (GWSWC(ii, jj) / ThickoVZ(ii, jj) - GridW0(ii, jj) / ThickoVZ(ii, jj)) * Sf / F);

            if isnan(Kg)
                Kg = 0.2;
            end
            if isnan(Ki)
                Ki = 0.5;
            end
            Grfc = RFC(ii, jj);
            Dp = DRMfc(ii, jj);
            ZUpper = AlUpper * ThickoVZ(ii, jj);
            ZLower = AlLower * ThickoVZ(ii, jj);
            ZDeeper = AlDeeper * ThickoVZ(ii, jj);
            RP = GridRP(ii, jj);

            if ZUpper > 300
                GWUM = (STFC(SType030(ii, jj)) - STWP(SType030(ii, jj))) * 300. + (STFC(SType30100(ii, jj)) - STWP(SType30100(ii, jj))) * (ZUpper - 300.);
                GWLM = (STFC(SType30100(ii, jj)) - STWP(SType30100(ii, jj))) * ZLower;
            else
                GWUM = (STFC(SType030(ii, jj)) - STWP(SType030(ii, jj))) * ZUpper;
                if ZUpper + ZLower > 300.
                    GWLM = (STFC(SType030(ii, jj)) - STWP(SType030(ii, jj))) * (300. - ZUpper) + (STFC(SType30100(ii, jj)) - STWP(SType30100(ii, jj))) * (ZLower - 300. + ZUpper);
                else
                    GWLM = (STFC(SType030(ii, jj)) - STWP(SType030(ii, jj))) * ZLower;
                end
            end
            GWDM = GWM - GWUM - GWLM;

            if Pe <= 0

                if (GS + Pe > 0.0)
                    GS = GS + Pe;
                    GEU = Ek;
                    GEL = 0;
                    GED = 0;
                    GE = Ek;

                else
                    Pe = GS + Pe;
                    GS = 0;

                    if (GWU + Pe < 0.0)
                        GEU = GWU + Ek + Pe;
                        GWU = 0.0;
                        if (GWLM < 0)
                            GWLM = 0.01;
                        end

                        GEL = (Ek - GEU) * GWL / GWLM;
                        if isnan(GEL)
                            GEL = 0;
                        end
                        if (GWL < DeeperC * GWLM)
                            GEL = DeeperC * (Ek - GEU);
                        end
                        if ((GWL - GEL) < 0.0)
                            GED = GEL - GWL;
                            GEL = GWL;
                            GWL = 0.0;
                            GWD = GWD - GED;
                        else
                            GED = 0.0;
                            GWL = GWL - GEL;
                        end
                    else
                        GEU = Ek;
                        GEL = 0.0;
                        GED = 0.0;
                        GWU = GWU + Pe;
                    end
                    GE = GEU + GEL + GED;
                    GW = GWU + GWL + GWD;
                    if (GW < 0)
                        GEU = 0.0;
                        GEL = 0.0;
                        GED = GridW(ii, jj);
                        GE = GED;
                        GWU = 0.0;
                        GWL = 0.0;
                        GWD = 0.0;
                        GW = 0.0;
                        RP =2;
                    end
                end

                SumDRMIca = Ek - GE;
                if  (SumDRMIca < 0.0)
                    SumDRMIca = 0.0;
                end
                if GW < GWM
                    Grfc = 0.0;
                else
                    Grfc = 1.0;
                end

                if  GW > GWM
                    RP =1;
                end

                GRs = 0.0;
                GRg = GS * Kg;
                GRi = GS * Ki;
                GS = GS * (1.0 - Kg - Ki);

            elseif Pe > 0
                GEU = Ek;
                GEL = 0;
                GED = 0;
                GE = Ek;

                if GW < GWM
                    if Pe > f
                        RP =2;
                    end
                else
                    RP =1;
                end
                if RP == 2
                    if Pe > f
                        GRs = Pe - f;
                    else
                        GRs = 0;
                    end
                    GRg = GS * Kg;
                    GRi = GS * Ki;
                else

                    nd = ceil(Pe / Div);
                    PPe = zeros(1, nd);
                    for m = 1 : nd - 1
                        PPe(m) = Div;
                    end

                    PPe(nd) = Pe - (nd - 1) * Div;
                    GRs = 0;
                    GRg = 0;
                    GRi = 0;
                    KKi = (1 - (1 - (Kg + Ki))^ (1 / nd)) / (Kg + Ki);
                    KKg = KKi * Kg;
                    KKi = KKi * Ki;
                    for m = 1:nd
                        if PPe(m) + GW < GWM
                            Grfc = 0.0;
                            if PPe(m) + GWU < GWUM
                                GWU = PPe(m) + GWU;
                            elseif GWL + GWU - GWUM + PPe(m) < GWLM
                                GWL = GWL + GWU - GWUM + PPe(m);
                                GWU = GWUM;
                            else
                                GWD = GWD + PPe(m) - (GWUM - GWU) - (GWLM - GWL);
                                GWU = GWUM;
                                GWL = GWLM;
                            end
                            GW = GWU + GWL + GWD;
                            GS = 0;
                            GRs = GRs;
                            GRg = GRg + GS * KKg;
                            GRi = GRi + GS * KKi;
                            GS = GS * (1.0 - KKg - KKi);
                        else
                            GR = PPe(m) + GW - GWM;
                            GWU = GWUM;
                            GWL = GWLM;
                            GWD = GWDM;
                            GW = GWM;
                            Grfc = 1.0;
                            if GR + GS <= GSM
                                GRs = GRs;
                                GS = GS + GR;
                                GRg = GRg + GS * KKg;
                                GRi = GRi + GS * KKi;
                                GS = GS * (1.0 - KKg - KKi);
                            else
                                GRs = GRs + GR + GS - GSM;
                                GS = GSM;
                                GRg = GRg + GSM * KKg;
                                GRi = GRi + GSM * KKi;
                                GS = GS * (1.0 - KKg - KKi);
                            end
                        end
                    end
                end
            end

            Total_GRs = Total_GRs +  GRs;
            Total_GRi = Total_GRi +  GRi;
            Total_GRg = Total_GRg +  GRg;
            E_act(ii,jj) = GE;
            if GridWi(ii,jj) >0
                E_act(ii,jj) = 0;
            end

            Ct = Garea / DT / 3.6;
            Cg = CCg ^ (DT / 24);
            Ci = CCi ^ (DT / 24);

            if Inrunmethod == 1
                GQi = GridQi(ii, jj);
                GQg = GridQg(ii, jj);
                GQs = GridQs(ii, jj);

                GQs = GRs * Ct;
                Gqch = 0;
                GQi = GQi * Ci + GRi * Ct * (1.0 - Ci);
                GQg = GQg * Cg + GRg * Ct * (1.0 - Cg);

                GridQi(ii, jj) = GQi;
                GridQg(ii, jj) = GQg;
                GridQs(ii, jj) = GQs;

                GW = GWU + GWL + GWD;
                GridW(ii, jj) = GW;
                GridS(ii, jj) = GS;
                if GW >= GWM
                    GridWM(ii, jj) = 1;
                    NoWM(i) = NoWM(i) + 1;
                else
                    GridWM(ii, jj) = 0;
                end

                if GS / (1.0 - KKg - KKi) >= GSM
                    GridSM(ii, jj) = 1;
                    NoSM(i) = NoSM(i) + 1;
                else
                    GridSM(ii, jj) = 0;
                end

                GridWU(ii, jj) = GWU;
                GridWL(ii, jj) = GWL;
                GridWD(ii, jj) = GWD;
                RFC(ii, jj) = Grfc;

            elseif Inrunmethod == 2
                GQi = GridQi(ii, jj);
                GQg = GridQg(ii, jj);
                GQs = GridQs(ii, jj);

                GQs = GRs * Ct * (1.0 - Dp);
                Gqch = GRs * Ct * Dp;
                GQi = GQi * Ci + GRi * Ct * (1.0 - Ci);
                GQg = GQg * Cg + GRg * Ct * (1.0 - Cg);

                GridQi(ii, jj) = GQi;
                GridQg(ii, jj) = GQg;
                GridQs(ii, jj) = GQs;

                GW = GWU + GWL + GWD;
                GridW(ii, jj) = GW;
                GridS(ii, jj) = GS;
                if GW >= GWM
                    GridWM(ii, jj) = 1;
                    NoWM(i) = NoWM(i) + 1;
                else
                    GridWM(ii, jj) = 0;
                end

                if GS / (1.0 - KKg - KKi) >= GSM
                    GridSM(ii, jj) = 1;
                    NoSM(i) = NoSM(i) + 1;
                else
                    GridSM(ii, jj) = 0;
                end

                GridWU(ii, jj) = GWU;
                GridWL(ii, jj) = GWL;
                GridWD(ii, jj) = GWD;
                RFC(ii, jj) = Grfc;
            end

            if GW - GridW0(ii, jj) <0
                GridF(ii, jj) = 0.00001;
            else
                GridF(ii, jj) = GW - GridW0(ii, jj);
            end
            GridRP(ii, jj) =  RP;

            switch FlowDirection(ii, jj)
                case 0
                    i1 = ii - 1;
                    j1 = jj + 1;
                case 1
                    i1 = ii;
                    j1 = jj + 1;
                case 2
                    i1 = ii + 1;
                    j1 = jj + 1;
                case 3
                    i1 = ii + 1;
                    j1 = jj;
                case 4
                    i1 = ii + 1;
                    j1 = jj - 1;
                case 5
                    i1 = ii;
                    j1 = jj - 1;
                case 6
                    i1 = ii - 1;
                    j1 = jj - 1;
                case 7
                    i1 = ii - 1;
                    j1 = jj;
                case 8
                    i1 = ii;
                    j1 = jj;
            end

            InflowQs(2, NextNoIJ(ii, jj)) = GQs + SumQs(NextNoIJ(ii, jj));

            if (i == 1)
                OutflowQs(2, NextNoIJ(ii, jj)) = InflowQs(2, NextNoIJ(ii, jj));
            end

            if (i > 1)
                MDt = DT;
                MC1 = (MXs * MKs + 0.5 * MDt) / ((1 - MXs) * MKs + 0.5 * MDt);
                MC2 = (0.5 * MDt - MXs * MKs) / ((1 - MXs) * MKs + 0.5 * MDt);
                MC3 = ((1 - MXs) * MKs - 0.5 * MDt) / ((1 - MXs) * MKs + 0.5 * MDt);
                OutflowQs(2, NextNoIJ(ii, jj)) = MC1 * InflowQs(1, NextNoIJ(ii, jj)) + MC2 * InflowQs(2, NextNoIJ(ii, jj)) + MC3 * OutflowQs(1, NextNoIJ(ii, jj));

                if (OutflowQs(2, NextNoIJ(ii, jj)) > 0)
                    OutflowQs(2, NextNoIJ(ii, jj)) = OutflowQs(2, NextNoIJ(ii, jj));
                else
                    OutflowQs(2, NextNoIJ(ii, jj)) = 0;
                end

                InflowQs(1, NextNoIJ(ii, jj)) = InflowQs(2, NextNoIJ(ii, jj));
                OutflowQs(1, NextNoIJ(ii, jj)) = OutflowQs(2, NextNoIJ(ii, jj));
            end

            if (FlowDirection(ii, jj) == 8)
                SumQs(NextNoIJ(i1, j1)) = OutflowQs(2, NextNoIJ(ii, jj));
            else
                SumQs(NextNoIJ(i1, j1)) = SumQs(NextNoIJ(i1, j1)) + OutflowQs(2, NextNoIJ(ii, jj));
            end

            InflowQi(2, NextNoIJ(ii, jj)) = GQi + SumQi(NextNoIJ(ii, jj));

            if i == 1
                OutflowQi(2, NextNoIJ(ii, jj)) = InflowQi(2, NextNoIJ(ii, jj));
            end

            if i > 1
                MDt = DT;
                MC1 = (MXi * MKi + 0.5 * MDt) / ((1 - MXi) * MKi + 0.5 * MDt);
                MC2 = (0.5 * MDt - MXi * MKi) / ((1 - MXi) * MKi + 0.5 * MDt);
                MC3 = ((1 - MXi) * MKi - 0.5 * MDt) / ((1 - MXi) * MKi + 0.5 * MDt);
                OutflowQi(2, NextNoIJ(ii, jj)) = MC1 * InflowQi(1, NextNoIJ(ii, jj)) + MC2 * InflowQi(2, NextNoIJ(ii, jj)) + MC3 * OutflowQi(1, NextNoIJ(ii, jj));

                if OutflowQi(2, NextNoIJ(ii, jj)) > 0.0
                    OutflowQi(2, NextNoIJ(ii, jj)) = OutflowQi(2, NextNoIJ(ii, jj));
                else
                    OutflowQi(2, NextNoIJ(ii, jj)) = 0.0;
                end

                InflowQi(1, NextNoIJ(ii, jj)) = InflowQi(2, NextNoIJ(ii, jj));
                OutflowQi(1, NextNoIJ(ii, jj)) = OutflowQi(2, NextNoIJ(ii, jj));
            end

            if FlowDirection(ii, jj) == 8
                SumQi(NextNoIJ(i1, j1)) = OutflowQi(2, NextNoIJ(ii, jj));
            else
                SumQi(NextNoIJ(i1, j1)) = SumQi(NextNoIJ(i1, j1)) + OutflowQi(2, NextNoIJ(ii, jj));
            end

            InflowQg(2, NextNoIJ(ii, jj)) = GQg + SumQg(NextNoIJ(ii, jj));

            if i == 1
                OutflowQg(2, NextNoIJ(ii, jj)) = InflowQg(2, NextNoIJ(ii, jj));
            end

            if i > 1
                MDt = DT;
                MC1 = (MXg * MKg + 0.5 * MDt) / ((1 - MXg) * MKg + 0.5 * MDt);
                MC2 = (0.5 * MDt - MXg * MKg) / ((1 - MXg) * MKg + 0.5 * MDt);
                MC3 = ((1 - MXg) * MKg - 0.5 * MDt) / ((1 - MXg) * MKg + 0.5 * MDt);
                OutflowQg(2, NextNoIJ(ii, jj)) = MC1 * InflowQg(1, NextNoIJ(ii, jj)) + MC2 * InflowQg(2, NextNoIJ(ii, jj)) + MC3 * OutflowQg(1, NextNoIJ(ii, jj));

                if OutflowQg(2, NextNoIJ(ii, jj)) > 0.0
                    OutflowQg(2, NextNoIJ(ii, jj)) = OutflowQg(2, NextNoIJ(ii, jj));
                else
                    OutflowQg(2, NextNoIJ(ii, jj)) = 0.0;
                end

                InflowQg(1, NextNoIJ(ii, jj)) = InflowQg(2, NextNoIJ(ii, jj));
                OutflowQg(1, NextNoIJ(ii, jj)) = OutflowQg(2, NextNoIJ(ii, jj));
            end

            if FlowDirection(ii, jj) == 8
                SumQg(NextNoIJ(i1, j1)) = OutflowQg(2, NextNoIJ(ii, jj));
            else
                SumQg(NextNoIJ(i1, j1)) = SumQg(NextNoIJ(i1, j1)) + OutflowQg(2, NextNoIJ(ii, jj));
            end

            InflowQch(2, NextNoIJ(ii, jj)) = Gqch + SumQch(NextNoIJ(ii, jj));

            if i == 1
                OutflowQch(2, NextNoIJ(ii, jj)) = InflowQch(2, NextNoIJ(ii, jj));
            end

            if i > 1
                MDt = DT;
                MC1 = (MXch * MKch + 0.5 * MDt) / ((1 - MXch) * MKch + 0.5 * MDt);
                MC2 = (0.5 * MDt - MXch * MKch) / ((1 - MXch) * MKch + 0.5 * MDt);
                MC3 = ((1 - MXch) * MKch - 0.5 * MDt) / ((1 - MXch) * MKch + 0.5 * MDt);
                OutflowQch(2, NextNoIJ(ii, jj)) = MC1 * InflowQch(1, NextNoIJ(ii, jj)) + MC2 * InflowQch(2, NextNoIJ(ii, jj)) + MC3 * OutflowQch(1, NextNoIJ(ii, jj));

                if OutflowQch(2, NextNoIJ(ii, jj)) > 0.0
                    OutflowQch(2, NextNoIJ(ii, jj)) = OutflowQch(2, NextNoIJ(ii, jj));
                else
                    OutflowQch(2, NextNoIJ(ii, jj)) = 0.0;
                end

                InflowQch(1, NextNoIJ(ii, jj)) = InflowQch(2, NextNoIJ(ii, jj));
                OutflowQch(1, NextNoIJ(ii, jj)) = OutflowQch(2, NextNoIJ(ii, jj));
            end


            if RiverPoint(i1, j1) == 1
                if FlowDirection(ii, jj) == 8
                    SumQch(NextNoIJ(i1, j1)) = OutflowQch(2, NextNoIJ(ii, jj)) + SumQi(NextNoIJ(i1, j1)) + SumQg(NextNoIJ(i1, j1));
                else
                    SumQch(NextNoIJ(i1, j1)) = SumQch(NextNoIJ(i1, j1)) + OutflowQch(2, NextNoIJ(ii, jj)) + SumQi(NextNoIJ(i1, j1)) + SumQg(NextNoIJ(i1, j1));
                end
                SumQi(NextNoIJ(i1, j1)) = 0.0;
                SumQg(NextNoIJ(i1, j1)) = 0.0;
            else
                if FlowDirection(ii, jj) == 8
                    SumQch(NextNoIJ(i1, j1)) = OutflowQch(2, NextNoIJ(ii, jj));
                else
                    SumQch(NextNoIJ(i1, j1)) = SumQch(NextNoIJ(i1, j1)) + OutflowQch(2, NextNoIJ(ii, jj));
                end
            end

            Qoutch = SumQch(NextNoIJ(i1, j1));
            Qouts = SumQs(NextNoIJ(i1, j1));
            Qouti = SumQi(NextNoIJ(i1, j1));
            Qoutg = SumQg(NextNoIJ(i1, j1));

            QSim(k, i + LagTime) = QSim(k, i + LagTime - 1) * CCS + (Qoutch + Qouts + Qouti + Qoutg) * (1 - CCS);

        end
        TotalPre_R(i, 1) = TotalPre_rain / DRMGCNo;
        TotalPre_R(i, 2) = TotalPre_snow / DRMGCNo;

        TotalPre_R(i, 3) = Total_GRs / DRMGCNo;
        TotalPre_R(i, 4) = Total_GRi / DRMGCNo;
        TotalPre_R(i, 5) = Total_GRg / DRMGCNo;

        totalSum_Wi = 0;
        totalSum_SFD = 0;
        totalSum_Eact = 0;

        for iii = 1:Nx
            for jjj = 1:Ny

                if GridWi(iii, jjj) >= 0
                    totalSum_Wi = totalSum_Wi + GridWi(iii, jjj);
                end
                if E_act(iii, jjj) >= 0
                    totalSum_Eact = totalSum_Eact + E_act(iii, jjj);
                end
            end
        end
        TotalPre_R(i, 6) = totalSum_Wi / DRMGCNo;
        TotalPre_R(i, 7) = totalSum_SFD / DRMGCNo;
        TotalPre_R(i, 8) = totalSum_Eact / DRMGCNo;

        for j = 1:DRMGCNo
            GridQSim(SortingRow(j), SortingCol(j), i) = QSim(j, i);
        end

    end

    ASC_Info = [Ny, Nx, XllCorner, YllCorner, DDem, Nodata];
    if save_state == 1
        SumQ = zeros(DRMGCNo, 4);
        AllFilpath.SumQFilepath = "../../Data/Initial_state/SumQ.txt";
        write_asc(AllFilpath.SumQFilepath, ASC_Info, [SumQch, SumQs, SumQi, SumQg]);
    end

    if save_state == 1
        AllFilpath.GridWFilepath = "../../Data/Initial_state/GridW.asc";
        write_asc(AllFilpath.GridWFilepath, ASC_Info, GridW);

        AllFilpath.GridSFilepath = "../../Data/Initial_state/GridS.asc";
        write_asc(AllFilpath.GridSFilepath, ASC_Info, GridS);

        AllFilpath.GridWUFilepath = "../../Data/Initial_state/GridWU.asc";
        write_asc(AllFilpath.GridWUFilepath, ASC_Info, GridWU);

        AllFilpath.GridWLFilepath = "../../Data/Initial_state/GridWL.asc";
        write_asc(AllFilpath.GridWLFilepath, ASC_Info, GridWL);

        AllFilpath.GridQiFilepath = "../../Data/Initial_state/GridQi.asc";
        write_asc(AllFilpath.GridQiFilepath, ASC_Info, GridQi);

        AllFilpath.GridQgFilepath = "../../Data/Initial_state/GridQg.asc";
        write_asc(AllFilpath.GridQgFilepath, ASC_Info, GridQg);

        AllFilpath.GridQsFilepath = "../../Data/Initial_state/GridQs.asc";
        write_asc(AllFilpath.GridQsFilepath, ASC_Info, GridQs);
    end

    if save_state == 1
        AllFilpath.GridDsFilepath = strcat("../../Data/Initial_state/", STCD, "_Snowdepth.asc");
        write_asc(AllFilpath.GridDsFilepath, ASC_Info, GridDs);

        AllFilpath.SimGridQFilepath = "../../Data/Initial_state/SimGridQ.asc";
        write_asc(AllFilpath.SimGridQFilepath, ASC_Info, GridQSim(:, :, ntime1));
    end
    toc
end

TotalQ_tmp = zeros(TSteps, 2*NoStation+2+3);
StaQSim = zeros(TSteps, NoStation);

for kk = 1 : NoStation
    if lonlat(kk,3) ~= Nodata
        ii = Nx-ceil((lonlat(kk,3)-YllCorner)/DDem);
        jj = ceil((lonlat(kk,2)-XllCorner)/DDem);

        for i = 1:TSteps
            StaQSim(i, kk) = GridQSim(ii, jj, i);
            TotalQ_tmp(i, 1) = TotalPre_R(i, 1);
            TotalQ_tmp(i, 2) = TotalPre_R(i, 2);
            TotalQ_tmp(i, 3) = TotalPre_R(i, 3);
            TotalQ_tmp(i, 4) = TotalPre_R(i, 4);
            TotalQ_tmp(i, 5) = TotalPre_R(i, 5);
            TotalQ_tmp(i, 6) = TotalPre_R(i, 6);
            TotalQ_tmp(i, 7) = TotalPre_R(i, 7);
            TotalQ_tmp(i, 8) = TotalPre_R(i, 8);

            TotalQ_tmp(i, 2*kk-1+3+3+2) = TotalQobs(i, kk);
            TotalQ_tmp(i, 2*kk+3+3+2) = GridQSim(ii, jj, i);

        end
    end
end

observation = TotalQ_tmp(:,9);
simulation = TotalQ_tmp(:, 10);
RE = relative_error(observation, simulation);
NSE = nash_sutcliffe_efficiency(observation, simulation);
start_date = datetime(2000, 1, 1);
end_date = datetime(2018, 12, 31);
date_series = (start_date : end_date)';

data_loc = ismember(date_series.Month, [3, 4, 5, 6]);
lowsim = data_loc.*simulation;
lowobs = data_loc.*observation;
low_obs = observation(data_loc);
low_sim = simulation(data_loc);
low_nse = nash_sutcliffe_efficiency(low_obs, low_sim);
low_re = relative_error(low_obs, low_sim);

disp("Saving results")
dates = datetime(2014,1,1) + caldays(0:ntime1);
AllFilpath.OutFilepath = "../output/StaQSim.txt";
% result=[IY2,IM2,ID2,TotalQ_tmp];
% write_result(AllFilpath.OutFilepath, result);
subplot(211)
plot(1:TSteps, observation, 1:TSteps, simulation)
xticks([1, 366, 731,1096,1461,1826,2191,2556,2921,3286,3651,4016]);
xticklabels({"200001", "200101","200201", "200301","200401", "200501","200601", "200701","200801", "200901","201001","201012"});
legend(["Obs", "Simu"])
text(400, 3100, ['RE=', num2str(RE)]);
text(400, 2900, ['NSE=', num2str(NSE)]);
xlabel('Date');
ylabel('Streamflow（m3/s）');
% saveas(gcf,"../output/StaQSim.png");
subplot(212)
plot(1:length(low_sim), low_obs, 1:length(low_sim), low_sim)
legend(["Obs", "Simu"])
text(400, 3100, ['lowre=', num2str(low_re)]);
text(400, 2900, ['lownse=', num2str(low_nse)]);