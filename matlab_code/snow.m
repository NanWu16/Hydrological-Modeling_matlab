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

[IY1, IM1, ID1, IH1, IY2, IM2, ID2, IH2] = data_length(IYB, IMB, IDB, 0, IYE, IME, IDE, 23);
ntime = int32(length(IY1));
ntime1 = int32(length(IY2));

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
    AllFilpath.SnowDepth_Filepath = strcat("../Data/ASC_file/", catid, '积雪深度.000');
else
    AllFilpath.SnowDepth_Filepath = strcat("../Data/Initial_state/", catid, '积雪深度.asc');
end
if ~exist(AllFilpath.SnowDepth_Filepath, 'file'), disp('The snow depth file is missing'); return; end
ctemp = readtable(AllFilpath.SnowDepth_Filepath, 'headerlines', 6, 'FileType', 'text');
GridDs = table2array(ctemp);

AllFilpath.OutletInfoFilepath = "../Data/Qobs_Station_of_Subbasin.txt";
if ~exist(AllFilpath.OutletInfoFilepath, 'file'), disp('The outlet info file is missing'); return; end
lonlat = textread(AllFilpath.OutletInfoFilepath, '', 'headerlines', 1);

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
    TSteps = ntime1;
   
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

        end 

        totalSum_Wi = 0;
        totalSum_SFD = 0;
       
        for iii = 1:Nx
            for jjj = 1:Ny
                
                if GridWi(iii, jjj) >= 0
                    totalSum_Wi = totalSum_Wi + GridWi(iii, jjj);
                end
            end
        end
        TotalPre_R(i, 6) = totalSum_Wi / DRMGCNo;
        TotalPre_R(i, 7) = totalSum_SFD / DRMGCNo;
        
       
    end 
end 

