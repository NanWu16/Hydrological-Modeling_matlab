function CON = read_conf()

    AllFilpath.CatchID_Filepath = '../../Data/CatID.dat';
    if ~exist(AllFilpath.CatchID_Filepath, 'file'), disp('The basin information file is missing'); return; end

    [CON.catid, CON.STCD, CON.NoBasin, CON.NoStation] = read_catchment(AllFilpath.CatchID_Filepath);
    disp(strcat("Catchment name:", CON.STCD));

    AllFilpath.Configure_Filepath = "../../Data/Configure.dat";
    if ~exist(AllFilpath.Configure_Filepath, 'file'), disp('The global congifure file is missing'); return; end
    [CON.IYB, CON.IMB, CON.IDB, CON.IYE, CON.IME, CON.IDE, CON.Nwarmday, CON.Nsm, CON.SCF, CON.SnowCoef, CON.FreeWaterCoef, CON.TensionWaterCoef, CON.Div, CON.save_state, CON.load_state, CON.Inrunmethod, CON.DT, CON.TimeStepCmp] = read_configure(AllFilpath.Configure_Filepath);
end