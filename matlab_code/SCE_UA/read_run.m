function RUN = read_run(CON)
    ntime1      = CON.ntime1;
    catid       = CON.catid;
    NoStation	= CON.NoStation;
    IYB         = CON.IYB;
    IMB         = CON.IMB;
    IDB         = CON.IDB;
    IYE         = CON.IYE;
    IME         = CON.IME;
    IDE         = CON.IDE;
    
    
    AllFilpath.Discharge_Filepath = strcat("../../Data/QH_", catid, ".xlsx");
    if ~exist(AllFilpath.Discharge_Filepath, 'file'), disp('The observation runoff file is missing'); return; end
    
    TotalQ = xlsread(AllFilpath.Discharge_Filepath);
    
    start_loc = find(((TotalQ(:, 1) == IYB) & (TotalQ(:, 2) == IMB) & (TotalQ(:, 3) == IDB)) == 1);
    endt_loc  = find(((TotalQ(:, 1) == IYE) & (TotalQ(:, 2) == IME) & (TotalQ(:, 3) == IDE)) == 1);
    
    RUN = TotalQ(start_loc : endt_loc, 4);
end