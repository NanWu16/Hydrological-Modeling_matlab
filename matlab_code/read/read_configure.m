function [IYB, IMB, IDB, IYE, IME, IDE, Nwarmday, Nsm, SCF, SnowCoef, FreeWaterCoef, TensionWaterCoef, Div, save_state, load_state, Inrunmethod, DT, TimeStepCmp] = read_configure(Configure_Filepath)
    fid = fopen(Configure_Filepath, 'r');

  
    ctemp = fgetl(fid);
    ctemp = fscanf(fid, '%d  %d  %d');
    IYB = ctemp(1); IMB = ctemp(2); IDB = ctemp(3);
    
   
    ctemp = fgetl(fid);
    ctemp = fscanf(fid, '%d  %d  %d');
    IYE = ctemp(1); IME = ctemp(2); IDE = ctemp(3);
    
   
    ctemp = fgetl(fid);
    Nwarmday = fscanf(fid, '%f');
    
   
    ctemp = fgetl(fid);
    Nsm = fscanf(fid, '%f');
    
   
    ctemp = fgetl(fid);
    SCF = fscanf(fid, '%f');
    
   
    ctemp = fgetl(fid);
    SnowCoef = fscanf(fid, '%f');
    
   
    ctemp = fgetl(fid);
    FreeWaterCoef = fscanf(fid, '%f');
    
   
    ctemp = fgetl(fid);
    TensionWaterCoef = fscanf(fid, '%f');
    
   
    ctemp = fgetl(fid);
    Div = fscanf(fid, '%f');
    
  
    ctemp = fgetl(fid);
    save_state = fscanf(fid, '%f');
    
  
    ctemp = fgetl(fid);
    load_state = fscanf(fid, '%f');
    
    
    ctemp = fgetl(fid);
    Inrunmethod = fscanf(fid, '%f');
    
   
    ctemp = fgetl(fid);
    DT = fscanf(fid, '%f');
    
   
    ctemp = fgetl(fid);
    TimeStepCmp = fscanf(fid, '%f');
    fclose(fid);
end