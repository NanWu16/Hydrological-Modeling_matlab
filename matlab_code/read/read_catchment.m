function [catid, STCD, SCF, NoStation] = read_catchment(CatchID_Filepath)
    fid = fopen(CatchID_Filepath, 'r');

   
    catid = fgetl(fid);
    
   
    STCD = fgetl(fid);
    
   
    ctemp = fgetl(fid);
    SCF = fscanf(fid, '%f');
    
    
    ctemp = fgetl(fid);
    NoStation = fscanf(fid, '%f');
    
    fclose(fid);
end