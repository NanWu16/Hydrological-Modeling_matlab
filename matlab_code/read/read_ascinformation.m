function [Ny, Nx, XllCorner, YllCorner, DDem, Nodata] = read_ascinformation(ASC_Filepath)
    fid = fopen(ASC_Filepath, 'r');
    pattern = '[0-9.]'; 
    
    Ny = fgetl(fid);
    Ny = str2num(cell2mat(regexp(Ny, pattern, 'match'))); 
    
    Nx = fgetl(fid);
    Nx = str2num(cell2mat(regexp(Nx, pattern, 'match')));
    
    XllCorner = fgetl(fid);
    XllCorner = str2num(cell2mat(regexp(XllCorner, pattern, 'match')));
    
    YllCorner = fgetl(fid);
    YllCorner = str2num(cell2mat(regexp(YllCorner, pattern, 'match')));

    DDem = fgetl(fid);
    DDem = str2num(cell2mat(regexp(DDem, pattern, 'match')));

    Nodata = fgetl(fid);
    Nodata = -1 * str2num(cell2mat(regexp(Nodata, pattern, 'match')));
    
    fclose(fid);
end