function [] = write_asc(Filepath, ASC_Info, Data)
    fid = fopen(Filepath, 'w');

    Ny = ASC_Info(1);
    Nx = ASC_Info(2);
    XllCorner = ASC_Info(3);
    YllCorner = ASC_Info(4);
    DDem = ASC_Info(5);
    Nodata = ASC_Info(6);

    fprintf(fid, 'ncols      %d\n', Ny);
    fprintf(fid, 'nrows      %d\n', Nx);
    fprintf(fid, 'xllcorner  %f\n', XllCorner);
    fprintf(fid, 'yllcorner  %f\n', YllCorner);
    fprintf(fid, 'cellsize   %f\n', DDem);
    fprintf(fid, 'NODATA_value %f\n', Nodata);

    [r, c] = size(Data);
    for ii = 1 : r
        for jj = 1 : c
            fprintf(fid, '%.3f\t', Data(ii, jj));
        end
        fprintf(fid, '\n');
    end

    fclose(fid);
end