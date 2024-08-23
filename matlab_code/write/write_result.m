function write_result(Filepath, Data)
    fid = fopen(Filepath, 'w');


    [r, c] = size(Data); 
    for ii = 1 : r
        for jj = 1 : c
            if jj <= 3
                fprintf(fid, '%d\t', Data(ii, jj));
            else
                fprintf(fid, '%.2f\t', Data(ii, jj));
            end
        end
        fprintf(fid, '\n');
    end

    fclose(fid);
end