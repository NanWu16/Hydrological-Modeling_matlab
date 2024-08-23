function CLM = read_climatic(CON)

    sim_startday = datetime([num2str(CON.IYB),'-',num2str(CON.IMB),'-',num2str(CON.IDB)], 'Inputformat', 'yyyy-MM-dd');
    sim_endday = datetime([num2str(CON.IYE),'-',num2str(CON.IME),'-',num2str(CON.IDE)], 'Inputformat', 'yyyy-MM-dd');


    data_startday = datetime('2000-01-01', 'InputFormat', 'yyyy-MM-dd');   
    loc_start = days(sim_startday - data_startday) + 1;
    loc_end = days(sim_endday - data_startday) + 1;
    load ../../Data/Prec/Precipitation.mat
    CLM.Prec = Precipitation(:, :, loc_start : loc_end);
    clear Precipitation;
    load ../../Data/PET/PET0.mat
    CLM.PET = PET0(:, :, loc_start : loc_end);
    clear PET0;
    load ../../Data/Temp/Temperature.mat
    CLM.Temp = Temperature(:, :, loc_start : loc_end);
    clear Temperature;
end