function value = objective_function(x, params, observation, UL, CON)

    simulation = GXAJ_Snow17_Frozensoil_func(x, params, UL, CON);
    
    % Whole series
    nse = nash_sutcliffe_efficiency(observation, simulation);
    % re  = abs(relative_error(observation, simulation));

    

    start_date = datetime(CON.IYB, CON.IMB, CON.IDB);
    end_date = datetime(CON.IYE, CON.IME, CON.IDE);
    date_series = (start_date : end_date)';
    
    data_loc = ismember(date_series.Month, [1, 2, 3, 4, 11, 12]);
    low_obs = observation(data_loc);
    low_sim = simulation(data_loc);
    low_nse = nash_sutcliffe_efficiency(low_obs, low_sim);

    value = (1 - nse) + (1 - low_nse);
end