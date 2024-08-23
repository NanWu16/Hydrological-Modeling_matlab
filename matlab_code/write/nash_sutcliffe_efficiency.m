function NSE = nash_sutcliffe_efficiency(observation, simulation)
    NSE = 1 - sum((observation - simulation) .^ 2) / sum((observation - mean(observation)) .^ 2);
end