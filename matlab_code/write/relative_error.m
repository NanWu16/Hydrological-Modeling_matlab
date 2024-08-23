function RE = relative_error(observation, simulation)
    mean_obs = nanmean(observation);
    mean_sim = nanmean(simulation);
    
    RE = (mean_sim - mean_obs) / mean_obs;
end