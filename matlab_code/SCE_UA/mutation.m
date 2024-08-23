function mutation_array = mutation(num_samples, upper_bound, lower_bound)
    num_params = numel(upper_bound);
    mutation_array = zeros(num_samples, num_params);


    for i = 1:num_params
        
        variance = (upper_bound(i) - lower_bound(i)) / 4; 
      
        mutation_array(:, i) = normrnd(0, sqrt(variance), num_samples, 1);
  
        mutation_array(:, i) = min(max(mutation_array(:, i), -0.25 * (upper_bound(i) - lower_bound(i))), 0.25 * (upper_bound(i) - lower_bound(i)));
    end
end