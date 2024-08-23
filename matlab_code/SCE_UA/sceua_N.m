function [best_params, calibration_results] = sceua_N(bounds, n_iterations, n_complexes, n_parameters, x_obs, y_obs, UL, CON)
    fn = @(x, params, observation, UL, CON) objective_function(x, params, observation, UL, CON);

    lower_bound = bounds(:, 1)';
    upper_bound = bounds(:, 2)';
    % population = rand(n_complexes, n_parameters) .* (upper_bound - lower_bound) + lower_bound;
    population = repmat(lower_bound, n_complexes, 1) + rand(n_complexes, n_parameters) .* repmat(upper_bound - lower_bound, n_complexes, 1);

    calibration_results = nan(n_iterations, 1);


    for iteration = 1:n_iterations
        

        obj_values = nan(size(population, 1), 1);

        for t = 1 : size(population, 1)
        
            obj_values(t) = objective_function(x_obs, population(t, :), y_obs, UL, CON);
        end
        % obj_values = arrayfun(@(params) objective_function(x_obs, params, y_obs), population);
        calibration_results(iteration) = min(obj_values);

        [sorted_obj_values, sorted_indices] = sort(obj_values);
        parents = population(sorted_indices(1:n_complexes/2), :);
        % CCEUA
        [snew, ~] = cceua_N(fn, parents, sorted_obj_values(1:n_complexes/2), lower_bound, upper_bound, x_obs, y_obs, UL, CON);
        parents(n_complexes/2, :) = snew;

        children = zeros(size(parents));
        for i = 1:n_complexes/2
            rand_indices = randperm(n_complexes/2, 2);
            parent1 = parents(rand_indices(1), :);
            parent2 = parents(rand_indices(2), :);
            beta = rand(1, n_parameters);
            children(i, :) = beta .* parent1 + (1 - beta) .* parent2;
        end

        mutation_scale = 0.1;  
        mutation_array = mutation(n_complexes/2, upper_bound, lower_bound);
        mutation_array = mutation_array * mutation_scale;
        children = children + mutation_array;
        
        population = [parents; children];
        
        population = max(population, repmat(lower_bound, n_complexes, 1));
        population = min(population, repmat(upper_bound, n_complexes, 1));
%         end_time = clock;
%         disp(strcat("Time for iteration ", num2str(iteration), " is ", num2str(round(etime(end_time, start_time), 2)), "s"))
    end

    [~, min_index] = min(obj_values);
    best_params = population(min_index, :);
end