function result = simulation_DelaySSA(sample_size, tmax, n_initial, t_initial, S_matrix, S_matrix_delay, k, reactant_matrix, delay_type, delaytime_list)

    if  ~isempty(delay_type) && ~isempty(delaytime_list) && ~isempty(S_matrix_delay)
        delay_effect_matrix = check_delay_relative(delay_type, delaytime_list, S_matrix, S_matrix_delay);
    else
        delay_effect_matrix = [];
    end
        result = arrayfun(@(x) delay_rejection(tmax, n_initial, t_initial, S_matrix, S_matrix_delay, k, reactant_matrix, delay_type, delaytime_list, delay_effect_matrix), 1:sample_size, 'UniformOutput', false);
end

