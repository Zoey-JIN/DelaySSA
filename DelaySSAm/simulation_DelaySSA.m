function result = simulation_DelaySSA(algorithm, sample_size, tmax, n_initial, t_initial, S_matrix, S_matrix_delay, k, reactant_matrix, delay_type, delaytime_list)
    % SIMULATION_DELAYSSA Simulates the reaction using the specified algorithm with delay handling.
    % Inputs:
    %   algorithm - The chosen simulation algorithm (string).
    %   sample_size - Number of samples to generate.
    %   tmax - Maximum time for the simulation.
    %   n_initial - Initial number of molecules.
    %   t_initial - Initial time.
    %   S_matrix - Stoichiometric matrix.
    %   S_matrix_delay - Stoichiometric matrix with delay (optional).
    %   k - Rate constants or function handle for rate constants.
    %   reactant_matrix - Matrix describing the reactants.
    %   delay_type - Vector indicating the type of delay.
    %   delaytime_list - List of delay times.
    
    % Check delay-related conditions if applicable
    if  ~isempty(delay_type) && ~isempty(delaytime_list) && ~isempty(S_matrix_delay)
        delay_effect_matrix = check_delay_relative(delay_type, delaytime_list, S_matrix, S_matrix_delay);
    else
        delay_effect_matrix = [];
    end
    
    % Choose the simulation function based on the algorithm
    switch algorithm
        case 'DelayRejection'
            result = arrayfun(@(x) delay_rejection(tmax, n_initial, t_initial, S_matrix, S_matrix_delay, k, reactant_matrix, delay_type, delaytime_list, delay_effect_matrix), 1:sample_size, 'UniformOutput', false);
        case 'DelayMNR'
            result = arrayfun(@(x) delay_modifiednextreaction(tmax, n_initial, t_initial, S_matrix, S_matrix_delay, k, @fun_fr, delay_type, delaytime_list, delay_effect_matrix), 1:sample_size, 'UniformOutput', false);
        case 'DelayDirect'
            result = arrayfun(@(x) delay_direct(tmax, n_initial, t_initial, S_matrix, S_matrix_delay, k, @fun_fr, delay_type, delaytime_list, delay_effect_matrix), 1:sample_size, 'UniformOutput', false);
        otherwise
            error('Error: No such algorithm!');
    end
end

