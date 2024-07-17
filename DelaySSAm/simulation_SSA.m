function result = simulation_SSA(algorithm, sample_size, tmax, n_initial, t_initial, S_matrix, k, reactant_matrix)
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

    

    % Choose the simulation function based on the algorithm
    switch algorithm
        case 'Direct'
            result = arrayfun(@(x) direct(tmax, n_initial, t_initial, S_matrix, k, reactant_matrix), 1:sample_size, 'UniformOutput', false);
        case 'MNR'
            result = arrayfun(@(x) modifiednextreaction(tmax, n_initial, t_initial, S_matrix, k, @fun_fr), 1:sample_size, 'UniformOutput', false);
        case 'NR'
            result = arrayfun(@(x) nextreaction(tmax, n_initial, t_initial, S_matrix, k, @fun_fr), 1:sample_size, 'UniformOutput', false);
        otherwise
            error('Error: No such algorithm!');
    end
end