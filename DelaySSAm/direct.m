function result = direct(tmax, n_initial, t_initial, S_matrix, k, reactant_matrix)
    % SIMULATE_REACTION Simulates the reaction process until the specified maximum time.
    % Inputs:
    %   tmax - Maximum simulation time.
    %   n_initial - Initial number of molecules.
    %   t_initial - Initial time.
    %   S_matrix - Stoichiometric matrix.
    %   k - Rate constants or function handle for rate constants.
    %   fun_fr - Function handle for propensity calculation.
    
    % Define the function for propensity calculation
    function result = fun_fr(k, n)
        % Check if k is a function handle
        if isa(k, 'function_handle')
            k_mask = k(n);
        else
            k_mask = k;
        end

        propensity_values = propensity_n(n, reactant_matrix);

        result = k_mask .* propensity_values;
    end

    

    function result = propensity_n(n, reactant_matrix)
        % PROPENSITY_N Computes the propensity values for the given reactant matrix.
        % Inputs:
        %   n - A vector of values.
        %   reactant_matrix - A matrix with reactant information.
        %
        % Outputs:
        %   result - A vector of propensity values for each column of the reactant_matrix.

        num_cols = size(reactant_matrix, 2);
        num_rows = size(reactant_matrix, 1); 

        result = ones(1, num_cols);

        for j = 1:num_cols

            prod_values = ones(1, num_rows);

            for i = 1:num_rows
                if reactant_matrix(i, j) == 0
                    prod_values(i) = 1;
                else

                    prod_values(i) = prod(n(i) - (0:(reactant_matrix(i, j) - 1)));
                end
            end

            result(j) = prod(prod_values);
        end
    end

    % Initialize time and molecule number vectors
    n_values = n_initial(:);  
    t_values = t_initial;
    n = n_initial;
    t = t_initial;

    % Simulation loop
    while t < tmax
        u1 = rand;  % Generate random number between 0 and 1
        u2 = rand;  % Generate random number between 0 and 1

        % Calculate propensity functions
        f_r = fun_fr(k, n);
        lambda_sum = sum(f_r);
        tau = -log(u1) / lambda_sum;

        % Find the reaction to occur
        cumulative_f_r = cumsum(f_r);
        r = find(cumulative_f_r > u2 * lambda_sum, 1);

        % Update the number of molecules and time
        n = n + S_matrix(:, r);
        t = t + tau;

        % Check if the current time is less than the last recorded time
        if t < t_values(end)
            break;
        end

        % Append new values to vectors
        t_values = [t_values; t];
        n_values = [n_values, n];
    end
    result = struct('t_values', t_values, 'n_values', n_values);
end
