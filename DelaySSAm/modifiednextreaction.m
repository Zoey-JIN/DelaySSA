function result = modifiednextreaction(tmax, n_initial, t_initial, S_matrix, k, reactant_matrix)
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

    n_values = n_initial(:);  
    t_values = t_initial;
    n = n_initial;
    t = t_initial;
    
    p_vec = zeros(1, size(S_matrix, 2));
    t_vec = zeros(1, size(S_matrix, 2));

    f_r = fun_fr(k, n);
    u2 = rand(1, size(S_matrix, 2));
    p_vec = log(1 ./ u2);

    while t < tmax
        tau_vec = (p_vec - t_vec) ./ f_r;
        [~, r] = min(tau_vec);
        tau = tau_vec(r);
        n = n + S_matrix(:, r);

        t_vec = t_vec + f_r * tau;
        u1 = rand;
        p_vec(r) = p_vec(r) + log(1 / u1);
        f_r = fun_fr(k, n);

        t = t + tau;
        if t < t_values(end)
            break;
        end

        % Append new values to vectors
        t_values = [t_values; t];
        n_values = [n_values, n];
    end
    result = struct('t_values', t_values, 'n_values', n_values);
end
