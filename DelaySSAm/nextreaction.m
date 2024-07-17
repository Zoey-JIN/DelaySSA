function result = nextreaction(tmax, n_initial, t_initial, S_matrix, k, reactant_matrix)
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

    f_r = fun_fr(k, n);

    u2 = rand(1, size(S_matrix, 2)); 
    tau_vec = -log(u2) ./ f_r;

    while t < tmax

        [tau, r] = min(tau_vec);

        n = n + S_matrix(:, r);

        f_r_update = fun_fr(k, n);

        if numel(tau_vec) > 1
            tau_vec(tau_vec ~= tau_vec(r)) = f_r(tau_vec ~= tau_vec(r)) ./ f_r_update(tau_vec ~= tau_vec(r)) .* (tau_vec(tau_vec ~= tau_vec(r)) - tau) + tau;
        end

        inf_idx = isinf(tau_vec) | isnan(tau_vec);
        if any(inf_idx)
            tau_vec(inf_idx) = -log(rand(1,sum(inf_idx))) ./ f_r_update(inf_idx) + tau;
        end

        u1 = rand;
        tau_vec(r) = 1 / f_r_update(r) * log(1 / u1) + tau;
        f_r = f_r_update;
        t = tau;

        if t < t_values(end)
            break;
        end

        % Append new values to vectors
        t_values = [t_values; t];
        n_values = [n_values, n];
    end
    result = struct('t_values', t_values, 'n_values', n_values);
end
