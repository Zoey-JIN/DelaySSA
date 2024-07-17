function result = delay_modifiednextreaction(tmax, n_initial, t_initial, S_matrix, S_matrix_delay, k, reactant_matrix, delay_type, delaytime_list, delay_effect_matrix)    
    
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
    
    p_vec = zeros(1, size(S_matrix, 2));
    t_vec = zeros(1, size(S_matrix, 2));
    Tstruct = []; 

    
    f_r = fun_fr(k, n);
    u2 = rand(1, size(S_matrix, 2));
    p_vec = log(1 ./ u2);

    while t < tmax
        tau_vec = (p_vec - t_vec) ./ f_r + t;

        [min_1,r_1] = min(tau_vec);
        if numel(Tstruct)==0
            min_2 = [];
            r_2 = [];
        else
            min_2 = Tstruct(1,1);
            r_2 = Tstruct(2,1);
        end
        if  isempty(min_2) || min_1 < min_2 || isnan(min_2) 
            r = r_1;
            tau = min_1 - t;
            t = min_1;
            
            if delay_type(r) == 0
                if ismember(r, delay_effect_matrix(1, :))
                    effect_r = delay_effect_matrix(2, delay_effect_matrix(1, :) == r);
                    drop_index = find(Tstruct(2) == effect_r);
                    drop_index = randsample(drop_index, 1);
                    Tstruct(:,drop_index) = [];
                end
                n = n + S_matrix(:, r);
            elseif delay_type(r) == 1
                tmp = tau_element(delaytime_list(r));
                add_tau = tmp{1} + t;
                tmp = [add_tau,r];
                if isempty(Tstruct)
                    Tstruct = tmp(:);
                else
                    index = find(Tstruct(1,:) >= add_tau, 1);
                    if isempty(index)
                        Tstruct = [Tstruct,tmp(:)];
                    else
                        Tstruct = [Tstruct(:,1:index-1), tmp(:), Tstruct(:,index:end)];
                    end
                end
            elseif delay_type(r) == 2
                n = n + S_matrix(:, r);
                tmp = tau_element(delaytime_list(r));
                add_tau = tmp{1} + t;
                tmp = [add_tau,r];
                if isempty(Tstruct)
                    Tstruct = tmp(:);
                else
                    index = find(Tstruct(1,:) >= add_tau, 1);
                    if isempty(index)
                        Tstruct = [Tstruct,tmp(:)];
                    else
                        Tstruct = [Tstruct(:,1:index-1), tmp(:), Tstruct(:,index:end)];
                    end
                end
            else
                error('Stop: wrong with the reaction type');
            end
            
            u1 = rand();
            p_vec(r) = p_vec(r) + log(1 / u1);
        else
            r = r_2;
            tau = min_2 - t;
            t = min_2;
            n = n + S_matrix_delay(:, r);
            Tstruct(:,1) = [];
        end

        t_vec = t_vec + f_r * tau;
        f_r = fun_fr(k, n);
        
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
