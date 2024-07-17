function result = delay_direct(tmax, n_initial, t_initial, S_matrix, S_matrix_delay, k, reactant_matrix, delay_type, delaytime_list, delay_effect_matrix)    
    
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
    Tstruct = [];

    while t < tmax
        u2 = rand;
        if numel(Tstruct)==0
            f_r = fun_fr(k, n);
            lambda_sum = sum(f_r);
            tau = -log(u2) / lambda_sum;
        else
            i_index = 1;
            f_r = fun_fr(k, n);
            lambda_sum = sum(f_r);
            a0 = 0;
            at = lambda_sum * (Tstruct(1,1) - t);
            F = 1 - exp(-at);
            
            while F < u2
                r = Tstruct(2,i_index);
                if delay_type(r) == 0
                    warning('Warning: delay_type is 0.');
                elseif delay_type(r) == 1
                    n = n + S_matrix_delay(:, r);
                elseif delay_type(r) == 2
                    n = n + S_matrix_delay(:, r);
                else
                    error('Stop: wrong with the reaction type.');
                end
                
                f_r = fun_fr(k, n);
                lambda_sum = sum(f_r);
                add = lambda_sum * (Tstruct(1,i_index + 1) - Tstruct(1,i_index));
                a0 = at;
                at = at + add;
                if isnan(add)
                    at = Inf;
                end
                F = 1 - exp(-at);
                t = Tstruct(1,i_index);
                i_index = i_index + 1;
                t_values = [t_values, t];
                n_values = [n_values, n];
            end
            
            i_index = i_index - 1;
            tau = -(log(1 - u2) + a0) / lambda_sum;
            if i_index < numel(Tstruct(1,:))
                Tstruct(1,:) = Tstruct(1,(i_index + 1):end);
                Tstruct(2,:) = Tstruct(2,(i_index + 1):end);
            else
                Tstruct = [];
            end
        end

        % 选择反应
        u1 = rand();
        [~, r] = max(cumsum(f_r) > u1 * lambda_sum);
        t = t + tau;

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
            error('Stop: wrong with the reaction type.');
        end
        
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
