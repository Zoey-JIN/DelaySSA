function result_matrix = check_delay_relative(delay_type, delaytime_list, S_matrix, S_matrix_delay)
    % CHECK_DELAY_RELATIVE Checks the delay type and delay time conditions.
    % If conditions are not met, it issues warnings.
    % Returns a matrix with matching column indices.
    %
    % Inputs:
    % delay_type - A vector indicating the delay type for each entry.
    % delaytime_list - A list of delay times corresponding to delay_type.
    % S_matrix - A matrix with the S matrix data.
    % S_matrix_delay - A matrix with the delayed S matrix data.
    %
    % Outputs:
    % result_matrix - A matrix of matched indices.

    % Check columns in S_matrix_delay for delay_type == 0
    zero_check = all(S_matrix_delay(:, delay_type == 0) == 0, 1);
    
    if ~all(zero_check)
        warning('Not all corresponding columns in S_matrix_delay are zeros for indices where delay_type is 0');
    end
    
    % Check delaytime_list for delay_type == 0
    tmp = [delaytime_list{delay_type == 0}];
    zero_check = tmp == 0;
    if ~all(zero_check)
        warning('Not all corresponding elements in delaytime_list are zeros for indices where delay_type is 0');
    end
    
    % Check columns in S_matrix for delay_type == 1
    zero_check = all(S_matrix(:, delay_type == 1) == 0, 1);
    if ~all(zero_check)
        warning('Not all corresponding columns in S_matrix are zeros for indices where delay_type is 1');
    end
    
    % Find non-zero delay types
    non_zero_indices = find(delay_type ~= 0);
    result_matrix = []; % Initialize an empty matrix
    
    for idx = non_zero_indices
        % Find matching columns
        match_index = find(all(S_matrix_delay(:, idx) == S_matrix));
        if ~isempty(match_index)
            result_matrix = [result_matrix, [match_index(1); repmat(idx, 1, length(match_index))]]; % Append to result_matrix
        end
    end
end