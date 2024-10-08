function result_matrix = check_delay_relative(delay_type, delaytime_list, S_matrix, S_matrix_delay)

    zero_check = all(S_matrix_delay(:, delay_type == 0) == 0, 1);
    
    if ~all(zero_check)
        warning('Not all corresponding columns in S_matrix_delay are zeros for indices where delay_type is 0');
    end
    
    tmp = [delaytime_list{delay_type == 0}];
    zero_check = tmp == 0;
    if ~all(zero_check)
        warning('Not all corresponding elements in delaytime_list are zeros for indices where delay_type is 0');
    end

    zero_check = all(S_matrix(:, delay_type == 1) == 0, 1);
    if ~all(zero_check)
        warning('Not all corresponding columns in S_matrix are zeros for indices where delay_type is 1');
    end

    non_zero_indices = find(delay_type ~= 0);

    for idx = non_zero_indices
        match_index = find(all(S_matrix_delay(:, idx) == S_matrix));
        if ~isempty(match_index)
            result_matrix = [result_matrix, [match_index(1); repmat(idx, 1, length(match_index))]]; % Append to result_matrix
        end
    end
end