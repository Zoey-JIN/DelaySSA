function n_output = picksample(result, i, t)
    if isscalar(t)
        n_output = get_value_at_time(result, i, t);
    else
        n_output = arrayfun(@(time_point) ...
            get_value_at_time(result, i, time_point), t);
    end

    function value = get_value_at_time(result, i, t)
        index = find(result.t_values <= t, 1, 'last');
        if isempty(index)
            error('No valid time value found.');
        end
        value = result.n_values(i, index);
    end
end