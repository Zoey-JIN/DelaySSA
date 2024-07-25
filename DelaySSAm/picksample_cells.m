function n_outputs = picksample_cells(num_cells, A, i, t)
    n_outputs = cell(length(num_cells), 1);  

    for n = 1:length(num_cells)
        idx = num_cells(n);
        result = A{idx};  
        n_outputs{idx} = picksample(result, i, t);  
    end
end