function n_outputs = picksample_cells(num_cells,A, i, t)
    if nargin == 3
        num_cells = numel(A);
    end
    n_outputs = cell(num_cells, 1);  % 初始化输出cell数组

    for idx = 1:num_cells
        result = A{idx};  % 获取当前cell的内容
        n_outputs{idx} = picksample(result, i, t);  % 调用picksample函数
    end
end