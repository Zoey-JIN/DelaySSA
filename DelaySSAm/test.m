tmax = 150;
n_initial = zeros(1, 1); 
t_initial = 0;
S_matrix = ones(1, 1); 
S_matrix_delay = -ones(1, 1); 
k = 10; 
reactant_matrix = zeros(1, 1); 
delay_type = 2; 

delaytime_list = {}; 

delaytime_list{1} = @fun_tau;

sample = 10;
A=simulation_DelaySSA(sample, tmax, n_initial, t_initial, S_matrix, S_matrix_delay, k, reactant_matrix, delay_type , delaytime_list)

T = [0:tmax];
plot_data = picksample_cells([1:1:10],A, 1, T);
mean_values = zeros(1, length(T));
% t_element = cellfun(@(x) x(2), plot_data);
% mean(t_element)
for i = 1:length(T)
    t_element = cellfun(@(x) x(i), plot_data);
    mean_values(i) = mean(t_element);
end

function delaytime_list = fun_tau()
    delaytime_list = gamrnd(7, 1);
end

