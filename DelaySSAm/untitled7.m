tmax = 150;
n_initial = zeros(1, 1); 
t_initial = 0;
S_matrix = ones(1, 1); 
S_matrix_delay = -ones(1, 1); 
k = 10; 
reactant_matrix = zeros(1, 1); 
delay_type = 2; 

delaytime_list = cell(1, 1); 
delaytime_list{1} = gamrnd(7, 1); 
sample = 10;
A=simulation_DelaySSA( "DelayRejection", sample, tmax, n_initial, t_initial, S_matrix, S_matrix_delay, k, reactant_matrix, delay_type , delaytime_list)

direct(tmax, n_initial, t_initial, S_matrix, k, fun_fr)
A=simulation_SSA( "Direct", sample, tmax, n_initial, t_initial, S_matrix, k, reactant_matrix)

picksample(A{1}, 1, 3)
