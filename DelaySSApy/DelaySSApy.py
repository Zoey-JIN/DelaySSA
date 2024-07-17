import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson
from numpy.random import uniform
from numba import njit
from scipy.special import perm
from sortedcontainers import SortedList
import bisect

def tau_element(element):
    if callable(element):
        return element()
    else:
        return element


def propensity_n(n, reactant_matrix):
    num_cols = reactant_matrix.shape[1]
    num_rows = reactant_matrix.shape[0]
    result = np.ones(num_cols)
    
    for j in range(num_cols):
        for i in range(num_rows):
            if reactant_matrix[i, j] != 0:
                prod = perm(n[i][0], reactant_matrix[i, j], exact=True)
            else:
                prod = 1.0

        result[j] = prod
        
    return result

def check_delay_relative(delay_type, delaytime_list, S_matrix, S_matrix_delay):
    delay_type = np.array(delay_type)
    delaytime_list = np.array(delaytime_list)
    S_matrix = np.array(S_matrix)
    S_matrix_delay = np.array(S_matrix_delay)
    if not set(delay_type.flatten()).issubset({0, 1, 2}):
        print("Stop: wrong with the reaction type")
        return False

    zero_check = np.all(S_matrix_delay[:, (delay_type == 0)[0]] == 0, axis=0)
    if not np.all(zero_check):
        print("Warning: Not all corresponding columns in S_matrix_delay are zeros for indices where delay_type is 0")


    zero_check = delaytime_list[(delay_type == 0)[0]] == 0
    if not np.all(zero_check):
        print("Warning: Not all corresponding elements in delaytime_list are zeros for indices where delay_type is 0")


    zero_check = np.all(S_matrix[:, (delay_type == 1)[0]] == 0, axis=0)
    if not np.all(zero_check):
        print("Warning: Not all corresponding elements in S_matrix are zeros for indices where delay_type is 1")

    non_zero_indices = np.where((delay_type != 0)[0])
    result_matrix = []

    for idx in non_zero_indices:
        match_index = np.where(np.all(S_matrix == S_matrix_delay[:, idx][:, np.newaxis], axis=0))[0]
        if match_index.size > 0:
            result_matrix.append([match_index[0], idx])

    return np.array(result_matrix).T  



class SortedList_withrecord:
    def __init__(self, value, record):
        self.value = value
        self.record = record
    def __lt__(self, other):
        return self.value < other.value


def simulate_reaction(tmax, n_initial, t_initial, S_matrix, reactant_matrix, k, fun_fr):
    n_values=[np.copy(n_initial)]
    t_values=[t_initial]
    n=np.copy(n_initial)
    t=t_initial
    while t<tmax:

        u1,u2 = np.random.rand(2)
        
        f_r = fun_fr(k, n, reactant_matrix)
        lambda_sum = np.sum(f_r[0])
        tau = -np.log(u1) / lambda_sum
        r = np.argmax(np.cumsum(f_r[0]) > u2 * lambda_sum)

        n += S_matrix[:,r].reshape(-1, 1)
        if tau<0:
            break
        t += tau
        if t < t_values[-1]:
            break

        t_values.append(t)
        n_values.append(np.copy(n))
        
    return t_values,n_values


def simulate_reaction_nextreaction(tmax, n_initial, t_initial, S_matrix, reactant_matrix, k, fun_fr):
    n_values=[np.copy(n_initial)]
    t_values=[t_initial]
    n=np.copy(n_initial)
    t=t_initial
    
    f_r = fun_fr(k, n, reactant_matrix)
    u2 = np.random.rand(S_matrix.shape[1])
    tau_vec = -np.log(u2) / f_r[0]
    
    while t < tmax:
        r = np.argmin(tau_vec)
        tau = tau_vec[r]
        n = n + S_matrix[:, r].reshape(-1, 1)
        f_r_update = fun_fr(k, n, reactant_matrix)
        
        if len(tau_vec) > 1:
            tau_vec = np.where(np.arange(len(tau_vec)) != r, 
                               f_r[0] / f_r_update[0]  * (tau_vec - tau) + tau, 
                               tau_vec)
        
        inf_nan_mask = np.isinf(tau_vec) | np.isnan(tau_vec)
        if np.any(inf_nan_mask):
            Inf_index = np.where(inf_nan_mask)
            tau_vec[Inf_index] = -np.log(np.random.rand(len(Inf_index))) / f_r_update[0][Inf_index] + tau
        
        u1 = np.random.rand()
        tau_vec[r] = 1 / f_r_update[0][r] * np.log(1 / u1) + tau
        
        f_r = f_r_update
        t = tau
        
        if t < t_values[-1]:
            break

        t_values.append(t)
        n_values.append(np.copy(n))
        
    return t_values,n_values

def simulate_reaction_modifiednextreaction(tmax, n_initial, t_initial, S_matrix, reactant_matrix, k, fun_fr):
    n_values=[np.copy(n_initial)]
    t_values=[t_initial]
    n=np.copy(n_initial)
    t=t_initial

    f_r = fun_fr(k, n, reactant_matrix)
    p_vec = np.zeros(S_matrix.shape[1])
    t_vec = np.zeros(S_matrix.shape[1])
    u2 = np.random.rand(S_matrix.shape[1])
    p_vec = np.log(1 / u2)

    while t < tmax:
        tau_vec = (p_vec - t_vec) / f_r[0]
        r = np.argmin(tau_vec)
        tau = tau_vec[r]
        n += S_matrix[:, r].reshape(-1, 1)

        t_vec += f_r[0] * tau
        u1 = np.random.rand()
        p_vec[r] += np.log(1 / u1)
        f_r = fun_fr(k, n, reactant_matrix)

        t += tau
        if t < t_values[-1]:
            break

        t_values.append(t)
        n_values.append(np.copy(n))
        
    return t_values,n_values



def simulate_reaction_delay_rejection(tmax, n_initial, t_initial, S_matrix, S_matrix_delay, reactant_matrix, k, fun_fr, delay_type, delaytime_list, delay_effect_matrix):
    n_values=[np.copy(n_initial)]
    t_values=[t_initial]
    n = np.copy(n_initial)
    t = t_initial

    Tstruct_time = SortedList()

    while t < tmax:
        u1 = np.random.rand()
        f_r = fun_fr(k, n, reactant_matrix)
        lambda_sum = np.sum(f_r[0])
        tau = -np.log(u1) / lambda_sum

        if len(Tstruct_time) > 0 and np.array(Tstruct_time[0].value) < t + tau:

            t = Tstruct_time[0].value
            r = Tstruct_time[0].record

            if delay_type[0][r] == 0:
                print("warning")
                n += S_matrix[:, r].reshape(-1, 1)
            elif delay_type[0][r] == 1:
                n += S_matrix_delay[:, r].reshape(-1, 1)
            elif delay_type[0][r] == 2:
                n += S_matrix_delay[:, r].reshape(-1, 1)
            
            Tstruct_time.pop(0)
        else:
            u2 = np.random.rand()
            r = np.argmax(np.cumsum(f_r[0]) > u2 * lambda_sum)
            t += tau

            if delay_type[0][r] == 0:
                if delay_effect_matrix == False:
                    break          
                if delay_effect_matrix.size > 0 and r in delay_effect_matrix[0, :]:
                    effect_r = delay_effect_matrix[1, np.where(delay_effect_matrix[0, :] == r)[0][0]]
                    drop_index = np.random.choice(np.where(np.array([item.record for item in Tstruct_time]) == effect_r)[0])
                    Tstruct_time.pop(drop_index)
                n += S_matrix[:, r].reshape(-1, 1)
            elif delay_type[0][r] == 1:
                add_tau = tau_element(delaytime_list[r]) + t
                add_Tmp = SortedList_withrecord(add_tau,r)
                Tstruct_time.add(add_Tmp)
            elif delay_type[0][r] == 2:
                n += S_matrix[:, r].reshape(-1, 1)
                add_tau = tau_element(delaytime_list[r]) + t
                add_Tmp = SortedList_withrecord(add_tau,r)
                Tstruct_time.add(add_Tmp)
        
        if t < t_values[-1]:
            break

        t_values.append(t)
        n_values.append(np.copy(n))
        
    return t_values,n_values



def simulate_reaction_delay_direct(tmax, n_initial, t_initial, S_matrix, S_matrix_delay, reactant_matrix, k, fun_fr, delay_type, delaytime_list, delay_effect_matrix):
    n_values=[np.copy(n_initial)]
    t_values=[t_initial]
    n = np.copy(n_initial)
    t = t_initial

    Tstruct_time = SortedList()

    while t < tmax:
        u2 = np.random.rand()
        if len(Tstruct_time) == 0:
            f_r = fun_fr(k, n, reactant_matrix)
            lambda_sum = np.sum(f_r[0])
            tau = -np.log(u2) / lambda_sum
        else:
            i = 1
            f_r = fun_fr(k, n, reactant_matrix)
            lambda_sum = np.sum(f_r[0])
            a0 = 0
            at = lambda_sum * (Tstruct_time[0].value - t)
            F = 1 - np.exp(-at)
            while F < u2:
                r = Tstruct_time[i - 1].record
                if delay_type[r] == 0:
                    print("warning")
                    n += S_matrix[:, r].reshape(-1, 1)
                elif delay_type[r] == 1:
                    n += S_matrix_delay[:, r].reshape(-1, 1)
                elif delay_type[r] == 2:
                    n += S_matrix_delay[:, r].reshape(-1, 1)
                f_r = fun_fr(k, n, reactant_matrix)
                lambda_sum = np.sum(f_r[0])
                add = lambda_sum * (Tstruct_time[i].value - Tstruct_time[i - 1].value)
                a0 = at
                at += add
                if np.isnan(add):
                    at = np.inf
                F = 1 - np.exp(-at)
                t = Tstruct_time[i - 1].value
                i += 1
                t_values.append(t)
                n_values.append(np.copy(n))
            i -= 1
            tau = -(np.log(1 - u2) + a0) / lambda_sum
            if i < len(Tstruct_time):
                Tstruct_time.remove_front(i)
            else:
                Tstruct_time = SortedList()
        
        u1 = np.random.rand()
        r = np.argmax(np.cumsum(f_r[0]) > u1 * lambda_sum)
        t += tau
        if delay_type[0][r] == 0:
            if delay_effect_matrix == False:
                break          
            if delay_effect_matrix.size > 0 and r in delay_effect_matrix[0, :]:
                effect_r = delay_effect_matrix[1, np.where(delay_effect_matrix[0, :] == r)[0][0]]
                drop_index = np.random.choice(np.where(np.array([item.record for item in Tstruct_time]) == effect_r)[0])
                Tstruct_time.pop(drop_index)
            n += S_matrix[:, r].reshape(-1, 1)
        elif delay_type[0][r] == 1:
            add_tau = tau_element(delaytime_list[r]) + t
            add_Tmp = SortedList_withrecord(add_tau,r)
            Tstruct_time.add(add_Tmp)
        elif delay_type[0][r] == 2:
            n += S_matrix[:, r].reshape(-1, 1)
            add_tau = tau_element(delaytime_list[r]) + t
            add_Tmp = SortedList_withrecord(add_tau,r)
            Tstruct_time.add(add_Tmp)

        
        if t < t_values[-1]:
            break

        t_values.append(t)
        n_values.append(np.copy(n))
        
    return t_values,n_values


def simulate_reaction_delay_modifiednextreaction(tmax, n_initial, t_initial, S_matrix, S_matrix_delay, reactant_matrix, k, fun_fr, delay_type, delaytime_list, delay_effect_matrix):
    n_values=[np.copy(n_initial) ]
    t_values=[t_initial]
    n = np.copy(n_initial) 
    t = t_initial
    p_vec = np.zeros(S_matrix.shape[1])
    t_vec = np.zeros(S_matrix.shape[1])
    Tstruct_time = SortedList()

    f_r = fun_fr(k, n, reactant_matrix)
    u2 = np.random.rand(S_matrix.shape[1])
    p_vec = -np.log(u2)

    while t < tmax:
        tau_vec = (p_vec - t_vec) / f_r[0] + t
        min_1 = np.min(tau_vec)
        r_1 = np.argmin(tau_vec)
        min_2 = Tstruct_time[0].value if Tstruct_time else np.nan
        r_2 = Tstruct_time[0].record if Tstruct_time else np.nan
        if np.isnan(min_2) or len(Tstruct_time) == 0 or min_1 < min_2:
            r = r_1
            tau = min_1 - t
            t = min_1

            if delay_type[0][r] == 0:
                if delay_effect_matrix == False:
                    break          
                if delay_effect_matrix.size > 0 and r in delay_effect_matrix[0, :]:
                    effect_r = delay_effect_matrix[1, np.where(delay_effect_matrix[0, :] == r)[0][0]]
                    drop_index = np.random.choice(np.where(np.array([item.record for item in Tstruct_time]) == effect_r)[0])
                    Tstruct_time.pop(drop_index)
                n += S_matrix[:, r].reshape(-1, 1)
            elif delay_type[0][r] == 1:
                add_tau = tau_element(delaytime_list[r]) + t
                add_Tmp = SortedList_withrecord(add_tau,r)
                Tstruct_time.add(add_Tmp)
            elif delay_type[0][r] == 2:
                n += S_matrix[:, r].reshape(-1, 1)
                add_tau = tau_element(delaytime_list[r]) + t
                add_Tmp = SortedList_withrecord(add_tau,r)
                Tstruct_time.add(add_Tmp)

            u1 = np.random.rand()
            p_vec[r] = p_vec[r] + np.log(1 / u1)
        else:
            r = r_2
            tau = min_2 - t
            t = min_2
            n = n + S_matrix_delay[:, r]
            Tstruct_time.pop(0)

        t_vec = t_vec + f_r[0] * tau

        f_r = fun_fr(k, n, reactant_matrix)

        if t < t_values[-1]:
            break

        t_values.append(t)
        n_values.append(np.copy(n))
        
    return t_values,n_values



def simulate_reaction_delay_rejection_matrix(tmax, n_initial, t_initial, S_matrix, S_matrix_delay, reactant_matrix, k, fun_fr, delay_type, delaytime_list, delay_effect_matrix):
    n_values=[np.copy(n_initial) ]
    t_values=[t_initial]
    n = np.copy(n_initial) 
    t = t_initial
    
    Tstruct = [[], []]

    while t < tmax:
        u1 = np.random.rand()
        f_r = fun_fr(k, n, reactant_matrix)
        lambda_sum = np.sum(f_r[0])
        tau = -np.log(u1) / lambda_sum

        if len(Tstruct[0]) > 0 and np.array(Tstruct[0][0]) < t + tau:

            t = Tstruct[0][0]
            r = Tstruct[1][0]
            if delay_type[0][r] == 0:
                print("warning")
                n += S_matrix[:, r].reshape(-1, 1)
            elif delay_type[0][r] == 1:
                n += S_matrix_delay[:, r].reshape(-1, 1)
            elif delay_type[0][r] == 2:
                n += S_matrix_delay[:, r].reshape(-1, 1)
            del Tstruct[0][0]
            del Tstruct[1][0] 
        else:
            u2 = np.random.rand()
            r = np.argmax(np.cumsum(f_r[0]) > u2 * lambda_sum)
            t += tau

            if delay_type[0][r] == 0:
                if delay_effect_matrix == False:
                    break          
                if delay_effect_matrix.size > 0 and r in delay_effect_matrix[0, :]:
                    effect_r = delay_effect_matrix[1, np.where(delay_effect_matrix[0, :] == r)[0][0]]
                    drop_index = np.random.choice(np.where(np.array(Tstruct[1]) == effect_r)[0])
                    del Tstruct[0][drop_index]
                    del Tstruct[1][drop_index]
                n += S_matrix[:, r].reshape(-1, 1)
            elif delay_type[0][r] == 1:
                add_tau = tau_element(delaytime_list[r]) + t
                index = bisect.bisect_left(Tstruct[0], add_tau)
                Tstruct[0].insert(index, add_tau)
                Tstruct[1].insert(index, r)
            elif delay_type[0][r] == 2:
                n += S_matrix[:, r].reshape(-1, 1)
                add_tau = tau_element(delaytime_list[r]) + t
                index = bisect.bisect_left(Tstruct[0], add_tau)
                Tstruct[0].insert(index, add_tau)
                Tstruct[1].insert(index, r)
        
        if t < t_values[-1]:
            break

        t_values.append(t)
        n_values.append(np.copy(n))
        
    return t_values,n_values



def simulate_reaction_delay_direct_matrix(tmax, n_initial, t_initial, S_matrix, S_matrix_delay, reactant_matrix, k, fun_fr, delay_type, delaytime_list, delay_effect_matrix):
    n_values=[np.copy(n_initial) ]
    t_values=[t_initial]
    n = np.copy(n_initial) 
    t = t_initial

    Tstruct = [[], []]

    while t < tmax:
        u2 = np.random.rand(5)
        if len(Tstruct[0]) == 0:
            f_r = fun_fr(k, n, reactant_matrix)
            lambda_sum = np.sum(f_r[0])
            tau = -np.log(u2) / lambda_sum
        else:
            i = 1
            f_r = fun_fr(k, n, reactant_matrix)
            lambda_sum = np.sum(f_r[0])
            a0 = 0
            at = lambda_sum * (Tstruct[0][0] - t)
            F = 1 - np.exp(-at)
            while F < u2:
                r = Tstruct[1][i - 1]
                if delay_type[r] == 0:
                    print("warning")
                    n += S_matrix[:, r].reshape(-1, 1)
                elif delay_type[r] == 1:
                    n += S_matrix_delay[:, r].reshape(-1, 1)
                elif delay_type[r] == 2:
                    n += S_matrix_delay[:, r].reshape(-1, 1)
                f_r = fun_fr(k, n, reactant_matrix)
                lambda_sum = np.sum(f_r[0])
                add = lambda_sum * (Tstruct[0][i] - Tstruct[0][i - 1])
                a0 = at
                at += add
                if np.isnan(add):
                    at = np.inf
                F = 1 - np.exp(-at)
                t = Tstruct[0][i - 1]
                i += 1
                t_values.append(t)
                n_values.append(np.copy(n))
            i -= 1
            tau = -(np.log(1 - u2) + a0) / lambda_sum
            if i < len(Tstruct[0]):
                Tstruct = [Tstruct[0][(i + 1):], Tstruct[1][(i + 1):]]
            else:
                Tstruct = [[], []]
        
        u1 = np.random.rand()
        r = np.argmax(np.cumsum(f_r[0]) > u1 * lambda_sum)
        t += tau

        if delay_type[0][r] == 0:
            if delay_effect_matrix == False:
                break          
            if delay_effect_matrix.size > 0 and r in delay_effect_matrix[0, :]:
                effect_r = delay_effect_matrix[1, np.where(delay_effect_matrix[0, :] == r)[0][0]]
                drop_index = np.random.choice(np.where(np.array(Tstruct[1]) == effect_r)[0])
                del Tstruct[0][drop_index]
                del Tstruct[1][drop_index]
            n += S_matrix[:, r].reshape(-1, 1)
        elif delay_type[0][r] == 1:
            add_tau = tau_element(delaytime_list[r]) + t
            index = bisect.bisect_left(Tstruct[0], add_tau)
            Tstruct[0].insert(index, add_tau)
            Tstruct[1].insert(index, r)
        elif delay_type[0][r] == 2:
            n += S_matrix[:, r].reshape(-1, 1)
            add_tau = tau_element(delaytime_list[r]) + t
            index = bisect.bisect_left(Tstruct[0], add_tau)
            Tstruct[0].insert(index, add_tau)
            Tstruct[1].insert(index, r)

        
        if t < t_values[-1]:
            break

        t_values.append(t)
        n_values.append(np.copy(n))
        
    return t_values,n_values


def simulate_reaction_delay_modifiednextreaction_matrix(tmax, n_initial, t_initial, S_matrix, S_matrix_delay, reactant_matrix, k, fun_fr, delay_type, delaytime_list, delay_effect_matrix):
    n_values=[np.copy(n_initial) ]
    t_values=[t_initial]
    n = np.copy(n_initial) 
    t = t_initial
    p_vec = np.zeros(S_matrix.shape[1])
    t_vec = np.zeros(S_matrix.shape[1])

    Tstruct = [[], []]

    f_r = fun_fr(k, n, reactant_matrix)
    u2 = np.random.rand(S_matrix.shape[1])
    p_vec = -np.log(u2)

    while t < tmax:
        tau_vec = (p_vec - t_vec) / f_r[0] + t

        min_1 = np.min(tau_vec)
        r_1 = np.argmin(tau_vec)
        min_2 = Tstruct[0][0] if Tstruct[0] else np.nan
        r_2 = Tstruct[1][0] if Tstruct[1] else np.nan
        if np.isnan(min_2) or len(Tstruct[0]) == 0 or min_1 < min_2:
            r = r_1
            tau = min_1 - t
            t = min_1

            if delay_type[0][r] == 0:
                if delay_effect_matrix == False:
                    break          
                if delay_effect_matrix.size > 0 and r in delay_effect_matrix[0, :]:
                    effect_r = delay_effect_matrix[1, np.where(delay_effect_matrix[0, :] == r)[0][0]]
                    drop_index = np.random.choice(np.where(np.array(Tstruct[1]) == effect_r)[0])
                    del Tstruct[0][drop_index]
                    del Tstruct[1][drop_index]
                n += S_matrix[:, r].reshape(-1, 1)
            elif delay_type[0][r] == 1:
                add_tau = tau_element(delaytime_list[r]) + t
                index = bisect.bisect_left(Tstruct[0], add_tau)
                Tstruct[0].insert(index, add_tau)
                Tstruct[1].insert(index, r)
            elif delay_type[0][r] == 2:
                n += S_matrix[:, r].reshape(-1, 1)
                add_tau = tau_element(delaytime_list[r]) + t
                index = bisect.bisect_left(Tstruct[0], add_tau)
                Tstruct[0].insert(index, add_tau)
                Tstruct[1].insert(index, r)

            u1 = np.random.rand()
            p_vec[r] = p_vec[r] + np.log(1 / u1)
        else:
            r = r_2
            tau = min_2 - t
            t = min_2
            n = n + S_matrix_delay[:, r].reshape(-1, 1)
            del Tstruct[0][0]
            del Tstruct[1][0] 

        t_vec = t_vec + f_r[0] * tau

        f_r = fun_fr(k, n, reactant_matrix)

        if t < t_values[-1]:
            break

        t_values.append(t)
        n_values.append(np.copy(n))
        
    return t_values,n_values





def simulation_delay_ssa(algorithm, sample_size, tmax, n_initial, t_initial, S_matrix, S_matrix_delay, reactant_matrix ,k, delay_type, delaytime_list):
    if all(v is not None for v in [delay_type, delaytime_list, S_matrix_delay]):
        delay_effect_matrix = check_delay_relative(delay_type, delaytime_list, S_matrix, S_matrix_delay)
    results = []
    def fun_fr(k, n, reactant_matrix):
        if callable(k):
            k_mask = k(n)
        else:
            k_mask = k
        return k_mask * propensity_n(n, reactant_matrix)
    for _ in range(sample_size):
        if algorithm == "DelayRejection":
            result = simulate_reaction_delay_rejection(tmax, n_initial, t_initial, S_matrix, S_matrix_delay,reactant_matrix, k, fun_fr, delay_type, delaytime_list, delay_effect_matrix)
        elif algorithm == "DelayMNR":
            result = simulate_reaction_delay_modifiednextreaction(tmax, n_initial, t_initial, S_matrix, S_matrix_delay,reactant_matrix, k, fun_fr, delay_type, delaytime_list, delay_effect_matrix)
        elif algorithm == "DelayDirect":
            result = simulate_reaction_delay_direct(tmax, n_initial, t_initial, S_matrix, S_matrix_delay,reactant_matrix, k, fun_fr, delay_type, delaytime_list, delay_effect_matrix)
        elif algorithm == "Direct":
            result = simulate_reaction(tmax, n_initial, t_initial, S_matrix,reactant_matrix, k, fun_fr)
        elif algorithm == "MNR":
            result = simulate_reaction_modifiednextreaction(tmax, n_initial, t_initial, S_matrix,reactant_matrix, k, fun_fr)
        elif algorithm == "NR":
            result = simulate_reaction_nextreaction(tmax, n_initial, t_initial, S_matrix,reactant_matrix, k, fun_fr)
        else:
            raise ValueError("Error: No such algorithm!")
        results.append(result)
    return results
    
def main():
    sample=400
    n_initial = np.array([0, 0]).reshape((2, 1)) 
    t_initial = 0

    S_matrix = np.array([1, 0, 0, -1]).reshape((2, 2)) 
    S_matrix_delay = np.array([-1, 0, 1, 0]).reshape((2, 2))

    
    def k(n):
        return np.array([1 / (1 + (n[1])**2), 1 / (1 + n[1])]).T
    
    reactant_matrix = np.array([[0, 0], [0, 1]])  

    delay_type = np.array([2, 0]).reshape((1, 2))  

    tau = 20
    delaytime_list = [tau, 0]  
    tmax = 400
    results = simulation_delay_ssa(algorithm="DelayRejection", sample_size=2, tmax=tmax, t_initial=t_initial, n_initial= n_initial, S_matrix=S_matrix, S_matrix_delay=S_matrix_delay, k=k, reactant_matrix=reactant_matrix, delay_type=delay_type, delaytime_list=delaytime_list)
    results[0][0]
    results[0]
    A=picksample_cells([0,1],results,1,range(10))
    A=picksample(results[0],1,range(10))
if __name__ == "__main__":
    main()




def picksample(result, i, t):
    def get_value_at_time(result, i, t):
        index = np.where(np.array(result[0]) <= t)[0]
        if len(index) == 0:
            raise ValueError('No valid time value found.')
        index = index[-1]  
        return result[1][index][i]

    if np.isscalar(t):
        n_output = get_value_at_time(result, i, t)
    else:
        n_output = np.array([get_value_at_time(result, i, time_point) for time_point in t])
        n_output = np.hstack(n_output)
    return n_output

def picksample_cells(num_cells, result, i, t):
    n_outputs = [None] * len(num_cells)
    for idx in range(len(num_cells)):
        result_cell = result[idx]  
        n_outputs[idx] = picksample(result_cell, i, t)  
    return n_outputs

# def analyze_simulation(num_cells, results, i, t):
#     A=picksample_cells(num_cells,results,i,t)
#     elements = [array[i, j] for array in A]

time_num=[]
mean_time_num=[]
for i in range(0,100):
    time_num.append([vector[i] for vector in num])
    mean_time_num.append(np.mean([vector[i] for vector in num]))
mean_time_num[0:0]=[0]

len(list(range(0,tmax+1)))
len(list(np.arange(0, 10, 0.1)))
len(mean_time_num)
# plt.plot(list(range(0,tmax+1)), mean_time_num, color='blue')
plt.plot(list(np.arange(0, 10.1, 0.1)), mean_time_num, color='blue')
plt.show()



vector = [vector[-1] for vector in return_n]
unique_values, frequencies = np.unique(vector, return_counts=True)
probabilities = frequencies / len(vector)
plt.plot(unique_values, probabilities, color='blue')
plt.show()






