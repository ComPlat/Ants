import numpy as np
import random

def pso(npar, atomic_numbers, lb, ub, npop, ngen, error_threshold, loss_fct):
    swarm = np.zeros((npop, npar))
    v = np.zeros((npop, npar))
    swarm_bests = np.zeros(npop)
    swarm_errors = np.zeros(npop)
    initial_cog = 2.5
    final_cog = 0.5
    initial_soc = 0.5
    final_soc = 2.5
    w = 0.5
    w_max = 0.9
    w_min = 0.4
    def generate_random_vec(n, lb, ub):
        v = np.zeros(n)
        for i in range(n):
            v[i] = random.uniform(lb, ub)
        return v

    for i in range(npop):
        swarm[i,:] = generate_random_vec(npar, lb, ub)
        swarm_errors[i] = loss_fct(swarm[i, :], atomic_numbers)

    swarm_bests = swarm_errors
    global_best = np.argmin(swarm_bests)
    global_best_vec = swarm[global_best,:]
    global_best_error = swarm_bests[global_best]
    swarm_bests_params = swarm
    k = 3

    def calc_neighberhood():
        neighberhood = np.zeros((npop, k), dtype = int)
        for i in range(npop):
            nneighbour = random.randint(0, k)
            neighbours = np.empty(k, dtype = int)
            neighbours[:] = -1
            for j in range(nneighbour):
                neighbours[j] = random.randint(0, npop - 1)
            neighberhood[i,:] = neighbours
        return neighberhood

    neighberhood = calc_neighberhood()
    no_improvement = 0
    iter = 1

    def correct(vec, lb, ub):
        for i in range(npar):
            if vec[i] < lb:
                vec[i] = lb
            if vec[i] > ub:
                vec[i] = ub
        return vec


    while iter < ngen:
        if (iter == 1) or (no_improvement!= 0):
            neighberhood = calc_neighberhood()
        w = w_max - iter * (w_max - w_min) / ngen
        cog = initial_cog - (initial_cog - final_cog) * (iter + 1) / ngen
        soc = initial_soc - (initial_soc - final_soc) * (iter + 1) / ngen
        for i in range(npop):
            current_neighberhood = neighberhood[i,:]
            current_neighberhood = current_neighberhood[current_neighberhood != -1]
            local_best_vec = swarm[i,:]
            if len(current_neighberhood) != 0:
                local_best = np.argmin(swarm_bests[current_neighberhood])
                local_best = current_neighberhood[local_best]
                local_best_vec = swarm[local_best,:]
            v[i,:] = w * v[i,:] + cog * random.uniform(0.0, 0.1) * (swarm_bests_params[i,:] - swarm[i,:]) + soc * random.uniform(0.0, 0.1) * (local_best_vec - swarm[i,:])
            swarm[i,:] = swarm[i,:] + v[i,:]
            swarm[i,:] = correct(swarm[i,:], lb, ub)
            error = loss_fct(swarm[i,:], atomic_numbers)
            if error < swarm_bests[i]:
                swarm_bests[i] = error
                swarm_bests_params[i,:] = swarm[i,:]
            if error < global_best_error:
                global_best = i
                global_best_vec = swarm[i,:]
                global_best_error = error
                no_improvement = 0
            else:
                no_improvement = no_improvement + 1
        iter = iter + 1
        print("Generation: " + str(iter) + " Best error: ", str(global_best_error))
        print(global_best_vec)
        if global_best_error < error_threshold:
            break
    print(global_best_vec)

    return global_best_error, global_best_vec


