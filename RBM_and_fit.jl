using BenchmarkTools
using DataFrames
using Statistics
using Plots
using SplitApplyCombine

"""
    Y_next_step(Y, a, Rc, beta, d, average_Y_other)

This function models `Y` at next step, which could represent the fruiting levels of an individual. The update depends on various biological parameters and random environmental cues.

# Arguments:
- `Y`: The current value of `Y` .
- `a`: A coefficient that controls the strength of the interaction between the resource and population.
- `Rc`: The cost of fruiting.
- `beta`: Pollen efficiency.
- `d`: The environmental inhibition threshold; set to 0 to ignore environmental cues.
- `average_Y_other`: The average of the other individuals' `Y` values in the population.

# Returns:
- The next value of `Y`, based on the current conditions, or 1 if environmental inhibition (`d`) is active.

# Notes:
- If `Y <= 0`, the value is incremented by 1.
- Two different versions of the equation are provided: the version in the paper and a more general version.
"""
function Y_next_step(Y, a, Rc, beta, d, average_Y_other)
    if Y <= 0
        return Y + 1
    else
        ## This is the version in the paper; it is equivalent.
        # return -(rand() > d) * (a * (Rc * (average_Y_other ^ beta) + 1) - 1) * Y + 1
        ## This is a more general version that allows for Y+1 instead of just 1.
        if rand() > d
            return -((a * (Rc * (average_Y_other ^ beta) + 1) - 1) * Y + 1)
        else
            return 1
        end
    end
end

"""
    simul_RBM(pop, Rc, a, beta, d, n_timestep, detail, ind=false)

Simulates a resource-based model (RBM) for masting behavior in a population over a number of timesteps.

# Arguments:
- `pop`: A vector representing the population, with each element corresponding to an individual's `Y` value.
- `Rc`: The cost of fruiting.
- `a`: Coefficient controlling the interaction strength between population resources.
- `beta`: Pollen efficiency.
- `d`: Threshold for environmental inhibition.
- `n_timestep`: Number of timesteps for the simulation.
- `detail`: If `true`, returns a `DataFrame` with detailed output for each timestep. If `false`, returns a simple array of results.
- `ind`: If `true`, tracks individual resource levels over time; otherwise, tracks the mean resource level.

# Returns:
- Depending on `detail` and `ind` flags, the function either returns:
    - A vector of mean `Y` values per timestep (`detail=false`, `ind=false`).
    - A 2D array of individual `Y` values per timestep (`ind=true`).
    - A `DataFrame` with detailed simulation results (`detail=true`).

# Notes:
- For better performance, intermediate population results are stored in the `pop` variable before saving them to the result array.
"""
function simul_RBM(pop, Rc, a, beta, d, n_timestep, detail, ind = false)
    ## Initialize output.
    if ind == false
        res = zeros(n_timestep)
        res[1] = mean(pop)
    else
        res = [zeros(length(pop)) for i in 1:n_timestep]
        res[1] = copy(pop)
    end
    ## For later use if all the parameters should be in one dictionary.
    # for key in keys(parameters) eval(:($(Symbol(key)) = $(parameters[key]))) end
    for i in 2:n_timestep
        # @ Faster to save as a pop and then output to res, than to access res multiple times.
        pop = Y_next_step.(pop, a, Rc, beta, d, max(sum((sum(pop) / (length(pop) - 1)) .- pop), 0.))
        if ind == false
            res[i] = mean(max.(pop, 0.0))
        else
            res[i] = max.(pop, 0.0)
        end
    end
    
    if detail == false
        if ind == false
            return res
        else
            return invert(res)
        end
    elseif detail == true && ind == false
        return DataFrame(a = fill(a, length(res)), beta = fill(beta, length(res)), Rc = fill(Rc, length(res)), d = fill(d, length(res)), t = 1:n_timestep, Y = res)
    else
        error("If ind is false, detail should be false (for now)")
    end
end

pop = rand(100)

## Example of a single run showing the mean Y for the population across time.
single_run = simul_RBM(pop, 2, 3, 2, 0.5, 500, true)
plot(single_run.t, single_run.Y)

## Example of a single run showing the Y of each individual across time.
single_run = simul_RBM(pop, 2, 2, 2, 0.0, 30, false, true)
plot(single_run, legend = false)

"""
    comparison(data_real::Vector{Float64}, Rc, a, beta, d, n_timestep::Integer, fun, skip=0)

Compares the real-world data series (`data_real`) to a simulated one, using a specified comparison function.

# Arguments:
- `data_real`: A vector of real-world data for comparison.
- `Rc`, `a`, `beta`, `d`: Model parameters for the simulation.
- `n_timestep`: Number of timesteps for the simulation.
- `fun`: A function to compare the real and simulated time series.
- `skip`: How many initial timesteps of the simulation to ignore in the comparison. If unspecified, it automatically keeps the same number of points as in the data.
.

# Returns:
- The result of the comparison function `fun` applied to the real data and the simulated data (ignoring the first `skip` timesteps).
"""
function comparison(data_real::Vector{Float64}, Rc, a, beta, d, n_timestep::Integer, fun, skip = 0)
    ## If skip is unspecified, it automatically keeps the same number of points as in the data.
    if skip == 0
        skip = n_timestep - length(data_real_test) + 1
    end
    if length(data_real) > n_timestep
        error("The number of time points in the simulation is smaller than in the real data")
    end
    data_simul = simul_RBM.(Ref(pop), Rc, a, beta, d, n_timestep, false)
    return fun(data_real, data_simul[skip+1:end])
end

## A function to compare time series. Use your favorite method.
function distance_mean_squared(x, y)
    return (mean(x) - mean(y))^2
end

## Example with fake data.
data_real_test = rand(pop, 500)
comparison(data_real_test, 2, 1.4, 2, 0, 500, distance_mean_squared)

"""
    find_best_parameters(data_real, range_Rc, range_a, range_beta, range_d, n_timestep, fun, method, n_replicates, skip=0)

Finds the best model parameters for the simulated RBM by comparing it to real data across multiple replicates and parameter sets.

# Arguments:
- `data_real`: The real-world time series data for comparison.
- `range_Rc`, `range_a`, `range_beta`, `range_d`: Parameter ranges for Rc, a, beta, and d to search over.
- `n_timestep`: Number of timesteps in the simulation.
- `fun`: A function to compare the real and simulated time series (e.g., `distance_mean_squared`).
- `method`: The search method (currently supports "brute_force").
- `n_replicates`: Number of replicates to run for each parameter set.
- `skip`: Number of initial timesteps to ignore in the comparison.

# Returns:
- A `DataFrame` with the best-fitting parameters and their respective scores.
"""
function find_best_parameters(data_real, range_Rc, range_a, range_beta, range_d, n_timestep, fun, method, n_replicates, skip = 0)
    if skip > length(data_real)
        error("skip is higher than the number of timesteps in the data")
    end
    res = DataFrame(i_replicate = [], Rc = [], a = [], beta = [], d = [], score = [])
    
    if method == "brute_force"
        for Rc in range_Rc
            for a in range_a
                for beta in range_beta
                    for d in range_d
                        for i_replicate in 1:n_replicates
                            push!(res, [i_replicate, Rc, a, beta, d, comparison(data_real, Rc, a, beta, d, n_timestep, fun, skip)])
                        end
                    end
                end
            end
        end
    else
        error("Unknown method")
    end
    
    return res
end

## Example
find_best_parameters(data_real_test, 1:0.1:2, 1, 2, 0:0.1:1, 500, distance_mean_squared, "brute_force", 2)



function coefficient_of_variation(data) m = mean(data); m == 0 ? NaN : std(data) / m end
