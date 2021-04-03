using Distributions
using DataFrames
using CSV
using Random
using ArgParse
using Distributed
using StatsBase
using LinearAlgebra

function update_function_Y(R_C, a, P, Y, theta)
    if Y <= 0.
        Y_after = Y + 1
    else
        Y_after = -theta*(a*(P*R_C +1) - 1)*Y + 1
    end
    return(Y_after) 
end

function null_if_negatif!(x)
    if x < 0
        x = 0
    end
    return(x)
end

R_C = 0.32
N_gen = 30
a=2.
N=10
beta=0.73
#Need to be drawn from uniform. If superior to threshold, 1, else 0
theta = 1

Y=zeros(N_gen,N)
Y[1,:] .= rand(N).-1


pollination_rate = zeros(N_gen,N)



for i in 1:(N_gen-1)
    if i > 0
        Y[i+1,:] .= update_function_Y.(R_C, a, pollination_rate[i,:], Y[i,:], theta)
    end
    Y_plus = null_if_negatif!.(Y[i,:])
    pollination_rate[i+1,:] .= (((1/(N-1)) *(sum(Y_plus).-Y_plus))).^1.1
end

df_res = DataFrame(gen = repeat(1:N_gen,inner=N), 
individual = repeat(1:N,outer=N_gen),
Y = dropdims(reshape(Y',(N*N_gen,1)),dims=2),
pollination_rate = dropdims(reshape(pollination_rate',(N*N_gen,1)),dims=2))

wd=pwd()*"/"
name_file = wd*"valentin_trees_res.csv"
CSV.write(name_file, df_res)