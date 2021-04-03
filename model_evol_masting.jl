using Distributions
using DataFrames
using CSV
using Random
using ArgParse
using Distributed
using StatsBase
using LinearAlgebra
using BenchmarkTools
#include("/mnt/c/Users/cedri/OneDrive/Research/B1-Codes/Utility.jl")


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--write", "-w"
            action = :store_true
            help = "Will the ouput written on a csv file?"
            dest_name = "write_file"
        "--split"
            action = :store_true
            dest_name = "split_simul"
            help = "write the output on a different file for each simul"
        "--nSimul", "-S"
            arg_type = Int64
            dest_name = "n_simul"
        "--nGen", "-G"
            arg_type = Int64
            help = "Total number of generation"
            dest_name = "n_gen"
        "--print"
            arg_type = Int64
            help = "Generation from which output is saved"
            dest_name = "n_print"
        "--jPrint"
            arg_type = Int64
            help = "Number of generation between each print"
            dest_name = "jump_print"
        "--de"
            arg_type = Int64
            help = "Level of detail of output (0 => by generation, 1 => by patch, 2 => by individual)"
            dest_name = "detail"
        "--nPop"
            arg_type = Int64
            help = "Total number of individuals"
            dest_name = "n_pop"
        "--mu"
            arg_type = Float64
            help = "Mutation rate"
            dest_name = "mu"
        "--nStep"
            arg_type = Int64
            help = "Number of reproduction by generation"
            dest_name = "n_step"
        "--meanK"
            arg_type = Float64
            help = "mean value of resources distribution"
            dest_name = "mean_K"
        "--sigK"
            arg_type = Float64
            help = "variance of resources distribution"
            dest_name = "sigma_K"    
        "--thr_swi"
            arg_type = Float64
            dest_name = "theta_switching"
        "--thr_sto"
            arg_type = Float64
            dest_name = "theta_storage"
        "--a_swi" 
            arg_type = Float64
            dest_name = "a_switching"
        "--a_sto" 
            arg_type = Float64
            dest_name = "a_storage"
        "-k"
            arg_type = Float64 
            help = "Strength of increasing return on benefits (pollen effiency hypothesis)"
            dest_name = "k"  
        "--N_mid"
            arg_type = Float64 
            help = "midpoint of the sigmoid of increase of fitness as a function of total number of seeds from other plants"
        end
    return parse_args(s)
end


#function calculate_alpha(strategy)
#    if strategy == "matching"
#        alpha = 1
#    elseif strategy == "switching"
#        alpha = 

#Could also be done by describing only the number of each strategy (dictionary) 
#but it needs that all individuals of same strategy are synchronised (start in the same time). For instance they all have the same alpha at a given time.
#We describe whole population in case we want to explore other cases

#Find the name of the model if in my directory
function get_name_model()
    #Check if we are on cedric directory
    if occursin("Research/A1-Projects",pwd())==true
        name_model = SubString(pwd(),findlast("_",pwd())[1]+1)
        name_model = SubString(name_model,1,findfirst("/",name_model)[1]-1)
    else
        name_model = "res"
    return(name_model)
    end
end

#To transform parameters into a string used for naming output file: name
#Or just give short name to parameter
function get_name_file(wd::String,parameters::Dict,parameters_to_omit::Array{String,1},format::String)
    parameters_copy=copy(parameters)
    delete!(parameters_copy,"write_file")
    delete!(parameters_copy,"split_simul")
    for i in parameters_to_omit
        delete!(parameters_copy,i)
    end
    #For having only the 3 letters after underscore
    #name_file=join(["-"*arg[1:(min(length(arg),findfirst("_",arg)[1]+3))] * "=" * string(val) for (arg,val) in parameters_copy],"")
    name_file=join(["-"*arg * "=" * string(val) for (arg,val) in sort(parameters_copy)],"")
    return(wd*name_file*format)
end

#When name_file is for each simulation
function get_name_file(wd::String,parameters::Dict,parameters_to_omit::Array{String,1},format::String,i_simul::Int64)
    parameters_copy=copy(parameters)
    delete!(parameters_copy,"write_file")
    delete!(parameters_copy,"split_simul")
    for i in parameters_to_omit
        delete!(parameters_copy,i)
    end
    name_file=join(["-"*arg * "=" * string(val) for (arg,val) in sort(parameters_copy)],"")
    return(wd*get_name_model()*name_file*"-S="*string(i_simul)*format)
end

function replicator(wd::String, parameters_to_omit)
    #If on my directory and to print in Res
    parameters = parse_commandline()
    df_res= DataFrame()
    #Could be done with list comprehension?
    #Maybe with @sync @distributed
    #Make the simulations
    #We could do array of data frame and then they have directly the good size
    for i_simul in 1:parameters["n_simul"]
        res=model(parameters,i_simul)
        if parameters["write_file"] == true && parameters["split_simul"] == true
            CSV.write(get_name_file(wd,parameters,parameters_to_omit,".csv",i_simul), res)
        else
            append!(df_res,res)
        end
    end
    if parameters["write_file"] == true && parameters["split_simul"] == false
        CSV.write(get_name_file(wd,parameters,parameters_to_omit,".csv"), df_res)
    elseif write_file == false
        return(df_res)
    end
end



function update_function_ecology(cumul_population,resources,stock, theta_switching, theta_storage, a_storage, a_switching)
    cumul_population
    alpha=fill(-1.,cumul_population[3])
    alpha[1:cumul_population[1]] .= 1.
    alpha[(cumul_population[1]+1):cumul_population[2]] .=  (resources > theta_switching)*a_switching
    alpha[(cumul_population[2]+1):cumul_population[3]] .=  (stock[(cumul_population[2]+1):cumul_population[3]] .> theta_storage).*a_storage
    n_seeds = alpha .* (stock .+ resources)
    remaining_stock = (1 .- alpha) .* (stock .+ resources)
    return(alpha, n_seeds, remaining_stock)
end


function model_ecology(n_step,population, distribution_resources, k, theta_switching, theta_storage, a_storage, a_switching, print_output::Bool)
    N = sum(population)
    resources = rand(distribution_resources, n_step)
    alpha = zeros(n_step, N)
    n_seeds = zeros(n_step, N)
    stock = zeros(n_step, N)
    fitness = zeros(n_step, N)
    cumul_population = cumsum(population)
    for i in 1:(n_step-1)
        #I think we can remove the assignement
        alpha[i+1,:], n_seeds[i+1,:],stock[i+1,:] = update_function_ecology(cumul_population, resources[i], stock[i,:], theta_switching, theta_storage, a_storage, a_switching)

        #Pollen efficiency
        #First part is total of other seeds produced, second is number of seeds produced by an individual
        #fitness[i+1,:] = ((sum(n_seeds[i+1,:]) .- n_seeds[i+1,:]) .* n_seeds[i+1,:]) .^ k
        fitness[i+1,:] = n_seeds[i+1,:] ./ (1 .+ exp.(-k .* ((sum(n_seeds[i+1,:]) .- n_seeds[i+1,:]) .- N_mid)))

    end

    #fitness_total = fitness ./ n_step
    if print_output == true
        return(DataFrame(n_step = repeat(1:n_step,inner=N),
        ind = repeat(1:N,outer=n_step),
        #; for vertical concatenation
        strategy = repeat([repeat(["Matching"],population[1]); repeat(["Switching"],population[2]); repeat(["Storage"],population[3])],n_step),
        alpha = dropdims(reshape(alpha',(1,N*n_step)),dims=1),
        n_seeds= dropdims(reshape(n_seeds',(1,N*n_step)),dims=1),
        stock= dropdims(reshape(stock',(1,N*n_step)),dims=1),
        fitness= dropdims(reshape(fitness',(1,N*n_step)),dims=1)))
    else
    #We want the relative fitness (If a tree has less generation, it will have less total number of seeds but more generations)
    #We measure fitness on a given period.
    #If mean, it was dimension 1
    return([sum(fitness[:,1:cumul_population[1]]),
        sum(fitness[:,(cumul_population[1]+1):cumul_population[2]]),
        sum(fitness[:,(cumul_population[2]+1):cumul_population[3]])])
    end
end



function reproduction_WG(n_pop,fitness,mu)
    relative_fitness = fitness / sum(fitness)
    new_population_complete = sample(1:3, Weights(relative_fitness),n_pop)
    new_population_complete = mutation.(mu,new_population_complete)
    return([count(x->x==1,new_population_complete),count(x->x==2,new_population_complete),count(x->x==3,new_population_complete)])
end

#function reproduction_replicator()
    #Weighted mean using matrix
#    mean_fitness = fitness' * prop_population
    #We can use the replicator equations (if we consider deterministic evolution)
    #We need to use the discrete version of replicator equations
    #Type 1 (Need to be careful of the value of k)
    #prop_population .+ prop_population .* k .* (fitness .- mean_fitness) 
    #Using Jorgen Weibul.Evolutionary Game Theory (Type 2 on wikipedia)
#    new_population = ((k .+ fitness)./((k .+ mean_fitness))) .* prop_population
    
    #Or we can simulate everything. Just do sampling with weigth.
#end

function mutation(mu, strategy)
    if rand() < mu
        return(rand(deleteat!([1,2,3],strategy)))
    else
        return(strategy)
    end
end


function model(parameters::Dict, i_simul::Int64)
    #Set parameters from dictionaries to local variable (only to improve readability)
    for key in keys(parameters) eval(:($(Symbol(key)) = $(parameters[key]))) end
    #Set seed
    Random.seed!(i_simul)
    n_gen_printed = floor(Int,(n_gen - n_print)/jump_print)
    distribution_resources=Truncated(Normal(mean_K,sigma_K),0,Inf)

    if detail == 1
        df_res = DataFrame(i_simul=repeat([i_simul],inner=n_gen_printed*3),
        gen = repeat(n_print:jump_print:(n_gen-1),inner=3), 
        strategy = repeat(["matching","switching","storage"],outer=n_gen_printed),
        n_ind = zeros(n_gen_printed*3),
        fitness_relative = zeros(n_gen_printed*3))
    end

    #Apparently it is not that simple to generate three random number that sums to N
    population = reproduction_WG(n_pop,[1.,1.,1.],0.)
    for i in 0:(n_gen-1)
        fitness = model_ecology(n_step,population, distribution_resources, k, theta_switching, theta_storage, a_storage, a_switching, false)
        #Write output
        if i%jump_print == 0
            if detail == 1
                df_res.n_ind[(3*(floor(Int,i/jump_print)-n_print)+1):(3*(1+floor(Int,i/jump_print)-n_print))] = population
                #We replace to avoid dividing by 0 when there are no individuals of a given strategy (faster than function that check)
                df_res.fitness_relative[(3*(floor(Int,i/jump_print)-n_print)+1):(3*(1+floor(Int,i/jump_print)-n_print))] = replace!(fitness./population, NaN=>0.)
            end
        end
        population = reproduction_WG(n_pop,fitness,mu)
    end
    return(df_res)
end    

#Make evolutionary simulation
wd=pwd()*"/"
replicator(wd,["write","jump_print","detail"])


#To make ecology simulation (Need to be updated)
#parameters = parse_commandline()
#for key in keys(parameters) eval(:($(Symbol(key)) = $(parameters[key]))) end
#df_res=model_ecology(n_step,[20,20,20], Truncated(Normal(mean_K,sigma_K),0,Inf), k, theta_switching, theta_storage, a_storage, a_switching, true)
#CSV.write(wd*"valentin_trees_res.csv", df_res)

