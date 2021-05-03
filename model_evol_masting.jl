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
        "--nYear"
            arg_type = Int64
            help = "Total number of yeareration"
            dest_name = "n_year"
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
        "--meanK"
            arg_type = Float64
            help = "mean value of resources distribution"
            dest_name = "mean_K"
        "--sigK"
            arg_type = Float64
            help = "variance of resources distribution"
            dest_name = "sigma_K" 
        "--N_y"
            arg_type = Float64
            help = "Interval for masting for alternate strategy"
        "--thr_swit"
            arg_type = Float64
            help = "Threshold for switching"
        "--thr_stor"
            arg_type = Float64
            help = "Threshold for storage"
        "-p" 
            arg_type = Float64
            help = "proportion of resources ALWAYS allocated to reproduction"
        "--coef_a" 
            arg_type = Float64
            help = "Coefficient for functional response"
            dest_name = "a"
        "--coef_h" 
            arg_type = Float64
            help = "handling rate for functional response"
            dest_name = "h"
        "-c"
            arg_type = Float64
            help = "Survival rate of predator"
        "--beta"
            arg_type = Float64 
            help = "Strength of increasing return on benefits (pollen effiency hypothesis)"
            dest_name = "beta"  
        "--init"
            arg_type = String
            help = "composition of the initial population either matching, switching, storage or random"
            dest_name = "pop_init"  
        "--D_zero"
            arg_type = Float64
            help = "CONSTANT proportion of dead individuals by year"
        "--D_mid"
            arg_type = Float64
            help = "Amount of resources for which half of the population die"
        "--D_inc"
            arg_type = Float64
            help = "shape of the sigmoid of mortality on resources"
        "--M_zero"
            arg_type = Float64
            help = "Initial number of predator when population is extinct"
        "--gamma_zero"
            arg_type = Float64
            help = "CONSTANT death rate of seeds"
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

function replicator(wd::String, name_model, parameters_to_omit)
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


function calculate_alpha(strategy, resources, stock,thr_swit, thr_stor, N_y, year, p)
    #Matching
    if strategy == 1
        return(1.)
    #Alternate
    elseif strategy == 2
        return(p + (1-p)*(year%N_y == 0))
    #Switching
    elseif strategy == 3
        return(p + (1-p)*(resources > thr_swit))
    #Reverse switching
    elseif strategy == 4
    return(p + (1-p)*(resources < thr_swit))
    #Storage
    elseif strategy == 5
        return(p + (1-p)*(stock > thr_stor))
    else
        print("Error: Unknown strategy")
    end
end

function mutation(mu, strategy)
    if rand() < mu
        return(rand(deleteat!([1,2,3,4,5],strategy)))
    else
        return(strategy)
    end
end

function calculate_fertilised_flowers(n_total_pollen, n_pop, beta)
    return(clamp(0.1*(n_total_pollen/(n_pop-1))^beta,0,1))
end

function functional_response(n_total_seeds,a,h)
    return((a*n_total_seeds)/(1+a*h*n_total_seeds))
end

function calculate_gamma(gamma_zero, n_total_seeds, n_predator, a, h)
    if gamma_zero != 0
        gamma = gamma_zero
    else
        gamma = 1 - exp(-(functional_response(n_total_seeds,a, h)*n_predator)/n_total_seeds)
    end
    return(gamma)
end

function model(parameters::Dict, i_simul::Int64)
    #Set parameters from dictionaries to local variable (only to improve readability)
    for key in keys(parameters) eval(:($(Symbol(key)) = $(parameters[key]))) end
    #Set seed
    Random.seed!(i_simul)
    n_year_printed = floor(Int,(n_year - n_print)/jump_print)
    distribution_resources=Truncated(Normal(mean_K,sigma_K),0,Inf)
    

    if detail == 0
        df_res = DataFrame(i_simul=repeat([i_simul],inner=n_year_printed*5),
        year = repeat(n_print:jump_print:(n_year-1),inner=5),
        strategy = repeat(["matching","alternate","switching","reversed","storage"],outer=n_year_printed),
        n_ind = zeros(n_year_printed*5),
        n_predator = zeros(n_year_printed*5))
    elseif detail == 1
        df_res = DataFrame(i_simul=repeat([i_simul],inner=n_year_printed*n_pop),
        year = repeat(n_print:jump_print:(n_year-1),inner=n_pop),
        population = zeros(n_year_printed * n_pop),
        resources = zeros(n_year_printed*n_pop),
        alpha = zeros(n_year_printed*n_pop),
        n_flowers = zeros(n_year_printed*n_pop),
        n_seeds = zeros(n_year_printed*n_pop),
        n_surviving_seeds = zeros(n_year_printed*n_pop),
        stock = zeros(n_year_printed*n_pop),
        n_predator = zeros(n_year_printed*n_pop))
    end


    if pop_init == "random"
        population = sample(1:5, Weights([1.,1.,1.,1.,1.]),n_pop)
    elseif pop_init == "matching"
        population = sample(1:5, Weights([1.,0,0,0,0]),n_pop)
    elseif pop_init == "alternate"
        population = sample(1:5, Weights([0,1.,0,0,0]),n_pop)
    elseif pop_init == "switching"
        population = sample(1:5, Weights([0,0,1.,0,0]),n_pop)
    elseif pop_init == "reversed"
        population = sample(1:5, Weights([0,0,0,1.,0]),n_pop)
    elseif pop_init == "storage"
        population = sample(1:5, Weights([0.,0.,0.,0.,1.]),n_pop)
    end

    alpha = zeros(n_pop)
    n_seeds = zeros(n_pop)
    n_surviving_seeds = zeros(n_pop)
    stock = zeros(n_pop)
    fitness = zeros(n_pop)
    n_predator = M_zero

    for i in 0:(n_year-1)
        #Allocation
        resources = rand(distribution_resources)
        alpha = calculate_alpha.(population, resources, stock,thr_swit, thr_stor, N_y, n_year, p)
        stock = (1 .- alpha) .* (stock .+ resources)
        n_flowers = alpha .* (stock .+ resources)

        #Fertilisation
        n_seeds = n_surviving_seeds .+ n_flowers .* calculate_fertilised_flowers.(sum(n_flowers).-n_flowers, n_pop, beta)

        #Predation
        gamma = calculate_gamma(gamma_zero, sum(n_seeds), n_predator, a, h)
        n_surviving_seeds = (1 .- gamma) .* n_seeds
        n_predator = c * sum(n_seeds) * gamma
        
        #Calculate number of dead adult individual
        if D_zero != 0
            n_dead = D_zero
        else
            n_dead = ceil(N/(1+exp(-D_inc*(resources-D_mid))))
        end

        #Reproduction
        for i in 1:n_dead
            parent = sample(1:n_pop, Weights(n_surviving_seeds))
            #We remove the seeds that grow from the bank of seeds
            n_surviving_seeds[parent] -= 1
            splice!(population,rand(1:n_pop),mutation(mu, population[parent] ))
        end

        #Write output
        if i%jump_print == 0
            if detail == 0
                df_res.n_ind[(5*(floor(Int,i/jump_print)-n_print)+1):(5*(1+floor(Int,i/jump_print)-n_print))] = [count(x->x==1,population),count(x->x==2,population),count(x->x==3,population),count(x->x==4,population),count(x->x==5,population)]
                df_res.n_predator[(5*(floor(Int,i/jump_print)-n_print)+1):(5*(1+floor(Int,i/jump_print)-n_print))] = repeat([n_predator],5)
            elseif detail == 1
                df_res.population[(n_pop*(floor(Int,i/jump_print)-n_print)+1):(n_pop*(1+floor(Int,i/jump_print)-n_print))] = population
                df_res.resources[(n_pop*(floor(Int,i/jump_print)-n_print)+1):(n_pop*(1+floor(Int,i/jump_print)-n_print))] = repeat([resources],n_pop)
                df_res.alpha[(n_pop*(floor(Int,i/jump_print)-n_print)+1):(n_pop*(1+floor(Int,i/jump_print)-n_print))] = alpha
                df_res.n_flowers[(n_pop*(floor(Int,i/jump_print)-n_print)+1):(n_pop*(1+floor(Int,i/jump_print)-n_print))] = n_flowers
                df_res.n_seeds[(n_pop*(floor(Int,i/jump_print)-n_print)+1):(n_pop*(1+floor(Int,i/jump_print)-n_print))] = n_seeds
                df_res.n_surviving_seeds[(n_pop*(floor(Int,i/jump_print)-n_print)+1):(n_pop*(1+floor(Int,i/jump_print)-n_print))] = n_surviving_seeds
                df_res.stock[(n_pop*(floor(Int,i/jump_print)-n_print)+1):(n_pop*(1+floor(Int,i/jump_print)-n_print))] = stock
                df_res.n_predator[(n_pop*(floor(Int,i/jump_print)-n_print)+1):(n_pop*(1+floor(Int,i/jump_print)-n_print))] = repeat([n_predator],n_pop)
            end
        end
    end
    return(df_res)
end    

#Make evolutionary simulation
wd=pwd()*"/"
replicator(wd,model,["write","jump_print","detail",])


#To make ecology simulation (Need to be updated)
#parameters = parse_commandline()
#for key in keys(parameters) eval(:($(Symbol(key)) = $(parameters[key]))) end
#df_res=model_ecology(n_step,[20,20,20], Truncated(Normal(mean_K,sigma_K),0,Inf), k, thr_swit, thr_stor, a_stor, a_swit, true)
#CSV.write(wd*"valentin_trees_res.csv", df_res)

