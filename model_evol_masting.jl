using Distributions
using DataFrames
using CSV
using Random
using ArgParse
using Distributed
using StatsBase
using LinearAlgebra
using BenchmarkTools
using Distances 
#include("/mnt/c/Users/cedri/OneDrive/Research/B1-Codes/Utility.jl")
include("generate_thr_storage.jl")

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
            default = -1.
        "--thr_rev"
            arg_type = Float64
            help = "Threshold for reversed switching"
            default = -1.
        "--thr_stor"
            arg_type = Float64
            help = "Threshold for storage"
            default = -1.
        "-p" 
            arg_type = Float64
            help = "proportion of resources ALWAYS allocated to reproduction"
        #Set up how quickly maximum predation is reached ()
        "--coef_a" 
            arg_type = Float64
            help = "Coefficient for functional response"
            dest_name = "a"
        #Maximum amount 1 - exp(-1/h). Set how maximum predation pressure 
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
        "--dis"
            arg_type = Float64
            help = "Parameter for dispersal. Average distance made with exponential function for dispersal"
            default = -1.
        "--str_sel"
            arg_type = Float64
            help = "Strength of selection"
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

function abbreviate_name(x)
    if typeof(x) == Float64
        return(round(x, digits=2))
    else
        return(x)
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
    name_file=join(["-"*arg * "=" * string(abbreviate_name(val)) for (arg,val) in sort(parameters_copy)],"")
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
    if parameters["thr_stor"] == -1
        parameters["thr_stor"] = estimate_thr_storage(parameters["N_y"],Truncated(Normal(parameters["mean_K"],parameters["sigma_K"]),0,Inf),"simulate")
    end
    if parameters["thr_swit"] == -1
        parameters["thr_swit"] = quantile(Truncated(Normal(parameters["mean_K"],parameters["sigma_K"]),0,Inf), 1 - (1 / parameters["N_y"]))
    end
    if parameters["thr_rev"] == -1
        parameters["thr_rev"] = quantile(Truncated(Normal(parameters["mean_K"],parameters["sigma_K"]),0,Inf), 1 / parameters["N_y"])
    end

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


function calculate_alpha(strategy, resources, stock,thr_swit, thr_rev, thr_stor, N_y, age, p)
    #Matching
    if strategy == 1
        return(1.)
    #Alternate
    elseif strategy == 2
        return(p + (1-p)*(age%N_y == 0))
    #Switching
    elseif strategy == 3
        return(p + (1-p)*(resources > thr_swit))
    #Reverse switching
    elseif strategy == 4
        return(p + (1-p)*(resources < thr_rev))
    #Storage
    elseif strategy == 5
        return(p + (1-p)*(stock > thr_stor))
    else
        println("Error: Unknown strategy")
    end
end

function mutation(mu, strategy)
    if rand() < mu
        return(rand(deleteat!([1,2,3,4,5],strategy)))
        #To consider only two strategies
        # if strategy == 3
        #     return(5)
        # else
        #     return(3)
        # end
    else
        return(strategy)
    end
end

function calculate_fertilised_flowers(n_total_pollen, n_pop, beta)
    return(clamp((0.1*n_total_pollen/(n_pop-1))^beta,0,1))
end

function functional_response(n_total_seeds,a,h)
    return((a*n_total_seeds)/(1+a*h*n_total_seeds))
end

function calculate_gamma(gamma_zero, n_total_seeds, n_predator, a, h)
    if n_total_seeds == 0
        gamma = 0.
    elseif gamma_zero != 0
        gamma = gamma_zero
    else
        gamma = 1 - exp(-(functional_response(n_total_seeds,a, h)*n_predator)/n_total_seeds)
    end
    return(gamma)
end

function create_map(N)
    # x=repeat([collect(1:2:(sqrt(N)*2)) ; collect(2:2:(sqrt(N)*2))],Integer(sqrt(N)/2))
    # y=repeat(collect(1:sqrt(N)),inner=Integer(sqrt(N)))

    x=repeat(collect(1:sqrt(N)),Integer(sqrt(N)))
    y=repeat(collect(1:sqrt(N)),inner=Integer(sqrt(N)))

    return(tuple.(x,y))
end


#We use a thin tailed for now. Note that it means that at long distance, almost no seeds from far. If we want that, we need thin tailed.
#Dispersal capacity = mean distance travelled with exponential
function dispersal_function(distance,dispersal_capacity)
    (1/(2*dispersal_capacity)) * exp(-abs(distance/dispersal_capacity))
end


function model(parameters::Dict, i_simul::Int64)
    #Set parameters from dictionaries to local variable (only to improve readability)
    for key in keys(parameters) eval(:($(Symbol(key)) = $(parameters[key]))) end
    #Set seed
    Random.seed!(i_simul)
    year_printed = n_print:jump_print:(n_year)
    n_year_printed = length(year_printed)
    distribution_resources=Truncated(Normal(mean_K,sigma_K),0,Inf)
    #distribution_resources=Truncated(Exponential(sigma_K),0,Inf)

    

    if detail == 0
        df_res = DataFrame(i_simul=repeat([i_simul],inner=n_year_printed*5),
        year = repeat(year_printed,inner=5),
        strategy = repeat(["matching","alternate","switching","reversed","storage"],outer=n_year_printed),
        resources = zeros(n_year_printed*5),
        n_ind = zeros(n_year_printed*5),
        n_predator = zeros(n_year_printed*5),
        n_dead = zeros(n_year_printed*5),
        gamma = zeros(n_year_printed*5),
        total_seeds = zeros(n_year_printed*5),
        total_flowers = zeros(n_year_printed*5),
        fertilisation_rate = zeros(n_year_printed*5))
    elseif detail == 1
        df_res = DataFrame(i_simul=repeat([i_simul],inner=n_year_printed*n_pop),
        year = repeat(year_printed,inner=n_pop),
        ID = repeat(1:1:n_pop,outer=n_year_printed),
        strategy = zeros(n_year_printed * n_pop),
        resources = zeros(n_year_printed*n_pop),
        age = zeros(n_year_printed*n_pop),
        alpha = zeros(n_year_printed*n_pop),
        n_flowers = zeros(n_year_printed*n_pop),
        fertilisation_rate = zeros(n_year_printed*n_pop),
        n_seeds = zeros(n_year_printed*n_pop),
        stock = zeros(n_year_printed*n_pop),
        n_predator = zeros(n_year_printed*n_pop),
        n_dead = zeros(n_year_printed*n_pop))
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
    fertilisation_rate = zeros(n_pop) 
    n_seeds = zeros(n_pop)
    n_surviving_seeds = zeros(n_pop)
    #If unlimited dispersal, we count only the total for each strategy
    if dis == -1
        bank_seeds = zeros(5)
    else
        bank_seeds = zeros(5, n_pop)
        coord_pop = create_map(n_pop)
        #Matrix R[i,j] gives distance between 
        distance_matrix = pairwise(Euclidean(), coord_pop, coord_pop)

        # foreach(i -> distance_matrix[i, i] = 100, 1:20)
        # distance_matrix=distance_matrix .- 2

        dispersal_matrix = dispersal_function.(distance_matrix,dis)

        #Proportion of seeds that each individual (row) send to patch (column)
        # normalised_dispersal_matrix= dispersal_matrix ./ sum(dispersal_matrix,dims = 2)
    end
    stock = zeros(n_pop)
    fitness = zeros(n_pop)  
    n_predator = M_zero
    age = zeros(Int64 ,n_pop)

    for i in 1:(n_year)
        #Allocation
        age = age .+ 1
        resources = rand(distribution_resources)
        alpha = calculate_alpha.(population, resources, stock,thr_swit, thr_rev,thr_stor, N_y, age, p)
        n_flowers = (alpha .* (stock .+ resources))
        stock = (1 .- alpha) .* (stock .+ resources)
          

        #Fertilisation
        fertilisation_rate = calculate_fertilised_flowers.(sum(n_flowers).-n_flowers, n_pop, beta)  
        n_seeds = (n_flowers .* fertilisation_rate )

        if dis == -1
            bank_seeds = bank_seeds .+ [sum((population .== 1) .* n_seeds),sum((population .== 2) .* n_seeds),sum((population .== 3) .* n_seeds),sum((population .== 4) .* n_seeds),sum((population .== 5) .* n_seeds)]
        else
            seeds_matrix = dispersal_matrix .* n_seeds
            bank_seeds = bank_seeds .+ [sum(seeds_matrix[population .== 1,:], dims = 1);
            sum(seeds_matrix[population .== 2,:], dims = 1);
            sum(seeds_matrix[population .== 3,:], dims = 1);
            sum(seeds_matrix[population .== 4,:], dims = 1);
            sum(seeds_matrix[population .== 5,:], dims = 1)]
        end

        #Predation
        gamma = calculate_gamma(gamma_zero, sum(bank_seeds), n_predator, a, h)
        n_predator = c * sum(bank_seeds) * gamma
        bank_seeds = (1 .- gamma) .* bank_seeds

        #Version with dispersal (We need to describe the seeds of each individual)
        #WARNING: This would work only without bank seeds and surviving seeds because it considers that seeds of dead disappear
        # n_seeds = n_flowers .* fertilisation_rate 
        # #Predation
        # gamma = calculate_gamma(gamma_zero, sum(n_seeds), n_predator, a, h)
        # n_predator = c * sum(n_seeds) * gamma
        # n_seeds = (1 .- gamma) .* n_seeds
        

        
        #Calculate number of dead adult individual
        if D_zero != 0
            n_dead = Int(ceil(D_zero*n_pop))
        else
            n_dead = Int(ceil(n_pop-n_pop/(1+exp(-D_inc*(resources-D_mid)))))
        end


        #Write output
        if i%jump_print == 0 && i >= n_print
            if detail == 0
                #interval = (5*(floor(Int,(i-n_print)/jump_print))+1):(5*(1+floor(Int,(i-n_print)/jump_print)))
                interval = (5*(floor(Int,(i-n_print)/jump_print))+1):5*(1+floor(Int,(i-n_print)/jump_print))
                df_res.n_ind[interval] = [count(x->x==1,population),count(x->x==2,population),count(x->x==3,population),count(x->x==4,population),count(x->x==5,population)]
                df_res.n_predator[interval] = repeat([n_predator],5)
                df_res.resources[interval] = repeat([resources],5)
                df_res.n_dead[interval] = repeat([n_dead],5)
                df_res.gamma[interval] = repeat([gamma],5)
                if dis == -1
                    df_res.total_seeds[interval] = bank_seeds
                else
                    df_res.total_seeds[interval] = sum(bank_seeds,dims=2)
                end
                df_res.total_flowers[interval] = [sum((population .== 1) .* n_flowers),sum((population .== 2) .* n_flowers),sum((population .== 3) .* n_flowers),sum((population .== 4) .* n_flowers),sum((population .== 5) .* n_flowers)]
                df_res.fertilisation_rate[interval] = [mean(fertilisation_rate[population .== 1]),mean(fertilisation_rate[population .== 2]),mean(fertilisation_rate[population .== 3]),mean(fertilisation_rate[population .== 4]),mean(fertilisation_rate[population .== 5])]
            elseif detail == 1
                interval = (n_pop*(floor(Int,(i-n_print)/jump_print))+1):n_pop*(1+floor(Int,(i-n_print)/jump_print))
                df_res.strategy[interval] = population
                df_res.resources[interval] = repeat([resources],n_pop)
                df_res.age[interval] = age
                df_res.alpha[interval] = alpha
                df_res.n_flowers[interval] = n_flowers
                df_res.fertilisation_rate[interval] = fertilisation_rate
                df_res.n_seeds[interval] = n_seeds
                df_res.stock[interval] = stock
                df_res.n_predator[interval] = repeat([n_predator],n_pop)
                df_res.n_dead[interval] = repeat([n_dead],n_pop)
            end
        end


        #Reproduction with dispersal
        #WARNING: This would work only without bank seeds and surviving seeds because it considers that seeds of dead disappear
        # deads = sample(1:n_pop,n_dead,replace=false)
        # for i_dead in deads
        #     #We draw random neighbours besides the dead individual
        #     distance_with_dead = euclidean.(Ref(coord_pop[i_dead]), coord_pop)
        #     n_seeds_dispersed_to_empty_patch = dispersal_function.(distance_with_dead,dis) .* n_seeds
        #     parent = sample(population, Weights(n_seeds_dispersed_to_empty_patch))
        #     splice!(population,i_dead,mutation(mu, parent ))
        #     stock[i_dead] = 0
        #     age[i_dead] = 0
        #     #n_surviving_seeds[dead] = 0 
        # end

        #Reproduction
        deads = sample(1:n_pop,n_dead,replace=false)
        for i_dead in deads
            if dis == -1
                #parent = sample(1:5, Weights(bank_seeds))
                parent = sample(1:5, Weights(bank_seeds .^ str_sel))

            else
                parent = sample(1:5, Weights(bank_seeds[:,i_dead] .^ str_sel))
            end
            splice!(population,i_dead,mutation(mu, parent ))
            stock[i_dead] = 0
            age[i_dead] = 0
            #n_surviving_seeds[dead] = 0 
        end



        #Migration for predators
        if n_predator == 0
            n_predator = M_zero
        end
    end
    return(df_res)
end    

#Make evolutionary simulation
wd=pwd()*"/"
replicator(wd,model,["write","jump_print","mu", "M_zero", "n_print","p","coef_a","coef_h","D_inc","thr_swi","thr_stor","mean_K","sigma_K"])


#To make ecology simulation (Need to be updated)
#parameters = parse_commandline()
#for key in keys(parameters) eval(:($(Symbol(key)) = $(parameters[key]))) end
#df_res=model_ecology(n_step,[20,20,20], Truncated(Normal(mean_K,sigma_K),0,Inf), k, thr_swit, thr_stor, a_stor, a_swit, true)
#CSV.write(wd*"valentin_trees_res.csv", df_res)

