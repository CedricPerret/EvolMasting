# using Distributions
# using DataFrames
# using CSV
# using Random
# using StatsBase
# using LinearAlgebra
# using BenchmarkTools

function mean_section_distribution(distribution, min, max)
    mean(collect(min:(max-min)/1000:max), Weights(pdf.(distribution_resources,collect(min:(max-min)/1000:max))))
end

function calculate_mean_stock(distribution,threshold, year)
    stock = 0
    for i in 1:(year-1)
        stock=stock + mean_section_distribution(distribution_resources,0,threshold-stock)
    end
    return(stock)
end



function calculate_proba_year(distribution_resources,threshold, year)
    proba_year = 0
    stock_for_each_year = calculate_mean_stock.(distribution_resources,threshold,1:year)
    proba_masting_at_given_year = 1 .- cdf.(distribution_resources,threshold .- stock_for_each_year)
    if year == 1
        proba_year = proba_masting_at_given_year[1]
    else
        proba_year = proba_masting_at_given_year[year] * prod((1 .- proba_masting_at_given_year)[1:(year-1)])
    end
    return(proba_year)
end

function calculate_mean_interval(distribution_resources, threshold)
    years_considered = collect(1:(ceil(Int64,(threshold/mean(distribution_resources))*3)))
    proba_for_each_year = calculate_proba_year.(distribution_resources,threshold,years_considered)
    return(sum(proba_for_each_year .* years_considered)/sum(proba_for_each_year))
end


function simulate_mean_interval(distribution_resources, threshold, n_simul)
    average_year_masting = 0
    for i in 1:n_simul
        year_masting = 0
        year_count = 1
        stock = 0
        while year_masting == 0
            stock = stock + rand(distribution_resources)
            if stock > threshold
                year_masting = year_count
            else
                year_count = year_count + 1
            end
        end
        average_year_masting = average_year_masting + year_masting/n_simul
    end 
    return(average_year_masting)
end

function estimate_thr_storage(mean_interval, distribution_resources, method)
    range = collect(0.01:0.01:(mean_interval*mean(distribution_resources)*1.1))
    if method == "calculate"
        res = calculate_mean_interval.(distribution_resources,range)
    elseif method == "simulate"
        res = simulate_mean_interval.(distribution_resources,range,100000)
    else
        println("Unknown method")
        return()
    end
    return(findmin(abs.(res .- mean_interval))[2] * 0.01)
end

# estimate_thr_storage(1.5,distribution_resources,"simulate")

# df = DataFrame(N_y = collect(2:0.5:5),thr_stor = estimate_thr_storage.(collect(2:0.5:5),distribution_resources,"simulate"))

# CSV.write("/mnt/c/Users/cedri/OneDrive/Research/A1-Projects/2021_EvolMasting/Code/thr_storage.csv",df)