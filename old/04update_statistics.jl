# using Pkg
# Pkg.add("DataFrames")
mutable struct AggregateStats
    n::Int32
    M1::Float32
    M2::Float32
    M3::Float32
    M4::Float32
    AggregateStats(n, m1, m2, m3, m4) = new(reinterpret(Int32, n), m1, m2, m3, m4)
end

function to_float32(x::AggregateStats)::Array{Float32}
    [reinterpret(Float32, x.n), x.M1, x.M2, x.M3, x.M4]
end

function from_float32(x)::AggregateStats
    AggregateStats(x[1], x[2], x[3], x[4], x[5])
end


function update_stats(agg::AggregateStats, x)
    n1 = agg.n
    agg.n += 1
    delta = x - agg.M1
    delta_n = delta / agg.n
    delta_n2 = delta_n * delta_n
    term1 = delta * delta_n * n1
    agg.M1 += delta_n
    agg.M4 += term1 * delta_n2 * (agg.n * agg.n - 3 * agg.n + 3) + 6 * delta_n2 * agg.M2 - 4 * delta_n * agg.M3
    agg.M3 += term1 * delta_n * (agg.n - 2) - 3 * delta_n * agg.M2
    agg.M2 += term1
    return agg
end

function combine_stats(agg1::AggregateStats, agg2::AggregateStats)
    combined = AggregateStats()
    combined.n = agg1.n + agg2.n
    
    delta = agg2.M1 - agg1.M1
    delta2 = delta * delta
    delta3 = delta * delta2
    delta4 = delta2 * delta2
    
    combined.M1 = (agg1.n * agg1.M1 + agg2.n * agg2.M1) / combined.n
    
    combined.M2 = agg1.M2 + agg2.M2 + 
                  delta2 * agg1.n * agg2.n / combined.n
    
    combined.M3 = agg1.M3 + agg2.M3 + 
                  delta3 * agg1.n * agg2.n * (agg1.n - agg2.n) / (combined.n * combined.n)
    combined.M3 += 3.0 * delta * (agg1.n * agg2.M2 - agg2.n * agg1.M2) / combined.n
    
    combined.M4 = agg1.M4 + agg2.M4 + delta4 * agg1.n * agg2.n * (agg1.n * agg1.n - agg1.n * agg2.n + agg2.n * agg2.n) / 
                  (combined.n * combined.n * combined.n)
    combined.M4 += 6.0 * delta2 * (agg1.n * agg1.n * agg2.M2 + agg2.n * agg2.n * agg1.M2) / (combined.n * combined.n) + 
                  4.0 * delta * (agg1.n * agg2.M3 - agg2.n * agg1.M3) / combined.n
    
    return combined
end

function agg_mean(agg::AggregateStats) 
    agg.M1
end

function agg_variance(agg::AggregateStats) 
    agg.M2 / (agg.n - 1.0);
end

function agg_sd(agg::AggregateStats)
    sqrt(agg_variance(agg))
end

function agg_skew(agg::AggregateStats)
    g = (sqrt(agg.n) * agg.M3) / (agg.M2^1.5)
    
    sqrt((agg.n * (agg.n - 1))) * g / (agg.n - 2)
end

function agg_kur(agg::AggregateStats)
    g = (agg.n * agg.M4) / (agg.M2 * agg.M2) - 3.0
    ((agg.n - 1) / ((agg.n - 2) * (agg.n - 3))) * ((agg.n + 1) * g + 6)
end

agg_n(agg::AggregateStats) = agg.n
