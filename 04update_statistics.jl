# using Pkg
# Pkg.add("DataFrames")
using RCall
using DataFrames
using HDF5

mutable struct AggregateStats
    n::UInt32
    M1::Float32
    M2::Float32
    M3::Float32
    M4::Float32
    AggregateStats() = new(0, 0, 0, 0, 0)
end

datatype(Complex(AggregateStats))
HDF5.h5t_create(H5T_COMPOUND, sizeof(AggregateStats))

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

agg_skew(agg)
agg_kur(agg)
agg_sd(agg)

arr = [
    [1, 3.4],
    [1, 4.6],
    [2, 6.6],
    [2, 4.2],
    [2, 1],
    [3, 5],
    [3, 4],
    [3, 3],
    [3, 1],
    [3, 6],
    [3, 4]
]

R"library(data.table)"
df=rcopy(R"data.table(indices = as.integer(c(1,1,2,2,2,3,3,3,3,3,3)), vals=c(3.4, 4, 6, 4, 1, 6, 5, 3, 1, 6, 4))")





using HDF5

l2b = h5open("E:/Documentos/Downloads/GEDI02_B_2019110232351_O02005_T03190_02_001_01.h5")
groups = keys(l2b)
groups = groups[groups .!= "METADATA"]


R"library(rGEDI)"
R"l2b=readLevel2B('E:/Documentos/Downloads/GEDI02_B_2019110232351_O02005_T03190_02_001_01.h5')"
R"l2b_df = getLevel2BVPM(l2b)"
R"l2b_pavd = getLevel2BPAVDProfile(l2b)"
R"cols = grep('pavd*',colnames(l2b_pavd), value=T)"
R"mask = l2b_pavd$l2b_quality_flag == 1"
pavd = rcopy(R"l2b_pavd[mask, ..cols]")
l2b_df = rcopy(R"l2b_df")

pavd = sum.(eachrow(pavd))

deg_met_factor = 1/111320
l2b_df = l2b_df[l2b_df[:,"l2b_quality_flag"] .== 1,:]
insertcols!(l2b_df, 1, :pavd => pavd)

resolution = 3000
res_factor = resolution = deg_met_factor
xsize = Int64(360/res_factor)
ysize = Int64(180/res_factor)

l2b_df = l2b_df[:,[:latitude_bin0, :longitude_bin0, :rh100, :fhd_normal, :pai, :cover, :pavd]]
l2b_df[:,:x_index] = Int32.(floor.((l2b_df[:, :longitude_bin0] .+ 180) ./ res_factor))
l2b_df[:,:y_index] = Int32.(floor.((90 .- l2b_df[:, :latitude_bin0]) ./ res_factor))
l2b_df[:,:index] = (l2b_df[:,:y_index] .* xsize) .+ l2b_df[:,:x_index]

output_h5 = h5open("new.h5", "cw")
rh100_ds = create_dataset(output_h5, "rh100", datatype(Float32), dataspace(xsize*ysize))
agg_ds = create_dataset(output_h5, "agg", , dataspace(xsize*ysize))
output_h5["agg"] = [AggregateStats(),AggregateStats()]
rh100_ds
Int64((180/res_factor) * (90/res_factor))

keys(l2b[groups[1]])[startswith.(keys(l2b[groups[1]]), "height")]