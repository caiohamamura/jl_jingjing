# using Pkg
# Pkg.add("ArgParse")
# Pkg.add("HDF5")
# Pkg.add("ArchGDAL")
using HDF5
using ArchGDAL
using Printf
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "base_rast"
            help = "Template raster to find the xsize and ysize of the matrix"
            required = true
        "h5_file"
            help = "HDF File where data is stored"
            required = true
    end

    return parse_args(s)
end
args = parse_commandline()

function list_recursive(obj, parent="/", groups=true, datasets=true)
    result = []
    try
        items = keys(obj)
        for item in items
            try
                g = obj[item]
                group_info(g)
                if (groups)
                    push!(result, string(parent, item))
                end
                append!(result, list_recursive(g, string(parent, item, "/"), groups, datasets))
            catch
                if (datasets)
                    push!(result, string(parent, item))
                end
            end
        end
        return result
    catch
        return result
    end
end


dataset = ArchGDAL.read(args["base_rast"], flags=ArchGDAL.OF_Raster | ArchGDAL.OF_ReadOnly)
ysize = ArchGDAL.height(dataset)
xsize = ArchGDAL.width(dataset)

cerrado_h5 = h5open(args["h5_file"], "r+")
groups = list_recursive(cerrado_h5, "/", true, false)

lines_to_read = 100

counter = 0

mean(sum::Float32, count::UInt8)::Float32 = count > 0 ? sum/count : -9999.0

function sd(sumsq::Float32, sum::Float32, count::UInt8)::Float32
    result = count > 1 ? (sumsq - ((sum^2)/count))/(count - 1) :
        count == 1 ? 0 :
        -9999.0

    result > 0 ? sqrt(result) :
    result == -9999.0 ? -9999.0 : 0
end




# group_name = groups[1]
for group_name in groups
    global counter += 1
    ds_count = cerrado_h5[@sprintf("%s/%s", group_name, "count")]
    ds_sum = cerrado_h5[@sprintf("%s/%s", group_name, "sum")]
    ds_sumsq = cerrado_h5[@sprintf("%s/%s", group_name, "sumsq")]

    ds_mean_name = @sprintf("%s/%s", group_name, "mean")
    ds_sd_name = @sprintf("%s/%s", group_name, "sd")
    if HDF5.haskey(cerrado_h5,ds_mean_name)
        HDF5.delete_object(cerrado_h5, ds_mean_name)
    end
    if HDF5.haskey(cerrado_h5,ds_sd_name)
        HDF5.delete_object(cerrado_h5, ds_sd_name)
    end
    
    ds_mean = create_dataset(cerrado_h5, ds_mean_name, datatype(Float32), HDF5.dataspace(ds_count))
    ds_sd = create_dataset(cerrado_h5, ds_sd_name, datatype(Float32), HDF5.dataspace(ds_count))

    println("Processing group: ", group_name, " (", counter, " of ", size(groups)[1], ")")
    # yy = 0
    for yy in 0:lines_to_read:ysize
        perc = 100 * yy / ysize
        @printf("\rCalculating... %.2f%%",perc)
        lines_read = lines_to_read
        
        if (yy + lines_read > ysize)
            lines_read = ysize - yy
        end

        
        
        lower_bound = yy*xsize+1
        upper_bound = (yy+lines_read)*xsize

        mat_count = ds_count[lower_bound:upper_bound]
        mat_sum = ds_sum[lower_bound:upper_bound]
        mat_sumsq = ds_sumsq[lower_bound:upper_bound]

        mat_mean = mean.(mat_sum, mat_count)
        mat_sd = sd.(mat_sumsq, mat_sum, mat_count)
        
        ds_mean[lower_bound:upper_bound] = mat_mean
        ds_sd[lower_bound:upper_bound] = mat_sd
    end
    @printf("\rCalculating... %.2f%%\n",100)
end
# n = M1 = M2 = M3 = M4 = 0

# n, M1, M2, M3, M4 = update_stats(n, M1, M2, M3, M4, 3)
# n, M1, M2, M3, M4 = update_stats(n, M1, M2, M3, M4, 2)

# a_n = a_M1 = a_M2 = a_M3 = a_M4 = 0

# a_n, a_M1, a_M2, a_M3, a_M4 = update_stats(a_n, a_M1, a_M2, a_M3, a_M4, 5)
# a_n, a_M1, a_M2, a_M3, a_M4 = update_stats(a_n, a_M1, a_M2, a_M3, a_M4, 8)

# n, M1, M2, M3, M4 = combine_stats(n, M1, M2, M3, M4, a_n, a_M1, a_M2, a_M3, a_M4)
# # Kurtosis
# g = (n*M4) / (M2*M2) - 3.0
# (n-1) / ( (n-2)*(n-3) ) * ( (n+1)*g + 6 )
# # Variance
# M2/(n-1.0)
# # Skewness
# g=(sqrt(n) * M3) / (M2^1.5)
# sqrt((n*(n-1)))*g / (n-2)
