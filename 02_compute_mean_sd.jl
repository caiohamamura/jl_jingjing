using Pkg
# Pkg.add("ThreadTools")
# using ThreadTools
# Pkg.add("HDF5")
using HDF5
# Pkg.add("ArchGDAL")
using ArchGDAL
using Printf

dataset = ArchGDAL.read("C:/Users/caioh/lsrc/r/prevFogo/cerrado_1000.tif", flags=ArchGDAL.OF_Raster | ArchGDAL.OF_ReadOnly)
ysize = ArchGDAL.height(dataset)
xsize = ArchGDAL.width(dataset)


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

cerrado_h5 = h5open("E:/Documentos/Downloads/cerrado_1000 - Copia.h5", "r+")
groups = list_recursive(cerrado_h5, "/", true, false)

lines_to_read = 100

counter = 0

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
        mask = mat_count .> 0
        
        mat_mean = mat_sum ./ mat_count
        mat_mean[.!mask] = repeat([-9999.0], sum(.!mask))
        
        mask = mat_count .> 1
        mat_var = repeat([-9999.0], size(mat_count)[1])
        mat_var[mask] = (mat_sumsq[mask] - (mat_sum[mask] .^ 2 ./ mat_count[mask])) ./ (mat_count[mask] .- 1)

        mask_no_data = mat_var .== -9999.0
        mask = mat_var .< 0
        mat_var[mask] = repeat([0], sum(mat_var .< 0))
        mat_sd = sqrt.(mat_var)
        mat_sd[mask_no_data] = repeat([-9999.0], sum(mask_no_data))
        
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
