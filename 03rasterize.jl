# using Pkg
# Pkg.add("ThreadTools")
# Pkg.add("HDF5")
# Pkg.add("ArchGDAL")
using HDF5
using ArchGDAL
using Printf

RES = 1000

base_byte_path = "../../r/prevFogo/base_byte.tif"
base_float_path = "../../r/prevFogo/base_float.tif"
cerrado_h5 = h5open("E:/Documentos/Downloads/cerrado_1000 - Copia.h5", "r")

tif = ArchGDAL.read(@sprintf("../../r/prevFogo/cerrado_%d.tif", RES), flags=ArchGDAL.OF_Raster | ArchGDAL.OF_ReadOnly)
ysize = ArchGDAL.height(tif)
xsize = ArchGDAL.width(tif)
ArchGDAL.destroy(tif)

function list_recursive(obj; parent="/", groups=true, datasets=true)
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

groups = list_recursive(cerrado_h5, groups=true, datasets=false)

lines_to_read = 100
half_ysize = Int32(ceil(ysize / 2))

counter = 0

# group_name = groups[1]
# ds_name = "mean"
for group_name in groups
    for ds_name in ["count", "mean", "sd"]
        global counter += 1
        base_rast_path = base_float_path
        if endswith(ds_name, "count")
            base_rast_path = base_byte_path
        end
        
        out_rast_name = @sprintf("%d_cerrado_%s_%s.tif", RES, group_name[2:end], ds_name)

        try
            rm(out_rast_name)
        catch
        end
        cp(base_rast_path, out_rast_name)
        base_rast = ArchGDAL.read(out_rast_name, flags=ArchGDAL.OF_Update)
        rasterband = ArchGDAL.getband(base_rast, 1)
        if endswith(ds_name, "count")
            ArchGDAL.setnodatavalue!(rasterband, 0)
        else
            ArchGDAL.setnodatavalue!(rasterband, -9999.0)
        end
    

        ds = cerrado_h5[@sprintf("%s/%s", group_name, ds_name)]
        println("Processing dataset: ", ds_name, " (", counter, " of ", 3 * size(groups), ")")
    # yy = 900
        for yy in 0:lines_to_read:ysize
            perc = 100 * yy / ysize
            @printf("\rWriting raster... %.2f%%",perc)
            lines_read = lines_to_read
        
            if (yy + lines_read > ysize)
                lines_read = ysize - yy
            end

            lower_bound = yy * xsize + 1
            upper_bound = (yy + lines_read) * xsize

            mat = ds[lower_bound:upper_bound]
            mat = reshape(mat, (xsize, lines_read))

            y_min = yy + 1
            y_max = yy + lines_read

            rasterband[:, y_min:y_max] = mat
        end
    
        HDF5.close(ds)
        ArchGDAL.destroy(rasterband)
        ArchGDAL.destroy(base_rast)
        @printf("\rWriting raster... %.2f%%\n",100)
    end 
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
