# using Pkg
# Pkg.add("HDF5")
# Pkg.add("ArchGDAL")
# Pkg.add("ArgParse")
using HDF5
using ArchGDAL
using Printf
using ArgParse
export main

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
    
cerrado_h5 = h5open(ARGS[1], "r")
base_byte_path = ARGS[2]
base_float_path = ARGS[3]
out_root = ARGS[4]


ArchGDAL.read(base_byte_path, flags=ArchGDAL.OF_Raster | ArchGDAL.OF_ReadOnly) do tif
    ysize = ArchGDAL.height(tif)
    xsize = ArchGDAL.width(tif)
    ArchGDAL.destroy(tif)
end

groups = list_recursive(cerrado_h5, groups=true, datasets=false)

lines_to_read = 100
half_ysize = Int32(ceil(ysize / 2))

counter = 0

for group_name in groups
    for ds_name in ["count", "mean", "sd"]
        global counter += 1
        base_rast_path = base_float_path
        if endswith(ds_name, "count")
            base_rast_path = base_byte_path
        end
            
        out_rast_name = @sprintf("%s_%s_%s.tif", out_root, group_name[2:end], ds_name)

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
