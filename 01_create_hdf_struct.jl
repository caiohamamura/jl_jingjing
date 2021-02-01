# using Pkg
# Pkg.add("HDF5")
# Pkg.add("ArchGDAL")
# Pkg.add("ArgParse")
using HDF5
using Printf
using Glob

include("old/04update_statistics.jl")

function slice_set(ds, index, value) 
    ds[index,:] = value
end

function slice_get(ds, index)
    ds[index,:]
end

# args = ["E:/Documentos/Downloads/", "E:/Documentos/out.h5"]
function main(args)
    out_h5 = HDF5.h5open(args[2], "w")
    
    dcpl = HDF5.create_property(HDF5.H5P_DATASET_CREATE)
    HDF5.h5p_set_deflate(dcpl, 5)
    
    lat_to_met_factor = 1 / 110540
    lon_to_met_factor = 1 / 111320
    resolution = 3000 # meters
    x_res = resolution * lon_to_met_factor
    y_res = resolution * lat_to_met_factor
    xmax = Int64(floor(360/x_res))
    ymax = Int64(floor(180/y_res))
    
    cols = ["rh100", "pai", "pavd", "cover", "fhd_normal"]
    out_ds = Dict(col => HDF5.create_dataset(
        out_h5, 
        @sprintf("agg_%s", col), 
        datatype(Float32), 
        dataspace((xmax*ymax,5));
        dcpl=dcpl
        ) for col in cols)
        
        
    
    list_h5 = glob("*02_B*.h5",args[1])
    total_counter = 0

    n_files = size(list_h5)[1]
    n_total = 8 * n_files * 5
    file_count = 0
    for h5_file_path in list_h5
        file_count += 1
        @printf("Processing %s (%d of %d)...\n", h5_file_path, file_count, n_files)
        in_h5 = h5open(args[1], "r")
        # HDF5.delete_object(ds_pavd)
        groups = keys(in_h5)
        groups = groups[startswith.(groups, "BEAM")]


        n_groups = size(groups)[1]
        group_count = 0
        #group = groups[1]
        for group in groups
            group_count += 1
            @printf("...processing %s (%d of %d)\n", group, group_count, n_groups)
            g = in_h5[group]
            n = HDF5.size(g["shot_number"])[1]
            
            mask = g["l2b_quality_flag"][:] .== 1
            seq = [1:1:n;]
            seq_mask = seq[mask]
            seq_max = maximum(seq_mask)
            seq_min = minimum(seq_mask)
            slice = seq_min:seq_max

            lats = g["geolocation/latitude_bin0"][slice][mask[slice]]
            y_ind = Int32.(floor.((90 .- lats) ./ y_res))
            lons = g["geolocation/longitude_bin0"][slice][mask[slice]]
            x_ind = Int32.(floor.((lons .+ 180) ./ x_res))
            indices = (y_ind .* xmax .+ x_ind)
            
            #col = "pavd"
            n_cols = size(cols)[1]
            
            for col in ["rh100", "pai", "cover", "fhd_normal", "pavd"]
                total_counter += 1
                @printf("\rProgress...%.2f%%", 100.0*total_counter/n_total)
                this_data = col == "pavd" ? 
                    pavd = dropdims(sum(g["pavd_z"][:,slice][:,mask[slice]], dims=1); dims=1) : 
                    g[col][slice][mask[slice]]
                this_agg = from_float32.(slice_get.((out_ds[col],), indices))
                update_stats.(this_agg, this_data)
                slice_set.((out_ds[col],), indices, to_float32.(this_agg))
            end #cols
            println()

        end #groups
    end #files
    
end #main

main(ARGS)
