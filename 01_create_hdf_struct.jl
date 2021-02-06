# using Pkg
# Pkg.add("HDF5")
# Pkg.add("ArchGDAL")
# Pkg.add("ArgParse")
using HDF5
using Printf
using Glob
using ArchGDAL
using DataFrames
using Printf

include("old/04update_statistics.jl")


function create_stat_tif(aggs, masked_inds, stat_agg, xsize, ysize, inds_x, inds_y, base_raster, out_root, col, stat_name)
    stat = [stat_agg(aggs[(x, y)]) for (x,y) in masked_inds]
    vals = fill(Float32(-9999), (xsize, ysize))
    setindex!.(Ref(vals), stat, inds_x, inds_y)
    out_tif = Printf.@sprintf("%s_%s_%s.tif", out_root, col, stat_name)

    cp(base_raster, out_tif)
    tif = ArchGDAL.read(out_tif, flags=ArchGDAL.OF_Update)
    band = ArchGDAL.getband(tif, 1)
    band[:,:] = vals
    ArchGDAL.setnodatavalue!(band, Float32(-9999.0))
    ArchGDAL.destroy(tif)
end

function slice_set(ds, x, y, value) 
    ds[x, y] = value
end

function slice_get(ds, pairs)
    ds[pairs[1], pairs[2], :]
end

function xy_to_lonlat(x, y, xres, yres, lon_min, lat_min)
    lon = x*xres + lon_min
    lat = y*(yres*-1) + lat_min
    [lon, lat]
end

# "E:/Documentos/Downloads/" "C:/Users/caioh/Desktop/mundo/gedi_mask.tif" "C:/Users/caioh/Desktop/mundo/base_float.tif" "C:/Users/caioh/Desktop/mundo/gedi_out"
# args = ["E:/Documentos/Downloads/", "C:/Users/caioh/Desktop/mundo/gedi_mask.tif", "C:/Users/caioh/Desktop/mundo/base_float.tif", "C:/Users/caioh/Desktop/mundo/gedi_out"]
function main(args)
    n_lines_read = 100
    lat_to_met_factor = 1 / 110540
    lon_to_met_factor = 1 / 111320
    resolution = 3000 # meters
    xres = resolution * lon_to_met_factor
    yres = resolution * lat_to_met_factor
    base_float_path = args[3]

    out_root = args[4]
    tif = ArchGDAL.read(args[2])
    band = ArchGDAL.getband(tif, 1)
    vals = band[:,:]
    mask = vals .== 1
    xsize, ysize = ArchGDAL.size(band)
    lon_min,xres,_,lat_min,_,yres = ArchGDAL.getgeotransform(tif)
    ArchGDAL.destroy(tif)
    the_inds = collect(Iterators.product(1:xsize, 1:ysize))[mask]
    
    dcpl = HDF5.create_property(HDF5.H5P_DATASET_CREATE)
    HDF5.h5p_set_deflate(dcpl, 5)
    HDF5.h5p_set_chunk(dcpl, 3, [Int64(floor(3.6 / xres)), Int64(floor(1.8 / yres*-1)), Int64(5)])
    
    
    cols = ["rh100", "pai", "pavd", "cover", "fhd_normal"]
        
    
    list_h5 = glob("*02_B*.h5", args[1])
    total_counter = 0

    n_files = size(list_h5)[1]
    n_total = 8 * n_files * 5
    file_count = 0
    col_count = 0

    # col = "rh100"
    for col in ["rh100", "pai", "cover", "fhd_normal", "pavd"]
        col_count += 1
        @printf("Processing %s (%d of %d)\n", col, col_count, 5)
        aggs = Dict(i => AggregateStats(0f32,0,0,0,0) for i in the_inds)
                
        # h5_file_path = list_h5[1]
        for h5_file_path in list_h5
            file_count += 1
            @printf("... Processing %s (%d of %d)\n", basename(h5_file_path), file_count, n_files)
            in_h5 = h5open(h5_file_path, "r")
            # HDF5.delete_object(ds_pavd)
            groups = keys(in_h5)
            groups = groups[startswith.(groups, "BEAM")]


            n_groups = size(groups)[1]
            x_inds = Array{Int32,1}()
            y_inds = Array{Int32,1}()
            vals = Array{Float32,1}()

            # group = groups[1]
            for group in groups
                total_counter += 1
                @printf("\r....... Progress: %.2f%%", 100.0 * total_counter / n_total)
                g = in_h5[group]
                n = HDF5.size(g["shot_number"])[1]
                
                mask = g["l2b_quality_flag"][:] .== 1
                seq = [1:n;]
                seq_mask = seq[mask]
                seq_max = maximum(seq_mask)
                seq_min = minimum(seq_mask)
                slice = seq_min:seq_max
                mask_slice = mask[slice]

                lats = g["geolocation/latitude_bin0"][slice][mask_slice]
                lons = g["geolocation/longitude_bin0"][slice][mask_slice]
                
                y_ind = Int32.(floor.((lats .- lat_min) ./ yres))
                append!(y_inds, y_ind)
                x_ind = Int32.(floor.((lons .- lon_min) ./ xres))
                append!(x_inds, x_ind)
                g["geolocation"]["elev_highestreturn"][1:10]
                this_vals = 
                    col == "pavd" ? 
                    dropdims(sum(g["pavd_z"][:,slice][:,mask_slice], dims=1); dims=1) : 
                    g[col][slice][mask_slice]
                
                append!(vals, this_vals)    
            end # groups
            println("")

            
            this_inds = collect(zip(x_inds, y_inds))
            mask = haskey.(Ref(aggs), this_inds).==1
            masked_inds = this_inds[mask]
            masked_vals = vals[mask]
        

            for (ind, val) in zip(masked_inds, masked_vals)
                update_stats(aggs[ind], val)
            end
            close(in_h5)
        end # files

        mask = agg_n.(values(aggs)).>0
        masked_inds = the_inds[mask]
        inds_x = getindex.(masked_inds, 1)
        inds_y = getindex.(masked_inds, 2)

        create_stat_tif(aggs, masked_inds, agg_n, xsize, ysize, inds_x, inds_y, base_float_path, out_root, col, "count")
        create_stat_tif(aggs, masked_inds, agg_mean, xsize, ysize, inds_x, inds_y, base_float_path, out_root, col, "mean")
        create_stat_tif(aggs, masked_inds, agg_sd, xsize, ysize, inds_x, inds_y, base_float_path, out_root, col, "sd")
        create_stat_tif(aggs, masked_inds, agg_kur, xsize, ysize, inds_x, inds_y, base_float_path, out_root, col, "kur")
        create_stat_tif(aggs, masked_inds, agg_skew, xsize, ysize, inds_x, inds_y, base_float_path, out_root, col, "skew")
    end # cols
    
end # main

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
