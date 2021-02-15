# using Pkg
# Pkg.add("HDF5")
# Pkg.add("ArchGDAL")
# Pkg.add("ArgParse")
using HDF5
using Printf
using Glob
using ArchGDAL

include("old/04update_statistics.jl")
LATITUDE_BIN0 = "geolocation/latitude_bin0"
LONGITUDE_BIN0 = "geolocation/longitude_bin0"
L2B_QUALITY_FLAG = "l2b_quality_flag"

function create_raster(output::AbstractString, xres::Float64, yres::Float64; ul_lat=51.7, ul_lon=-180, lr_lat=-51.7, lr_lon=180, dtype=Float32, tile_size=1024)::ArchGDAL.Dataset
    height = Int32(ceil((ul_lat - lr_lat) / yres))
    width = Int32(ceil((lr_lon - ul_lon) / xres))
    tif = ArchGDAL.unsafe_create(
        output;
        driver=ArchGDAL.getdriver("GTiff"),
        width=width,
        height=height,
        nbands=1,
        dtype=dtype,
        options=["COMPRESS=DEFLATE", "BIGTIFF=YES", "TILED=YES", "BLOCKXSIZE=$tile_size", "BLOCKYSIZE=$tile_size"])
    try
        ArchGDAL.setgeotransform!(tif, [ul_lon, xres, 0.0, ul_lat, 0.0, -yres])
        if !occursin(r"^U",string(dtype))
            ArchGDAL.setnodatavalue!(ArchGDAL.getband(tif,1), dtype(-9999))
        else
            ArchGDAL.setnodatavalue!(ArchGDAL.getband(tif, 1), dtype(0))
        end 
    catch
        ArchGDAL.destroy(tif)
    end
    return tif
end


# julia 02_create_hdf_raster_notmemory.jl "G:\01_GEDI\02_Level2B" "I:\01_GlobalRasters_L2B\gedi_out" 51.7 -180 -51.7 180 25 1024
# args = ["E:/Documentos/Downloads/", "C:/Users/caioh/Desktop/mundo/gedi_out", "51.7", "-180", "-51.7", "180", "1000", "1024"]
function main(args)
    n_lines_read = 100
    lat_to_met_factor = 1 / 110540
    lon_to_met_factor = 1 / 111320
    resolution = parse(Float32, args[7]) # meters
    xres = resolution * lon_to_met_factor
    yres = -resolution * lat_to_met_factor
    ul_lat, ul_lon, lr_lat, lr_lon = parse.(Float32, args[3:6])

    # tif = ArchGDAL.read(args[2])
    # band = ArchGDAL.getband(tif, 1)
    out_root = args[2]
    tile_size = parse(Int32, args[8])
    
    lon_min = ul_lon
    lat_min = ul_lat

    xsize = Int32(ceil((lr_lon  - ul_lon) / xres))
    ysize = Int32(ceil((lr_lat - ul_lat) / yres))
    
    cols = ["rh100", "pai", "pavd", "cover", "fhd_normal"]
    
    list_h5 = glob("*02_B*.h5", args[1])
    total_counter = 0

    n_files = size(list_h5)[1]
    n_total = 8 * n_files * 5
    col_count = 0

    # col = "rh100"
    for col in ["rh100", "pai", "cover", "fhd_normal", "pavd"]
        file_count = 0
        col_count += 1
        @printf("\x1b[2K\rProcessing %s (%d of %d)\n", col, col_count, 5)
        tif_m1_path = Printf.@sprintf("%s_%s_%s.tif", out_root, col, "m1")
        tif_m2_path = Printf.@sprintf("%s_%s_%s.tif", out_root, col, "m2")
        tif_m3_path = Printf.@sprintf("%s_%s_%s.tif", out_root, col, "m3")
        tif_m4_path = Printf.@sprintf("%s_%s_%s.tif", out_root, col, "m4")
        tif_count_path = Printf.@sprintf("%s_%s_%s.tif", out_root, col, "count")
        tif_m1 = create_raster(tif_m1_path, xres, -yres, ul_lat=ul_lat, ul_lon=ul_lon, lr_lat=lr_lat, lr_lon=lr_lon, dtype=Float32, tile_size=tile_size)
        tif_m2 = create_raster(tif_m2_path, xres, -yres, ul_lat=ul_lat, ul_lon=ul_lon, lr_lat=lr_lat, lr_lon=lr_lon, dtype=Float32, tile_size=tile_size)
        tif_m3 = create_raster(tif_m3_path, xres, -yres, ul_lat=ul_lat, ul_lon=ul_lon, lr_lat=lr_lat, lr_lon=lr_lon, dtype=Float32, tile_size=tile_size)
        tif_m4 = create_raster(tif_m4_path, xres, -yres, ul_lat=ul_lat, ul_lon=ul_lon, lr_lat=lr_lat, lr_lon=lr_lon, dtype=Float32, tile_size=tile_size)
        tif_count = create_raster(tif_count_path, xres, -yres, ul_lat=ul_lat, ul_lon=ul_lon, lr_lat=lr_lat, lr_lon=lr_lon, dtype=UInt32, tile_size=tile_size)
        band_m1 = ArchGDAL.getband(tif_m1, 1)
        band_m2 = ArchGDAL.getband(tif_m2, 1)
        band_m3 = ArchGDAL.getband(tif_m3, 1)
        band_m4 = ArchGDAL.getband(tif_m4, 1)
        band_count = ArchGDAL.getband(tif_count, 1)
                
        # h5_file_path = list_h5[1]
        for h5_file_path in list_h5
            file_count += 1
            @printf("\x1b[2K\r... Processing %s (%d of %d)\n", basename(h5_file_path), file_count, n_files)
            in_h5 = try
                h5open(h5_file_path, "r")
            catch e
                @warn "... Error reading: $h5_file_path"
                continue
            end
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
                dataset_name = col == "pavd" ? "pavd_z" : col
                if !(HDF5.haskey(g, L2B_QUALITY_FLAG) &&
                    HDF5.haskey(g, dataset_name) &&
                    HDF5.haskey(g, LATITUDE_BIN0) &&
                    HDF5.haskey(g, LONGITUDE_BIN0)
                    ) 
                    @warn "\nThe H5 file $h5_file_path is missing required columns!\n\n"
                    continue 
                end
                n = HDF5.size(g[L2B_QUALITY_FLAG])[1]
                
                mask = g[L2B_QUALITY_FLAG][:] .== 1
                seq = [1:n;]
                seq_mask = seq[mask]
                if size(seq_mask)[1] == 0 continue end

                seq_max = maximum(seq_mask)
                seq_min = minimum(seq_mask)
                slice = seq_min:seq_max
                mask_slice = mask[slice]

                lats = g[LATITUDE_BIN0][slice][mask_slice]
                lons = g[LONGITUDE_BIN0][slice][mask_slice]
                
                y_ind = Int32.(floor.((lats .- lat_min) ./ yres))
                append!(y_inds, y_ind)
                x_ind = Int32.(floor.((lons .- lon_min) ./ xres))
                append!(x_inds, x_ind)                

                this_vals = 
                    col == "pavd" ? 
                    dropdims(sum(g["pavd_z"][:,slice][:,mask_slice], dims=1); dims=1) : 
                    g[col][slice][mask_slice]
                
                append!(vals, this_vals)    
            end # groups
            

            blocks_x = Int32.(floor.(x_inds / tile_size))
            blocks_y = Int32.(floor.(y_inds / tile_size))

            # xblock, yblock = 9, 1
            

            for (xblock, yblock) in unique(zip(blocks_x, blocks_y))
                mask = (blocks_x .== xblock) .& (blocks_y .== yblock)
                m1_buffer = Array{Float32}(undef, tile_size, tile_size)
                m2_buffer = Array{Float32}(undef, tile_size, tile_size)
                m3_buffer = Array{Float32}(undef, tile_size, tile_size)
                m4_buffer = Array{Float32}(undef, tile_size, tile_size)
                count_buffer = Array{UInt32}(undef, tile_size, tile_size)
                ArchGDAL.readblock!(band_m1, xblock, yblock, m1_buffer)
                ArchGDAL.readblock!(band_m2, xblock, yblock, m2_buffer)
                ArchGDAL.readblock!(band_m3, xblock, yblock, m3_buffer)
                ArchGDAL.readblock!(band_m4, xblock, yblock, m4_buffer)
                ArchGDAL.readblock!(band_count, xblock, yblock, count_buffer)
                aggs = AggregateStats.(count_buffer, m1_buffer, m2_buffer, m3_buffer, m4_buffer)
                update_stats.(getindex.(Ref(aggs), x_inds[mask].-(xblock*tile_size).+1, y_inds[mask].-(yblock*tile_size).+1), vals[mask])
                
                ArchGDAL.writeblock!(band_m1, xblock, yblock, getfield.(aggs, :M1))
                ArchGDAL.writeblock!(band_m2, xblock, yblock, getfield.(aggs, :M2))
                ArchGDAL.writeblock!(band_m3, xblock, yblock, getfield.(aggs, :M3))
                ArchGDAL.writeblock!(band_m4, xblock, yblock, getfield.(aggs, :M4))
                ArchGDAL.writeblock!(band_count, xblock, yblock, getfield.(aggs, :n))
            end

            
            close(in_h5)
            print("\x1b[1A\r")
        end # files

        ArchGDAL.destroy(tif_m1)
        ArchGDAL.destroy(tif_m2)
        ArchGDAL.destroy(tif_m3)
        ArchGDAL.destroy(tif_m4)
        ArchGDAL.destroy(tif_count)
        print("\x1b[1A")
    end # cols
end # main

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end


