# using Pkg
# Pkg.add("HDF5")
# Pkg.add("ArchGDAL")
# Pkg.add("ArgParse")
using HDF5
using Printf
using Glob

LATITUDE_BIN0 = "geolocation/latitude_bin0"
LONGITUDE_BIN0 = "geolocation/longitude_bin0"
L2B_QUALITY_FLAG = "l2b_quality_flag"


# "E:/Documentos/Downloads/" "C:/Users/caioh/Desktop/mundo/gedi_mask.tif" "C:/Users/caioh/Desktop/mundo/base_float.tif" "C:/Users/caioh/Desktop/mundo/gedi_out"
# args = ["E:/Documentos/Downloads/", "C:/Users/caioh/Desktop/mundo/gedi_mask.tif", "C:/Users/caioh/Desktop/mundo/base_float.tif", "C:/Users/caioh/Desktop/mundo/gedi_out"]
function main(args)
    cols = ["rh100", "pai", "pavd_z", "cover", "fhd_normal"]
        
    
    list_h5 = glob("*02_B*.h5", args[1])
    total_counter = 0

    n_files = size(list_h5)[1]
    n_total = 8 * n_files * 5
    col_count = 0

    # col = "rh100"
    file_count = 0
                
        # h5_file_path = list_h5[1]
        for h5_file_path in list_h5
            file_count += 1
            @printf("\x1b[2K\r... Processing %s (%d of %d)\n", basename(h5_file_path), file_count, n_files)
            in_h5 = try
                h5open(h5_file_path, "r")
            catch e
                @warn "... Error reading: $h5_file_path\n\n"
                continue
            end
            # HDF5.delete_object(ds_pavd)
            groups = keys(in_h5)
            groups = groups[startswith.(groups, "BEAM")]

            check_datasets = deepcopy(cols)
            append!(check_datasets, [L2B_QUALITY_FLAG, LATITUDE_BIN0, LONGITUDE_BIN0])
            check_all_datasets = ["$i/$j" for (i,j) in  Iterators.product(groups, check_datasets)]
            if !(all(HDF5.haskey.(Ref(in_h5), check_datasets))) 
                @warn "\nThe H5 file $h5_file_path is missing required columns!\n\n"
                continue 
            end
            close(in_h5)
            print("\x1b[1A\r")
        end # files
end # main

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
