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

function update_stats(n, M1, M2, M3, M4, x)
    n1 = n
    n += 1
    delta = x - M1
    delta_n = delta / n
    delta_n2 = delta_n * delta_n
    term1 = delta * delta_n * n1
    M1 += delta_n
    M4 += term1 * delta_n2 * (n*n - 3*n + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3
    M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2
    M2 += term1
    return n, M1, M2, M3, M4
end

function combine_stats(a_n, a_M1, a_M2, a_M3, a_M4, b_n, b_M1, b_M2, b_M3, b_M4)
    combined_n = a_n + b_n
    
    delta = b_M1 - a_M1
    delta2 = delta*delta
    delta3 = delta*delta2
    delta4 = delta2*delta2
    
    combined_M1 = (a_n*a_M1 + b_n*b_M1) / combined_n
    
    combined_M2 = a_M2 + b_M2 + 
                  delta2 * a_n * b_n / combined_n
    
    combined_M3 = a_M3 + b_M3 + 
                  delta3 * a_n * b_n * (a_n - b_n)/(combined_n*combined_n)
    combined_M3 += 3.0*delta * (a_n*b_M2 - b_n*a_M2) / combined_n
    
    combined_M4 = a_M4 + b_M4 + delta4*a_n*b_n * (a_n*a_n - a_n*b_n + b_n*b_n) / 
                  (combined_n*combined_n*combined_n)
    combined_M4 += 6.0*delta2 * (a_n*a_n*b_M2 + b_n*b_n*a_M2)/(combined_n*combined_n) + 
                  4.0*delta*(a_n*b_M3 - b_n*a_M3) / combined_n
    
    return combined_n, combined_M1, combined_M2, combined_M3, combined_M4
end

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

datasets = list_recursive(cerrado_h5, "/", false, true)

lines_to_read = 100
half_ysize = Int32(ceil(ysize/2))

counter = 0

for ds_name in datasets
    global counter += 1
    ds = cerrado_h5[ds_name]
    println("Processing dataset: ", ds_name, " (", counter, " of ", size(datasets)[1], ")")
    for yy in 0:lines_to_read:half_ysize
        perc = 100 * yy / half_ysize
        @printf("\rReversing... %.2f%%",perc)
        lines_read = lines_to_read
        
        if (yy + lines_read > half_ysize)
            lines_read = half_ysize - yy
        end

        rev_yy = ysize - yy - lines_read
        
        
        lower_bound = yy*xsize+1
        upper_bound = (yy+lines_read)*xsize

        mat = ds[lower_bound:upper_bound]
        mat = reshape(mat, (xsize, lines_read))
        mat = mat[:,lines_read:-1:1]

        rev_lower_bound = rev_yy*xsize+1
        rev_upper_bound = (rev_yy+lines_read)*xsize
        rev_mat = ds[rev_lower_bound:rev_upper_bound]
        rev_mat = reshape(rev_mat, (xsize, lines_read))
        rev_mat = rev_mat[:,lines_read:-1:1]
        
        ds[lower_bound:upper_bound] = rev_mat
        ds[rev_lower_bound:rev_upper_bound] = mat
    end
    @printf("\rReversing... %.2f%%\n",100)
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
