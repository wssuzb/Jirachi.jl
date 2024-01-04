# push!(LOAD_PATH, "../../lib_julia_external/Jirachi")
using Jirachi
using StatsBase
band = ["g", "r", "i", "z"]

band_pair = []

for i=1:lastindex(band)
   for j=2:lastindex(band)
       if i<j
           push!(band_pair, [band[i], band[j]])
       end
   end
end

tunits = 24 * 3600
sf_bin = 0.05
cv_bin = 0.1
sf_bin_edges=0:sf_bin:5
cv_bin_edges=0:cv_bin:5
nsigma = 3#parse(Int, ARGS[1])
erron=true
nsim=1000
mode="both" #fr/rss
# fi_np::String="./run_result/run_all.h5"
lower_bounds = [0, 0, 0, 0.001]
upper_bounds = [10, 2e4, 2, 0.1]

t_cad = 103.68 #mean(diff(lc.time))
t_start = 0 - (t_cad / 2) #lc.time[1] - (t_cad / 2)
t_end = 25000 + (t_cad /2 )
lc_edges = t_start:t_cad:t_end

sf_noise_sigma = 2

println("================= input paramerters =================")

println("light curve bin size: ", string(t_cad))
println("lc bin used: ", lc_edges)

println("sf bin used: ", sf_bin_edges)
println("cv bin used: ", cv_bin_edges)
println("nsigma used in cv: ", string(nsigma))
println("erron used in cv: ", string(erron))
println("number of simul.: ", string(nsim))
println("mode used in bootstrap: ", mode)


for i = 1: lastindex(band_pair)
    println("=================", band_pair[i], "=================")

    lc1 = load_data("../test/ngc4395/Montano22_n1/" * band_pair[i][1] * "_4395.txt", [1, 2, 3];  band=band_pair[i][1])
    
    lc2 = load_data("../test/ngc4395/Montano22_n1/" * band_pair[i][2] * "_4395.txt", [1, 2, 3]; band=band_pair[i][2])
    
    fi_np="./run_bin_" * string(t_cad)  * "_sf_noise_cut_" * string(round(sf_noise_sigma, digits=2)) * "_sf_bin_" * string(sf_bin) * "_cv_bin_" * string(cv_bin) *  "_montano22_n1_nsim_" * string(nsim) * "_mode_" * string(mode) *  "_nsigma_" * string(nsigma) * "_band_" * band_pair[i][1] * "_" * band_pair[i][2] * "_sun14.h5"
    
    runall(lc1, lc2, lc_edges; tunits=tunits, label="montano22_n1", sf_bin_edges=sf_bin_edges, cv_bin_edges=cv_bin_edges, nsigma=nsigma, erron=erron, nsim=nsim, fi_np = fi_np , lower_bounds=lower_bounds, upper_bounds=upper_bounds, p0=[1, 1e3, 1, 0.005], mode=mode, sf_noise_sigma = sf_noise_sigma)
    
    println("================= end =================")

end
