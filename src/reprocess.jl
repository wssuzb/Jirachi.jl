export thin_disk, thin_disk_temperature, blackbody_cgs, reprocessing

function thin_disk_temperature(rmid_zones, d::disk_parameter_update; three="no")
    radius = reshape(transpose(rmid_zones), (1, length(rmid_zones)))
    radius_cm = radius * d.rgravcm
    # tmp = 3 * 6.67384E-8 * d.bhmass * d.disk_mdot ./ (8 * pi * (radius_cm .^ 3) * 5.670373E-05)
    tmp = @. 3 * 6.67384E-8 * d.bhmass * d.disk_mdot / (8 * pi * (radius_cm ^ 3) * 5.670373E-05)
    if three == "yes"
        tmp = tmp * 27
    end
    temperature = (tmp .* (1 .- sqrt.(d.rin ./ radius))) .^0.25
    # temperature = @. (tmp * (1 - sqrt(d.rin / radius))) ^ 0.25
    return temperature
end


function blackbody_cgs(wavelength_in_bb::Any, d::disk_parameter_update, temp::Any)
    wavelength_in_bb = wavelength_in_bb .* 1.e-8
    tmp = d.con_hp * d.con_c / (d.con_kb * temp)
    sed = 2 * d.con_hp * d.con_c * d.con_c ./ (reshape(wavelength_in_bb.^5, (length(wavelength_in_bb),1))) ./ (expm1.(tmp ./wavelength_in_bb))
    # sed = @. 2 * d.con_hp * d.con_c * d.con_c / (reshape(wavelength_in_bb^5, (length(wavelength_in_bb),1))) / (expm1.(tmp /wavelength_in_bb))
    # sed = 2 * con_hp * con_c*con_c
    return sed * 1.E-8

end


function thin_disk(band::Vector{String}, in_time::Vector{Float64}, in_flux::Vector{Float64}, in_fluxe::Vector{Float64}, wave_filter::Matrix{Float64}, d::disk_parameter_update)
    ## 定义吸积盘相关变量
    # 吸积盘区域划分
    wavelength = wave_filter[:, 1]
    disk_zones = zeros((Int(d.numring), Int(d.fracring)))
    # 吸积盘内部半径
    rin_zones = copy(disk_zones)
    rin_zones[1, :] = fill!(rin_zones[1, :], d.rin)

    @inbounds @simd for i=1: size(rin_zones, 1)-1
        rin_zones[i+1, :] = fill!(rin_zones[i+1, :], rin_zones[i] * d.rratio)
    end

    # 吸积盘外部半径
    rout_zones = copy(disk_zones) #2
    rout_zones[1, :] = fill!(rout_zones[1,:], d.rin * d.rratio)
    for i=1: size(rout_zones, 1)-1
        rout_zones[i+1, :] = fill!(rout_zones[i+1, :], rout_zones[i] * d.rratio)
    end
    # 吸积盘mid半径以及面积
    rmid_zones = sqrt.(rin_zones.*rout_zones) #3
    area_zones = (pi * (rout_zones .^2 - rin_zones .^2)) / d.fracring #5
    delta_ang = 360. / d.fracring
    angle_zones = copy(disk_zones) #4
    @inbounds for i=1: size(angle_zones, 2)
        angle_zones[:, i] = fill!(angle_zones[:, i], (Int(i-1) * delta_ang + delta_ang / 2.) / 360. * 2 * pi)
    end
    temperature_zones_original = thin_disk_temperature(rmid_zones, d; three="no") #6
    temperature_zones_illuminate = copy(temperature_zones_original) #8

    # println(temperature_zones_original)
    # 定义温度框架，计算能谱，计算每个透过率曲线卷积的结果。tmp_temperature, tmp_sed

    tmp_temperature = 1:0.001:6
    # tmp_flux = copy(tmp_temperature)
    # wavelength_in_bb = wavelength
    tmp_sed_zones = zeros((length(tmp_temperature), length(wavelength)))
    @inbounds for i=1:length(tmp_temperature)
        tmp_sed_zones[i, :] = blackbody_cgs(wavelength, d, 10 ^ tmp_temperature[i])
    end

    tmp_band_flux = Dict()
    for i in 1:length(band)
        tmp_band_flux[band[i]] = collect(copy(tmp_temperature))
    end

    area_ = reshape(transpose(area_zones), (1, Int(d.numring)*Int(d.fracring)))
    tem_ori = reshape(transpose(temperature_zones_original), (1, Int(d.numring)*Int(d.fracring)))
    sed_ori = zeros((length(tem_ori),length(wavelength)))

    @inbounds for i = 1: length(tem_ori)
        sed_ori[i, :] = blackbody_cgs(wavelength,d, tem_ori[i]) * area_[i] * (d.rgravcm .^ 2)
    end
    sed_sum = collect(sum(transpose(sed_ori), dims=2))
    for j = 1:length(band)
        for i=1:length(tmp_temperature)
            tmp_band_flux[band[j]][i] = sum(wave_filter[:, j+1] .* tmp_sed_zones[i, :]) / sum(wave_filter[:, j+1])
        end
    end
    ori_flux = Dict()
    @inbounds for j = 1:length(band)
        ori_flux[band[j]] = sum(collect(sum(transpose(sed_ori), dims=2)) .* wave_filter[:, j+1]) / sum(wave_filter[:, j+1])
    end 
    # ori_flux

    rmid = reshape(transpose(rmid_zones), (1, length(rmid_zones)))
    angle = reshape(transpose(angle_zones), (1, length(angle_zones)))
    lighttravel_time = sqrt.(d.hcor ^2 .+ rmid .^ 2) .- rmid .* cos.(angle * pi / 180) * sin(d.disk_ia) .+ d.hcor * cos(d.disk_ia)
    lighttravel_time_zones = transpose(reshape(lighttravel_time, (Int(d.fracring), Int(d.numring)))) #7

    # in_time = readdlm(lcpath)[:, 1]
    # in_flux = readdlm(lcpath)[:, 2]
    # in_fluxe = readdlm(lcpath)[:, 3]

    avg = mean(in_flux)
    sig = std(in_flux, corrected=false)
    lc_sigma = 0.3
    lc_mean = 1.
    in_flux = (in_flux .- avg) ./ sig .* lc_sigma .+ lc_mean
    in_flux = in_flux * d.lumin.val

    obs_per_day = Int(round(1 / (in_time[2] - in_time[1])))
    total_days = Int(round(in_time[end] - in_time[1]))
    star_i = 0
    star_j = star_i + 50 *obs_per_day
    end_index = length(in_flux) - 50 * obs_per_day
    idx = range(Int(star_j), stop=Int(end_index))


    delay_time = lighttravel_time_zones * d.rgravcm / d.con_c
    delay_time = delay_time / (86400. / obs_per_day) # 9


    return (band=band,disk_zones = disk_zones, rin_zones = rin_zones, rout_zones = rout_zones, rmid_zones = rmid_zones, area_zones = area_zones,
    delta_ang = delta_ang, angle_zones = angle_zones, temperature_zones_original = temperature_zones_original, temperature_zones_illuminate = temperature_zones_illuminate,
    tmp_temperature = tmp_temperature, tmp_sed_zones = tmp_sed_zones, lighttravel_time = lighttravel_time,
    in_time = in_time, in_flux = in_flux, idx = idx, delay_time = delay_time, tmp_band_flux = tmp_band_flux,
    sed_ori = sed_sum, ori_flux = ori_flux)
end


function reprocessing(td::Any, items::Any, d::disk_parameter_update) # td: thin_disk return parameter
    band = td.band

    index = items[2]

    rcor = copy(td.rmid_zones)
    rcorcm = sqrt.(d.hcor ^2 .+ rcor .^2) * d.rgravcm
    s_index = Int.(round.(index .- td.delay_time) .+ 1)
    lcor_ = copy(td.rmid_zones)
    for i=1:Int(d.numring)
        for j=1:Int(d.fracring)
            lcor_[i, j] = td.in_flux[s_index[i, j]]
        end
    end

    temp = transpose(reshape(td.temperature_zones_original, (24,38)))
    albedo = 0.5
    fenmu = rcorcm .^ 3 .* 4 .* pi .* 5.670373E-05
    temp4 = temp .^ 4 .+ (1-albedo) .* lcor_ .* d.hcorcm ./ fenmu
    temperature_zones_illuminate = temp4 .^ 0.25
    # area_ =reshape(transpose(td.area_zones), (1, 912))

    flux_itp = Dict()
    band_flux = Dict()
    for i=1:length(band)
        flux_itp[band[i]] = LinearInterpolation(td.tmp_temperature, td.tmp_band_flux[band[i]], extrapolation_bc=Flat())
        band_flux[band[i]] = flux_itp[band[i]].(log10.(temperature_zones_illuminate)).* td.area_zones .*  (d.rgravcm .^ 2)./(4*pi * d.disk_ld ^2)
    end

    # for i=1:length(band)
    #     band_flux[band[i]] = flux_itp[band[i]].(log10.(temperature_zones_illuminate)).* td.area_zones .*  (d.rgravcm .^ 2)./(4*pi * d.disk_ld ^2)
    # end

    # return (idx = items[1], ary=[td.in_time[index], sum(band_flux["uw2"]),sum(band_flux["um2"]),sum(band_flux["uw1"]),sum(band_flux["uuu"]),sum(band_flux["ubb"]),
	# 	sum(band_flux["uvv"])])
    return (idx = items[1], ary=[td.in_time[index], sum.(band_flux[i] for i in band)])
end