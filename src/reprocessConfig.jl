export disk_parameter, disk_parameter_update, disk_systerm, trans_curve_info, trans_curve, wavelength_filter

Base.@kwdef mutable struct disk_parameter
    disk_ia::Float64 = 45
    lumin::Float64 = 2.82e+44
    redshift::Float64 = 0.017
    bhmass::Float64 = 32000000.0
    disk_η::Float64 = 0.083
    fracring::Float64 = 24 
    illuminating_height::Float64 = 10 
    numrings::Float64 = 38
    inner_radius::Float64 = 6.0
    radius_ratio::Float64 = 1.2
end

mutable struct disk_parameter_update
    rin::Float64
    fracring::Float64
    numring::Float64
    rratio::Float64
    hcor::Float64
    bhmass::Float64
    disk_z::Float64
    lumin::Any
    rgravcm::Float64
    hcorcm::Float64
    disk_eta::Float64
    disk_ld::Float64
    disk_mdot::Float64

    con_hp::Float64
    con_kb::Float64
    con_c::Float64
    disk_ia::Float64
end


function Base.show(disk::disk_parameter)
    println("==================Disk===================")
    println(" ")
    
    println("disk inclination angle: ", disk.disk_ia  * pi /180.)
    println("disk luminosity: ", disk.lumin *  1u"erg" / 1u"s")
    println("disk redshift: ", disk.redshift)
    println("disk black hole: ", uconvert(u"g", disk.bhmass * (1u"Msun"|>upreferred)).val)
    println("disk accretion rate eta: ", disk.disk_η)
    println("disk gravity radius in units of cm: ", uconvert(u"cm/g", 1u"G" * disk.bhmass / (1u"c" * 1u"c")|> upreferred).val)
    println("disk hcorcm: ", 10 * uconvert(u"cm/g", 1u"G" * disk.bhmass / (1u"c" * 1u"c")|> upreferred).val) # 10 * rgravcm
    println(" ")
    println("===============Simulaiton================")
    println(" ")
    println("pieces of a ring used: ", disk.fracring)
    println("number of a ring used: ", disk.numrings)
    println("inner radius: ", disk.inner_radius)
    println("radius ratio: ", disk.radius_ratio)
    println("illuminating height: ", disk.illuminating_height)
    println(" ")
    println("===============Constants================")
    println(" ")
    println("disk luminosity distance: ", 2.2760023267972407e+26) ## need change
    println("disk M dot: ", 3.780329106109885e+24) ## need change
    println("plank constant hp: ", 6.62606957e-27)
    println("boltzman constant kb:", 1.3806488e-16)
    println("light speed: ", Float64(uconvert(u"cm/s", (1u"c" |> upreferred)).val))
    println(" ")
    println("==================END====================")
end


function disk_systerm(d::disk_parameter)
    show(d)
    rin = d.inner_radius
    fracring = d.fracring
    numring = d.numrings
    rratio = d.radius_ratio
    hcor = d.illuminating_height
    bhmass = uconvert(u"g", d.bhmass * (1u"Msun"|>upreferred)).val
    disk_z = d.redshift
    lumin = d.lumin * 1u"erg" / 1u"s"
    rgravcm = uconvert(u"cm/g", 1u"G" * bhmass / (1u"c" * 1u"c")|> upreferred).val
    hcorcm = 10 * rgravcm
    disk_eta = d.disk_η
    disk_ld = 2.2760023267972407e+26
    disk_mdot = 3.780329106109885e+24

    con_hp = 6.62606957e-27 
    con_kb =  1.3806488e-16
    con_c = Float64(uconvert(u"cm/s", (1u"c" |> upreferred)).val) #29979245800
    disk_ia = d.disk_ia * pi /180.

    disk_update = disk_parameter_update(rin, fracring, numring, rratio, hcor, bhmass, disk_z, lumin, rgravcm, hcorcm, disk_eta, disk_ld, disk_mdot, con_hp, con_kb, con_c, disk_ia)
    return disk_update
end


mutable struct trans_curve_info
    band::Vector{String}
    weff::Vector{Float64}
end


function Base.show(trans_info::trans_curve_info)
    println("Band used in reprocessing model: ", trans_info.band)
    println("Band W_eff: ", trans_info.weff)
end


function trans_curve(filepath::String, band::Vector{String}, weff::Vector{Float64})#; plot_trans_curve::Bool=true)

    println("The transmission file located in: ", filepath)
    println("Transmission curve file name: ", readdir(filepath))
    
    trans_info = trans_curve_info(
        band, weff
    )

    show(trans_info)

    # if plot_trans_curve
    #     p = plot(
    #         legend = true, 
    #         framestyle = :box,
    #         # aspect_ratio = 1, 
    #         size = (400, 200),
    #         format = :svg
    #     )        
    #     for i in 1:length(trans_info.band)
    #         t = readdlm(filepath * "/total_" * band[i] * ".txt")
    #         plot!(p, t[:, 1], t[:, 2], label=trans_info.band[i])
    #     end
    #     display(p)
    #     trans = Dict()
    #     for i in 1:length(trans_info.band)
    #         t = readdlm(filepath * "/total_" * band[i] * ".txt")
    #         trans[trans_info.band[i]] = (wavelength=t[:, 1], strength=t[:, 2])
    #     end
    # else
    trans = Dict()
    for i in 1:length(trans_info.band)
        t = readdlm(filepath * "/total_" * band[i] * ".txt")
        trans[trans_info.band[i]] = (wavelength=t[:, 1], strength=t[:, 2])
    end
    # end
    return (trans=trans, band=band, weff=weff)
end


function wavelength_filter(w_start::Float64, w_end::Float64,transmission, band::Vector{String})
    wavelength = range(w_start, stop=w_end,length=Int(w_end - w_start) + 1)
    six_filter = zeros((Int(w_end - w_start) + 1, length(transmission)+1))
    six_filter[:, 1] = wavelength

    for i=1:length(transmission)
        tran_w = transmission[band[i]].wavelength # wavelength
        resp = transmission[band[i]].strength # specresp
        itp = LinearInterpolation(tran_w, resp, extrapolation_bc=Flat())
        six_filter[:, i+1] = itp(wavelength)
    end
    return six_filter
end