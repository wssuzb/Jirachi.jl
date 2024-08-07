export corsig, xcor, peakcent, xcor_mc, interpolate_with_max_gap


function corsig(r, v)
    tst = r * sqrt(v / (1 - r ^ 2))
    pvalue = 1 - cdf( TDist(v), tst)
    return pvalue
end

function xcor(lc1::lightcurve, lc2::lightcurve, trange::Tuple{Float64, Float64}, tunit::Float64, imode::Int64, itp_gap::Tuple{Symbol, Float64}=(:no, 60.0))

    tlagmin, tlagmax = trange
    
    t1, y1 = lc1.time, lc1.flux
    t2, y2 = lc2.time, lc2.flux

    n1 = length(lc1.time)
    n2 = length(lc2.time)

    (n1<2) | (n2 <2) && throw(DomainError("Length of lc1 is $n1, of lc2 is $n2, should be larger than 2!!!"))

    safe = tunit * 0.1
    taulist12, taulist21 = [], []
    npts12, npts21 = [], []
    ccf12, ccf21 = [], []

    tau_max = tlagmax + safe

    imode != 1 ? tau = tlagmin : tau = tau_max
    
    while tau < tau_max
        t2new = @. t1 + tau
        selin = minimum(t2) .<= t2new .<= maximum(t2)
        knot = sum(selin)
        if knot >0
            if (itp_gap[1] == :no)
                # y2new = interpolate_with_max_gap(t2, y2, t2new[selin], itp_gap[2])
                y2new = do_with_gap(t2, y2, t2new[selin], itp_gap[2])
                idx = all.(isfinite, y2new) # find values without nan
                y2new = y2new[idx]

                new_itp = Interpolations.linear_interpolation(t2new, y1, extrapolation_bc = Line())
                y1new = new_itp.(t2new[selin][idx])

                y1sum = sum(y1new)
                y1sqsum = sum(y1new .* y1new)

                y2sum = sum(y2new)
                y2sqsum = sum(y2new .* y2new)
                y1y2sum = sum(y1new .* y2new)
                
                knot = length(y1new)

            else
                itp = Interpolations.linear_interpolation(t2, y2, extrapolation_bc=Line())
                y2new = itp.(t2new[selin])
                y1new = y1[selin]

                y1sum = sum(y1new)
                y1sqsum = sum(y1new .* y1new)
                y2sum = sum(y2new)
                y2sqsum = sum(y2new .* y2new)
                y1y2sum = sum(y1new .* y2new)
            end

            
            fn = float(knot)
            rd1_sq = @. fn * y2sqsum - y2sum * y2sum
            rd2_sq = @. fn * y1sqsum - y1sum * y1sum

            (rd1_sq > 0) ? (rd1 = sqrt(rd1_sq)) : (rd1 = 0.0)
            (rd2_sq > 0) ? (rd2 = sqrt(rd2_sq)) : (rd2 = 0.0)

            (rd1 * rd2 == 0) ? (r = 0.0) : (r = (fn * y1y2sum - y2sum * y1sum) / (rd1 * rd2))

            push!(ccf12, Float64(r))
            push!(taulist12, Float64(tau))
            push!(npts12, Int64(knot))

        end
        tau += tunit
    end

    (imode != 2) ? (tau = tlagmin) : (tau = tau_max)

    while tau < tau_max

        t1new = @. t2 - tau
        selin = minimum(t1) .<= t1new .<= maximum(t1)
        knot = sum(selin)

        if knot > 0
            if (itp_gap[1] == :no)
                # y1new = interpolate_with_max_gap(t1, y1, t1new[selin], itp_gap[2])
                y1new = do_with_gap(t1, y1, t1new[selin], itp_gap[2])
                idx = all.(isfinite, y1new) # find values without nan
                y1new = y1new[idx]
                new_itp = Interpolations.linear_interpolation(t1new, y2, extrapolation_bc = Line())
                y2new = new_itp.(t1new[selin][idx])
                y1sum = sum(y1new)
                y1sqsum = sum(y1new .* y1new)
                y2sum = sum(y2new)
                y2sqsum = sum(y2new .* y2new)
                y1y2sum = sum(y1new .* y2new)
                knot = length(y1new)
            else
                itp = Interpolations.linear_interpolation(t1, y1,extrapolation_bc=Line())
                y1new = itp.(t1new[selin])
                y2new = y2[selin]
                y2sum = sum(y2new)
                y2sqsum = sum(y2new .* y2new)
                y1sum = sum(y1new)
                y1sqsum = sum(y1new .* y1new)
                y1y2sum = sum(y1new .* y2new)
            end

            fn = float(knot)
            rd1_sq = @. fn * y2sqsum - y2sum * y2sum
            rd2_sq = @. fn * y1sqsum - y1sum * y1sum

            (rd1_sq > 0) ? (rd1 = sqrt(rd1_sq)) : (rd1 = 0.0)
            (rd2_sq > 0) ? (rd2 = sqrt(rd2_sq)) : (rd2 = 0.0)

            (rd1 * rd2 == 0) ? (r = 0.0) : (r = (fn * y1y2sum - y2sum * y1sum) / (rd1 * rd2))

            push!(ccf21, r)
            push!(taulist21, tau)
            push!(npts21, knot)
        end
        tau += tunit
    end

    if imode == 0
        # make sure taulist12 and taulist21 have the same size!!!
        if isequal(taulist12, taulist21)
            ccf = @. (ccf12 + ccf21) * 0.5
            taulist = taulist12
            npts = npts12
        else
            taulist = intersect(taulist12, taulist21)
            sel_cb12 = findall(in(taulist21), taulist12)
            sel_cb21 = findall(in(view(taulist12, sel_cb12)), taulist21)

            ccf = @. (ccf12[sel_cb12] + ccf21[sel_cb21]) * 0.5
            npts = @. (npts12[sel_cb12] + npts21[sel_cb21]) * 0.5
        end
    elseif imode == 1
        ccf = ccf21
        taulist = taulist21
        npts = npts21
    else
        ccf = ccf12
        taulist = taulist12
        npts = npts12
    end
    return (ccf=Float64.(ccf), taulist=Float64.(taulist), npts=Float64.(npts))
end

function peakcent(
    lc1::lightcurve, lc2::lightcurve,
    trange::Tuple{Float64, Float64}, tunit::Float64, 
    thres::Float64, siglevel::Float64, imode::Int, sigmode::Float64, itp_gap::Tuple{Symbol, Float64}=(:no, 60.0)
    )

    tlagmin, tlagmax = trange

    alpha = 1.0 - siglevel

    ccf_pack = xcor(lc1, lc2, trange, tunit, imode, itp_gap)
    
    max_indx = argmax(ccf_pack.ccf)
    max_rval = ccf_pack.ccf[max_indx]
    
    safe_sqrt = -0.0001

    ccf_pack.npts[max_indx] > 2.0 ? peak_pvalue = corsig(ccf_pack.ccf[max_indx] + safe_sqrt, ccf_pack.npts[max_indx] - 2.0) : peak_pvalue = 1.0

    if sigmode > 0
        if (max_rval >= sigmode) && (ccf_pack.taulist[max_indx] > tlagmin) && (ccf_pack.taulist[max_indx] < tlagmax)
            tlag_peak = ccf_pack.taulist[max_indx]
            max_rval = max_rval
            status_peak = 1
            status_rval = 1
            status_centroid = 0
            tlag_centroid = -9999.0
        else

            max_rval = -9999.0
            tlag_peak = -9999.0
            tlag_centroid = -9999.0
            status_peak = 0
            status_rval = 0
            status_centroid = 0

        end
    else
        if (peak_pvalue < alpha) && (ccf_pack.taulist[max_indx] > tlagmin) && (ccf_pack.taulist[max_indx] < tlagmax)
            tlag_peak = ccf_pack.taulist[max_indx]
            max_rval = max_rval
            status_peak = 1
            status_rval = 1
            status_centroid = 0
            tlag_centroid = -9999.0
        else
            max_rval = -9999.0
            tlag_peak = -9999.0
            tlag_centroid = -9999.0
            status_peak = 0
            status_rval = 0
            status_centroid = 0
        end


    end

    if status_peak == 1

        rcent = thres * max_rval
        rdif_neg = @. (ccf_pack.ccf - rcent < 0.0)
        tlag_rneg = @. ccf_pack.taulist[rdif_neg] - tlag_peak
        tlag_leftall = @. abs(tlag_rneg[tlag_rneg .< 0])
        tlag_rightall = @. abs(tlag_rneg[tlag_rneg .> 0])

        if (length(tlag_leftall) > 0) && (length(tlag_rightall) > 0)
            rdif_pos = @. (ccf_pack.ccf - rcent >= 0.0)
            
            if sum(rdif_pos) > 0
                tlag_centroid = sum(ccf_pack.ccf[rdif_pos] .* ccf_pack.taulist[rdif_pos]) / sum(ccf_pack.ccf[rdif_pos])
                if isa(tlag_centroid, Number) #&& (typeof(tlag_centroid) != Int)
                    status_centroid = 1
                    tlag_centroid = tlag_centroid
                else
                    status_centroid = 0
                    tlag_centroid = -9999.0
                end
            end
        end
    end

    if status_centroid == 0
        status_peak = 0
        tlag_peak = -9999.0
        max_rval = -9999.0
        status_rval = 0
    end

    return (tlag_peak = tlag_peak, status_peak = status_peak, tlag_centroid = tlag_centroid, status_centroid = status_centroid, ccf_pack = ccf_pack, max_rval = max_rval, status_rval = status_rval, peak_pvalue = peak_pvalue)
    
end

function xcor_mc(lc1::lightcurve, lc2::lightcurve, trange::Tuple{Float64, Float64}, tunit::Float64, thres::Float64, siglevel::Float64, imode::Int64, nsim::Int64, mcmode::String, sigmode::Float64, itp_gap::Tuple{Symbol, Float64}=(:no, 60.0))

    numt1, numt2 = length(lc1.time), length(lc2.time)
    
    (numt1 < 2) | (numt2 < 2) && throw(DomainError("Length of lc1 is $numt1, of lc2 is $numt2, should be larger than 2!!!"))

    tlags_peak = zeros(nsim)
    tlags_centroid = zeros(nsim)
    pvals = zeros(nsim)
    nsuccess_peak = 0
    nsuccess_rvals = 0
    nfail_peak = 0
    nsuccess_centroid = 0
    nfail_centroid = 0
    nfail_rvals = 0
    max_rvals = zeros(nsim)

    for i=1: nsim

        lc1_new = lc_bootstrapped(lc1; seed = i, mode = mcmode)
        lc2_new = lc_bootstrapped(lc2; seed = i, mode = mcmode)

        (length(lc1_new.time) < 2) | (length(lc2_new.time) < 2) && continue

        pc_pack = peakcent(lc1_new, lc2_new, trange, tunit, thres, siglevel, imode, sigmode, itp_gap)

        if (pc_pack.status_peak == 1)
            tlags_peak[i] = pc_pack.tlag_peak
            pvals[i]  = pc_pack.peak_pvalue
            nsuccess_peak += 1
        elseif pc_pack.status_peak == 0
            nfail_peak += 1
        end

        if (pc_pack.status_centroid == 1)
            tlags_centroid[i] = pc_pack.tlag_centroid
            nsuccess_centroid += 1
        else
            nfail_centroid += 1
        end

        if pc_pack.status_rval == 1
            max_rvals[i] = pc_pack.max_rval
            nsuccess_rvals += 1
        else
            nfail_rvals += 1
        end
    end

    return (tlags_peak=tlags_peak, tlags_centroid=tlags_centroid, nsuccess_peak=nsuccess_peak, nfail_peak=nfail_peak, nsuccess_centroid=nsuccess_centroid, nfail_centroid=nfail_centroid, max_rvals=max_rvals, nfail_rvals=nfail_rvals, pvals=pvals)
end


function interpolate_with_max_gap(orig_x, orig_y, target_x, max_gap=9999, orig_x_is_sorted=false, target_x_is_sorted=false)

    if !orig_x_is_sorted
        idx = sortperm(orig_x)
        orig_x = orig_x[idx]
        orig_y = orig_y[idx]
    end

    if !target_x_is_sorted
        target_idx = sortperm(target_x)
        target_idx_for_reverse = sortperm(target_idx)
        target_x = target_x[target_idx]
    end

    target_y = zeros(length(target_x))
    idx_orig = 1 # julia index start from 1, while python from 0.
    orig_gone_through = false

    for (idx_target, x_new) in enumerate(target_x)
        while !orig_gone_through
            if idx_orig +1 >= length(orig_x)
                orig_gone_through = true
            elseif x_new > orig_x[idx_orig + 1]
                idx_orig += 1
            else
                break
            end
        end
        
        if orig_gone_through
            target_y[idx_target] = NaN
            continue
        end

        x1 = orig_x[idx_orig]
        y1 = orig_y[idx_orig]
        x2 = orig_x[idx_orig + 1]
        y2 = orig_y[idx_orig + 1]

        if x_new < 1
            target_y[idx_target] = NaN
            continue
        end

        Δx = x2 - x1
        if Δx > max_gap
            target_y[idx_target] = NaN
            continue
        end

        Δy = y2 - y1
        if Δx == 0
            target_y[idx_target] = NaN
            continue
        end

        k = Δy / Δx
        Δx_new = x_new - x1
        Δy_new = k * Δx_new
        y_new = y1 + Δy_new
        target_y[idx_target] = y_new
    end

    if !target_x_is_sorted
        res = target_y[target_idx_for_reverse]
    else
        res = target_y
    end

    return res
end


function find_gap(lc1::lightcurve, lc2::lightcurve, max_gap::T=60) where {T}
    # get gap in lightcurve 1
    idx1 = findall(vcat(0, diff(lc1.time)) .>= max_gap)
    gap1 = (lc1.time[idx1] .+ lc1.time[idx1.-1]) * 0.5
    # get gap in lightcurve 2
    idx2 = findall(vcat(0, diff(lc2.time)) .>= max_gap)
    gap2 = (lc2.time[idx2] .+ lc2.time[idx2.-1]) * 0.5
    
    length(gap1) > length(gap2) ? gap = gap1 : gap = gap2

    # add start and end to the gap
    gap = vcat(minimum([minimum(lc1.time), minimum(lc2.time)]) - max_gap, gap, maximum([maximum(lc1.time), maximum(lc2.time)]) + max_gap)
    return gap
end



function do_with_gap(x0, y0, xint, maxgap)
    
    # some problems... need to debug

    itp = LinearInterpolation(x0, y0, extrapolation_bc=Line())#Interpolations.linear_interpolation
    yint = itp.(xint)

    x_index = searchsortedfirst.(Ref(x0), xint)
    x_index = clamp!(x_index, 0, length(x0) - 1)

    dx = vcat(0, diff(x0))

    index = @. (dx[x_index] > maxgap)

    x_index = searchsortedfirst.(Ref(x0), xint)
    x_index = clamp!(x_index, 0, length(x0) - 1)
    dx = vcat(0, diff(x0))
    index = @. (index) & (dx[x_index] > maxgap)

    yint[index] .= NaN
    return yint
end