export stochastic_process

function stochastic_process(t, tau, sigma, m, beta)
    time_2d = repeat(t', length(t))
    dt_all = time_2d - transpose(time_2d)
    cov = @. sigma ^ 2 * exp(-(abs(dt_all/tau)) ^ beta)
    u, sv, v = svd(cov)
    s = diagm(sv)
    # Random.seed!(2)
    n = rand(Normal(0, 1), length(t))
    y = m .+ (u * sqrt(s)) * n
    return y
end