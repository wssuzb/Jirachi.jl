export sigma_G, sigma_Gz, test

sigma_G(mag::T) where{T} = 0.741 * iqr(mag)

function sigma_Gz(mag::T, err::T) where{T}
    
    # ...

end

function test()
    println("hello")
end
