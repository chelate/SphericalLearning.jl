
# Lets check normalization
using GeometryTypes
using Makie
using Colors
using JuliennedArrays

export spherical_plot, complex_spherical_plot

function sampleSphere(n::Int)
    map(normalize, Slices(randn(3,n),1))
end

function complexImage(vals::Array{Complex{Float64}}) 
    max = maximum(abs.(vals))
    map(x -> HSV(angle(x)*360/(2π),(abs(x)/max)^.5,1),vals) 
end


function spherical_plot(fun::Function, n::Int)
    # n sets the resolution
    θ = [0;(0.5:n-0.5)/n;1]
    φ = [(0:2n-2)*2/(2n-1);2]
    x = [cospi(φ)*sinpi(θ) for θ in θ, φ in φ]
    y = [sinpi(φ)*sinpi(θ) for θ in θ, φ in φ]
    z = [cospi(θ) for θ in θ, φ in φ]
    vals = map(((x,y,z) -> fun( [x,y,z])), x,y,z);
    s = Makie.surface(x, y, z, color = complexImage(vals))
end

function comptosphere(z::Complex)
    # solve for the 
    x = real(z)
    y = imag(z)
    t = 2/(1+abs(z)^2)
    [t*x,t*y,1-t]
end

function complex_spherical_plot(fun::Function, n::Int)
    # plots a function of the complex plane on the riemann sphere
    function spheretocomp(r::Array)
        # solve for the complex plane point
        x = r[1]
        y = r[2]
        t = 1-r[3]
        z = (x + im*y)/t
    end
    spherical_plot((x->fun(spheretocomp(x))), n)
end
