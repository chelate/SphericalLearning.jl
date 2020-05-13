include("SphericalTensor.jl")

using LinearAlgebra
using SpecialFunctions



export helmholtzG, helmholtzGapprox, hologram

function sphericalbesselj(nu, x::T) where {T}
    besselj_nuhalf_x = besselj(nu + one(nu)/2, x)
    if abs(x) ≤ sqrt(eps(real(zero(besselj_nuhalf_x))))
        nu == 0 ? one(besselj_nuhalf_x) : zero(besselj_nuhalf_x)
    else
        √((float(T))(π)/2x) * besselj_nuhalf_x
    end
end


sphericalbessely(nu, x::T) where {T} = √((float(T))(π)/2x) * bessely(nu + one(nu)/2, x)

sphericalhankel1(nu, x::T) where {T} = sphericalbesselj(nu, x) + im*sphericalbessely(nu, x)


function helmholtzG(x,y;k=1)
    d = norm(x .- y)
    exp(im*k*d)/(4π*d)
end

function sphericalbesselj(nu, x::T) where {T}
    besselj_nuhalf_x = besselj(nu + one(nu)/2, x)
    if abs(x) ≤ sqrt(eps(real(zero(besselj_nuhalf_x))))
        nu == 0 ? one(besselj_nuhalf_x) : zero(besselj_nuhalf_x)
    else
        √((float(T))(π)/2x) * besselj_nuhalf_x
    end
end

function helmR(l,x;k=1) # inside solution
    conj(sphericalbesselj(l, k*norm(x)))*
        harmonic(l, x)
end

function helmS(l,x;k=1) # outside solution
    sphericalhankel1(l, k*norm(x))*
        harmonic(l, x)
end

function hologram(l,x;k=1,R=1) # projection on sphere 
    # outputs a spherical tensor a, order l
    # can reconstruct it's portion of the image with harmonic(a,x)
    if norm(x) < R
        im*k*conj(helmR(l,x;k=1))*sphericalhankel1(l,k*R)
    else
        im*k*conj(sphericalbesselj(l, k*R)*helmS(l,x;k=1))
    end
end

function helmholtzGapprox(x,y;k=1,lmax=10)
    if norm(x) < norm(y)
        sum(im*k*dot(helmR(l,x;k=k),helmS(y,l;k=k)) for l in 0:lmax)
    else
        sum(im*k*dot(helmR(l,y;k=k),helmS(x,l;k=k)) for l in 0:lmax)
    end
end

