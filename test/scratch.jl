using SphericalLearning
using Gadfly
##


##
spherical_plot((x->harmonic(hologram(4,[0,0,.9]),x)),50)
##
spherical_plot((x->helmholtzG([0,0,.9],x;k=5)),50)

##

x = randn(3)
y = randn(3)


## these should show that the functions
Gadfly.plot([
    (z->imag(helmholtzG(z*x,y,k=1+0.2im))  ),
    (z->imag(helmholtzGapprox(z*x,y,k=1+0.2im))),
    (z->real(helmholtzG(z*x,y,k=1+0.2im))  ),
    (z->real(helmholtzGapprox(z*x,y,k=1+0.2im)) )
    ],0.2,7)
##
Gadfly.plot([
    (x->imag(sum(harmonic(hologram(l,[0,0,.9]),[0,1,x]) for l in 0:12))),
    (x->imag(helmholtzG([0,0,.9],normalize([0,1,x])))+0.01),
    (x->real(sum(harmonic(hologram(l,[0,0,.9]),[0,1,x]) for l in 0:12))),
    (x->real(helmholtzG([0,0,.9],normalize([0,1,x]))))
    ], -1,1)