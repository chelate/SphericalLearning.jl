using SphericalLearning
using Gadfly
##
x = randn(3)
y = randn(3)

##
spherical_plot((x->harmonic(hologram(4,[0,0,.9]),x)),50)
##
spherical_plot((x->helmholtzG([0,0,.9],x;k=5)),50)
##
Gadfly.plot([
    (z->imag(helmholtzG(z*x,y,k=1+0.2im))  ),
    (z->imag(helmholtzGapprox(z*x,y,k=1+0.2im)) ),
    (z->real(helmholtzG(z*x,y,k=1+0.2im))  ),
    (z->real(helmholtzGapprox(z*x,y,k=1+0.2im)) )
    ],0.2,7)
