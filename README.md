# LTB-solution-for-perturbation-scale-factor-equation
This software, will solve a second order differential equation using Explicit Euler method, with a variable time lapse.
The reason that lead to use a variable time lapse, is the divergent behaviour of b(t) while in the limit t -> 0.
By using a variable time lapse, the max variation on b(t) first derivative may be set to a fixed precision, in our case 10^-4.
The equation, was derived using a LTB approach to a matter + cosmological constant universe, where a(r,t) = a(t) + b(t)f(r), on which a(t) is the standard FRW scale factor for a matter + Cosmological constant universe, while b(t)f(r) is a small perturbation term that will allow the resulting universe to have small radial dishomogeneities.
