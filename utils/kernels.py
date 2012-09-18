#!/opt/local/bin/python2.7
from __future__ import division
from sympy import *
x, q = symbols('x q')
f = symbols('f',cls=Function)
term1 = sympify((2-q)**4)
term2 = -5*(sympify(6)/5 - q)**4
term3 = 10*(sympify(2)/5 - q)**4
f = Piecewise((term1 + term2 + term3,q < sympify(2)/5), (term1 + term2, q < sympify(6)/5), (term1, q < 2), (0, q > 2), (0, True))
print "\nW:"
print f
print "\nFirst derivative:"
g = diff(f,q)
print g
print "\n2nd derivative:"
h = diff(g,q)
print h
print "\n1D normalisation:"
norm1D = sympify(1)/(2*integrate(f,(q,0,2)))
print norm1D
print "\n2D normalisation:"
norm2D = sympify(1)/(integrate(2*pi*q*f,(q,0,2)))
print norm2D
print "\n3D normalisation:"
norm3D = sympify(1)/(integrate(4*pi*q*q*f,(q,0,2)))
print norm3D
