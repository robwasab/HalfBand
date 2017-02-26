import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import pdb
pi = np.pi
sqrt = np.sqrt
tan = np.tan
sin = np.sin
cos = np.cos
power = np.power

# ELLIPTIC FILTER DESIGN FOR A CLASS OF GENERALIZED
#               HALFBAND FILTERS
#
#               BY RASHID ANSARI
#
#                INTEPRETED BY
#                  R. J. TONG
# 
# CITATION:
#    R. Ansari, "Elliptic filter design for a class of generalized halfband filters," 
#    in IEEE Transactions on Acoustics, Speech, and Signal Processing, vol. 33, no. 5, pp. 1146-1150, Oct 1985.
#    doi: 10.1109/TASSP.1985.1164709
#    URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1164709&isnumber=26193
#
# DATE:
#    2/26/2017
# 
# INPUTS:
#    fp = pass band freq. 0<fp<0.4999 i.e 0.44
#    dp = pass band ripple i.e. 0.0001
#    ds = stop band ripple i.e. 0.1

fp = 0.495
dp = 0.001
ds = 0.001

wp = fp * pi
ws = pi - wp

print 'wp: %.3e*pi'%(wp/pi)
print 'ws: %.3e*pi'%(ws/pi)

d1 = ds
d2 = sqrt(2.0*dp - power(dp,2.0))
d = min([d1, d2])/2.0
d_prime = 1.0 - sqrt(1.0 - power(d,2.0))

print 'd1: %f'%d1
print 'd2: %f'%d2
print 'd : %f'%d
print 'd\': %e'%d_prime


k = power(tan(wp/2.0),2.0)
k_prime = sqrt(1.0 - power(k,2.0))

print 'k : %.3e'%k
print 'k\': %.3e'%k_prime

delta = (1.0 - power(d,2.0))/power(d,2.0)
p = 0.5 * (1.0 - sqrt(k_prime)) / (1.0 + sqrt(k_prime))
q = p + 2.0*power(p,5.0) + 15.0*power(p,9.0) + 150.0*power(p,13.0)

print 'p : %.6e'%p
print 'q : %.6e'%q
N = int(np.ceil(2.0*np.log(4.0*delta)/-np.log(q)))
if N % 2 == 0:
	print 'N is even, adding +1'
	N += 1

L = (N - 1)/2

print 'N: %d'%N
print 'L: %d'%L
print 'delta: %.3e'%delta

omegas = np.zeros(L)
m_max = 10 
for n in range(L):
	i = n + 1.0
	num_sigma = 0.0
	#print 'Calculating numerator...'
	for m in range(0, m_max):
		A = (-1.0)**m
		B = power(q,m*(m+1.0))
		C = sin((2.0*m + 1.0)*pi*i/float(N))
		num_sigma += A*B*C
		pad = ' '*m
		#print '%s->%e'%(pad,num_sigma)

	den_sigma = 0.0
	#print 'Calculating denominator...'
	for m in range(1, m_max):
		A = (-1.0)**m
		B = power(q, float(m*m))
		C = cos(2.0*m*pi*i/float(N))
		den_sigma += A*B*C
		pad = ' '*(m-1)
		#print '%s->%e'%(pad,den_sigma)

	omega_i = 2.0*power(q,0.25) * num_sigma / (1.0 + 2.0*den_sigma)
	omegas[n] = omega_i
	print 'omega[%d]: '%n, omega_i

omegas_2 = omegas * omegas
ris = sqrt((1.0-k*omegas_2)*(1.0-omegas_2/k))
cos_thetas = (-1.0)**((np.arange(L)+1.0)+1.0)*ris/(1.0+omegas**2.0)
ais = (1.0-cos_thetas)/(1.0+cos_thetas)
print 'ai: %s'%(str(ais))

h0_ai = ais[ais >= 1.0]
h1_ai = ais[ais < 1.0]

print 'H0(z): %s'%(str(h0_ai))
print 'H1(z): %s'%(str(h1_ai))

pts = 1000
F = np.arange(pts)/float(pts) * 0.5
z = np.exp(2.0j*pi*F)
z_2 = power(z, 2.0)

h0s = [(z_2 + ai) / (ai*z_2 + 1.0) for ai in h0_ai]
h1s = [(ai*z_2 + 1.0) / (z_2 + ai) for ai in h1_ai]

h0_total = np.ones(pts, dtype=np.complex64)
for h0 in h0s:
	h0_total *= h0

h1_total = np.ones(pts, dtype=np.complex64)
for h1 in h1s:
	h1_total *= h1

hs = 1.0/z*h0_total + h1_total

plt.plot(F, 20.0*np.log10(np.abs(hs)))
plt.show()
