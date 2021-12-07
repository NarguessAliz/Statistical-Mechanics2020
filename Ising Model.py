import numpy as np
import matplotlib.pyplot as plt

#a n*n random square lattice:
n= 50
a= np.ones ([n,n])
for i in range (n):
	for j in range (n):
		x= np.random.rand()
		if x>0.5:
			a[i][j]= -1		
#print (a)

#hamiltonian of system
J=1                                     #normalization coef
B=0                                     #magnetic field 
mu=0.33                                 #
def H(a):				#first neighbor
	h=0
	hh=0
	for i in range (n):
		for j in range (n):
			h+= -J*a[i,j]*((a[(i+1)%n,j])+(a[(i-1)%n,j])+(a[i,(j-1)%n])+(a[i,(j+1)%n]))
			hh+= a[i,j]*B*mu
                
	return h + hh

#variation in energy
def dE(a,i,j):
	d=-2*J*a[i,j]*((a[(i+1)%n,j])+(a[(i-1)%n,j])+(a[i,(j-1)%n])+(a[i,(j+1)%n]))
	return d

#magetization:
def m(a):
	m= 0
	for i in range (n):
		for j in range (n):
			m+=a[i,j]
	return m
	
c=0                                             #counter of accepted rotatese
k= 1
T= 0.2					        #temprature
ham0= H(a)					#energy before rotation
for i in range (10000):
	p= int(n*(np.random.rand()))
	q= int(n*(np.random.rand()))
	(a[p][q])= (a[p][q])*(-1)
	#ham1= H(a)					#energy after rotation
	delta= dE(a,p,q)
	if delta >0.0:
		p1= np.random.rand()
		prob= np.exp(-delta/(k*T))		#probability of acceptation dependent with Temparature
		if p1>prob:
			a[p][q]= a[p][q]*(-1)
	else:
		c+=1
		ham0= ham0+delta
	#print (c,"	",i, "	",ham0)
plt.imshow(a,  interpolation="nearest")
plt.suptitle('Ising Model', fontsize=14, fontweight='bold')
plt.title(f'Temparature= {T}K')
plt.show()


