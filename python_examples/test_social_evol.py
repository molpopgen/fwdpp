#Simulate a population of N=500 diploids for 10N generations
#with theta = rho = 50.  Take a sample of size n = 20 chromosomes,
#and get the mean SFS
import social_evol
import random
import sys

SEED=int(sys.argv[1])
gsl_rng=social_evol.GSLrng(SEED)
NREPS=100
NSAM=20
N=500
GENS=10*N
THETA=50
RHO=THETA
mu=THETA/(4.*N)
mu_del = 0.1*mu
s=10/(2*N)
h=2
r=RHO/(4.*N)
SFS=[0]*(NSAM-1)
for i in range(0,NREPS):
    XX=social_evol.evolve(gsl_rng,N,GENS,mu,mu_del,s,h,r,0.1,0.1,0.005,0.0025)
    XX.sfs=social_evol.sfs(gsl_rng,XX,NSAM)
    for i in range(0,(NSAM-1)):
        SFS[i]=SFS[i]+XX.sfs[i]


#get mean of SFS
for i in range(0,(NSAM-1)):
    SFS[i]=float(SFS[i])/float(NREPS)
print(SFS)
