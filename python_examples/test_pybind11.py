#Simulate a population of N=500 diploids for 10N generations
#with theta = rho = 50.  Take a sample of size n = 20 chromosomes,
#and get the mean SFS
import fwdpp_pybind11
import random
import sys

SEED=int(sys.argv[1])
gsl_rng=fwdpp_pybind11.GSLrng(SEED)
NREPS=100
NSAM=20
N=500
GENS=10*N
THETA=50
RHO=THETA
mu=THETA/(4.*N)
r=RHO/(4.*N)
SFS=[0]*(NSAM-1)
for i in range(0,NREPS):
    XX=fwdpp_pybind11.evolve(gsl_rng,N,GENS,mu,r)
    XXsfs=fwdpp_pybind11.sfs(gsl_rng,XX,NSAM)
    for i in range(0,(NSAM-1)):
        SFS[i]=SFS[i]+XXsfs[i]
    #Get the mutations that are not extinct
    #This shows how we are able to use pybind11
    #to expose C++ types as dicts, taking advantage
    #of pybind11's ability to auto-convert
    #C++'s vector to Python's list.
    muts=[m.as_dict() for m,n in zip(XX.mutations,XX.mcounts) if n > 0]

    #The next line is commented out,
    #but if you printed it, you'd see
    #that mutation containers in fwdpp
    #contain both extant and extinct mutations.
    #print len(muts)," ",len(XX.mutations)," ",len(XX.fixations)
    
    #Similarly, let's get the neutral mutations in the first gamete in each diploid
    for dip in range(len(XX.diploids)):
        #Again, use pybind11 to provide an easy means to implement as_dict() for a mutation behind the scenes
        muts_dip_gam1 = [XX.mutations[m].as_dict() for m in XX.gametes[XX.diploids[dip][0]].mutations]

#get mean of SFS
for i in range(0,(NSAM-1)):
    SFS[i]=float(SFS[i])/float(NREPS)
print(SFS)
