# A Revealed Preference Approach to Identification and Inference in Producer-consumer Models

using DelimitedFiles
using DataFrames
using CSV
using LinearAlgebra
using Statistics
using Random
using Distributions


#Set up the directory
rootdir=@__DIR__
dirdata=rootdir*"/Data"
dirresults=rootdir*"/Results"

####################################################################
#NielsenIQ data
####################################################################
subsample = "singles"

#1645 singles, 6364 couples, 3539 households size 3+
global const N = 1645     # Number of households
global const K = 4        # Number of goods
global T = 6              # Number of periods
#Load data
dum0=CSV.read(dir*"/paper_data.csv",DataFrame,datarow=2)
dum1=reshape(dum0[:,1],K,T,N)
dum2=reshape(dum0[:,2],K,T,N)
dum3=reshape(dum0[:,3],K,T,N)
@eval purch=$dum1[:,:,:]
@eval p=$dum2[:,:,:]
@eval trips=$dum3[:,:,:]

#Include procedures to get the MCMC
include(rootdir*"/MCMC_generator.jl")

#Number of trials
nsim = 200000000
#Number of draws and burns for the MCMC
npdraws = 2500
burns = 500
#length of the MCMC
cl = npdraws - burns

Random.seed!(123)
#Generate MCMC
genMCMC(nsim,npdraws,N,K,T,purch,trips,p)


##################################
#OPTIMIZATION
##################################
#Load the MCMC
cprices1=CSV.read(dirresults*"/MCMC_singles_cprices.csv",DataFrame,datarow=1)
cprices1=Matrix(cprices1);
cprices = reshape(cprices1,K,T,cl,N);

omega1=CSV.read(dirresults*"/MCMC_singles_search.csv",DataFrame,datarow=1)
omega1=Matrix(omega1);
omega = reshape(omega1,K,T,cl,N);

alpha1=CSV.read(dirresults*"/MCMC_singles_alpha.csv",DataFrame,datarow=1)
alpha1=Matrix(alpha1)
alpha = reshape(alpha1,K,cl,N);

#Include procedures for computing moments
include(rootdir*"/Moments.jl")

#Compute moments for each draw of the MCMC
chainM = logmoments30search(N,K,T,p,cprices,omega,collect(1:N))

#Include procedures for the optimization
include(rootdir*"/Optimization.jl")

#Initial guess of gamma
dg = 30
guessgamma = zeros(dg)

#Test the model
Random.seed!(123)
@time solve = Optim.optimize(guessgamma->objMCcu2(N,dg,guessgamma,chainM,cl), guessgamma, ParticleSwarm());
TSMC=2*Optim.minimum(solve)*N
solvegamma=Optim.minimizer(solve)

#Refinement
Random.seed!(123)
@time solve = Optim.optimize(solvegamma->objMCcu2(N,dg,solvegamma,chainM,cl), solvegamma, NelderMead());
TSMC=2*Optim.minimum(solve)*N




########################################
#Confidence sets on expected parameters
########################################

#Number of moments
dg = K*T + T + 1
#Subsample (default is all)
group = collect(1:N)

#Initial guess of gamma
guessgamma = zeros(dg)
#Lower bound, upper bound, and step size for inference over expected elasticity
lb,ub,step = -0.1,-0.2,0.025
#Indicate parameter to estimate
param = 1
#Recover confidence set on the expected elasticity w.r.t. search intensity
recovertheta(lb,ub,step,N,K,T,p,cprices,omega,group,alpha,param)


#Initial guess of gamma
guessgamma = zeros(dg)
#Lower bound, upper bound, and step size for inference over expected search cost
lb,ub,step = -10,-50,-10
#Indicate parameter to estimate
param = 2
#Recover confidence set on the expected search cost
recovertheta(lb,ub,step,N,K,T,p,cprices,omega,group,alpha,param)
