#This file contains the procedures needed for the optimization step

using CUDA



function preobjMCcu2(n,nmem,gamma,chainM,cl,valf,geta,gtry,dvecM,logunif)
    #This function computes maximum-entropy moments


    #Distribute work on the GPU
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x

    #Compute maximum-entropy moments
    for i=index:stride:n
        for j=1:cl
            valf[i]=0.0
            for t=1:nmem
                gtry[i,t] = chainM[i,t,j]
                valf[i] += gtry[i,t]*gamma[t]-geta[i,t]*gamma[t]
            end
            for t=1:nmem
                geta[i,t] = logunif[i,j] < valf[i] ? gtry[i,t] : geta[i,t]
                dvecM[i,t] += geta[i,t]/cl
            end
        end
    end

    return nothing
end


function objMCcu2(n,nmem,gamma0,chainM,cl)
    #This function computes the test statistic


    #Define variables on the GPU
    chainMcu = cu(chainM)
    valf=cu(zeros(n))
    geta=cu(zeros(n,nmem))
    gtry=cu(zeros(n,nmem))
    @inbounds geta[:,:]=chainMcu[:,:,1]
    dvecM=cu(zeros(n,nmem))
    logunif=log.(CuArray(rand(n,cl)))
    gamma=cu(gamma0)

    #Define GPU architecture and compute maximum-entropy moments
    numblocks = ceil(Int, n/64)
    @cuda threads=64 blocks=numblocks preobjMCcu2(n,nmem,gamma,chainMcu,cl,valf,geta,gtry,dvecM,logunif)

    #Transfer maximum-entropy moments to the CPU and compute the average
    dvecM=Array(dvecM)*1.0
    dvec=sum(dvecM,dims=1)'/n

    #Compute covariance
    numvar=zeros(dg,dg)
    @simd for i=1:n
        BLAS.syr!('U',1.0/n,dvecM[i,:],numvar)
    end

    var=numvar+numvar'- Diagonal(diag(numvar))-dvec*dvec'

    (Lambda,QM)=eigen(var)
    inddummy=Lambda.>0
    An=QM[:,inddummy]

    dvecdum2=An'*(dvec)
    vardum3=An'*var*An
    Omega2=inv(vardum3)

    #Compute objective function
    Qn2=1/2*dvecdum2'*Omega2*dvecdum2

    return Qn2[1]
end


function recovertheta(lb,ub,step,N,K,T,p,cprices,omega,group,alpha,param)
    #This function computes the test statistic at various value of the expected elasticity
    #param is in {1,2} and determines which expected parameter is recovered

    for theta = lb:step:ub
        #Initialize random seed.
        Random.seed!(123)

        if param == 1
            chainM2 = logmoments30searchavgelasticity(N,K,T,p,cprices,omega,group,alpha,theta)
        elseif param == 2
            chainM2 = logmoments30searchcost(N,K,T,p,cprices,omega,group,alpha,theta)
        end

        #First Optimization
        Random.seed!(123)
        @time solve = Optim.optimize(guessgamma->objMCcu2(size(group,1),dg,guessgamma,chainM2,cl), guessgamma, ParticleSwarm());
        TSMC=2*Optim.minimum(solve)*size(group,1)
        solvegamma=Optim.minimizer(solve)

        #Second Optimization for refinement
        Random.seed!(123)
        @time solve = Optim.optimize(solvegamma->objMCcu2(size(group,1),dg,solvegamma,chainM2,cl), solvegamma, NelderMead());
        TSMC=2*Optim.minimum(solve)*size(group,1)
        guessgamma = Optim.minimizer(solve)

        #Write ouput as it goes for safety
        if param == 1
            open(dirresults*"/avgelasticity.csv", "a",) do io
                writedlm(io, [theta TSMC], ',')
                GC.gc()
            end
        elseif param == 2
            open(dirresults*"/avgsearchcost.csv", "a",) do io
                writedlm(io, [theta TSMC], ',')
                GC.gc()
            end
        end
    end

    return nothing
end
