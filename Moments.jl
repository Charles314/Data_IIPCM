#This file contains procedures to re-order the MCMC for the optimization step

using Optim



function logmoments30search(n,K,T,p,pc,omega,group)
    chainM = ones(size(group,1),T*K + T,cl)
    counter = 0

    #Moments on measurement error
    for k=1:K,t=1:T
        counter += 1
        for s=1:cl
            @inbounds chainM[:,counter,s] = (log.(p[k,t,group]) - log.(pc[k,t,s,:]))[:]
        end
    end
    #Moments on productivity shocks
    for t=1:T, s=1:cl
        @inbounds chainM[:,T*K + t,s] = (1/K)*sum(omega[:,t,s,:],dims=1)
    end
    return(chainM)
end


function logmoments30searchavgelasticity(n,K,T,p,pc,omega,group,alpha,theta0)
    chainM = ones(size(group,1),T*K + T + 1, cl)
    counter = 0

    #Moments on measurement error
    for k=1:K,t=1:T
        counter += 1
        for s=1:cl
            @inbounds chainM[:,counter,s] = (log.(p[k,t,group]) - log.(pc[k,t,s,group]))[:]
        end
    end
    #Moments on productivity shocks
    for t=1:T, s=1:cl
        @inbounds chainM[:,T*K + t,s] = (1/K)*sum(omega[:,t,s,group],dims=1)
    end
    #Moment on the expected average elasticity of price with respect to shopping intensity
    for s=1:cl
        @inbounds chainM[:,end,s] = (1/K)*(1/K)*sum(alpha[:,s,group],dims=1)[:] - theta0*ones(size(group,1))
    end
    return(chainM)
end


function logmoments30searchcost(n,K,T,p,pc,omega,group,alpha,theta0)
    chainM = ones(size(group,1),T*K + T + 1, cl)
    counter = 0

    #Moments on measurement error
    for k=1:K,t=1:T
        counter += 1
        for s=1:cl
            @inbounds chainM[:,counter,s] = (log.(p[k,t,group]) - log.(pc[k,t,s,group]))[:]
        end
    end
    #Moments on productivity shocks
    for t=1:T, s=1:cl
        @inbounds chainM[:,T*K + t,s] = (1/K)*sum(omega[:,t,s,group],dims=1)
    end
    #Moment on expected search cost
    for s=1:cl
        searchcostt = zeros(size(group,1))
        for t=1:T
            @inbounds searchcostt .+= (1/T)*(1/K)*sum(  ( (alpha[:,s,group] .* pc[:,t,s,group] .* purch[:,t,group]) ) ,dims=1  )[:]
        end
        @inbounds chainM[:,end,s] = searchcostt - theta0*ones(size(group,1))
    end
    return(chainM)
end
