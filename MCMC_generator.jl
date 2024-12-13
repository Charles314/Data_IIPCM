#This file contains procedures to sample solutions consistent with the model with maximum-entropy

using CUDA
using JuMP
using Clp



#Initialize optimal solution status
dum = Model(Clp.Optimizer);
@objective(dum,Min,0.0);
optimize!(dum);
optimal = termination_status(dum);

function constraintsQL(model,u,counterkr,K,T,purch,trips,pc,rhoc)
    #This function imposes the restrictions of the model in the optimization problem

    for i=1:T
        for j=1:T
            if i != j
                @constraint(model,u[j] - u[i] + (pc[:,j,counterkr]'*(purch[:,i,counterkr] - purch[:,j,counterkr]) + (rhoc[:,j,counterkr] .* purch[:,j,counterkr])'*(trips[:,i,counterkr] - trips[:,j,counterkr])) >= 0)
            end
        end
    end
end

function QL(counterkr,K,T,purch,trips,pc,rhoc)
    #This function sets up the optimization problem

    #Create model
    model = Model(Clp.Optimizer)
    set_silent(model::Model)

    #Create constrained variables
    @variable(model, u[i=1:T] >= 0.0)
    #Define the objective to optimize
    @objective(model,Min,0.0)
    #Add QL constraints
    constraintsQL(model,u,counterkr,K,T,purch,trips,pc,rhoc)
    #Solve the problem
    optimize!(model)

    return(termination_status(model))
end

function logpricessearch30(T,pc,wc,ind,counterkr,logpd,wd)
    #This function computes the acceptance ratio for the Metropolis-Hastings algorithm
    #and determines whether to accept or reject the candidate draw of latent variables

     logtrydens = (-sum(sum( (log.(p[:,1:T,ind]) - log.(pc[:,:,counterkr])).^2 ,dims=1),dims=2) + sum(sum( (log.(p[:,1:T,ind]) - logpd[:,:,counterkr]).^2 ,dims=1),dims=2) - sum( wc[:,:,counterkr].^2,dims=1:2) + sum( wd[:,:,counterkr].^2 ,dims=1:2) )[1,1,1]
     alpha = log(rand()) < min(0, logtrydens)
     return(alpha)
end

function update_MCMC_chain(alpha,indexind,keeprun,pck,rhock,wck,elasticityk,ind,counterkr,  logpd,rhod,wd,elasticityd,logpd0,rhod0,wd0,elasticity0,soldet0,truerho,truep,truew,trueelasticity,R)
    #Update default draw and add realization to the MCMC chain

    if alpha == true
        #Update default
        logpd0[:,:,ind] = log.(pck[:,:,counterkr])
        elasticity0[:,ind] = elasticityk[:,counterkr]
        wd0[:,:,ind] = wck[:,:,counterkr]
        #Accept new draw
        push!(truep[ind], pck[:,:,counterkr])
        push!(trueelasticity[ind],elasticityk[:,counterkr])
        push!(truew[ind],wck[:,:,counterkr])

        #Keep track of
        soldet0[ind, findlast(soldet0[ind,:]) + 1] = true
    end
    if alpha == false
        #Accept default draw
        push!(truep[ind],exp.(logpd[:,:,counterkr]))
        push!(trueelasticity[ind],elasticityd[:,counterkr])
        push!(truew[ind],wd[:,:,counterkr])

        soldet0[ind, findlast(soldet0[ind,:]) + 1] = true
    end
    #Remove individual if chain is complete
    if sum(soldet0[ind,:]) > R
        keeprun = keeprun[1:end .!= indexind]
    end

    return(keeprun)
end

function metropolis_hastings(n, keeprun, K, T, pck, rhock, wck, elasticityk,ind, counterkr,  logpd,rhod,wd,elasticityd,logpd0,rhod0,wd0,elasticity0,soldet0,truerho,truep,truew,trueelasticity,R)
    #This function runs the Metropolis-Hastings algorithm

    #Get the index of ind in keeprun
    indexind = findall(keeprun .== ind)[1]
    #Metropolis-Hastings
    alpha = logpricessearch30(T,pck,wck,ind,counterkr,logpd,wd)
    keeprun = update_MCMC_chain(alpha,indexind,keeprun,pck,rhock,wck,elasticityk,ind,counterkr,  logpd,rhod,wd,elasticityd,logpd0,rhod0,wd0,elasticity0,soldet0,truerho,truep,truew,trueelasticity,R)

    return(keeprun)
end

function cpu_power_QL(n, keeprun, K, T, purch, trips, pck, rhock, wck, elasticityk,  logpd,rhod,wd,elasticityd,logpd0,rhod0,wd0,elasticity0,soldet0,truerho,truep,truew,trueelasticity,R)
    #This function checks whether the candidate draw of latent variables is consistnet with the model and decides whether to include it in the MCMC

    #Position in vector of nleft individals
    counterkr = 0

    for ind=1:n
        #Check if QL passes for some draw
        if ind in keeprun && sum(soldet0[ind,:]) <= R
            counterkr += 1
            #QL
            sol = QL(counterkr,K,T,purch,trips,pck,rhock)
            s = sol == optimal

            #If so, no more draws to check
            if s == true
                keeprun = metropolis_hastings(n, keeprun, K, T, pck, rhock, wck, elasticityk, ind, counterkr,  logpd,rhod,wd,elasticityd,logpd0,rhod0,wd0,elasticity0,soldet0,truerho,truep,truew,trueelasticity,R)
                break
            end
        end
    end
    return(keeprun)
end

function MCMC_QL_search(nsim, R, n, K, T, purch, trips, p, group)
    #This function generates the MCMC


    #Store draws satisfying the model and went through the Metropolis-Hastings algorithm
    truep = [Array{Float64}[] for i=1:n]
    truerho = [Array{Float64}[] for i=1:n]
    trueelasticity = [Array{Float64}[] for i=1:n]
    truew = [Array{Float64}[] for i=1:n]

    #Determines whether the MCMC has reached its target length
    soldet0 = [trues(n,1) falses(n,R+1)]
    #Determines the set of individuals that hasn't reached the targeted length of the MCMC
    keeprun = group

    #Initial guess of default draw
    logpd0 = zeros(K,T,n)
    logpd0[1,:,:],logpd0[2,:,:],logpd0[3,:,:],logpd0[4,:,:] = log.(p[1,:,:]), log.(p[2,:,:]), log.(p[3,:,:]), log.(p[4,:,:])
    rhod0 = -ones(K,T,n)
    wd0 = ones(K,T,n)
    elasticity0 = -ones(K,n)

    #Generate the MCMC by rejection sampling; try new candidates latent variables until chain is of size R
    for sim=1:nsim

        #If number of total draws not attained for everyone
        if size(keeprun,1) > 0

            #Set number of individuals with incomplete MCMC chain
            nleft = size(keeprun)[1]
            #Update default draw
            logpd = logpd0[:,1:T,keeprun]
            rhod = rhod0[:,1:T,keeprun]
            wd = wd0[:,1:T,keeprun]
            elasticityd = elasticity0[:,keeprun]

            ######################
            #Guess true variables
            ######################
            #Draw elasticity of price with respect to shopping time. Support (0,1) ensures bound on rho is satisfied.
            elasticityc = -reshape( 1.0*rand.(K*nleft), (K,nleft))
            #Draw log prices
            logpc = zeros(K,T,nleft)
            logpc[1,:,:],logpc[2,:,:],logpc[3,:,:],logpc[4,:,:] = reshape(rand.(Normal.(logpd[1,:,:][:],0.35)),(T,nleft)), reshape(rand.(Normal.(logpd[2,:,:][:],0.35)),(T,nleft)), reshape(rand.(Normal.(logpd[3,:,:][:],0.35)),(T,nleft)), reshape(rand.(Normal.(logpd[4,:,:][:],0.35)),(T,nleft))
            pc = exp.(logpc)
            #Draw intercept
            intercept = zeros(K,nleft)
            intercept1,intercept2,intercept3,intercept4 = mean(logpc[1,:,:],dims=1)[:] .- (elasticityc[1,:]).*(mean(log.(trips[1,:,nleft]),dims=1)[:]), mean(logpc[2,:,:],dims=1)[:] - (elasticityc[2,:]).*(mean(log.(trips[2,:,nleft]),dims=1)[:]), mean(logpc[3,:,:],dims=1)[:] - (elasticityc[3,:]).*(mean(log.(trips[3,:,nleft]),dims=1)[:]), mean(logpc[4,:,:],dims=1)[:] - (elasticityc[4,:]).*(mean(log.(trips[4,:,nleft]),dims=1)[:])
            intercept[1,:], intercept[2,:], intercept[3,:], intercept[4,:] = intercept1, intercept2, intercept3, intercept4
            #Get implied search ability and dpda
            rhoc = zeros(K,T,nleft)
            wc = zeros(K,T,nleft)
            for t=1:T
                wc[:,t,:] = intercept .+ elasticityc.*log.(trips[:,t,keeprun]) .-  logpc[:,t,:]
                rhoc[:,t,:] =(elasticityc.*exp.(intercept)).*((trips[:,t,keeprun]).^(elasticityc .-1)).*(exp.(-wc[:,t,:]))
            end

            #Sequentially check the restrictions of the model and apply Metropolis-Hastings
            keeprun = cpu_power_QL(n, keeprun, K, T, purch[:,:,keeprun], trips[:,:,keeprun], pc, rhoc, wc, elasticityc,  logpd,rhod,wd,elasticityd,logpd0,rhod0,wd0,elasticity0,soldet0,truerho,truep,truew,trueelasticity,R)

        end
        #If nobody left to check, stop
        if size(keeprun,1) <= 0
            break
        end
    end

    return(truep,trueelasticity,truew)
end

function genMCMC(nsim, npdraws, N, K, T, purch, trips, p)
    #This function constructs the MCMC and saves the results as it goes
    #Beware, obtaining the MCMC may take some time

    #Generate the MCMC
    for i=1:N
        trueprices1 = MCMC_QL_search(nsim, npdraws, N, K, T, purch, trips, p, [i]);

        cprices = Float64[]
        alpha = Float64[]
        search = Float64[]

        if size(trueprices1[1][i],1) >= npdraws
            for j=1:npdraws-burns
                append!(cprices,vec(trueprices1[1][i][end+1-j]))
                append!(alpha,vec(trueprices1[2][i][end+1-j]))
                append!(search,vec(trueprices1[3][i][end+1-j]))
            end
        end

        open(dirresults*"/MCMC_"*subsample*"_cprices1.csv", "a",) do io
            writedlm(io, cprices, ',')
            GC.gc()
        end;
        open(dirresults*"/MCMC_"*subsample*"_alpha1.csv", "a",) do io
            writedlm(io, alpha, ',')
            GC.gc()
        end;
        open(dirresults*"/MCMC_"*subsample*"_search1.csv", "a",) do io
            writedlm(io, search, ',')
            GC.gc()
        end;

        #Progress tracker
        if mod(i,25) == 0
            print( round((i/N)*100, digits=2), "%  ")
        end
    end
    println("")

    print("MCMC complete!")

    return nothing
end
