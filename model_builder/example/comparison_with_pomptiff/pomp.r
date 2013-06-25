require('pomp')

DELTA_STO <- 14

logit <- function(x){
    return(log(x/(1.0-x)))
}

inv.logit <- function(x){
    return (1.0/(1.0+exp(-x)))
}

sir.proc.sim <- function (x, t, params, delta.t, ...) {
    ## unpack the parameters
    N <- params["N"] # population size
    mu <- params["mu"] # birth rate = death rate
    v <- 1/(0.00001 + params["v"]/7) # recovery rate (user parametrize v as a duration in days
    beta <- params["r0"]*v # contact rate
    foi <- beta*x["I"]/N # the force of infection
    trans <- c(
               rpois(n=1,lambda=mu*N*delta.t), # births are assumed to be Poisson
               reulermultinom(n=1,size=x["S"],rate=c(foi,mu),dt=delta.t), # exits from S
               reulermultinom(n=1,size=x["I"],rate=c(v,mu),dt=delta.t), # exits from I
               reulermultinom(n=1,size=x["R"],rate=c(mu),dt=delta.t) # exits from R
               )
    ## now connect the compartments
    x[c("S","I","R","cases")]+c(
                                trans[1]-trans[2]-trans[3],
                                trans[2]-trans[4]-trans[5],
                                trans[4]-trans[6],
                                trans[4] # accumulate the recoveries
                                )
}

data <- read.csv('data/data.csv', header=TRUE)[,2]
N_DATA <- length(data)

sir <- pomp(
            data=data.frame(
            time= 1:N_DATA,
            reports=data
            ),
            times="time",
            t0=0,
            rprocess=euler.sim(
            step.fun=sir.proc.sim,
            delta.t=1/DELTA_STO
            ),
            parameter.transform=function(params,...){

                params[c('r0', 'v')] <- exp(params[c('r0', 'v')])
                params[c('I0')] <- inv.logit(params[c('I0')])

                params
            },
            parameter.inv.transform=function(params,...){

                params[c('r0', 'v')] <- log(params[c('r0', 'v')])
                params[c('I0')] <- logit(params[c('I0')])

                params
            },
            rmeasure =  function(x,t,params,...){

                yobs = rnorm(n=1, mean=params['rep2']*x['cases'], sd=sqrt((1-params['rep2'])*params['rep2']*x['cases'] + (params['phi']*params['rep2']*x['cases'])^2))
                if(yobs>0){
                    return(yobs)
                } else {
                    return(0.0)
                }

            },
            dmeasure = function(y,x,t,params,log,...){

                my.mean=params['rep2']*x['cases']
                my.sd=sqrt((1-params['rep2'])*params['rep2']*x['cases'] + (params['phi']*params['rep2']*x['cases'])^2)

                if(y>0.0){
                    p = pnorm(y+0.5, mean = my.mean, sd = my.sd, log.p = FALSE)-pnorm(y-0.5, mean = my.mean, sd = my.sd, log.p = FALSE)
                } else{
                    p = pnorm(y+0.5, mean = my.mean, sd = my.sd, log.p = FALSE)
                }

                if(log){
                    return(log(p))
                } else {
                    return(p)
                }

            },
            initializer=function(params, t0, ic.pars, comp.names, ...){
                x0 <- c(S=0,I=0,R=0,cases=0)
                N <- params["N"]
                fracs <- params[ic.pars]
                x0[comp.names] <- round(N*fracs/sum(fracs))

                x0
            },
         
            zeronames=c("cases"), # 'cases' is an accumulator variable
            ic.pars=c("S0","I0","R0"), # initial condition parameters
            comp.names=c("S","I","R") # names of the compartments
            )


theta.guess <- c(
                 N=1000000,
                 r0=20,
                 v=11,
                 mu=0.00027,
                 rep2=0.6,
                 phi=0.1,
                 S0=0.07,
                 I0=1e-05,
                 R0=1-0.07-1e-05
                 )

rw.sd <- c(
           r0=0.02,
           v=0.02,
           I0=0.01
           )/sqrt(N_DATA)

estpars <- c('r0', 'v')
estic <- c('I0')

simulate(
         sir,
         params=theta.guess,
         seed=677573454L
          ) -> simul

plot(simul)

pf <- pfilter(sir, params=theta.guess, Np=100, max.fail=1000, tol=1e-17)
logLik(pf)

mif <- mif(
           tol=1e-17,
           sir,
           Nmif=50,
           start=theta.guess,
           transform=TRUE,
           pars=estpars,
           ivps=estic,
           rw.sd=rw.sd,
           Np=100,
           var.factor=2,
           ic.lag=39,
           method='mif',
           cooling.factor=0.975,
           max.fail=2000,
           verbose=FALSE,
           seed=677573454L
           )

plot(mif)

logLik(mif)
res <- conv.rec(mif)

layout(matrix(1:6,3,2))
for (p in estpars){
    plot(exp(res[,p]), type='l', ylab=p)
}

for (p in estic){
    plot(inv.logit(res[,p]), type='l', ylab=p)
}

plot(res[,'loglik'], type='l', ylab='loglik')
