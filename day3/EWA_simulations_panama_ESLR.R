library(rethinking)
library(truncnorm)
#require(RColorBrewer)
#col.pal <- brewer.pal(3, "Dark2")
col.pal=c("#1B9E77", "#D95F02", "#7570B3") #graphing color pallette

##create a softmax function to simply code
Softmax <- function(x){
exp(x)/sum(exp(x))
} 

#data sims
n <- 50                                     #number of individuals/pop size
nbouts <- 100                               #timesteps

#simulate values for options
techmeans <- c( 4 , 4.2 , 5)                  #mean efficiency of techniques
techvar <- c( 1, 1 , 1)                         #variance of techniques

#plot to visualize overlap of payoffs
dens(rnorm( 10000 , mean=techmeans[1] , sd=techvar[1] ) ,col=col.pal[1] , xlim=c(0,10) )
dens(rnorm( 10000 , mean=techmeans[2] , sd=techvar[2] ) ,col=col.pal[2] , xlim=c(0,10) , add=TRUE )
dens(rnorm( 10000 , mean=techmeans[3] , sd=techvar[3] ) ,col=col.pal[3] , xlim=c(0,10) , add=TRUE )

#parameter sims
fc.sim <- log(1.5)				            ##frequency dependecy parameter on log scale
phi.sim <- logit(0.18)                       ## attraction updating parameter on log-odds scale
gamma.sim <- logit(0.20)                     ## weight of social info parameter on log-odds scale
k.lambda <- 1                             ##sensitivity to individual payoffs


#varying effects offsets for individuals
gamma.sim_i <- rnorm( n , mean=0 , sd=0.5 ) 
phi.sim_i <- rnorm( n , mean=0 , sd=0.5 ) 
fc.sim_i <- rnorm( n , mean=0 , sd=0.5 ) 

#unique parameters for each individual, to visualize heterogeneity and for plotting individual predictions
gamma.sim_id <- round( logistic(gamma.sim + gamma.sim_i), digits=2) ##simulated gammas for all n individuals
phi.sim_id <- round(logistic(phi.sim +phi.sim_i), digits=2)  ##simulated phis for all n individuals
fc.sim_id <- round(exp(fc.sim + fc.sim_i), digits=2)  ##simulated strength of frequency dependent learning for all n individuals
gamma.sim_id
phi.sim_id
fc.sim_id

#begin to simulate data
dsim_s <- data.frame( id=0 , bout=0 , tech=0 , y1=0 , y2=0, y3=0 , s1=0 , s2=0 , s3=0 , A1=0 , A2=0 , A3=0 , Pr1=0 , Pr2=0 , Pr3=0 )
therow <- 1

AC <- matrix(1,ncol=3,nrow=n) #attraction scores for each tech
AC[,3] <- 0.2 #change to make high payoff option less likely
Softmax(AC[1,]) #run to see initial prob of choosing a behavior

S1 <- S2 <- S3 <- rep(0,n+1) # num of individuals choosing each tech in previous bout
PS1 <- PS2 <- PS3 <- rep(0,nbouts+1) # empty vector for mean observed in previous rounds
s_temp <-  rep(0,3)
# S1[1] <- 0.5*nbouts
# S2[1] <- 0.25*nbouts
# S3[1] <- 0.25*nbouts

for ( r in 1:nbouts ) {
    for ( i in 1:n ) {  

		prtech_i <-  Softmax(k.lambda*AC[i,]) 
        my.gam <- logistic( gamma.sim + gamma.sim_i[i] ) #social info weight for individual i
        my.phi <- logistic( phi.sim + phi.sim_i[i] ) #social info weight for individual i
        my.fconf <- exp( fc.sim + fc.sim_i[i])  #strength of frequency dependence for individual i
        prtech_su <- c(S1[r],S2[r],S3[r]) #attraction score for individual i

        #//conformity aspect below
        if ( r > 1 ) {
            if (sum( prtech_su ) > 0 ) {

                #// compute frequency cue
                for ( j in 1:3 ){ s_temp[j] <- prtech_su[j]^my.fconf}

                prtech_s <- s_temp/sum(s_temp)
                prtech <- (1-my.gam)*prtech_i + my.gam*prtech_s

            } else {
                prtech <- prtech_i
            }
        } else {
            prtech <- prtech_i
         }
# choose tech
        tech <- sample( 1:3 , size=1 , prob=prtech)
        yield <- rtruncnorm( 1 , a=0 , b=Inf, mean=techmeans[tech] , sd=techvar[tech] )
# update attractions
        yields <- rep(0,3)
        yields[tech] <- yield#makes payoff yield
        for (k in 1:3){
        	AC[i,k] <- (1-my.phi)*AC[i,k] + my.phi*yields[k]
        }

        dsim_s[therow,] <- c( i , r , tech , yields[1] , yields[2] , yields[3] , S1[r] , S2[r] , S3[r] , AC[i,1] , AC[i,2] , AC[i,3] ,  prtech[1] , prtech[2], prtech[3] )
        therow <- therow + 1
    } #i
    S1[r+1] <- length( dsim_s$tech[dsim_s$tech==1 & dsim_s$bout==r] )
    S2[r+1] <- length( dsim_s$tech[dsim_s$tech==2 & dsim_s$bout==r] )
    S3[r+1] <- length( dsim_s$tech[dsim_s$tech==3 & dsim_s$bout==r] )

 }

o <- order( dsim_s$i )
dsim <- dsim_s[o,]

#plot raw data of group level effects
plot(s1/n ~ bout, data=dsim, col=col.pal[1] , ylim=c(0,1) , xlim=c(2,nbouts+1), pch=19 , xlab="Time" , ylab="Proportion of Individuals Choosing Option" , main="Population Mean ")
points(s2/n ~ bout, data=dsim , col=col.pal[2], pch=19)
points(s3/n ~ bout, data=dsim , col=col.pal[3], pch=19)
legend("topleft", cex=1 , as.character(techmeans), pch=19 ,col=col.pal, horiz=TRUE , bty="y")


# ds <- list(
# N = nrow(dsim),
# J = length( unique(dsim$i)),
# K=max(dsim$tech),
# tech = dsim$tech,
# y = cbind( dsim$y1 , dsim$y2 , dsim$y3 ),
# s = cbind(dsim$s1 , dsim$s2 , dsim$s3),
# id = dsim$i,
# bout = dsim$bout ,
# N_effects=3
# )
# parlistcombo=c("lambda" ,"a_id" , "mu" , "fconf", "dev" , "log_lik")

# fit_combo <- stan( file = 'PN_social_combo.stan', data = ds , 
#     iter = 2000, warmup=1000, chains=2, cores=2, pars=parlistcombo, 
#     control=list( adapt_delta=0.98 ) )

# post <- extract(fit_combo)
# par(mfrow=c(2, 2), oma = c(0, 0, 0, 0) , mar=c(3,.5,2,.5))
# par(cex = 0.6)
# par(tcl = -0.25)
# par(mgp = c(2, 0.6, 0))

# dens(logistic(post$mu[,1]) , main=expression(paste(phi)) , xlim=c(0,1), ylab='' , xlab= "a. weight of new experience" , col="white", yaxt='n' , cex.lab=1.5)##phi
# abline( v=logistic(phi.sim) , col="cornflowerblue" ) 
# shade( density(logistic(post$mu[,1])) , lim= as.vector(HPDI(logistic(post$mu[,1]), prob=0.9999)) , col = col.alpha("cornflowerblue", 0.25))

# dens(logistic(post$mu[,2]) , main=expression(paste(gamma))  ,  xlim=c(0,1) , ylab='', xlab= "b. weight of social information" , col="white", yaxt='n', , cex.lab=1.5) #gamma
# abline(v=logistic(gamma.sim) , col="cornflowerblue") 
# shade( density(logistic(post$mu[,2])) , lim= as.vector(HPDI(logistic(post$mu[,2]), prob=0.9999)) , col = col.alpha("cornflowerblue", 0.25))

# dens(exp(post$mu[,3]) , main=expression(paste( "\u0192"[c])),  xlim=c(0,4) , xlab="c. strength of frequency dependence" , col="white", ylab='', yaxt='n',  cex.lab=1.5)##fconf
# abline(v=exp(fc.sim) , col="red" ) #fconf
# shade( density(exp(post$mu[,3])) , lim= as.vector(HPDI(exp(post$mu[,3]), prob=0.9999)) , col = col.alpha("red", 0.25))

# dens(post$mu[,4]  ,main=expression(paste(beta)[pay]) ,  xlim=c(-4,4), xlab="d. strength of payoff bias" , ylab='',col="white", yaxt='n', cex.lab=1.5)##fpay
# abline( v=beta.p , col=col.pal[1]  ) 
# shade( density(post$mu[,4]) , lim= as.vector(HPDI((post$mu[,4]), prob=0.9999)) , col = col.alpha(col.pal[1], 0.25))

# dens(as.vector(post$lambda) , main=expression(paste(lambda)) , ylab='', xlab="i. sensitivity to individual payoff" , col="white" , yaxt='n', cex.lab=1.5)
# abline( v=k.lambda , col="black" ) 
# shade( density(post$lambda) , lim= HPDI(as.vector(post$lambda), prob=0.9999) , col = col.alpha("black", 0.25) )

# dens(post$lambda)

# gam.k <- logistic(gamma.sim + gamma.sim_i)
# gam.pred <- rep(0,n)
# for(i in 1:n){gam.pred[i] <- mean(logistic(post$mu[,2] + post$a_id[,i,2]))}
# ###plot for 2 options probability frequencyxtime(evolutionary dynamics)
# plot(gam.k,gam.pred , pch=19 , col=col.pal[1] , xlim=c(0,0.7) , ylim=c(0,0.7) )
# abline(a = 0, b = 1)


# phi.k<- logistic(phi.sim + phi.sim_i)
# phi.pred <- rep(0,n)
# for(i in 1:n){phi.pred[i] <- median(logistic(post$mu[,1] + post$a_id[,i,1]))}
# plot(phi.k,phi.pred , pch=19 , col="slateblue" , xlim=c(0,0.6) , ylim=c(0,0.6) )
# abline(a = 0, b = 1)
pdf("freq_sims.pdf",width=9,height=11)
par(mfrow=c(5,2))##set up the plot
par( mar=c(4,5,0.6,0.6) , oma=c(1,1,.1,.1) )
par(cex = 0.5)

###main plot
plot(s1/n ~ (bout-1), data=dsim, col=col.pal[1] , ylim=c(0,1.1) , pch=19 , xlab="Time" , ylab="Proportion of Individuals Choosing Option" , xlim=c(2,nbouts) )
points(s2/n ~ (bout-1), data=dsim , col=col.pal[2], pch=19)
points(s3/n ~ (bout-1), data=dsim , col=col.pal[3], pch=19)
title(main = paste("Pop. Mean: lambda=",k.lambda ,", gamma=",round(logistic(gamma.sim), digits=2),", phi=",round(logistic(phi.sim),digits=2),", f=", exp(round(fc.sim,digits=2 ))) , line = -1.1, outer = FALSE)
legend("top", inset=.8, cex=1 , as.character(techmeans), pch=19 ,col=col.pal, horiz=TRUE , bty="n")

for(i in 1:n){
    plot(Pr1 ~ (bout-1), data=dsim[dsim$id==i,] , col=col.pal[1] , ylim=c(0,1.1) , pch=19 , xlab="Time" , ylab="Proportion of Individuals Choosing Option" , xlim=c(2,nbouts) )
    points(Pr2 ~ (bout-1), data=dsim[dsim$id==i,] , col=col.pal[2], pch=19)
    points(Pr3 ~ (bout-1), data=dsim[dsim$id==i,] , col=col.pal[3], pch=19)
    title(main = paste("id=",i ,", lambda=",k.lambda ,", gamma=",gamma.sim_id[i],", phi=",phi.sim_id[i],", f=", fc.sim_id[i] ) , line = -1.1, outer = FALSE)
    legend("top", inset=.8, cex=1 , as.character(techmeans), pch=19 ,col=col.pal, horiz=TRUE , bty="n")

}

dev.off()