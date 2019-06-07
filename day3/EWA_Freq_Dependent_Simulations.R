#######################################################
######set color pallete and load functions#############
#######################################################

col.pal=c("#1B9E77", "#D95F02", "#7570B3") #graphing color pallette

#########create a softmax function to simply code, add in logit and logistic fcns
Softmax <- function(x){
    exp(x)/sum(exp(x))
} 

logit <- function(p){ 
    (log(p/(1-p)))
}

logistic <- function(x){ 
    (1/(1+exp(-x)))
}

#######add convenient density plooting function, code lifted from rethinking by McElreath
dens <- function (x, adj = 0.5, norm.comp = FALSE, main = "", show.HPDI = FALSE, 
    show.zero = FALSE, rm.na = TRUE, add = FALSE, ...) 
{
    if (inherits(x, "data.frame")) {
        n <- ncol(x)
        cnames <- colnames(x)
        set_nice_margins()
        par(mfrow = make.grid(n))
        for (i in 1:n) {
            dens(x[, i], adj = adj, norm.comp = norm.comp, show.HPDI = show.HPDI, 
                show.zero = TRUE, xlab = cnames[i], ...)
        }
    }
    else {
        if (rm.na == TRUE) 
            x <- x[!is.na(x)]
        thed <- density(x, adjust = adj)
        if (add == FALSE) {
            set_nice_margins()
            plot(thed, main = main, ...)
        }
        else lines(thed$x, thed$y, ...)
        if (show.HPDI != FALSE) {
            hpd <- HPDI(x, prob = show.HPDI)
            shade(thed, hpd)
        }
        if (norm.comp == TRUE) {
            mu <- mean(x)
            sigma <- sd(x)
            curve(dnorm(x, mu, sigma), col = "white", lwd = 2, 
                add = TRUE)
            curve(dnorm(x, mu, sigma), add = TRUE)
        }
        if (show.zero == TRUE) {
            lines(c(0, 0), c(0, max(thed$y) * 2), lty = 2)
        }
    }
}

set_nice_margins <- function () {
    par_mf <- par("mfrow", "mfcol")
    if (all(unlist(par_mf) == 1)) {
        par(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1, 
            tck = -0.02)
    }
}

#######################################################
######begin data simulations###########################
#######################################################

#data sims
n <- 50                                  #number of individuals/pop size
nbouts <- 75                             #timesteps

#simulate values for options
techmeans <- c(8 , 8 , 8)               #mean efficiency of techniques
techvar <- c(0, 0 , 0)                  #variance of techniques

#parameter sims
fc.sim <- log(2)				          #frequency dependecy parameter on log scale
phi.sim <- logit(0.18)                    #memory/attraction updating parameter on log-odds scale
gamma.sim <- 1                            #weight given to social info parameter on log-odds scale
k.lambda <- 1                             #sensitivity to individual payoffs

#varying effects offsets for individuals
gamma.sim_i <- rnorm( n , mean=0 , sd=0) ##frequency dependent offsets per individual
phi.sim_i <- rnorm( n , mean=0 , sd=0)   
fc.sim_i <- rnorm( n , mean=0 , sd=0) 

#plot to visualize overlap of payoffs
dens(rnorm( 10000 , mean=techmeans[1] , sd=techvar[1] ) ,col=col.pal[1] , xlim=c(0,20) )
dens(rnorm( 10000 , mean=techmeans[2] , sd=techvar[2] ) ,col=col.pal[2] , xlim=c(0,20) , add=TRUE )
dens(rnorm( 10000 , mean=techmeans[3] , sd=techvar[3] ) ,col=col.pal[3] , xlim=c(0,20) , add=TRUE )

#unique parameters for each individual, to visualize heterogeneity and for plotting individual predictions
gamma.sim_id <- round( logistic(gamma.sim + gamma.sim_i), digits=2) 		##simulated gammas for all n individuals
phi.sim_id <- round(logistic(phi.sim +phi.sim_i), digits=2) 				##simulated phis for all n individuals
fc.sim_id <- round(exp(fc.sim + fc.sim_i), digits=2)  						##simulated strength of frequency dependent learning for all n individuals
gamma.sim_id
phi.sim_id
fc.sim_id

#begin to simulate data
dsim_s <- data.frame( id=0 , bout=0 , tech=0 , y1=0 , y2=0, y3=0 , s1=0 , s2=0 , s3=0 , A1=0 , A2=0 , A3=0 , Pr1=0 , Pr2=0 , Pr3=0 )
therow <- 1

AC <- matrix(1,ncol=3,nrow=n) 												#attraction scores for each behavior
AC[,3] <- 0.2 																#change to make high payoff option less likely
Softmax(AC[1,]) 															#run to see initial prob of choosing a behavior

S1 <- S2 <- S3 <- rep(0,n+1) 												# num of individuals choosing each tech in previous bout
PS1 <- PS2 <- PS3 <- rep(0,nbouts+1) 										# empty vector for mean observed in previous rounds
s_temp <-  rep(0,3)

# S1[1] <- 0.48*nbouts
# S2[1] <- 0.32*nbouts
# S3[1] <- 0.2*nbouts

for ( r in 1:nbouts ) {
    for ( i in 1:n ) {  

		prtech_i <-  Softmax(k.lambda*AC[i,]) 
        my.gam <- logistic( gamma.sim + gamma.sim_i[i] ) 	#social info weight for individual i
        my.phi <- logistic( phi.sim + phi.sim_i[i] )		#weight given to past experience for individual i
        my.fconf <- exp( fc.sim + fc.sim_i[i])  			#strength of frequency dependence for individual i
        prtech_su <- c(S1[r],S2[r],S3[r]) 					#social info individual i observed

        #//frequency dependent social learning aspect below
        if ( r >= 1 ) { 
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
        yield <- rnorm( 1 , mean=techmeans[tech] , sd=techvar[tech] )


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
plot(s1/n ~ bout, data=dsim[dsim$bout>1,], col=col.pal[1] , ylim=c(0,1.1) , xlim=c(2,nbouts+1), pch=19 , xlab="Time" , ylab="Proportion of Individuals Choosing Option" , main="Population Mean ")
points(s2/n ~ bout, data=dsim[dsim$bout>1,] , col=col.pal[2], pch=19)
points(s3/n ~ bout, data=dsim[dsim$bout>1,] , col=col.pal[3], pch=19)
legend("topleft", cex=1 , as.character(techmeans), pch=19 ,col=col.pal, horiz=TRUE , bty="n", title="Payoffs")

pdf("freq_sims.pdf",width=9,height=11)
par(mfrow=c(5,2))##set up the plot
par( mar=c(4,5,0.6,0.6) , oma=c(1,1,.1,.1) )
par(cex = 0.5)

###main plot
plot(s1/n ~ bout, data=dsim[dsim$bout>1,], col=col.pal[1] , ylim=c(0,1.2) , pch=19 , xlab="Time" , ylab="Proportion of Individuals Choosing Option" , xlim=c(2,nbouts) )
points(s2/n ~ bout, data=dsim[dsim$bout>1,] , col=col.pal[2], pch=19)
points(s3/n ~ bout, data=dsim[dsim$bout>1,] , col=col.pal[3], pch=19)
title(main = paste("Pop. Mean: lambda=",k.lambda ,", gamma=",round(logistic(gamma.sim), digits=2),", phi=",round(logistic(phi.sim),digits=2),", f=", exp(round(fc.sim,digits=2 ))) , line = -1.2, outer = FALSE)
legend("top", inset=.05, cex=1 , as.character(techmeans), pch=19 ,col=col.pal, horiz=TRUE , bty="n")

for(i in 1:n){
    plot(Pr1 ~ (bout-1), data=dsim[dsim$id==i & dsim$bout>1,] , col=col.pal[1] , ylim=c(0,1.2) , pch=19 , xlab="Time" , ylab="Proportion of Individuals Choosing Option" , xlim=c(2,nbouts) )
    points(Pr2 ~ (bout-1), data=dsim[dsim$id==i & dsim$bout>1,] , col=col.pal[2], pch=19)
    points(Pr3 ~ (bout-1), data=dsim[dsim$id==i & dsim$bout>1,] , col=col.pal[3], pch=19)
    title(main = paste("id=",i ,", lambda=",k.lambda ,", gamma=",gamma.sim_id[i],", phi=",phi.sim_id[i],", f=", fc.sim_id[i] ) , line = -1.2, outer = FALSE)
    legend("top", inset=.05, cex=1 , as.character(techmeans), pch=19 ,col=col.pal, horiz=TRUE , bty="n")

}

dev.off()

