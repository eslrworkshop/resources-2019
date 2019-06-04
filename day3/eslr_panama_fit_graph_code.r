##################################################################
#######prepare required packages and data files###################
##################################################################

library(rethinking)
require("repmis") #to read rdata files from github
require("RCurl") #to read csv's from github

#load simplified posterior samples from GitHub
source_data("https://github.com/bjbarrett/resources-2019/blob/barretteslr2019/eslr_pn_posterior.Rdata?raw=True")

#load data from github
d  <- read.csv(text=getURL("https://raw.githubusercontent.com/bjbarrett/resources-2019/master/panama_data_14days.csv"), header=T)
mono_index_all  <- read.csv(text=getURL("https://raw.githubusercontent.com/bjbarrett/resources-2019/master/mono_indexing_all.csv"), header=T)
mono_index_subset  <- read.csv(text=getURL("https://raw.githubusercontent.com/bjbarrett/resources-2019/master/mono_indexing.csv"), header=T)


###Below is code to load locally if internet ist kaputt
# load("/Users/brendanbarrett/eslr_pn_posterior.Rdata")#change directory to local pathway
# d <- read.csv("~/Dropbox/eslr/ESLR stuff/panama_data_14days.csv" , header=TRUE) #load data file from in github
# mono_index_all <- read.csv("~/Dropbox/eslr/ESLR stuff/mono_indexing_all.csv") #monkey indexing data with HF and WZ (non-foragers)

#modify data internally for graphing purposes

#creates index variable that includes all observers
mono_index_all$mono <- mono_index_all$monos
mono_index_subset$mono <- mono_index_subset$monos

mono_index <- mono_index_all

ages <- subset(mono_index_subset, select=c(yob,mono) )
ages$age <- 2015-ages$yob 
ages$age.c <-ages$age - mean(ages$age)

#softmax function to simplify code
Softmax <- function(x){
exp(x)/sum(exp(x))
} 


##use same data list that was used when fitting stan model
d_global <- list(
    N = nrow(d),                                                                        #length of dataset
    J = length( unique(d$mono_index) ),                                                 #number of individuals
    K=max(d$tech_index),                                                                #number of processing techniques
    tech = d$tech_index,                                                                #technique index
    y = cbind( d$y1 , d$y2 , d$y3 , d$y4 , d$y5 , d$y6 , d$y7 ),                        #individual processing times for all techniques at each bout N (individual payoff)
    s = cbind(d$s1 , d$s2 , d$s3 , d$s4 ,d$s5 ,d$s6 , d$s7 ),                           #observed counts of all K techniques to individual J (frequency-dependence)
    ps = cbind(d$ps1 , d$ps2 , d$ps3 , d$ps4 ,d$ps5 ,d$ps6 , d$ps7 ),                   #observed mean payoffs of all K techniques to individual J (payoff bias)
    ks = cbind(d$ks1 , d$ks2 , d$ks3 , d$ks4 ,d$ks5 ,d$ks6 , d$ks7 ),                   #observed matrilineal kin cues of all K techniques to individual J (matrilineal kin-bias)
    cohos = cbind(d$cos1 , d$cos2 , d$cos3 , d$cos4 ,d$cos5 ,d$cos6 , d$cos7 ),         #observed cohort cues of all K techniques to individual J (age-cohort/similarity bias)
    yobs = cbind(d$Yobs1 , d$Yobs2 , d$Yobs3 , d$Yobs4 ,d$Yobs5 ,d$Yobs6 , d$Yobs7 ),   #observed age cues of all K techniques to individual J (age-bias)
    press = cbind(d$prs1 , d$prs2 , d$prs3 , d$prs4 ,d$prs5 ,d$prs6 , d$prs7 ),         #observed rank cues of all K techniques to individual J (age-bias)
    bout = d$forg_bout,#bout is forg index here                                         #processing bout unique to individual J
    id = as.integer(as.factor(d$mono_index)),                                           #individual ID
    N_effects=8,                                                                        #number of parameters to estimates
    age=d$age.c                                                                         #centered age of individuals
)

#scale payoffs/cues by dividing by max value
d1 <- d_global
d1$yobs <- d1$yobs / max(d1$yobs)
d1$cohos <- d1$cohos / max(d1$cohos)
d1$ks <- d1$ks / max(d1$ks)
d1$ps <- d1$ps / max(d1$ps)
d1$press <- d1$press / max(d1$press)

##################################################################
#######Code to fit Model which we will comment out################
##################################################################

# #parameter list to save in posterior extractions
# parlistglobalage=c("lambda" ,"a_id" , "mu" ,"Bpay","Bkin","Bpres","Bcoho","Byob", "fconf", "dev" , "log_lik" , "b_age" , "Rho" , "sigma" )

# #global model assuming file is in working directory
# fit_global_age <- stan( file = 'PN_social_global_age.stan', data = d1 , 
#     iter = 2000, warmup=1000 , chains=3, cores=3,
#     control=list( adapt_delta=0.98 ) ,pars=parlistglobalage )


##################################################################
#######plots of all main effect posteriors against priors#########
##################################################################

###########main effects posterior graphs##################
pdf("Main_effects_postgraphs.pdf",height=8,width=8)
par(mfrow=c(3,3), oma = c(0, 0, 0, 0) , mar=c(3,.5,2,.5))
par(cex = 0.6)
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

dens(logistic(post$mu[,1]) , main=expression(a.~ weight~of~new~experience~(paste(phi)) ) , xlim=c(0,1), ylab='' , xlab= '' , col="white", yaxt='n')##phi
abline(v=median(logistic(post$mu[,1])) , col="cornflowerblue" ) 
shade( density(logistic(post$mu[,1])) , lim= as.vector(HPDI(logistic(post$mu[,1]), prob=0.9999)) , col = col.alpha("cornflowerblue", 0.25))
dens( logistic(rnorm( 10000 , 0 , 1)) , lty=2  , col="cornflowerblue" , adj=1 , add=TRUE)

dens(logistic(post$mu[,2]) , main=expression(b.~weight~of~social~information~(paste(gamma)))  ,  xlim=c(0,1) , ylab='', xlab= '' , col="white", yaxt='n') #gamma
abline(v=median(logistic(post$mu[,2])) , col="cornflowerblue") 
shade( density(logistic(post$mu[,2])) , lim= as.vector(HPDI(logistic(post$mu[,2]), prob=0.9999)) , col = col.alpha("cornflowerblue", 0.25))
dens( logistic(rnorm( 10000 , 0 , 1)) , lty=2  , col="cornflowerblue" , adj=1 , add=TRUE)

dens(exp(post$mu[,3]) , main=expression(c.~strength~of~frequency~dependence~(paste( "\u0192"[c]))),  xlim=c(0,4) , xlab="c. strength of frequency dependence" , col="white", ylab='', yaxt='n')##fconf
abline(v=median(exp(post$mu[,3])) , col="red" ) #fconf
shade( density(exp(post$mu[,3])) , lim= as.vector(HPDI(exp(post$mu[,3]), prob=0.9999)) , col = col.alpha("red", 0.25))
dens( exp(rnorm( 10000 , 0 , 1)) , lty=2  , col="red" , adj=1 , add=TRUE)

dens(post$mu[,4]  ,main=expression(d.~strength~of~payoff~bias~(paste(beta)[pay])) ,  xlim=c(-4,4), xlab='' , ylab='',col="white", yaxt='n')##fpay
abline( v=median(post$mu[,4]) , col="orange"  )
shade( density(post$mu[,4]) , lim= as.vector(HPDI((post$mu[,4]), prob=0.9999)) , col = col.alpha("orange", 0.25))
curve( dnorm( x , 0 , 1) , lty=2 , add=TRUE , col="orange" )

dens(post$mu[,5]  ,main=expression(e.~strength~of~kin~bias~(paste(beta)[kin])) ,  xlim=c(-4,4), ylab='',xlab='', col="white", yaxt='n')##fpay
abline( v=median(post$mu[,5]) , col="orange" ) 
shade( density(post$mu[,5]) , lim= as.vector(HPDI((post$mu[,5]), prob=0.9999)) , col = col.alpha("orange", 0.25))
curve( dnorm( x , 0 , 1) , lty=2 , add=TRUE , col="orange" )

dens(post$mu[,6] ,main=expression(f.~strength~of~rank~bias~(paste(beta)[rank])) ,  xlim=c(-4,4), ylab='', col="white", xlab='', yaxt='n')##fpay
abline( v=median(post$mu[,6]) , col="orange" )
shade( density(post$mu[,6]) , lim= as.vector(HPDI((post$mu[,6]), prob=0.9999)) , col = col.alpha("orange", 0.25))
curve( dnorm( x , 0 , 1) , lty=2 , add=TRUE , col="orange" )

dens(post$mu[,7]  ,main=expression(g.~strength~of~age-cohort~bias~(paste(beta)[coho])) ,  xlim=c(-4,4), ylab='', xlab='' , col="white", yaxt='n')##fpay
abline( v=median(post$mu[,7]) , col="orange" ) 
shade( density(post$mu[,7]) , lim= as.vector(HPDI((post$mu[,7]), prob=0.9999)) , col = col.alpha("orange", 0.25))
curve( dnorm( x , 0 , 1) , lty=2 , add=TRUE , col="orange" )

dens(post$mu[,8]  ,main=expression(h.~strength~of~age~bias~(paste(beta)[age])) , col="white", xlim=c(-4,4), ylab='',xlab='', yaxt='n')
abline( v=median(post$mu[,8]) , col="orange" ) 
shade( density(post$mu[,8]) , lim= as.vector(HPDI((post$mu[,8]), prob=0.9999)) , col = col.alpha("orange", 0.25))
#curve( dnorm( x , 0 , 1) , lty=2 , add=TRUE , col="orange" )
curve( dnorm( x , 0 , 1) , lty=2 , add=TRUE , col="orange" )

dens(as.vector(post$lambda) , main=expression(i.~sensitivity~to~individual~payoff~(paste(lambda))) , ylab='', xlab='' , col="white" , yaxt='n', xlim=c(0,25) )
abline( v=median(post$lambda) , col="black" ) 
shade( density(post$lambda) , lim= HPDI(as.vector(post$lambda), prob=0.9999) , col = col.alpha("black", 0.25) )
curve( dexp( x , rate=1) , lty=2 , add=TRUE , col="black" )

dev.off()

##########################phi age effects graph################

vfphi <- matrix(0,nrow=length(post$lambda),ncol=23)
varefphi <- rep(0,23)
for(i in 1:23){vfphi[,i] <- logistic(post$mu[,1] + post$a_id[,i,1] + post$b_age[,1]*ages$age.c[i] ) }
for(i in 1:23){varefphi[i] <- mean(vfphi[,i])}

rando.samps <- sample(1:nrow(post$lambda), size=100, replace = FALSE, prob = NULL)
     
cairo_pdf("phi_age_varef.pdf", width=8 , height=8)
par(mar=c(5,5,0.5,0.5))
plot(varefphi~ ages$age.c , ylab="attraction toward new experience (\u0278)" , xlab="age (years)" , pch=19 , col="orange" ,xlim=c( (min(ages$age.c)-1) , (max(ages$age.c) + 2)), ylim=c(0,.6) , cex=1.5 , cex.lab=2.4  , xaxt="n" , yaxt="n" , axes=FALSE)

axis(1, at = seq(from=(min(ages$age.c)-1) , to=(max(ages$age.c) + 2), by = 5) , labels=seq(from=0 , to=25 , by=5), tck=-0.02 , cex.axis=1.5)
axis(2, at = seq(from=0 , to=0.6, by = 0.2) , tck=-0.02 , cex.axis=1.5)
axis(1, at = seq(from=(min(ages$age.c)-1) , to=(max(ages$age.c) + 2), by = 1 ), labels=F  , tck=-0.01, cex.axis=1.5)
axis(2, at = seq(from=0 , to=0.6, by = 0.1) ,tck=-0.01 , labels=F , cex.axis=1.5)

age.seq <- seq(from=min(ages$age.c) , to=max(ages$age.c) , length=25 )
pred.mean <- sapply(age.seq , function(z)
     mean(logistic(post$mu[,1] + post$b_age[,1]*z) ))

pred.ci<- sapply(age.seq , function(z)
     HPDI(logistic(post$mu[,1] + post$b_age[,1]*z) ))

pred.sims <- sapply(age.seq , function(z)
     logistic(post$mu[,1] + post$b_age[,1]*z) )

for (i in rando.samps){lines(age.seq , pred.sims[i,] , lw=2 , col=col.alpha( "orange" , alpha = 0.1 ))}
text(ages$age.c, varefphi, mono_index$mono, cex=.4, col="black")
lines(age.seq , pred.mean , lw=2)
dev.off()


##############gamma age effects graph###################

vfgamma <- matrix(0,nrow=length(post$lambda),ncol=23)
varefgamma <- rep(0,23)
for(i in 1:23){vfgamma[,i] <- logistic(post$mu[,2] + post$a_id[,i,2] + post$b_age[,2]*ages$age.c[i] ) }
for(i in 1:23){varefgamma[i] <- mean(vfgamma[,i])}
lines(age.seq , pred.mean , lw=2)

cairo_pdf("gamma_age_varef.pdf", width=8 , height=8)
par(mar=c(5,5,0.5,0.5))
plot(varefgamma~ ages$age.c , ylab="weight given to social information (\u03B3)" , xlab="age (years)" , pch=19 , col="cornflowerblue" ,xlim=c( (min(ages$age.c)-1) , (max(ages$age.c) + 2)), ylim=c(0,.6) , cex=1.7 , cex.lab=2.4  , xaxt="n" , yaxt="n" , axes=FALSE )
axis(1, at = seq(from=(min(ages$age.c)-1) , to=(max(ages$age.c) + 2), by = 5) , labels=seq(from=0 , to=25 , by=5), tck=-0.02, cex.axis=1.5)
axis(2, at = seq(from=0 , to=0.6, by = 0.2) , tck=-0.02 , cex.axis=1.5)
axis(1, at = seq(from=(min(ages$age.c)-1) , to=(max(ages$age.c) + 2), by = 1 ), labels=F  , tck=-0.01, cex.axis=1.5)
axis(2, at = seq(from=0 , to=0.6, by = 0.1) ,tck=-0.01 , labels=F , cex.axis=1.5)
age.seq <- seq(from=min(ages$age.c) , to=max(ages$age.c) , length=25 )
pred.mean <- sapply(age.seq , function(z)
     mean(logistic(post$mu[,2] + post$b_age[,2]*z) ))
pred.ci<- sapply(age.seq , function(z)
     HPDI(logistic(post$mu[,2] + post$b_age[,2]*z) ))

pred.sims <- sapply(age.seq , function(z)
     logistic(post$mu[,2] + post$b_age[,2]*z) )
for (i in rando.samps){lines(age.seq , pred.sims[i,] , lw=2 , col=col.alpha( "cornflowerblue" , alpha = 0.1 ))}
text(ages$age.c, varefgamma, mono_index$mono, cex=.4, col="black")
lines(age.seq , pred.mean , lw=2)

dev.off()


##################################################################
#######code to make individual level plots########################
##################################################################

d1$fruit_index <- d$fruit_index
d1$date_index <- d$date_index
Preds = array(0,dim=c(nrow(d),7,23)) #predictions for all individuals, all techniques, across all timesteps
Preds2 = array(0,dim=c(nrow(d),7)) ##predictions for all individuals, all techniques, at times when they foraged

lambda = mean(post$lambda)
AC=PrS=PrA=lin_mod=s_temp=rep(0,7) #stroage slots for calculating predictions

for ( i in 1:max(d1$N) ) {
    #if ( d1$bout[i]==1 ) {
        #// calculate new individual's parameter values
        phi = mean(logistic(post$mu[,1] + post$a_id[,d1$id[i],1] + post$b_age[,1]*ages$age.c[d1$id[i]] ))
        fconf = mean(exp( post$mu[,3] + post$a_id[,d1$id[i],3] ))
        fpay = mean( post$mu[,4] + post$a_id[,d1$id[i],4] )
        gamma = mean(logistic( post$mu[,2] + post$a_id[,d1$id[i],2] + post$b_age[,2]*ages$age.c[d1$id[i]]))
        ##PICK UP HERE
        fkin = mean(( post$mu[,5] + post$a_id[,d1$id[i],5] ))
        fpres = mean(( post$mu[,6] + post$a_id[,d1$id[i],6] ))
        fcoho = mean(( post$mu[,7] + post$a_id[,d1$id[i],7] ))
        fyob = mean(( post$mu[,8] + post$a_id[,d1$id[i],8] ))
        #}
        #//update attractions
        for ( j in 1:max(d1$tech) ) {
            if ( d1$bout[i] > 1 ) {
                AC[j] = (1-phi)*AC[j] + phi*d1$y[i-1,j]
            } else {
                AC[j] = 0;
            }
        }#//j

        for (j in 1:max(d1$tech)){PrA[j] = exp(lambda*AC[j])/sum(exp(lambda*AC))}
        
        #//conformity aspect below
        if ( d1$bout[i] > 1 ) {
            if (sum( d1$s[i,] ) > 0 ) {

                #// compute non-frequency cues as log-linear model
                for ( j in 2:max(d1$tech) ) {
                    lin_mod[j] = exp( fpay*d1$ps[i,j] + fkin*d1$ks[i,j] + fpres*d1$press[i,j] + fcoho*d1$cohos[i,j] + fyob*d1$yobs[i,j])
                }
                lin_mod[1] = 1; #// aliased outcome

                #// compute frequency cue
                for ( j in 1:max(d1$tech) ){ s_temp[j] = d1$s[i,j]^fconf }
                for ( j in 1:max(d1$tech) ){ lin_mod[j] = lin_mod[j] * s_temp[j] }

                for (j in 1:7){PrS[j] = lin_mod[j]/sum(lin_mod)}
                
                for(j in 1:7){ Preds[d1$fruit_index[i],j,d1$id[i]] = (1-gamma)*PrA[j] + gamma*PrS[j] 
                               Preds2[i,j] = (1-gamma)*PrA[j] + gamma*PrS[j]
                }

            } else {
                for(j in 1:7){ Preds[d1$fruit_index[i],j,d1$id[i]]= PrA[j] 
                               Preds2[i,j] = PrA[j]}
            }
        } else {
            for(j in 1:7){ Preds[d1$fruit_index[i],j,d1$id[i]]= PrA[j]
                           Preds2[i,j] = PrA[j] }
         }
     }#//i  



#col_index= c("darkgreen","blue","red","gold","grey","violet","orange")
col_index= c("#1B9E77", "#D95F02" ,"#7570B3", "#E7298A" ,"#66A61E", "#E6AB02" ,"#A6761D")
d$date_index <- as.integer(as.factor(d$date_index))
matF <- matrix(, nrow = 75, ncol = 8)
matP <- matrix(, nrow = 75, ncol = 8)
mat <- matrix(, nrow = 75, ncol = 8)

for (i in 1:75){
    for( j in 1:7) {
        matF[i,j] <- length(unique(d$fruit_index[d$date_index==i & d$tech_index==j]))/length(unique(d$fruit_index[d$date_index==i]))
        matP[i,1] <- mean(d$ps1[d$date_index==i])
        matP[i,2] <- mean(d$ps2[d$date_index==i])
        matP[i,3] <- mean(d$ps3[d$date_index==i])
        matP[i,4] <- mean(d$ps4[d$date_index==i])
        matP[i,5] <- mean(d$ps5[d$date_index==i])
        matP[i,6] <- mean(d$ps6[d$date_index==i])
        matP[i,7] <- mean(d$ps7[d$date_index==i])
    }
}

tech_id <- as.vector(sort(unique(d$tech)))
mat[,8] <- c(1:75)

    for(k in 1:23){
        for (i in 1:1440){
        for (j in 1:7){
            Preds[i+1,j,k] <- ifelse(Preds[i+1,j,k]==0 & Preds[i,j,k]!=0,Preds[i,j,k],Preds[i+1,j,k] )
        }
    }
}




##plots where probability is constant between fruits



##below plot was not in paper but plots time series of behavioral change against all foraging behaviors, useful for analyst

pdf("individual_panama_predictions_fruit_pred_all.pdf",width=11,height=8.5) 

par(mfrow=c(4, 1) , oma = c(5,3,1,0) , mar=c(3,3,1.5,0))
par(cex = 0.6)
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

for(k in 1:23){
    plot(d1$date_index~d1$fruit_index, ylim=c(0,1.2) , col="white" , xlab="" , ylab="" , cex.lab=.8 ,main=mono_index$mono[k], axes=FALSE, cex.lab=1.1)
           # legend("top", inset=0.01, c(tech_id) , fill=col_index, horiz=TRUE,cex=0.6,bty = "n")
    axis(1, at = seq(from=0 , to=1450, by = 100) , tck=-0.03)
    axis(2, at = seq(from=0 , to=1 , by = 0.2) , tck=-0.03 )
    axis(1, at = seq(from=0 , to=1450, by = 25) , tck=-0.015 , labels=FALSE )
    axis(2, at = seq(from=0 , to=1 , by = 0.1) , tck=-0.015 , labels=FALSE)
    for (j in 1:7){    
        for (i in 1:1441){
           points(Preds[i,j,k]~i, col=col.alpha(col_index[j], alpha=0.99) , pch=20 , cex=0.5 )
        }
        lls <- d$timestep[d$fruit_index<1500 & d$tech_index==j & d$open==1 & d$mono_index==k]
        llf <- d$timestep[d$fruit_index<1500 & d$tech_index==j & d$open==0 & d$mono_index==k]
        llo <-  d$timestep[grepl(mono_index$monos[k],d[,"RCI"])==TRUE & d$tech_index==j]

            #points( mono_index$yob_order[j] ~ d$timestep[i]  , pch=15 , cex=2 , col=col.alpha("black" , .2) )
        points(llo,rep(1.15,length(llo))  , pch="*" , cex=1.2 , col=col_index[j] )
        points(lls,rep(1.05,length(lls)), pch=23 , cex=0.75 , col=col_index[j] , bg=col.alpha(col_index[j] , 0.5) )
        points(llf,rep(1.1,length(llf)), pch=4 , cex=0.75 , col=col_index[j] , bg=col.alpha(col_index[j] , 0.5) )    
    
    }
}

mtext("fruit processing event time", side = 1, outer = TRUE, cex = 1.5, line = .1)
mtext("probability of choosing technique", side = 2, outer = TRUE, cex = 1.5, line = .1)

legend(x=200 , y=-.5, inset= -.1, tech_id , fill=col_index, border=col_index ,horiz=TRUE,cex=1, bty = "n", xpd=NA)
legend(x=500 , y=-.65, inset= -.1, c(" tech succesful", "tech failure", "observed tech") , col=1, horiz=TRUE, cex=1, pch=c(23,4,8), xpd=NA , bty="n" , bg=col.alpha("black" , 0.5) )

dev.off()


PredsPop = array(0,dim=c(75,7))
PredsPop2 = array(0,dim=c(75,7))
PredsMono = array(0,dim=c(75,7,23))

for(i in 1:75){
    for(k in 1:23){
        for (j in 1:7){
            PredsMono[i,j,k] <- mean(Preds[d1$fruit_index[d1$date_index==i & d1$tech==j],j,k])
        }
    }
}

for(i in 1:75){
        for (j in 1:7){
            PredsPop[i,j] <- mean(Preds[d1$fruit_index[d1$date_index==i & d1$tech==j],j,])
        }
    }


pdf("individual_panama_predictions2.pdf",width=8.5,height=11) 
par( mfrow=c(12, 2) , mar=c(1,1,1,1) , oma=c(4,4,.5,.5) )
par(cex = 0.5)
par(tcl = -0.2)
par(mgp = c(2, 0.6, 0))

for(k in 1:23){
        plot(mat[,1]~mat[,8], ylim=c(0,1) , col="white" , xlab="" , ylab="" ,main=""  , cex.axes=0.4)
        #legend("top", inset=0.01, c(tech_id) , fill=col_index,border=col_index, horiz=TRUE,cex=0.75,bty = "n")
#
    for(i in 1:75){
        for (j in 1:7){
            points(PredsMono[i,j,k]~i ,col=col.alpha(col_index[j], alpha=0.99) , pch=18 )
            lines (PredsMono[i,j,k]~i ,col=col.alpha(col_index[j], alpha=0.99) , pch=18 )
 
        }
mtext(mono_index$mono[k], side = 3, line = -2, adj = 0.01, cex = .8) 


    }
}
 plot(mat[,1]~mat[,8], ylim=c(0,1) , col="white" , xlab="" , ylab="" ,main="" , xaxt="n" , yaxt="n" , axes=FALSE)
 legend("top", inset=0.1, c(tech_id[1:3]) , fill=col_index[1:3],border=col_index[1:3], horiz=TRUE,cex=1.2,bty = "n")
 legend("bottom", inset=0.1, c(tech_id[4:7]) , fill=col_index[4:7],border=col_index[4:7], horiz=TRUE,cex=1.2,bty = "n")
mtext("Experimental Days (N=75)", side = 1, line = 1.2, cex = 1, outer=TRUE) 
mtext("daily average probability of choosing technique", side = 2, line = 1.2, cex = 1.2 , outer=TRUE) 

#axis(1, at = seq(from=min(mat[,8]) , to=max(mat[,8]) , by = 5 ), labels=F  , tck=-0.01)
#axis(2, at = seq(from=0 , to=1, by = 0.25) ,tck=-0.01 , labels=TRUE )
dev.off()


plot(mat[,1]~mat[,8], ylim=c(0,1.1) , col="white" , xlab="Experimental Days (N=75)" , ylab="probability of choosing technique" , cex.lab=1.5 )
    for(i in 1:75){
        for (j in 1:7){
            points(PredsPop[i,j]~i ,col=col.alpha(col_index[j], alpha=0.99) , pch=18 )
            lines( smooth.spline(date_rep, y=PredsPop[,j] , spar=.9) , col=col_index[j] , lw=1)
        }
    }


###old
for(i in 1:75){
        for (j in 1:7){
            PredsPop2[i,j] <- mean(Preds2[d1$fruit_index[d1$date_index==i & d1$tech==j],j])
        }
    }
###trial
d1$row.seq <- seq(1:d1$N)
d$row.seq <- seq(1:nrow(d))

for(i in 1:75){
        for (j in 1:7){
            PredsPop2[i,j] <- mean(Preds2[d$row.seq[d$date_index==i & d$tech_index==j],j])
        }
    }

d[d$date_index==1 & d$tech_index==1,]


PredsPop2[is.nan(PredsPop2)] = 0
date_rep <- seq(1:75)


#raw data
mat[,8] <- c(1:75)

matlo <- mat
for (i in 1:7){ 
    matlo[,i] <- ifelse(matlo[,i]==0, 0.001, matlo[,i])
    matlo[,i] <- ifelse(matlo[,i]==1, 0.999, matlo[,i])
}


#####code for figure s3, individual predictions by day

cairo_pdf("individual_panama_predictions_with_pop.pdf",width=8.5,height=11)
m <- matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,25,13,14,15,16,17,18,19,20,21,22,23,24,25), nrow = 13, ncol = 2)
##set up the plot
layout(m)
par( mar=c(1,1,0.6,0.6) , oma=c(2,4,.1,.1) )
par(cex = 0.5)
par(tcl = -0.2)
par(mgp = c(2, 0.6, 0))
        plot(mat[,1]~mat[,8], ylim=c(0,1) , col="white" , xlab="" , ylab="" ,main=""  , cex.axes=0.4)
    for(i in 1:75){
        for (j in 1:7){
            points(PredsPop2[i,j]~i ,col=col.alpha(col_index[j], alpha=0.99) , pch=18 ,  )
            #ss <- smooth.spline(x=date_rep, y=logit(PredsPop2[,j]) , spar=0.94  ) 
            #lines( ss$x, logistic(ss$y), col=col_index[j] , lw=4)
        }
    }

mtext("Population Mean", side = 3, line = -2, adj = 0.01, cex = .8) 

for(k in 1:23){
        plot(mat[,1]~mat[,8], ylim=c(0,1) , col="white" , xlab="" , ylab="" ,main=""  , cex.axes=0.4)
        #legend("top", inset=0.01, c(tech_id) , fill=col_index,border=col_index, horiz=TRUE,cex=0.75,bty = "n")
#
    for(i in 1:75){
        for (j in 1:7){
            points(PredsMono[i,j,k]~i ,col=col.alpha(col_index[j], alpha=0.99) , pch=18 )
            lines (PredsMono[i,j,k]~i ,col=col.alpha(col_index[j], alpha=0.99) , pch=18 )
 
        }
mtext(mono_index$mono[k], side = 3, line = -2, adj = 0.01, cex = .8) 
mtext(mono_index$yob[k], side = 3, line = -3, adj = 0.01, cex = .5) 


    }
}

 plot(mat[,1]~mat[,8], ylim=c(0,0.8) , col="white" , xlab="" , ylab="" ,main="" , xaxt="n" , yaxt="n" , axes=FALSE)
 legend("top", inset=0.1, c(tech_id[1:3]) , fill=col_index[1:3],border=col_index[1:3], horiz=TRUE,cex=1.3,bty = "n")
 legend("bottom", inset=0.1, c(tech_id[4:7]) , fill=col_index[4:7],border=col_index[4:7], horiz=TRUE,cex=1.3,bty = "n")
mtext("Experimental Days (N=75)", side = 1, line = 0.01, cex = 1.2, outer=TRUE) 
mtext("daily average probability of choosing technique", side = 2, line = 1.2, cex = 1.2 , outer=TRUE) 

#axis(1, at = seq(from=min(mat[,8]) , to=max(mat[,8]) , by = 5 ), labels=F  , tck=-0.01)
#axis(2, at = seq(from=0 , to=1, by = 0.25) ,tck=-0.01 , labels=TRUE )
dev.off()

