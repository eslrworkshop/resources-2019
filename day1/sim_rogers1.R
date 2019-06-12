# basic evolution social learning
# Rogers' model
# Q: what is the benefit of social learning to a species?

# "Roughgarden method" - after Joan Roughgarden
# three kinds of assumptions
# (1) what is the population structure?
# (2) individual states
# (3) life cycle

# (1) population structure
# finite population, well-mixed (everyone can interact with everyone else)
# population all experiences same environment
# environment changes periodically

# (2) individual states
# heritable states: learning strategies - 
# (a) individual learning: pay a personal cost to figure out what to do - c is cost, s is chance you figure out the right thing
# (b) social learning: pay no cost! copy a random adult and maybe get a good behavior
# b is benefit of doing the right thing

# (3) life cycle
# (a) birth
# (b) environmental change
# (c) babies learn
# (d) all adults die

sim_rogers <- function( 
		num_generations = 1e3,
		population_size = 200,
		c = 1,                 # cost of individual learning
		s = 1,                 # chance learn right thing
		b = 2,                 # benefit of right thing
		w0 = 10,               # baseline fitness
		u = 0.1,               # prob environment changes
		mu = 0.001			   # mutation rate
 ) 
{

	# initialize the population
	population_strategies <- sample( 1:2 , size=population_size , replace=TRUE )
	population_behavior <- rep( 0 , population_size )

	# initialize history
	history <- matrix( NA , nrow=num_generations , ncol=4 )

	# loop over generations
	for ( t in 1:num_generations ) {

		# record history
		# freq of individual learners
		history[t,1] <- sum(population_strategies==1)/population_size
		# freq of social learners
		history[t,2] <- sum(population_strategies==2)/population_size
		# freq of right behavior
		history[t,3] <- sum(population_behavior==1)/population_size

		# reproduction / birth
		# compute fitness of each adult
		fitness <- w0 + (population_behavior==1)*b - (population_strategies==1)*c
		# record fitness
		history[t,4] <- mean(fitness)
		# produce babies in proportion to fitness
		babies_strategies <- sample( population_strategies , size=population_size , replace=TRUE , prob=fitness )

		# mutation
		for ( i in 1:population_size ) {
			if ( runif(1) < mu )
				babies_strategies[i] <- 3 - babies_strategies[i]
		}#i

		# environment changes?
		if ( runif(1) < u ) population_behavior <- rep( 0 , population_size )

		# learning
		baby_behavior <- rep( 0 , population_size )
		for ( i in 1:population_size ) {
			if ( babies_strategies[i]==1 )
				if ( runif(1) < s ) baby_behavior[i] <- 1
			if ( babies_strategies[i]==2 )
				baby_behavior[i] <- sample( population_behavior , size=1 )
		}#i

		# death and growth
		population_strategies <- babies_strategies
		population_behavior <- baby_behavior


	}#t

	return(history)

}# end of sim_rogers

h <- sim_rogers( u=0.2 , b=3 , c=1 , s=1 )

plot( h[,1] , ylim=c(0,1)  , xlab="generation" , ylab="frequency", type="l" )
#lines( 1:nrow(h) , h[,2] , col="purple" )
lines( 1:nrow(h) , h[,3] , col="indianred" )

plot( h[,4] , type="l" )

x <- 500:700
plot( h[x,1] , ylim=c(0,1)  , xlab="generation" , ylab="frequency", type="l" )
#lines( 1:length(x) , h[x,2] , col="purple" )
lines( 1:length(x) , h[x,3] , col="indianred" )

