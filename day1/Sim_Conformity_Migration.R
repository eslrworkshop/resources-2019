# simulation model of evolution of social learning
# includes both spatial and temporal environmental variation
# includes individual, social, and conformist learning

# INSTRUCTIONS
# Paste the contents of this file into your R terminal.
# Then the simulation can be run with:
#       result <- sim_rogers2()
# Add arguments to change the simulation parameters. For example:
#   result <- sim_rogers2( migration_rate=0 )

# The result contains an array with the history of the population.
# The structure is:
#   result[ generation , group , trait ]
# where 'trait' is one of:
# 1: frequency of individual learners
# 2: frequency of social learners
# 3: frequency of conformist learners
# 4: frequency of adaptive behavior

# the function sim_plot (at bottom of this script) will plot these frequencies, averaging over groups.

sim_rogers2 <- function( 
	number_of_generations=1e3 , # length of simulation
	n_groups=2 ,                # number of sub-populations
	pop_size=200 ,              # size of each sub-population
	rate_of_env_change=0.01 ,   # prob environment changes each generation
	migration_rate=0.01 ,       # rate of migration among groups
	cost_of_indiv_learning=1 ,  # fitness cost of individual learning
	prob_learn_success=1 ,      # prob individual learning produces adaptive behavior
	benefit_of_adapted=2 ,      # fitness benefit of adaptive behavior
	baseline_fitness=10 ,       # baseline fitness
	mutation_rate=0.001 )       # rate of mutation among learning strategies
{
	# initialize population
	# 1 means individual learner
	# 2 means social learner
	# 3 means conformist - sample 3 adults, take majority
	pop <- array( NA , c( n_groups , pop_size , 2 ) )
	for ( g in 1:n_groups ) {
		pop[g,,1] <- sample( 1:3 , size=pop_size , replace=TRUE )
		pop[g,,2] <- rep( 0 , pop_size ) # none adapted at start
	}

	# initial envrionment states
	env <- sample( 1:100 , size=n_groups )

	# initialize record keeping
	# 1-3: proportions of strategies
	# 4: proportion adapted
	history <- array( NA , c( number_of_generations , n_groups , 4 ) )

	# loop over generations
	for ( t in 1:number_of_generations ) {

		# record history
		for ( g in 1:n_groups ) {
			for ( strat in 1:3 )
				history[t,g,strat] <- length(which( pop[g,,1]==strat )) / pop_size
			history[t,g,4] <- mean( pop[g,,2] )
		}

		# babies
		# calculate fitness of each adult
		babies <- array( NA , c( n_groups , pop_size , 2 ) )
		for ( g in 1:n_groups ) {
			fitness <- pop[g,,2]*benefit_of_adapted + baseline_fitness - cost_of_indiv_learning*(pop[g,,1]==1)
			# vector of babies
			babies[g,,1] <- sample( pop[g,,1] , size=pop_size , prob=fitness , replace=TRUE )
			# mutation
			for ( i in 1:pop_size ) {
				if ( runif(1) < mutation_rate )
					babies[g,i,1] <- sample( 1:3 , size=1 )
			}
			# init non-adapted
			babies[g,,2] <- 0
		}

		# environmental change
		for ( g in 1:n_groups ) {
			if ( runif(1) < rate_of_env_change ) {
				pop[g,,2] <- rep( 0 , pop_size )
			}
		}#g

		# learning
		# each social learning baby samples a random adult
		# each indiv learning baby learns for itself
		for ( g in 1:n_groups ) {
			for ( i in 1:pop_size ) {
				if (babies[g,i,1]==2) { # 
					babies[g,i,2] <- sample( pop[g,,2] ,size=1)
				}
				if (babies[g,i,1]==1) { # individual learning
					babies[g,i,2] <- rbinom( 1 , size=1 , prob=prob_learn_success )
				}
				if (babies[g,i,1]==3) { # conformity
					sample_of_3_adults <- sample( pop[g,,2] ,size=3)
					babies[g,i,2] <- ifelse(sum(sample_of_3_adults)>1,1,0)
				}
			}#i
		}#g

		# replace adults with babies
		pop <- babies

		# migration
		for ( g in 1:n_groups ) {
			for ( i in 1:pop_size ) {
				if ( runif(1) < migration_rate ) {
					# get another group
					destination <- g + sample( c(-1,1) , size=1 )
					if ( destination < 1 ) destination <- n_groups
					if ( destination > n_groups ) destination <- 1
					# swap with same i individual in destination group
					swap <- pop[destination,i,]
					pop[destination,i,] <- pop[g,i,]
					pop[g,i,] <- swap
					pop[g,i,2] <- 0
					pop[destination,i,2] <- 0
				}
			}#i
		}#g

	}#t

	return(
		history
	)
}

sim_plot <- function( x ) {
  x <- apply( x , c(1,3) , mean )
  n <- dim(x)[1]
  par(mar=c(5, 4, 2, 8.5), xpd=TRUE)
  plot( x[,1] , type="l" , ylim=c(0,1) , ylab="Frequency" , xlab="generation" )
  lines( 1:n , x[,2] , col="orange" )
  lines( 1:n , x[,3] , col="green" )
  lines( 1:n , x[,4] , col="blue" , lty=2 )
  legend("right", inset=c(-0.22,0), legend=c("Individual Learn.","Social Learn.", "Conformist Learn.", "Adapted agents"), 
         col=c("black", "orange", "green", "blue"),lty=c(1,1,1,2), title = "Proportion", cex = 0.8)
}
