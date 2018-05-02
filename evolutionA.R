
# where is the proteo-data?
file_path = "proteomics_dataset.csv"


### Load libraries
  library(class)
  library(ggplot2)
  library(reshape2)

# Sub functions
  
  generate_initial_population = function(length_of_chromosome, 
                                         size_of_population)
    { # creates initial random population
    init.mat = matrix(data=NA,ncol=length_of_chromosome,nrow= size_of_population)
    init.mat = apply(init.mat,2,function(x)round(runif(n=size_of_population)))
    return(init.mat)
  }
  
  calculate_fitness = function(individuals = temp.mat)
    { # one-max fitness function, the sum of the string
    fitness =  apply(individuals,1,sum)     
    return(fitness)
  }
  
  mutate = function(individuals = temp.mat, 
                    mutation_rate = 0.01)
    { # changes a certain % of cells from 0->1 or 1->0
    
    mut_loci = sample(x = 1:length(individuals),
                      size = ceiling(length(individuals)*mutation_rate),replace=F)
  individuals[mut_loci] = abs(individuals[mut_loci]-1)
  return(individuals)
  }
  
  crossover = function(parent1, parent2,length_of_chromosome)
    {  # combines chromosomes of two parents and randomly and creates a child chromosome
    crossover_point = round(sample(0:length_of_chromosome,1))
    if(crossover_point==0){
      child = parent2
    } else {
      if(crossover_point==length_of_chromosome){
        child=parent1
      } else {
        child = c(parent1[1:crossover_point], parent2[(crossover_point+1):length_of_chromosome])
      }
    }
    return(child)
  }
  
  select_parent = function(scores)
    { # randomly selects parents, whereby prob. is dependend on fitness (return the index of the parent!)
    parent = sample(1:length(scores),1,prob=scores)
    return(parent)
  }
  

############################################################
# evolution MAXONE  ########################################
############################################################


simpleGA = function(string.length = 10, # the length of the chromosome to be evolved
                    number_of_generations = 100, # # generations the evolution will last
                    size_of_population = 20,
                    data = df,
                    mutation_rate = 0.01){
      string.length = string.length
        population = generate_initial_population(length_of_chromosome = string.length, 
                                                 size_of_population = size_of_population)
        
        new_pop = matrix(data=NA,ncol=string.length,nrow= size_of_population)
        fitness = matrix(nrow=number_of_generations,ncol= size_of_population)
        
        for(gen in 1:number_of_generations){
          fitness[gen,] = calculate_fitness(population)
          population = mutate(individuals = population, mutation_rate = mutation_rate)
          for(j in 1:size_of_population){
            p1 = select_parent(fitness[gen,])
            p2 = select_parent(fitness[gen,])
            new_pop[j,] = crossover(population[p1,],population[p2,],string.length)
          }
          population = new_pop 
        }
        
        results = max(fitness[number_of_generations,])
        mean.evolution = rowMeans(fitness)
        max.evolution = apply(fitness,1,max)
        min.evolution = apply(fitness,1,min)
        plot(mean.evolution,type="l",ylim=c(0,string.length))
        lines(max.evolution,col="gray")
        lines(min.evolution,col="gray")
        return(results)
        }



simpleGA(string.length = 10, # the length of the chromosome to be evolved
         number_of_generations = 100, # # generations the evolution will last
         size_of_population = 20,
         data = df,
         mutation_rate = 0.01)



############################################################
# Feature selection for knn evolution ######################
############################################################

# load the proteo data set
  proteo = read.csv(file=file_path)

# adjust the fittness function
  fit.proteo = function(individuals,size,dat,y,kk)
  { # adjustment of the fitness function for proteo data : 
    # fitness = absolute number of correctly classfieid Y in test data above chance!
    # % or normal aboslute value dont give good weights for parent selection?!
    indiv.fit = matrix(ncol= nrow(individuals),nrow=1)
    correct = matrix(ncol= nrow(individuals),nrow=1)
    for(p in 1:size){
      train.set = sample(1:nrow(dat),nrow(dat)*2/3)
       temp = sum(class::knn(train = dat[train.set,as.logical(individuals[p,])],
                                     test = dat[-train.set,as.logical(individuals[p,])],
                                     cl = y[train.set],
                                     k = kk) == y[-train.set]) # use the total number of correctly classfieid
       # for plotting later
       correct[,p] = temp/length(y[-train.set])
       
       # for fittness
       temp = temp - max(table(y[-train.set]))
       temp = ifelse(temp <=0,0,temp)
      indiv.fit[,p] = temp
    }
    
    fitness =  list(indiv.fit=indiv.fit,
                    correct=correct)
    return(fitness)
  }

#####  takes a few minutes to run !
# adjust the main function
  GA.proteo = function(length_of_chromosome = ncol(proteo)-1, # the length of the chromosome to be evolved
                      number_of_generations = 50, # # generations the evolution will last
                      size_of_population = 3,   # AT LEAST 3!
                      data = proteo[,-length(proteo)],
                      Y = proteo[,length(proteo)],
                      mutation_rate = 0.01,
                      k = 5){
    start.time = Sys.time()
    population = generate_initial_population(length_of_chromosome = length_of_chromosome, 
                                             size_of_population = size_of_population)
    
    new_pop = matrix(data=NA,ncol=length_of_chromosome,nrow= size_of_population)
    fitness = matrix(nrow=number_of_generations,ncol= size_of_population)
    classification.history = matrix(nrow=number_of_generations,ncol= size_of_population)
    
    
    for(gen in 1:number_of_generations){
      cat("\n Generation ",gen," of ",number_of_generations,": -- ",round(gen/number_of_generations,4)*100,"%",sep="")
      # assess fitness
      temp.list = fit.proteo(individuals=population,size=size_of_population,dat=data,y=Y,kk=k )
      classification.history[gen,] = temp.list$correct # correctly classfiied
      fitness[gen,] = temp.list$indiv.fit
      
      population = mutate(individuals = population, mutation_rate = mutation_rate)
      
      for(j in 1:size_of_population){
        p1 = select_parent(fitness[gen,])
        p2 = select_parent(fitness[gen,])
        new_pop[j,] = crossover(population[p1,],population[p2,],length_of_chromosome)
      }
      population = new_pop 
    }
    
    results = max(classification.history[number_of_generations,])
    mean.evolution = rowMeans(classification.history)
    max.evolution = apply(classification.history,1,max)
    min.evolution = apply(classification.history,1,min)
    best.guy = which(colMeans(classification.history) == max(colMeans(classification.history)))
    best.guy = classification.history[,best.guy]
    
    best.guy.plot = ggplot() + 
      geom_line(aes(x=1:number_of_generations,y=best.guy)) +
      ylab("% correctly classified in test set (1/3 of data)") +
      xlab("Generation") +
      ylim(0,1) +
      ggtitle("Best guy performance - (best mean classification)")
    
    plot.evolution = 
      ggplot() + 
          geom_line(aes(x=1:number_of_generations,y=mean.evolution)) +
          geom_line(aes(x=1:number_of_generations,y=max.evolution),col="gray",alpha=0.5) +
          geom_line(aes(x=1:number_of_generations,y=min.evolution),col="gray",alpha=0.5) +
          ylab("% correctly classified in test set (1/3 of data)") +
          xlab("Generation") +
          ylim(0,1) 
    
    
    individual.generation.process = melt(fitness)
    names(individual.generation.process) = c("Generation","Individual","value")
    individual.generation.process$Individual = as.factor(individual.generation.process$Individual)
    
    plot.individual = 
      ggplot(individual.generation.process) +
       geom_line(aes(x=Generation,y=value,col=Individual))
      
    result.list = list(final.pop = population,
                       mean.evolution = mean.evolution,
                       max.evolution = max.evolution,
                       fitness.process=fitness,
                       plot.individual = plot.individual,
                       best.guy.plot=best.guy.plot,
                       plot.evolution = plot.evolution)
    elapsed.time = Sys.time() - start.time
    cat("\n \n \n --> DONE \n Time elapsed:",round(elapsed.time,2)  ,attributes(elapsed.time)$units )
    return(result.list)
  }
  

####################################
########  TEST RUN   !    ##########
####################################  
  
  
  
test.run = GA.proteo(length_of_chromosome = ncol(proteo)-1, # length of string within indviduals = n of predictors
                      number_of_generations = 500,  # number of generations
                      size_of_population = 50, # number of individuals
                      data = proteo[,-length(proteo)], # predictors
                      Y = proteo[,length(proteo)],  # class to be predicted
                      mutation_rate = 0.005,  # how many values are mutated every round
                      k = 10) # assess fitness using k-nearest neighbours

# how did each individual evolve?
# test.run$plot.individual

# how did the best guy evolve? (the best, on average!)
test.run$best.guy.plot

# how did the group of individuals evolve, on average
test.run$plot.evolution




