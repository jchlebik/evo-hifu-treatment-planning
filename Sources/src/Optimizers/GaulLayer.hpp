#include <FitnessFunctions/FitnessFunction.hpp>

#include <Optimizers/Optimizer.hpp>

#include <DataModels/GeneralOptions.hpp>

#include <Logging/EvoLogGenerator.hpp>

extern "C"
{
  #include <gaul/gaul.h>
}

#include <fmt/format.h>

#include <vector>
#include <limits>
#include <cmath>

class DifferentialEvolution;
class GeneticAlgorithm;
class SimulatedAnnealing;
class TabuSearch;

class GaulLayer
{
  friend DifferentialEvolution;
  friend GeneticAlgorithm;
  friend SimulatedAnnealing;
  friend TabuSearch;

protected:

  static inline int mPopId;

  static bool gaulScoreCallback(population* pop, 
                                entity* entity)
  {
    Optimizer* optimizer = reinterpret_cast<Optimizer**>(pop->data)[0];

    std::vector<float> chromo(pop->len_chromosomes);
    for (int i = 0; i < pop->len_chromosomes; i++)
    {
      chromo.push_back(static_cast<float>(static_cast<double*>(entity->chromosome[0])[i]));
    }

    float fitness = optimizer->score(chromo);

    //GAUL is a maximization tool. The problems we are solving are minimization problems, we need to convert them.
    entity->fitness = -1.0f * fitness;   
    optimizer->mRunParams.mFitnessCallsPerformed++;
    return true;
  }

  static bool gaulUniformBoundedSeedCallback(population* pop, 
                                             entity* adam)
  {
    Optimizer* optimizer = reinterpret_cast<Optimizer**>(pop->data)[0];

    /* Checks. */
    if (!pop) 
    {
      die("Null pointer to population structure passed.");
    }
    if (!adam)
    {
       die("Null pointer to entity structure passed.");
    }

    std::vector<float> upperBounds = optimizer->mFitnessFunction.getUpperBounds();
    std::vector<float> lowerBounds = optimizer->mFitnessFunction.getLowerBounds();

    for (int chromo = 0; chromo < pop->num_chromosomes; chromo++)
    {
      // Seeding. 
      for (int gene = 0; gene < pop->len_chromosomes; gene++)
      {
        reinterpret_cast<double**>(adam->chromosome)[chromo][gene] =
            random_double_range(static_cast<double>(lowerBounds[gene]), static_cast<double>(upperBounds[gene]));
      }
    }
    return true;
  }


  static void mutateDoubleMultipoint(population*  pop, 
                                     entity*      father, 
                                     entity*      son)
  {
    Optimizer* optimizer      = reinterpret_cast<Optimizer**>(pop->data)[0];
    float*     mutationStdDev = reinterpret_cast<float**>(pop->data)[1];
    //float*     mutationDev  = reinterpret_cast<float**>(pop->data)[1];

    /* Checks */
    if (!father || !son)
    {
      die("Null pointer to entity structure passed");
    }

    /* Copy chromosomes of parent to offspring. */
    for (int i = 0; i < pop->num_chromosomes; i++)
    {
      memcpy(son->chromosome[i], father->chromosome[i], pop->len_chromosomes * sizeof(double));
    }
    
    std::vector<float> lowerBounds = optimizer->mFitnessFunction.getLowerBounds();
    std::vector<float> upperBounds = optimizer->mFitnessFunction.getUpperBounds();

    /*
    * Mutate by tweaking alleles.
    */
    for (int chromo = 0; chromo < pop->num_chromosomes; chromo++)
    {
      for (int gene = 0; gene < pop->len_chromosomes; gene++)
      {
        if (random_boolean_prob(pop->allele_mutation_prob))
        {
          float scale = (1.0f - random_unit_uniform() * 2.0f) * (*mutationStdDev);
          float mutation = std::abs(upperBounds[gene] - lowerBounds[gene]) * scale;
          reinterpret_cast<double**>(son->chromosome)[chromo][gene] += mutation;
        }
      }
    }
  }


  static bool checkTabuListRoundingFloatsCallback(population* pop,
				                                          entity*     putative,
				                                          entity*     tabu)
  {
    double* a = nullptr;
    double* b = nullptr;         /* Comparison double arrays. */

    /* Checks. */
    if ( !pop ) die("Null pointer to population structure passed.");
    if ( !putative || !tabu ) die("Null pointer to entity structure passed.");

    for (int i = 0; i < pop->num_chromosomes; i++)
    {
      a = static_cast<double*>(putative->chromosome[i]);
      b = static_cast<double*>(tabu->chromosome[i]);

      for (int j = 0; j < pop->len_chromosomes; j++)
      {
        if (static_cast<int>(std::round(a[j])) != static_cast<int>(std::round(b[j]))) 
        {
          return false;
        }
      }
    }

    return true;
  }

  static bool gaulGenerationCallback(const int currentGen,
                                     population* pop)
  {
    Optimizer* optimizer = reinterpret_cast<Optimizer**>(pop->data)[0];

    entity* bestEntity = ga_get_entity_from_rank(pop, 0);
    const float bestNow = -1.0f * bestEntity->fitness;

    if (bestNow < optimizer->mRunParams.mBestEver) 
    { 
      optimizer->mRunParams.mBestEver = bestNow; 
      for (int i = 0; i < optimizer->mRunParams.mProblemDimension; i++)
      {
        optimizer->mRunParams.mBestEverChromosome[i] =
            static_cast<float>(reinterpret_cast<double**>(bestEntity->chromosome)[0][i]);
      }
    }

    std::vector<float> popFitnesses(pop->size);
    for(int i = 0; i < pop->size; i++)
    {
      if (pop->entity_iarray[i] != nullptr 
          && pop->entity_iarray[i]->fitness != GA_MIN_FITNESS)
      {
        //GAUL maximazes but we want to minimize. For this goal, we invert the scores before 
        //passing them to GAUL. Here we convert them back for logging.
        popFitnesses.push_back(-1.0f * static_cast<float>(pop->entity_iarray[i]->fitness));
      }
    }

    std::string genString =
        optimizer->mLogGenerator.createEvolutionLogString(currentGen, popFitnesses, optimizer->mRunParams.mBestEver);
    fmt::print(optimizer->mRunParams.mLogFile, "{:s}", genString);

    optimizer->mRunParams.mGenerationsPerformed++;
    return !optimizer->stoppingCriteriaMet(optimizer->mRunParams.mBestEver);
  }


  static bool gaulIterationCallback_SA(const int currentGen,
                                       entity* bestEntity)
  {
    population* pop = ga_get_population_from_id(mPopId);
    Optimizer* optimizer = reinterpret_cast<Optimizer**>(pop->data)[0];

    const float bestNow = -1.0f * bestEntity->fitness;

    if (bestNow < optimizer->mRunParams.mBestEver) 
    { 
      optimizer->mRunParams.mBestEver = bestNow; 
      for (int i = 0; i < optimizer->mRunParams.mProblemDimension; i++)
      {
        optimizer->mRunParams.mBestEverChromosome[i] =
            static_cast<float>(reinterpret_cast<double**>(bestEntity->chromosome)[0][i]);
      }
    }

    //hacky solution to find the candidate solution of current generation for SA
    //FOR SA
    //somewhere inside the array are 3 entity structures
    //1) 'initial'  
    //2) 'best'       - we get this one in parameters
    //3) 'putative'
    //with the way the SA is implemented in GAUL, three possible combinations exist.
    //a) the putative solution was the best found so far
    //---> putative is new best and is copied to initial. 
    //---> putative solution is in the array thrice. we GET it in params
    //b) the putative solution was worse than initial but accepted by SA. 
    //---> putative solution is new best but NOT copied to initial. 
    //---> putative solution is in the array twice while being the smaller fitness. we DIDNT get it in params
    //c) the candidate solution was worse and wasnt accepted by SA.
    //---> putative solution is the smallest value in array
    //---> WE DIDNT get it in params

    double bestFit = bestEntity->fitness;
    double currentFit = std::numeric_limits<double>::infinity();

    //sort_population(pop);
    int bestCounter     = 0;
    int otherCounter    = 0;
    float smallestValue = bestFit;

    for(int i = 0; i < pop->size; i++)
    {
      if (pop->entity_iarray[i] != nullptr                      // GAUL reserves more space than it uses
          && pop->entity_iarray[i]->fitness != GA_MIN_FITNESS)  // GAUL preallocates. GA_MIN_FITNESS is default value
      {
        pop->entity_iarray[i]->fitness == bestFit ? bestCounter++ : otherCounter++;

        if (pop->entity_iarray[i]->fitness < smallestValue)
        {
          smallestValue = pop->entity_iarray[i]->fitness;
        }
      }
    }

    if (bestCounter > otherCounter) //a)
    {
      currentFit = bestFit;
    }
    else //b) and c)
    {
      currentFit = smallestValue;
    }

    std::string genString =
        optimizer->mLogGenerator.createEvolutionLogString(currentGen, -1.0f * currentFit, optimizer->mRunParams.mBestEver);
    fmt::print(optimizer->mRunParams.mLogFile, "{:s}", genString);

    optimizer->mRunParams.mGenerationsPerformed++;
    return !optimizer->stoppingCriteriaMet(optimizer->mRunParams.mBestEver);
  }



  static bool gaulIterationCallback_Tabu(const int currentGen,
                                         entity* bestEntity)
  {
    population* pop = ga_get_population_from_id(mPopId);
    Optimizer* optimizer = reinterpret_cast<Optimizer**>(pop->data)[0];

    const float bestNow = -1.0f * bestEntity->fitness;

    if (bestNow < optimizer->mRunParams.mBestEver) 
    { 
      optimizer->mRunParams.mBestEver = bestNow; 
      for (int i = 0; i < optimizer->mRunParams.mProblemDimension; i++)
      {
        optimizer->mRunParams.mBestEverChromosome[i] =
            static_cast<float>(reinterpret_cast<double**>(bestEntity->chromosome)[0][i]);
      }
    }

    double bestFit = bestEntity->fitness;
    double currentFit = std::numeric_limits<double>::infinity();

    std::vector<float> putativeList;

    for(int i = 0; i < pop->size; i++)
    {
      if (pop->entity_iarray[i] != nullptr                      // GAUL reserves more space than it uses
          && pop->entity_iarray[i]->fitness != GA_MIN_FITNESS)  // GAUL preallocates. GA_MIN_FITNESS is default value
      {
        putativeList.push_back(-1.0f * pop->entity_iarray[i]->fitness);
      }
    }

    std::string genString =
        optimizer->mLogGenerator.createEvolutionLogString(currentGen, putativeList, optimizer->mRunParams.mBestEver);
    fmt::print(optimizer->mRunParams.mLogFile, "{:s}", genString);

    optimizer->mRunParams.mGenerationsPerformed++;
    return !optimizer->stoppingCriteriaMet(optimizer->mRunParams.mBestEver);
  }

  static bool gaulTournamentSelect(population *pop, 
                                   entity **mother, 
                                   entity **father)
  {
    Optimizer* optimizer = reinterpret_cast<Optimizer**>(pop->data)[0];
    entity	*challenger1, *challenger2;	/* Random competitors. */

    if (!pop) 
    {
      die("Null pointer to population structure passed.");
    }

    if (pop->orig_size < 2)
    {
      *mother = nullptr;
      *father = nullptr;
      return true;
    }

    *mother = pop->entity_iarray[random_int(pop->orig_size)];
    *father = pop->entity_iarray[random_int(pop->orig_size)];

    while (*father == *mother)
    {
      *father = pop->entity_iarray[random_int(pop->orig_size)];
    }

    for (int i = 0; i < optimizer->mRunParams.mGATournamentSelectSize; i++)
    {
      challenger1 = pop->entity_iarray[random_int(pop->orig_size)];
      challenger2 = pop->entity_iarray[random_int(pop->orig_size)];

      if (challenger1->fitness > (*mother)->fitness)
      {
        *mother = challenger1;
      }

      if (challenger2 != *mother && challenger2->fitness > (*father)->fitness)
      {
        *father = challenger2;
      }
    }
    
    pop->select_state++;

    return pop->select_state > (pop->orig_size * pop->crossover_ratio);
  }
};