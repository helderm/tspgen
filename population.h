#ifndef POPULATION_H_
#define POPULATION_H_

#include "common.h"
#include "map.h"

typedef struct{
    int *chromosome;
    int fitness;
  	int index;
}tspsIndividual_t;

typedef struct{
    tspsIndividual_t *individuals;
    int numIndividuals;
}tspsPopulation_t;

int generatePopulation(tspsPopulation_t *pop, tspsConfig_t *config);
int *generateRandomChromosome(int chSize);
void swap(int *a, int *b);
int calculateFitnessChromosome(int *chromosome, tspsMap_t *map);
int generateNewPopulation(tspsPopulation_t *pop, tspsConfig_t *config);
int compare (const void *a, const void *b);
int sortPopulation(tspsPopulation_t *pop);
int crossoverPopulation(tspsPopulation_t *pop, tspsConfig_t *config);
int mutatePopulation(tspsPopulation_t *pop, tspsConfig_t *config);
int calculateFitnessPopulation(tspsPopulation_t *pop, tspsMap_t *map);
int migrateIndividuals(tspsPopulation_t *pop, int mpiId, int numProcs);
#endif /* POPULATION_H_ */

