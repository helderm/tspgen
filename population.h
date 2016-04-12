#ifndef POPULATION_H_
#define POPULATION_H_

#include "common.h"

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
int calculateFitnessChromosome(int *chromosome);
int generateNewPopulation(tspsPopulation_t *pop, tspsConfig_t *config);
int compare (const void *a, const void *b);
#endif /* POPULATION_H_ */

