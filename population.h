#ifndef POPULATION_H_
#define POPULATION_H_

#include "common.h"

typedef struct{
    int *chromosome;
    int fitness;
}tspsIndividual_t;

typedef struct{
    tspsIndividual_t *individuals;
    int numIndividuals;
}tspsPopulation_t;

int generatePopulation(tspsPopulation_t *pop, tspsConfig_t *config);
int *generateRandomChromosome(int chSize);
void swap(int *a, int *b);
int calculateFitnessChromosome(int *chromosome);

#endif /* POPULATION_H_ */

