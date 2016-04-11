#ifndef POPULATION_H_
#define POPULATION_H_

#include "common.h"

typedef struct{
    int chromossome[NUM_NODES];
    int fitness;
}tspsIndividual_t;

typedef struct{
    tspsIndividual_t *individuals;
    int numIndividuals;
}tspsPopulation_t;

int generatePopulation(tspsPopulation_t *pop, tspsConfig_t *config);

#endif /* POPULATION_H_ */

