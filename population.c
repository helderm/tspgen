#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include "common.h"
#include "population.h"

int generatePopulation(tspsPopulation_t *pop, tspsConfig_t *config){

    pop->numIndividuals = config->populationSize;
    pop->individuals = (tspsIndividual_t*)malloc(pop->numIndividuals * sizeof(tspsIndividual_t));

    int i;
    for(i=0; i<pop->numIndividuals; i++){
        //pop->individuals[i].chromossome = <random vector of unique nodes>
    }

    return TSPS_RC_SUCCESS;
}