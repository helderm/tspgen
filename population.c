#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include "common.h"
#include "population.h"
#include "map.h"

int generatePopulation(tspsPopulation_t *pop, tspsConfig_t *config){

    	pop->numIndividuals = config->populationSize;
    	pop->individuals = (tspsIndividual_t*)malloc(pop->numIndividuals * sizeof(tspsIndividual_t));

    	int i;
	
    	for(i=0; i<10; i++){//pop->numIndividuals; i++){
		pop->individuals[i].chromosome = generateRandomChromosome(NUM_NODES); //<random vector of unique nodes>	
		pop->individuals[i].fitness = calculateFitnessChromosome(pop->individuals[i].chromosome);
	}
	

    	return TSPS_RC_SUCCESS;
}

int *generateRandomChromosome(int chSize){	
	int *arr = (int	*)malloc(chSize*sizeof(int)); 
	int i; 
	for(i=0; i<chSize; i++){
		arr[i]=i;   // city index starts from zero
	}
	
	struct timeval tv;
    	gettimeofday(&tv, NULL);
    	int usec = tv.tv_usec;
    	srand48(usec);
	
	int n =chSize;
    	if (n > 1) {
        	size_t i;
        	for (i = n - 1; i > 0; i--) {
            		size_t j = (unsigned int) (drand48()*(i+1));
            		int t = arr[j];
            		arr[j] = arr[i];
            		arr[i] = t;
        	}
    	}	
	return arr;
}

int calculateFitnessChromosome( int *chromosome){
	int fitnessValue=0;	
	tspsMap_t map;
	if(parseMap(&map) != TSPS_RC_SUCCESS){
        	printf("Error! Unable to read map 'maps/brazil58.tsp'!\n");
        	return TSPS_RC_FAILURE;
    	}
	int i, firstCity, secondCity; 
	for (i=0; i<NUM_NODES-1; i++){
		firstCity = chromosome[i];
		secondCity = chromosome[i+1];
		fitnessValue = fitnessValue +  map.weights[firstCity][secondCity];		
	}
	printf("%d ", fitnessValue);	
	return fitnessValue;
}

void swap (int *a, int *b){
    	int temp = *a;
    	*a = *b;
    	*b = temp;
}
