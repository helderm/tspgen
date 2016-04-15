#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <mpi.h>
#include "common.h"
#include "population.h"
#include "map.h"

int generatePopulation(tspsPopulation_t *pop, tspsConfig_t *config){

    	pop->numIndividuals = config->populationSize;
    	pop->individuals = (tspsIndividual_t*)malloc(pop->numIndividuals * sizeof(tspsIndividual_t));

    	int i;
    	for(i=0; i<pop->numIndividuals; i++){
    		pop->individuals[i].chromosome = generateRandomChromosome(NUM_NODES); //<random vector of unique nodes>

            //pop->individuals[i].fitness = calculateFitnessChromosome(pop->individuals[i].chromosome);
    		pop->individuals[i].index = i;
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

int calculateFitnessChromosome(int *chromosome, tspsMap_t *map){
	int fitnessValue=0;
	int i, firstCity, secondCity;
	for (i=0; i<NUM_NODES-1; i++){
		firstCity = chromosome[i];
		secondCity = chromosome[i+1];
		fitnessValue = fitnessValue +  map->weights[firstCity][secondCity];
	}
	//printf("%d ", fitnessValue);
	return fitnessValue;
}

void swap (int *a, int *b){
    	int temp = *a;
    	*a = *b;
    	*b = temp;
}

int generateNewPopulation(tspsPopulation_t *pop, tspsConfig_t *config){

	qsort(pop->individuals, config->populationSize, sizeof(tspsIndividual_t), compare);
	int i;
	for (i=0; i < config->populationSize; i++){
		printf("%d ", pop->individuals[i].fitness);
	}
	return TSPS_RC_SUCCESS;
}

int compare (const void *a, const void *b)
{

  tspsIndividual_t * popA = (tspsIndividual_t *)a;
  tspsIndividual_t * popB = (tspsIndividual_t *)b;

  return ( popB->fitness - popA->fitness );
}

int mutatePopulation(tspsPopulation_t *pop, tspsConfig_t *config){

    tspsIndividual_t *ind = NULL;
    int alreadySwaped[NUM_NODES];
    int mutationRate = config->mutationRate * 100;
    int index1, index2;     //the index of the nodes to be swapped
    int aux;
    int i, j;

    for(i=config->numElitism; i<pop->numIndividuals; i++){
        if(rand()%100 > mutationRate)
            continue;

        memset(alreadySwaped, 0, sizeof(alreadySwaped));
        ind = &pop->individuals[i];

        //mutate!
        //swap mutationSize nodes in the chromosome
        for(j=0; j<config->mutationSize; j++){
            index1 = rand() % NUM_NODES;

            //if already swaped, jump to the next of the list
            while(alreadySwaped[index1] !=0){
                if(index1 + 1 < NUM_NODES)
                    index1++;
                else
                    index1 = 0;
            }
            alreadySwaped[index1] = 1;

            index2 = rand() % NUM_NODES;

            //if already swaped, jump to the next of the list
            while(alreadySwaped[index2] !=0){
                if(index2 + 1 < NUM_NODES)
                    index2++;
                else
                    index2 = 0;
            }
            alreadySwaped[index2] = 1;

            //swap the nodes
            aux = ind->chromosome[index1];
            ind->chromosome[index1] = ind->chromosome[index2];
            ind->chromosome[index2] = aux;
        }

        //recalculate the fitness
        //ind->fitness = calculateFitnessChromosome(ind->chromosome, map);
    }

    return TSPS_RC_SUCCESS;

}

int sortPopulation(tspsPopulation_t *pop){
    qsort(pop->individuals, pop->numIndividuals, sizeof(tspsIndividual_t), compare);
    return TSPS_RC_SUCCESS;
}

int crossoverPopulation(tspsPopulation_t *pop, tspsConfig_t *config){
    return 0;
}

int calculateFitnessPopulation(tspsPopulation_t *pop, tspsMap_t *map){
    int i;

    for(i=0; i<pop->numIndividuals; i++){
        pop->individuals[i].fitness = calculateFitnessChromosome(pop->individuals[i].chromosome, map);
    }

    return TSPS_RC_SUCCESS;
}

int migrateIndividuals(tspsPopulation_t *pop, int mpiId, int numProcs){

    int *emigrant1 = NULL, *emigrant2 = NULL;
    int *imigrant1 = NULL, *imigrant2 = NULL;
    int i=0, j=-1;
    MPI_Status status;
    int chromSize = sizeof(int) * NUM_NODES;

    //choose randomly for now the individual to migrate
    i = rand() % pop->numIndividuals;
    emigrant1 = pop->individuals[i].chromosome;
    imigrant1 = (int*)malloc(chromSize);

    // choose another emigrant if they are in the middle of the process chain
    if(mpiId > 0 && mpiId < numProcs-1){
        do{
            j = rand() % pop->numIndividuals;
        } while(j == i);

        emigrant2 = pop->individuals[j].chromosome;
        imigrant2 = (int*)malloc(chromSize);
    }

    /* red/black communication*/
    if (mpiId % 2 != 0){
        /* i'm red */

        if(mpiId < numProcs - 1){
            /* send emigrant to black at my right */
            //printf("+ Pr[%d]: Sending to [%d]... \n", mpiId, mpiId+1);
            MPI_Send(emigrant1, NUM_NODES, MPI_INT, mpiId+1, MPI_MIGRATION_TAG, MPI_COMM_WORLD);

            /* receive imigrant from black at my right */
            //printf("+ Pr[%d]: Receiving from [%d]... \n", mpiId, mpiId+1);
            MPI_Recv(imigrant1, NUM_NODES, MPI_INT, mpiId+1, MPI_MIGRATION_TAG, MPI_COMM_WORLD, &status );
        }

        if(mpiId > 0){
            /* send emigrant to black at my left */
            //printf("+ Pr[%d]: Sending to [%d]... \n", mpiId, mpiId-1);
            MPI_Send((mpiId < numProcs - 1? emigrant2 : emigrant1), NUM_NODES, MPI_INT, mpiId-1, MPI_MIGRATION_TAG, MPI_COMM_WORLD);

            /* receive imigrant from black at my left */
            //printf("+ Pr[%d]: Receiving from [%d]... \n", mpiId, mpiId-1);
            MPI_Recv((mpiId < numProcs - 1? imigrant2 : imigrant1), NUM_NODES, MPI_INT, mpiId-1, MPI_MIGRATION_TAG, MPI_COMM_WORLD, &status );
        }

    } else {
        /* i'm black */
        if(mpiId > 0){
            /* receive imigrant from red at my left */
            //printf("- Pr[%d]: Receiving from [%d]... \n", mpiId, mpiId-1);
            MPI_Recv(imigrant1, NUM_NODES, MPI_INT, mpiId-1, MPI_MIGRATION_TAG, MPI_COMM_WORLD, &status );

            /* send emigrant to red at my left */
            //printf("- Pr[%d]: Sending to [%d]... \n", mpiId, mpiId-1);
            MPI_Send(emigrant1, NUM_NODES, MPI_INT, mpiId-1, MPI_MIGRATION_TAG, MPI_COMM_WORLD);
        }

        if(mpiId < numProcs-1){
            /* receive imigrant from red at my right */
            //printf("- Pr[%d]: Receiving from [%d]... \n", mpiId, mpiId+1);
            MPI_Recv((mpiId > 0? imigrant2 : imigrant1), NUM_NODES, MPI_INT, mpiId+1, MPI_MIGRATION_TAG, MPI_COMM_WORLD, &status );

            /* send emigrant to red at my right */
            //printf("- Pr[%d]: Sending to [%d]... \n", mpiId, mpiId+1);
            MPI_Send((mpiId > 0? emigrant2 : emigrant1), NUM_NODES, MPI_INT, mpiId+1, MPI_MIGRATION_TAG, MPI_COMM_WORLD);
        }
    }

    logg("* Migration ended!\n");

    // move imigrants to our current population
    if(numProcs > 1){
        memcpy(pop->individuals[i].chromosome, imigrant1, chromSize);
        pop->individuals[i].fitness = INFINITY;
        logg("- imigrant = %d, %d, %d...\n", pop->individuals[i].chromosome[0], pop->individuals[i].chromosome[1], pop->individuals[i].chromosome[2]);
    }
    free(imigrant1);

    if(mpiId > 0 && mpiId < numProcs-1){
        memcpy(pop->individuals[j].chromosome, imigrant2, chromSize);
        pop->individuals[j].fitness = INFINITY;
        free(imigrant2);
    }
}
