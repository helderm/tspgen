#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include "mpi.h"
#include "common.h"
#include "map.h"
#include "population.h"

int mpiId = 0;

int main(int argc, char **argv){
	tspsMap_t map;
    tspsConfig_t config;
    tspsPopulation_t population;

    unsigned long int numGenerations = 0;
    int mpiNumProcs = 0;

    //starting MPI directives
    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD,&mpiNumProcs);
    MPI_Comm_rank(MPI_COMM_WORLD,&mpiId);

    logg("* Starting tspgen...\n");

    // parse the command line args
    if(readConfig(&config, argc, argv) != TSPS_RC_SUCCESS){
        return TSPS_RC_FAILURE;
    }

    // parse the map
    logg("* Parsing map...\n");
    if(parseMap(&map) != TSPS_RC_SUCCESS){
        logg("Error! Unable to read map 'maps/brazil58.tsp'!\n");
        return TSPS_RC_FAILURE;
    }

    // initialize random seed:
    srand ( time(NULL)*mpiId );

    logg("* Generating population...\n");
    if(generatePopulation(&population, &config) != TSPS_RC_SUCCESS){
        logg("Error! Failed to create a new random population!");
        return TSPS_RC_FAILURE;
    }

    // start a timer (mpi_barrier + mpi_wtime)
    while(1){

        numGenerations++;

        calculateFitnessPopulation(&population, &map);

        if(numGenerations % 1000 == 0){
            logg("- Generation %d...\n", numGenerations);
            printIndividual(&population.individuals[0], "Current Top Individual");
        }

        sortPopulation(&population);

        if(config.numGenerations > 0 && numGenerations == config.numGenerations){
            logg("* Max number of generations [%d] reached!\n", config.numGenerations);
            break;
        }
	/*	
	int i,j;
        for(i=0; i< population.numIndividuals; i++){
               	for(j=0; j<NUM_NODES; j++){
                       	printf("%d ", population.individuals[i].chrom[j]);
               	} printf("\n");
       	}
       	printf("\n\n AFTER CROSSOVER:   ....................\n");
	*/
        crossoverPopulation(&population,  &config);
/*
	for(i=0; i< population.numIndividuals; i++){
              	for(j=0; j<NUM_NODES; j++){
                       	printf("%d ", population.individuals[i].chrom[j]);
               	} printf("\n");
        }

        printf("\n\n");
*/	
        mutatePopulation(&population, &config);

        // migrate population at every n generation
        if(numGenerations % config.migrationRate == 0){
            migrateIndividuals(&population, mpiId, mpiNumProcs);
        }
    }

    // join all the populations
    joinPopulations(&population, mpiId, mpiNumProcs);
    sortPopulation(&population);

    // get the best inidividual and print it
    printIndividual(&population.individuals[0], "Top Individual");
    printIndividual(&population.individuals[config.populationSize-1], "Worst (of the top ones) Individual");

    // stop the timer

    logg("* tspgen finished!\n");
    free(population.individuals);
    MPI_Finalize();
    return TSPS_RC_SUCCESS;
}



