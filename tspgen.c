#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include "mpi.h"
#include "common.h"
#include "map.h"

int mpiId = 0;

int main(int argc, char **argv){
    tspsMap_t map;
    tspsConfig_t config;
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

    // initialize random seed:
    srand ( time(NULL)*mpiId );

    logg("* Parsing map...\n");

    // parse the map
    if(parseMap(&map) != TSPS_RC_SUCCESS){
        printf("Error! Unable to read map 'maps/brazil58.tsp'!\n");
        return TSPS_RC_FAILURE;
    }

    // start a timer (mpi_barrier + mpi_wtime)

    while(1){

        numGenerations++;

        // sort population by fitness

        // crossover population

        // mutate population

        // end execution
        if(config.numGenerations > 0 && numGenerations == config.numGenerations){
            logg("* Max number of generations [%d] reached!\n", config.numGenerations);
            break;
        }

        // migrate population at every n generation
    }

    // join all the populations

    // get the best inidividual and print it

    // stop the timer

    logg("* tspgen finished!\n");
    MPI_Finalize();
    return TSPS_RC_SUCCESS;

}
