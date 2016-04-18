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
    	for(i=0; i<pop->numIndividuals; i++){
    		pop->individuals[i].chromosome = generateRandomChromosome(NUM_NODES, i); //<random vector of unique nodes>	
    		//pop->individuals[i].fitness = calculateFitnessChromosome(pop->individuals[i].chromosome);
    		//pop->individuals[i].index = i;		
        }  	

	return TSPS_RC_SUCCESS;
}

int *generateRandomChromosome(int chSize, int index){	
	int *arr = (int	*)malloc(chSize*sizeof(int)); 
	//int arr[NUM_NODES]={0};
	int i; 
	for(i=0; i<chSize; i++){
		arr[i]=i;   // city index starts from zero
	}
	
	struct timeval tv;
    	gettimeofday(&tv, NULL);
    	int usec = tv.tv_usec+index;
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
	/*
	for(i=0; i<chSize; i++){
		printf("%d ", arr[i]);   // city index starts from zero
	}printf("\n");
	*/	
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

/*int generateNewPopulation(tspsPopulation_t *pop, tspsConfig_t *config){
	
	qsort(pop->individuals, config->populationSize, sizeof(tspsIndividual_t), compare);	
	int i;
	for (i=0; i < config->populationSize; i++){
		printf("%d ", pop->individuals[i].fitness);
	}
	return TSPS_RC_SUCCESS;
}*/
	
int compare (const void *a, const void *b)
{

  tspsIndividual_t * popA = (tspsIndividual_t *)a;
  tspsIndividual_t * popB = (tspsIndividual_t *)b;

  return ( popA->fitness - popB->fitness );
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

int crossoverPopulation(tspsPopulation_t *pop, tspsPopulation_t *pop_buffer, tspsConfig_t *config){

	int i, j, k, l;
	//pop_buffer->numIndividuals = config->populationSize;
        //pop_buffer->individuals = (tspsIndividual_t*)malloc(pop->numIndividuals * sizeof(tspsIndividual_t));
	
	double fitness_sum;
	for (i=0; i < config->populationSize; i++){
		fitness_sum = fitness_sum +  pop->individuals[i].fitness;
	}
	for (i=0; i < config->populationSize; i++){	
		pop->individuals[i].probability = pop->individuals[i].fitness/fitness_sum  ;
	//	printf("%lf\n", pop->individuals[i].probability);
	}
	
	/*for(i=0; i<config->numElitism; i++ ){
		pop_buffer->individuals[i] = pop->individuals[i];
	}*/	
	/*
	printf("\nTHE ELITES \n");
	for (i=0; i<config->numElitism; i++){
		for (j=0; j<NUM_NODES; j++){
			printf("%d ", pop->individuals[i].chromosome[j] );
		}
		printf("\n");		
	}
	printf("\n");
	*/	
	double rndNumber_one, rndNumber_two ;
	int cross_pone, cross_ptwo;
	
	double offset_one = 0.0;
	double offset_two = 0.0;
	
	int pick_one = 0, pick_two = 0;
	int count = config->numElitism;
	
	int temp;
	int max = NUM_NODES -1;
	int min = 1;
	// numElitism should be a even number	
	while(count < config->populationSize){
		rndNumber_one = rand() / (double) RAND_MAX;
		rndNumber_two = rand() / (double) RAND_MAX;
		count = count+2;
		for (i = 0; i < config->populationSize; i++) {
			offset_one += pop->individuals[i].probability;	
			if (rndNumber_one < offset_one) {
        			pick_one = i;
        			break;
    			}	
		}	
		for (i = 0; i < config->populationSize; i++) {
			offset_two += pop->individuals[i].probability;
    			
			if (rndNumber_two < offset_two) {
        			pick_two = i;
        			break;
    			}	
		}
		offset_one=0; offset_two=0;
		
		//pmx function
		cross_pone = (max - min + 1)*(double)rand()/RAND_MAX + min;		
		cross_ptwo = (max - min + 1)*(double)rand()/RAND_MAX + min; 
		
		if(cross_ptwo<cross_pone){
			temp= cross_ptwo;
			cross_ptwo = cross_pone;
			cross_pone = temp; 
		}
			
	//	printf("count: %d:%d :: %d %d", count-1, count, cross_pone, cross_ptwo);
	//		printf("\n ");
			
		int *child_1 = pop->individuals[count-2].chromosome; //parent(1) number 10
		int *child_2 = pop->individuals[count-1].chromosome; //parent(2) number 11
		/*	
		for (i=0; i<NUM_NODES; i++){
			printf("%d ",child_1[i] );
		}		
			printf("\n ");
		for (i=0; i<NUM_NODES; i++){
			printf("%d ",child_2[i] );
		}		
		*/
		int *vec_1 = (int*)malloc((cross_ptwo-cross_pone)*sizeof(int));
		int *vec_2 = (int*)malloc((cross_ptwo-cross_pone)*sizeof(int));
		int num=0;
		//printf("\n ");
		
		for (i = cross_pone; i < cross_ptwo; i++ ){
			int tem = child_1[i];
			child_1[i] = child_2[i]; //
			child_2[i] = tem; //		 
			
			vec_1[num] = child_1[i];
			vec_2[num] = child_2[i];
			num++;
		}
		
		for(i = 0; i<cross_pone; i++){
			child_1[i] = pop->individuals[count-2].chromosome[i];
			child_2[i] = pop->individuals[count-1].chromosome[i];		 
		}		
		for(i = cross_ptwo; i<NUM_NODES; i++){
			child_1[i] = pop->individuals[count-2].chromosome[i];
			child_2[i] = pop->individuals[count-1].chromosome[i];	
		}
		/*	
		printf("\n");
		for (i = 0; i < num; i++ ){
			printf("%d ", vec_1[i] );
		}
		printf("\n");
		
		for (i = 0; i < num; i++ ){
			printf("%d ", vec_2[i] );	
		}
		printf("\n\n");
		*/
		int buf, flag=0;
		for (i = 0; i < num; i++ ){
			if(flag==1){
				i=i-1;
				flag=0;
			}
			for (j = 0; j < num; j++ ){
				if((vec_1[i]==vec_2[j])&&(vec_1[i]!=vec_2[i] )){
					buf = vec_1[i];
					vec_1[i] = vec_1[j];
					vec_1[j] = buf; 	
					flag=1;
					continue;
				}
			}
		}
		/*	
		for (i = 0; i < num; i++ ){
			printf("%d ",vec_1[i] );
		}
		printf("\n");
		
		for (i = 0; i < num; i++ ){
			printf("%d ",vec_2[i] );	
		}
		printf("\nfinals........................");
		*/
		for (i=0; i<num; i++){
			for(j = 0; j<cross_pone; j++){
				if(vec_1[i]==child_1[j]){
					child_1[j] = vec_2[i]; 	
				}
				if(vec_2[i]==child_2[j]){
					child_2[j] = vec_1[i]; 	
				}
						 
			}		
			for(j = cross_ptwo; j<NUM_NODES; j++){
				if(vec_1[i]==child_1[j]){
					child_1[j] = vec_2[i]; 	
				}
				if(vec_2[i]==child_2[j]){
					child_2[j] = vec_1[i]; 	
				}		
			}	
		}
		/*	
		printf("\n");
		
		for (i = 0; i < NUM_NODES; i++ ){
			printf("%d ", child_1[i] );
		}
		printf("\n");
		
		for (i = 0; i < NUM_NODES; i++ ){
			printf("%d ", child_2[i] );	
		}
		printf("\n\n\n");
		*/
	}
	
	printf("%d \n", pop->individuals[0].fitness);
	for(j=0; j<NUM_NODES; j++){     
        	printf("%d ", pop->individuals[0].chromosome[j]);  
        } printf("\n");
	return 0;
}

int calculateFitnessPopulation(tspsPopulation_t *pop, tspsMap_t *map){
    	int i,j;
	
	/*
	for(i=0; i<pop->numIndividuals; i++){
		for(j=0; j<NUM_NODES; j++){	
			printf("%d ", pop->individuals[i].chromosome[j]);   // city index starts from zero
		} printf("\n");
	}
	
	printf("\n");
	*/
    	for(i=0; i<pop->numIndividuals; i++){
        	pop->individuals[i].fitness = calculateFitnessChromosome(pop->individuals[i].chromosome, map);
	}

    	return TSPS_RC_SUCCESS;
}

