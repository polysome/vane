#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include "QUAD16.c"

int 
main(int argc, char *argv[]) 
{

	const int numberofruns = 100;
	div_t divresult;
	printf("QUAD over GF(16)\n");
	printf("1 for cyclicQUAD or 2 for QUADLRS\n");
	choice = 1;
	printf("Number of variables: \n");
	nvar = 20;
	printf("Length of keystream in byte:\n");
	L = 32;
	divresult = div (2 * L,nvar);
	klength = divresult.quot;
	if (divresult.rem > 0)
		klength = klength + 1;
	printf("Repetitions: %d\n", numberofruns);
	printf("Keylength: %d\n", klength);
	D = (nvar + 1) * nvar / 2;
	
	state = (int *)malloc(nvar * sizeof(int));
	usedstate = malloc(numberofruns * klength * sizeof(int*));
	
	printf("Beginning state: \n");
	for (int i = 0; i < nvar; i++){
		state[i] = (rand() % 16);
		printf("%d ", state[i]);
	}
	printf("\n");

if (choice == 1)
{
	int i;
	for(i = 0; i < numberofruns; i++){
		generate_MPcyclic();
		gen_keystreamcyclic();
    }
}

if (choice == 2){
	for(int i = 0; i < numberofruns; i++){
		generate_Gamma();
		generate_MPLRS();
		gen_keystreamLRS(i);	
	}
	for (int j = 0; j < numberofruns*klength; j++){
		free (usedstate[j]);
	}
	free (Gamma);
	free (Delta);
}

	printf("Final state: \n");
	for (int i = 0; i < nvar; i++)
		printf("%d ", state[i]);
	printf("\n");
	free (state);
	free (usedstate);
	return 0;
}