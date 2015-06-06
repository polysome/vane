#include "QUAD16.h"

unsigned long long 
rdtsc(void)
{
  unsigned long long int x;
     __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
     return x;
}

double * 
ave_max_min(double time[], int numberofruns)
{
	static double result_time[3];
	result_time[0] = 0;
	result_time[1] = 0;
	result_time[2] = 10000000;

	for (int i = 0; i < numberofruns; i++){
		result_time[0] = result_time[0] + time[i];

		if (result_time[1] <= time[i])
			result_time[1] = time[i];


		if (result_time[2] >= time[i])
			result_time[2] = time[i];

	}

	result_time[0] = result_time[0] / (double)numberofruns;
	return result_time;
}

unsigned long long * 
ave_max_min_cpu(unsigned long long cpucycles[], int numberofruns)
{
	static unsigned long long result[3];
	result[0] = 0;
	result[1] = 0;
	result[2] = 1000000000;

	for (int i = 0; i < numberofruns; i++){
		result[0] = result[0] + cpucycles[i];

		if (result[1] <= cpucycles[i])
			result[1] = cpucycles[i];


		if (result[2] >= cpucycles[i])
			result[2] = cpucycles[i];

	}

	result[0] = result[0] / numberofruns;
	return result;
}

void 
free_arrays (int ***array, int count)
{
	for (int i = 0; i < count; i++){
		for (int j = 0; j < count + 1; j++){
			free (array[i][j]);
		}
		free (array[i]);
	}
	free (array);
}

void 
generate_MPK()
{
	PK = malloc (nvar * sizeof (int *));
	for (int i = 0; i < nvar; i++)
		PK[i] = malloc (D * sizeof (int));

	for (int i = 0; i < nvar; i++){
		for (int j = 0; j < D; j++)
			PK[i][j] = (rand () % 16);
	}

	QK = malloc(nvar * sizeof (int *));
	for (int i = 0; i < nvar; i++)
		QK[i] = malloc (D * sizeof (int));

	for (int i = 0; i < nvar; i++){
		for (int j = 0; j < D; j++)
			QK[i][j] = (rand () % 16);
	}
}

void 
generate_MP()
{
	P = malloc (nvar * sizeof (int **));
	for (int i = 0; i < nvar; i++) {
		P[i] = malloc ((nvar + 1) * sizeof (int *));
		for (int j = 0; j < (nvar + 1); j++)
			P[i][j] = malloc ((nvar + 1) * sizeof (int));
	}

	for (int i = 0; i < nvar; i++){
		for (int j = 0; j < nvar; j++){
			for (int k = j; k < nvar; k++){
				P[i][j][k] = (rand () % 16);
			}
		}
		for (int j = 0; j < nvar; j++){
			for (int k = 0; k < j; k++){
				P[i][j][k] = 0;
			}
		}
	}

	Q = malloc (nvar * sizeof (int **));
	for (int i = 0; i < nvar; i++) {
		Q[i] = malloc ((nvar + 1) * sizeof (int *));
		for (int j = 0; j < (nvar + 1); j++)
			Q[i][j] = malloc ((nvar + 1) * sizeof (int));
	}

	for (int i = 0; i < nvar; i++){
		for (int j = 0; j < nvar; j++){
			for (int k = j; k < nvar; k++){
				Q[i][j][k] = (rand () % 16);
			}
		}
		for (int j = 0; j < nvar; j++){
			for (int k = 0; k < j; k++){
				Q[i][j][k] = 0;
			}
		}
	}
}

void 
generate_Gamma()
{
	Gamma = malloc ((nvar + 1) * sizeof (int));
	Delta = malloc ((nvar + 1) * sizeof (int));

	for (int i = 0; i < nvar; i++){
		Gamma[i] = (rand () % 16);
		Delta[i] = (rand () % 16);
	}
}

void 
generate_MPLRS()
{
	P = malloc (nvar * sizeof(int **));
	for (int i = 0; i < nvar; i++) {
		P[i] = malloc ((nvar + 1) * sizeof (int *));
		for (int j = 0; j < (nvar + 1); j++)
  			P[i][j] = malloc ((nvar + 1) * sizeof (int));
	}

	for (int k = 0; k < nvar; k++){
		P[k][0][0] = 1;
		for (int j = 1; j < nvar; j++)
			P[k][0][j] = lookup[Gamma[k] * 16 + P[k][0][j - 1]];
		for (int i = 1; i < nvar; i++){
			P[k][i][i] = lookup[Gamma[k] * 16 + P[k][i - 1][nvar - 1]];
			for (int j = i + 1; j < nvar; j++)
				P[k][i][j] = lookup[Gamma[k] * 16 + P[k][i][j - 1]];
		}
		for (int i = 0; i < nvar; i++){
			for (int j = 0; j < i; j++){
				P[k][i][j] = 0;
			}
		}
	}

	Q = malloc (nvar * sizeof (int **));
	for (int i = 0; i < nvar; i++) {
		Q[i] = malloc ((nvar + 1) * sizeof (int *));
		for (int j = 0; j < (nvar + 1); j++)
			Q[i][j] = malloc ((nvar + 1) * sizeof (int));
	}

	for (int k = 0; k < nvar; k++){
		Q[k][0][0] = 1;
		for (int j = 1; j < nvar; j++)
			Q[k][0][j] = lookup[Delta[k] * 16 + Q[k][0][j - 1]];
		for (int i = 1; i < nvar; i++){
			Q[k][i][i] = lookup[Delta[k] * 16 + Q[k][i - 1][nvar - 1]];
			for (int j = i + 1; j < nvar; j++)
				Q[k][i][j] = lookup[Delta[k] * 16 + Q[k][i][j - 1]];
		}
		for (int i = 0; i < nvar; i++){
			for (int j = 0; j < i; j++){
				Q[k][i][j] = 0;
			}
		}
	}
}

void 
generate_MPcyclic()
{

	P = malloc (nvar * sizeof (int **));
	for (int i = 0; i < nvar; i++) {
    	P[i] = malloc ((nvar + 1) * sizeof (int *));
		for (int j = 0; j < (nvar + 1); j++)
      		P[i][j] = malloc ((nvar + 1) * sizeof (int));
	}

	//first polynomial
	for (int i = 0; i < nvar; i++){
		for (int j = i; j < nvar; j++){
			P[0][i][j] = (rand () % 16);
		}
	}

	for (int i = 0; i < nvar; i++){
		for (int j = 0; j < i; j++)
				P[0][i][j] = 0;
	}

	//polynomials 2 - n
	for (int k = 1; k < nvar; k++){
		P[k][0][0] = P[k - 1][nvar - 1][nvar - 1];
		for (int j = 1; j < nvar; j++)
			P[k][0][j] = P[k - 1][0][j - 1];
		for (int i = 1; i < nvar; i++) {
			P[k][i][i] = P[k - 1][i - 1][nvar - 1];
			for (int j = i + 1; j < nvar; j++)
				P[k][i][j] = P[k - 1][i][j - 1];
		}
		for (int i = 0; i < nvar; i++){
			for (int j = 0; j < i; j++){
				P[k][i][j] = 0;
			}
		}
	}

	Q = malloc (nvar * sizeof (int **));
	for (int i = 0; i < nvar; i++) {
    	Q[i] = malloc ((nvar + 1) * sizeof (int **));
		for (int j = 0; j < (nvar + 1); j++)
      		Q[i][j] = malloc ((nvar + 1) * sizeof (int **));
	}

	// first polynomial
	Q[0][0][0] = P[nvar - 1][nvar - 1][nvar - 1];
	for (int j = 1; j < nvar; j++)
		Q[0][0][j] = P[nvar - 1][0][j - 1];
	for (int i = 1; i < nvar; i++){
		Q[0][i][i] = P[nvar - 1][i - 1][nvar - 1];
		for (int j = i + 1; j < nvar; j++){
			Q[0][i][j] = P[nvar - 1][i][j - 1];
		}
	}

	for (int i = 0; i < nvar; i++){
		for (int j = 0; j < i; j++)
			Q[0][i][j] = 0;
	}

	//polynomials 2 - n
	for (int k = 1; k < nvar; k++){
		Q[k][0][0] = Q[k - 1][nvar - 1][nvar - 1];
		for (int j = 1; j < nvar; j++)
			Q[k][0][j] = Q[k - 1][0][j - 1];
		for (int i = 1; i < nvar; i++) {
			Q[k][i][i] = Q[k - 1][i - 1][nvar - 1];
			for (int j = i + 1; j < nvar; j++)
				Q[k][i][j] = Q[k - 1][i][j - 1];
		}
		for (int i = 0; i < nvar; i++){
			for (int j = 0; j < i; j++)
				Q[k][i][j] = 0;
		}
	}
}

int * 
evaluatePK(int **MPK, int *state)
{

	int *mon = malloc (D * sizeof (int));
	int *res = malloc (nvar * sizeof (int));
	int counter = 0;

	for (int j = 0; j < nvar; j++){
		for (int l = j; l < nvar; l++){
			mon[counter] = lookup[state[j] * 16 + state[l]];
			counter++;
		}
	}

	for (int i = 0; i < nvar; i++){
		res[i] = 0;
		for (int j = 0; j < D; j++)
			res[i] = res[i] ^ lookup[MPK[i][j] * 16 + mon[j]];
	}

	free (mon);
	return res;
}

int * 
evaluate(int ***MP, int *state)
{
	int *temp = malloc (nvar * sizeof (int));
	int *res = malloc (nvar * sizeof (int));

	for (int i = 0; i < nvar; i++){
	  	res[i] = 0;
		for (int j = 0; j < nvar; j++){
			temp[j] = 0;
			for (int k = 0; k <= j; k++){
				temp[j] = temp[j] ^ lookup[MP[i][k][j] * 16 + state[k]];
			}
		}
		for (int j = 0; j < nvar; j++){
			res[i] = res[i] ^ lookup[temp[j] * 16 + state[j]];
		}
	}

	free (temp);
	return res;
}

int * 
evaluateLRS(int ***MP, int *Beta, int *state)
{
	int *temp = malloc (nvar * sizeof (int));
	int *res = malloc (nvar * sizeof (int));

	for (int i = 0; i < nvar; i++){
	  	res[i] = 0;
	  	temp[0] = state[0];
		for (int j = 1; j < nvar; j++){
			temp[j] = lookup[Beta[i] * 16 + temp[j - 1]] ^ lookup[MP[i][j][j] * 16 + state[j]];
		}
		for (int j = 0; j < nvar; j++){
			res[i] = res[i] ^ lookup[temp[j] * 16 + state[j]];
		}
	}
	
	free (temp);
	return res;
}

int * 
evaluatecyclic(int ***MP, int ***MQ, int *state)
{
	int *temp = malloc (nvar * sizeof (int));
	int *res = malloc (2 * nvar * sizeof (int));

	// first polynomial
	res[0] = 0;
	for (int i = 0; i < nvar; i++){
		temp[i] = 0;
		for (int j = 0; j <= i; j++){
			temp[i] = temp[i] ^ lookup[MP[0][j][i] * 16 + state[j]];
		}
	}

	for (int i = 0; i < nvar; i++){
		res[0] = res[0] ^ lookup[temp[i] * 16 + state[i]];
	}

	//polynomials 2 - n
	for (int k = 1; k < nvar; k++){
		res[k] = 0;
		for (int i = (nvar - 1); i > 0; i--){
			temp[i] = temp[i - 1] ^ lookup[MP[k][i][i] * 16 + state[i]];
		}
		temp[0] = lookup[MP[k][0][0] * 16 + state[0]];
		for (int i = 0; i < nvar; i++){
			res[k] = res[k] ^ lookup[temp[i] * 16 + state[i]];
		}
	}

	//Q
	for (int k = nvar; k < 2 * nvar; k++){
		res[k] = 0;
		for (int i = (nvar - 1); i > 0; i--){
			temp[i] = temp[i - 1] ^ lookup[MQ[k - nvar][i][i] * 16 + state[i]];
		}
		temp[0] = lookup[MQ[k - nvar][0][0] * 16 + state[0]];
		for (int i = 0; i < nvar; i++){
			res[k] = res[k] ^ lookup[temp[i] * 16 + state[i]];
		}
	}

	free (temp);
	return res;
}

void 
gen_keystream()
{
	keystream = malloc (nvar * sizeof (int));
	for (int i = 0; i < klength; i++){
		keystream = evaluate (Q, state);
		state = evaluate (P, state);
	}
}

void 
gen_keystreamPK()
{
	keystream = malloc (nvar * sizeof (int));
	for (int i = 0; i < klength; i++){
		keystream = evaluatePK (QK, state);
		state = evaluatePK (PK, state);
	}
}

void 
gen_keystreamLRS(int n)
{
	keystream = malloc (nvar * sizeof (int));
	int **used = malloc (klength * sizeof (int *));

	for (int i = 0; i < klength; i++){
		keystream = evaluateLRS (Q, Delta, state);
		state = evaluateLRS (P, Gamma, state);
		used[i] = keystream;
		usedstate[n * nvar + i] = state;
		free (used[i]);
	}

	free (used);
	free_arrays (P, nvar);
	free_arrays (Q, nvar);
}

void 
gen_keystreamcyclic()
{
	int **used = malloc (klength * sizeof (int *));
	for (int i = 0; i < klength; i++){
		used[i] = evaluatecyclic (P, Q, state); 
		for (int k = 0; k < nvar; k++)
			state[k] = used[i][k];
		free (used[i]);
	}

	free (used);
	free_arrays (P, nvar);
	free_arrays (Q, nvar);

}