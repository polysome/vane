#include "abc.h"
#include <math.h>
#include <stdio.h>
#include <time.h>



WORD Logtable[FIELD] = 
{ 0, 0, 25, 1, 50, 2, 26, 198, 75, 199, 27, 104, 51, 238, 223, 3, 
100, 4, 224, 14, 52, 141, 129, 239, 76, 113, 8, 200, 248, 105, 28, 193, 
125, 194, 29, 181, 249, 185, 39, 106, 77, 228, 166, 114, 154, 201, 9, 120, 
101, 47, 138, 5, 33, 15, 225, 36, 18, 240, 130, 69, 53, 147, 218, 142, 
150, 143, 219, 189, 54, 208, 206, 148, 19, 92, 210, 241, 64, 70, 131, 56, 
102, 221, 253, 48, 191, 6, 139, 98, 179, 37, 226, 152, 34, 136, 145, 16, 
126, 110, 72, 195, 163, 182, 30, 66, 58, 107, 40, 84, 250, 133, 61, 186, 
43, 121, 10, 21, 155, 159, 94, 202, 78, 212, 172, 229, 243, 115, 167, 87, 
175, 88, 168, 80, 244, 234, 214, 116, 79, 174, 233, 213, 231, 230, 173, 232, 
44, 215, 117, 122, 235, 22, 11, 245, 89, 203, 95, 176, 156, 169, 81, 160, 
127, 12, 246, 111, 23, 196, 73, 236, 216, 67, 31, 45, 164, 118, 123, 183, 
204, 187, 62, 90, 251, 96, 177, 134, 59, 82, 161, 108, 170, 85, 41, 157, 
151, 178, 135, 144, 97, 190, 220, 252, 188, 149, 207, 205, 55, 63, 91, 209, 
83, 57, 132, 60, 65, 162, 109, 71, 20, 42, 158, 93, 86, 242, 211, 171, 
68, 17, 146, 217, 35, 32, 46, 137, 180, 124, 184, 38, 119, 153, 227, 165, 
103, 74, 237, 222, 197, 49, 254, 24, 13, 99, 140, 128, 192, 247, 112, 7};


WORD Alogtable[FIELD] = 
{ 1, 3, 5, 15, 17, 51, 85, 255, 26, 46, 114, 150, 161, 248, 19, 53, 
95, 225, 56, 72, 216, 115, 149, 164, 247, 2, 6, 10, 30, 34, 102, 170, 
229, 52, 92, 228, 55, 89, 235, 38, 106, 190, 217, 112, 144, 171, 230, 49, 
83, 245, 4, 12, 20, 60, 68, 204, 79, 209, 104, 184, 211, 110, 178, 205, 
76, 212, 103, 169, 224, 59, 77, 215, 98, 166, 241, 8, 24, 40, 120, 136, 
131, 158, 185, 208, 107, 189, 220, 127, 129, 152, 179, 206, 73, 219, 118, 154, 
181, 196, 87, 249, 16, 48, 80, 240, 11, 29, 39, 105, 187, 214, 97, 163, 
254, 25, 43, 125, 135, 146, 173, 236, 47, 113, 147, 174, 233, 32, 96, 160, 
251, 22, 58, 78, 210, 109, 183, 194, 93, 231, 50, 86, 250, 21, 63, 65, 
195, 94, 226, 61, 71, 201, 64, 192, 91, 237, 44, 116, 156, 191, 218, 117, 
159, 186, 213, 100, 172, 239, 42, 126, 130, 157, 188, 223, 122, 142, 137, 128, 
155, 182, 193, 88, 232, 35, 101, 175, 234, 37, 111, 177, 200, 67, 197, 84, 
252, 31, 33, 99, 165, 244, 7, 9, 27, 45, 119, 153, 176, 203, 70, 202, 
69, 207, 74, 222, 121, 139, 134, 145, 168, 227, 62, 66, 198, 81, 243, 14, 
18, 54, 90, 238, 41, 123, 141, 140, 143, 138, 133, 148, 167, 242, 13, 23, 
57, 75, 221, 124, 132, 151, 162, 253, 28, 36, 108, 180, 199, 82, 246, 1};


WORD 
add(WORD a, WORD b) 
{
	return a ^ b;
}

WORD 
sub(WORD a, WORD b) 
{
	return a ^ b;
}


WORD 
mul(WORD a, WORD b) 
{
	/* Multiply two elements of GF(2^m) */
	if (a && b)
		return Alogtable[(Logtable[a] + Logtable[b]) % (FIELD - 1)];
	else 
		return 0;
}


WORD 
divs(WORD a, WORD b) 
{
	int j;
	if (b == 0) {
		printf( "Division by zero\n" );
		abort();
	}
	if (a == 0) 
		return (0);

	if ((j = Logtable[a] - Logtable[b]) < 0) 
		j += (FIELD - 1);

	return (Alogtable[j]);
}


WORD 
inv(WORD in) 
{
	/* 0 is self inverting */
	if (in == 0) 
		return 0;
	else
        return Alogtable[((FIELD - 1) - Logtable[in])];
}

/* Compute  A^-1 */
int 
matrixinv(WORD a[], int n)
{ 
   	int *is, *js, i, j, k, l, u, v;
    WORD d, p;
    is = (int *)malloc(n * sizeof(int));
    js = (int *)malloc(n * sizeof(int));
    for (k = 0; k <= n - 1; k++){ 
		d=0;
		for (i = k; i <= n - 1; i++){
			for (j = k; j <= n - 1; j++){ 
				l = i * n + j; p = (a[l]);
				if (p > d){ 
					d = p; 
					is[k] = i; 
					js[k] = j;
				}
			}
		}
		if (d == 0){ 
			free(is); 
			free(js); 
			printf("err**not inv\n");
			return (0);
		}
		if (is[k] != k){
			for (j = 0; j <= n - 1; j++){ 
		      	u = k * n + j; v = is[k] * n + j;
				p = a[u]; 
				a[u] = a[v]; 
				a[v] = p;
			}
		}
		if (js[k]!=k){
			for (i=0; i<=n-1; i++){ 
				u = i * n + k; v = i * n + js[k];
				p = a[u]; 
				a[u] = a[v]; 
				a[v] = p;
			}
		}
		l= k * n + k;
		a[l] = inv(a[l]);
		for (j = 0; j <= n - 1; j++){
			if (j != k){ 
				u= k * n + j; 
				a[u] = mul(a[u], a[l]);
			}
		}
		for (i = 0; i <= n - 1; i++){
			if (i != k){
				for (j = 0; j <= n - 1; j++){
					if (j != k){ 
						u= i * n + j;
						a[u] = a[u] ^ (mul(a[i * n + k], a[k * n + j]));
					}
				}
			}
		}
		for (i = 0; i <= n - 1; i++){
			if (i != k){ 
				u = i * n + k; a[u] = mul(a[u], a[l]);
			}
		}
	}
	for (k = n - 1; k >= 0; k--){ 
		if (js[k] != k){
			for (j = 0; j <= n - 1; j++){ 
				u= k * n + j; 
				v = js[k] * n + j;
				p = a[u]; 
				a[u] = a[v]; 
				a[v] = p;
			}
		}
		if (is[k] != k){
			for (i = 0; i <= n - 1; i++){ 
				u = i * n + k; 
				v = i * n + is[k];
				p = a[u]; 
				a[u] = a[v]; 
				a[v] = p;
			}
		}
	}
    free(is); 
    free(js);
    return (1);
}

/* A is m*n matrix, B is n*k matrix, C=A*B */
void 
matrixmul(WORD a[], WORD b[], int m, int n, int k, WORD c[])
{ 
  	int i, j, l, u;
    for (i = 0; i <= m - 1; i++){
	    for (j = 0; j <= k - 1; j++){ 
			u = i * k + j; 
			c[u] = 0;
			for (l = 0; l <= n - 1; l++)
				c[u] = c[u] ^ (mul(a[i * n + l], b[l * k + j]));
		}
    }
    return;
}

/* Compute transpose(A) */
int 
matrixtranspose(WORD *a,int n,WORD *b)
{
	int i, j;
	for (i = 0; i <= n - 1; i++){
		for(j = 0; j <= n - 1; j++)
			b[i * n + j] = a[j * n + i];
	}
	return (1);
}


/* c = (a_1 * x_1 +...+ a_n * x_n) * (b_1 * x_1 +...+ b_n * x_n)
 * input: a = [a_1,..,a_n] and b = [b_1,...,b_n]
 * output: the coefficents of c = [c_1,...] */
int 
tensorproduct(WORD a[], WORD b[], WORD c[], int n, int r)
{
	int i, j, k, temp = 0;
	WORD t = 0;
	for (i = 0; i < n; i++){
		for (j = i; j < n; j++){
			if (i == j) {
				for (k = 0; k < r; k++){
					t = add(t, mul(a[k * n + i],b[k * n + i]));
				}
				c[temp++] = t;
				t = 0;
			}
			if (i != j) {
				for (k = 0; k < r; k++){
					t = add(t, add(mul(a[k * n + i], b[k * n + j]), 
						mul(a[k * n + j], b[k * n + i])));
				}
				c[temp++] = t;
				t = 0;
			}
		}
	}
	return (1);
}

/* Compute the reduce echelon form of a matrix */
int 
echelonform(WORD M[], int nrows, int ncols)
{
	int lead = 0, rix = 0, iix = 0, j, k;
	WORD temp, de, sb;
	for (rix = 0; rix < nrows; rix++){
		if ( lead >= ncols) {
			return 0;
		}
		iix = rix;
		while (M[iix*ncols+lead] == 0){
			iix++;
			if (iix == nrows){
				iix = rix;
				lead++;
				if (ncols == lead){
					return 0;
				}
			}
		}
		for (j = 0; j < ncols; j++){
			temp = M[rix * ncols + j];
			M[rix * ncols + j] = M[iix * ncols + j];
			M[iix * ncols + j] = temp;
		}
		de = M[rix * ncols + lead];
		if (de != 0) {
			for (j = 0 ; j < ncols; j++){
				M[rix*ncols + j] = divs(M[rix * ncols + j], de);
			}
		}
		for (j = 0; j < nrows; j++){
			if (j != rix){
				sb = M[j * ncols + lead];
				for (k = 0; k < ncols; k++){
					M[j * ncols + k] = sub(M[j * ncols + k], 
									   mul(sb, M[rix * ncols + k]));
				}
			}
		}
		lead++;
	}
	return 1;
}

/* Key pair generation */

int 
keypair(WORD *sk, WORD *pk) 
{
	WORD  S[VARIABLE*VARIABLE], INVS[VARIABLE*VARIABLE], ST[VARIABLE*VARIABLE]; 
	WORD  M[VARIABLE*VARIABLE], U[VARIABLE*VARIABLE], V[VARIABLE*VARIABLE];
    WORD  T[EQUATION*EQUATION], INVT[EQUATION*EQUATION];
	WORD  B[BROW*BCOL*VARIABLE], C[CROW*CCOL*VARIABLE], A[AROW*ACOL*VARIABLE];
	WORD  TEMPA[ACOL*VARIABLE], TEMPB[BROW*VARIABLE], TEMPC[CROW*VARIABLE], TEMPF[CENTRAL_MAP_SIZE];
	WORD *F;
	WORD *FB;
	int i, j, s, t, ij = 0, eof, count, flag = 0;
	WORD W =  0;
	F = (WORD *) malloc (((AROW * BCOL + AROW * CCOL) * CENTRAL_MAP_SIZE) * 
		sizeof (WORD));
	FB = (WORD *) malloc (((AROW * BCOL + AROW * CCOL) * CENTRAL_MAP_SIZE) * 
		 (sizeof (WORD)));

	printf( "************************************************* \n");
	printf( "Start ABC key pair generation \n");
	printf( "************************************************* \n");

	/* Create L1, check that L1 is not singular and
	 * initialize linear transformation  L1 in S */
	eof = 0 ;
	while (eof == 0) {
		for (i = 0 ; i < VARIABLE; i ++){
			for (j = 0 ; j < VARIABLE; j++){
				INVS[i * VARIABLE + j] = S[i * VARIABLE + j] = rand() % FIELD;
			}
		}
		eof =  matrixinv(INVS, VARIABLE);
	}

	/* Create L2, check that L2 is not singular and
	 * initialize linear transformation  L2 in T 
	 * and also put T^-1 on sk */
	eof = 0;
	while (eof == 0){
		for (i = 0; i < EQUATION; i ++){
			for (j = 0; j < EQUATION; j++){
				INVT[i * EQUATION + j] = T[i * EQUATION + j] = rand() % FIELD;
			}
		}
		eof =  matrixinv(INVT, EQUATION);
	}
	/* Store T^(-1) in secrete key sk */
	for (i = 0; i < EQUATION ; i++){
		for (j = 0; j < EQUATION; j++){
			sk[ij++] = INVT[i * EQUATION + j];
		}
	}
	/* Generate the central map */
	/* Generate A */
	for (i = 0; i < AROW * ACOL * VARIABLE; i++) {
		A[i] = rand() % FIELD ;
	}
	/* Generate B and stor B to sk */
	for (i = 0; i < BROW * BCOL * VARIABLE; i++) {
		sk[ij++] = B[i] = rand() % FIELD;
	}
	/* Generate C and stor C to sk */
	for (i = 0; i < CROW * CCOL * VARIABLE; i++) {
		sk[ij++] = C[i] = rand() % FIELD;
	}
	/* Compute A * B */
	count = 0;
	flag = 0;
	for (i = 0;  i < AROW; i++){
		for (j = 0; j < ACOL * VARIABLE; j++) {
			TEMPA[j] = A[i * ACOL * VARIABLE + j];
		}
		for (j = 0; j < BCOL; j++){
			for (s = 0 ; s < BROW ; s++) {
				for (t = 0; t < VARIABLE; t++) {
					TEMPB[s * VARIABLE + t] = B[j * VARIABLE + s * BCOL * 
											  VARIABLE + t];
				}
			}

			tensorproduct(TEMPA, TEMPB, TEMPF, VARIABLE, BROW);
			
			for (s=0; s<VARIABLE; s++){
				for (t=s; t<VARIABLE; t++){
						sk[ij++] = F[flag++] = TEMPF[count++];
				}
			}
			count = 0;
		}
	}
	/* Compute A * C */
	for (i = 0;  i < AROW; i++){
		for (s = 0;  s < ACOL; s++) {
			for (t = 0; t < VARIABLE; t++) {
				TEMPA[s * VARIABLE + t] = A[i * ACOL * VARIABLE + s * 
										  VARIABLE + t];
			}
		}
		for (j = 0; j < CCOL; j++){
			for (s = 0; s < CROW; s++) {
				for (t = 0; t < VARIABLE; t++) {
					TEMPC[s * VARIABLE + t] = C[j * VARIABLE + s * CCOL * 
											  VARIABLE + t];
				}
			}

			tensorproduct(TEMPA, TEMPC, TEMPF, VARIABLE, CROW);

			for (s = 0; s < VARIABLE; s++){
				for (t = s; t < VARIABLE; t++){
						sk[ij++] = F[flag++] = TEMPF[count++];
				}
			}
			count = 0;
		}	
	}	

	/* Compute FoS = Transpose(S) * F_i * S, i in 0...i < AROW * (BCOL + CCOL) - 1 */
	matrixtranspose(S, VARIABLE, ST);
	
	count = 0;
	flag = 0;
	int num = 0;
	for (i = 0; i < AROW * (BCOL + CCOL); i++) {
		for (s=0; s < VARIABLE * VARIABLE; s++){
			M[s] = (WORD)0;
		}
		for (s = 0; s < VARIABLE; s++){
			for (t = s; t < VARIABLE; t++){
					M[s * VARIABLE + t] = F[count++];
			}
		}

		matrixmul(ST, M, VARIABLE, VARIABLE, VARIABLE, U);
		matrixmul(U, S, VARIABLE, VARIABLE, VARIABLE, V);

		for (s = 0; s < VARIABLE; s++){
			for (t = s; t < VARIABLE; t++){
				if (t == s){
					FB[flag++] = V[s * VARIABLE + t];
				}
				if (t != s){
					FB[flag++] = add((WORD)V[s * VARIABLE + t], 
								 (WORD)V[t * VARIABLE + s]);
				}
			}
		}
	}
	/* Compute ToFoS */
	matrixmul(T, FB, EQUATION, EQUATION, CENTRAL_MAP_SIZE, pk);

	/* Store S^-1 to the secret key sk */
	for (i = 0; i < VARIABLE; i++) {
		for (j = 0; j < VARIABLE; j++){
			sk[ij++] = INVS[i * VARIABLE + j];
		}
	}
	free (F);
	free (FB);
	return (1);
}



int 
encryption(WORD *pk, WORD *plaintext, WORD *ciphertext, int plaintextlen) 
{
	int i, j, k, ij = 0, eof = 0 ;
	WORD temp, x[VARIABLE];
	if (plaintextlen != VARIABLE) 
		return (0);

	printf( "************************************************* \n");
	printf( "Start ABC encryption \n");
	printf( "************************************************* \n");

	for (k = 0; k < VARIABLE; k++) 
		x[k] = (WORD)plaintext[k];

	for (k = 0; k < EQUATION; k++) {
		temp = 0;
		for (i = 0; i < VARIABLE; i++) {
			for (j = i; j < VARIABLE; j++) {
				temp =temp ^ mul((WORD)pk[ij++], mul(x[i], x[j]));	
			}
		}
		ciphertext[k] = temp;
	}
	return (1);
}

int 
decryption(WORD *sk, WORD *decrypttext, WORD *ciphertext,  int ciphertextlen) 
{
	int i, j, k, ij = 0, eof = 0;
	WORD  x[AROW * ACOL + VARIABLE], y[EQUATION], y1[EQUATION], y2[EQUATION], 
		  X[VARIABLE];
	WORD coeff[ACOL * (BCOL + CCOL) * (ACOL * AROW + VARIABLE)];
	WORD TP, xn;
	if (ciphertextlen != EQUATION) 
		return (0) ;


	printf( "************************************************* \n");
	printf( "Start ABC decryption \n");
	printf( "************************************************* \n");


	for (k = 0; k < EQUATION; k++) 
		y[k] = (WORD)ciphertext[k];
	for (k = 0; k < AROW * ACOL + VARIABLE; k++) 
		x[k] = (WORD)0;
	for (i = 0; i < ACOL * (BCOL + CCOL) * (ACOL * AROW + VARIABLE); i++)
		coeff[i] = (WORD)0;

	WORD Sqrttable[FIELD];
	for (i = 0; i < FIELD; i++) {
		Sqrttable[i] = mul(i, i);
	}
	/* Compute first T^{-1} */
	k = 0;
	for (i = 0; i < EQUATION; i++){
		y1[i] = 0;
		for (j = 0; j < EQUATION; j++) {
			y1[i] = add(y1[i], mul( (WORD)sk[ij++], y[j]));
		}
	}
    /* Next find F^{-1} */
	for (k = 0; k < ACOL; k++){
		for (i = 0; i < BCOL; i++){
			for (j = 0; j < AROW; j++){
				coeff[(k * BCOL + i) * (ACOL * AROW + VARIABLE) + k * AROW + j] = 
				y1[j * BCOL + i];
			}
		}
		for (i = 0; i < CCOL; i++){
			for (j = 0; j < AROW; j++){
				coeff[ACOL * BCOL * (ACOL * AROW + VARIABLE) + (k * CCOL + i) * 
				(ACOL * AROW + VARIABLE) + k * AROW + j] = 
				y1[AROW * BCOL + j * CCOL + i];
			}
		}
	}
	for (i = 0; i < ACOL * BCOL; i++){
		for (j = AROW * ACOL; j < (ACOL * AROW + VARIABLE); j++){
			coeff[i * (AROW * ACOL + VARIABLE) + j] = (WORD)sk[ij++];
		}
	}
	for (i = 0; i < ACOL * CCOL; i++){
		for (j = AROW * ACOL; j < (ACOL * AROW + VARIABLE); j++){
			coeff[ACOL * BCOL * (AROW * ACOL + VARIABLE) + i * 
			(AROW * ACOL + VARIABLE) + j] = (WORD)sk[ij++];
		}
	}

	echelonform(coeff, ACOL * (BCOL + CCOL), AROW * ACOL + VARIABLE);
	
	WORD WM[VARIABLE];
	for (i = 0; i < VARIABLE - 1; i++){
		WM[i] = coeff[(i + AROW * ACOL) * (ACOL * AROW + VARIABLE) + 
				(ACOL * AROW + VARIABLE - 1)];
	}
	WM[VARIABLE - 1] = 1;
	for (k = 0; k < EQUATION; k++){
		TP = 0;
		for (i = 0; i < VARIABLE; i++){
			for (j = i; j < VARIABLE; j++){
				TP = add(TP, mul(sk[ij++], mul(WM[i], WM[j])));
			}
		}
		y2[k] = TP;
	}
	for (i = 0; i < EQUATION; i++){
		if (y2[i] != 0){
			break;
		}
	}
	TP = mul(inv(y2[i]), y1[i]);
	for (i = 0; i < FIELD; i++){
		if (Sqrttable[i] == TP){
			xn = (WORD)(i);
		}
	}
	X[VARIABLE-1] = xn;
	for (i = 0; i < VARIABLE - 1; i++){
		X[i] = mul(WM[i], X[VARIABLE - 1]);
	}

	/* Compute S^-1 */
	for (i = 0; i < VARIABLE; i++){
		decrypttext[i] = 0;
		for (j = 0; j < VARIABLE; j++) {
			decrypttext[i] = add(decrypttext[i], mul( (WORD)sk[ij++], X[j]));
		}
	}
	return (1);
}



