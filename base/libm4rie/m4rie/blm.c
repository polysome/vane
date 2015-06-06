#include <m4ri/djb.h>
#include <m4rie/blm.h>
#include <m4rie/mzd_ptr.h>


void _mzd_ptr_apply_blm_mzd(mzd_t **X, const mzd_t **A, const mzd_t **B, const blm_t *f) {
	assert((f->H!=NULL) & (f->F!=NULL) & (f->G!=NULL) &
			(f->H->ncols == f->F->nrows) &   
			(f->F->nrows == f->G->nrows));

	mzd_t *t0 = mzd_init(A[0]->nrows, B[0]->ncols);
	mzd_t *t1 = mzd_init(A[0]->nrows, A[0]->ncols);
	mzd_t *t2 = mzd_init(B[0]->nrows, B[0]->ncols);

	for(rci_t i=0; i < f->F->nrows; i++) {
		mzd_set_ui(t1, 0);
		for(rci_t j=0; j < f->F->ncols; j++) {
			if(mzd_read_bit(f->F, i, j)) {
				mzd_add(t1, t1, A[j]);
			}
		}

		mzd_set_ui(t2, 0);
		for(rci_t j=0; j < f->G->ncols; j++) {
			if(mzd_read_bit(f->G, i, j)) {
				mzd_add(t2, t2, B[j]);
			}
		}


		mzd_mul(t0, t1, t2, 0);

		for(rci_t j=0; j < f->H->nrows; j++)
			if(mzd_read_bit(f->H, j, i))
				_mzd_ptr_add_modred(NULL, t0, X, j);
	}

	mzd_free(t0);
	mzd_free(t1);
	mzd_free(t2);
}

blm_t *_blm_djb_compile(blm_t *f) {
	assert((f->f == NULL) && (f->g == NULL) && (f->h == NULL));
	assert((f->F != NULL) && (f->G != NULL) && (f->H != NULL));


	mzd_t *F = mzd_copy(NULL, f->F);
	f->f = djb_compile(F);
	mzd_free(F);
	if(mzd_equal(f->F, f->G))
		f->g = f->f;
	else {
		mzd_t *G = mzd_copy(NULL, f->G);
		f->g = djb_compile(G);
		mzd_free(G);
	}
	mzd_t *H = mzd_copy(NULL, f->H);
	f->h = djb_compile(H);
	mzd_free(H);

	return f;
}

void _mzd_ptr_apply_blm_djb(mzd_t **X, const mzd_t **A, const mzd_t **B, const blm_t *f) {
	assert((f->H!=NULL) & (f->F!=NULL) & (f->G!=NULL) &   \
			(f->H->ncols == f->F->nrows) &                 \
			(f->F->nrows == f->G->nrows));


	mzd_t **t0 = (mzd_t**)m4ri_mm_malloc(sizeof(mzd_t*)*f->F->nrows);
	mzd_t **t1 = (mzd_t**)m4ri_mm_malloc(sizeof(mzd_t*)*f->F->nrows);
	mzd_t **t2 = (mzd_t**)m4ri_mm_malloc(sizeof(mzd_t*)*f->F->nrows);

	for(rci_t i=0; i<f->F->nrows; i++) {
		t1[i] = mzd_init(A[0]->nrows, A[0]->ncols);
		t2[i] = mzd_init(B[0]->nrows, B[0]->ncols);
	}

	djb_apply_mzd_ptr(f->f, t1, A);
	djb_apply_mzd_ptr(f->g, t2, B);

	for(rci_t i=0; i<f->F->nrows; i++) {
		t0[i] = mzd_init(A[0]->nrows, B[0]->ncols);
		mzd_mul(t0[i], t1[i], t2[i], 0);
		mzd_free(t1[i]);
		mzd_free(t2[i]);
	}



	djb_apply_mzd_ptr(f->h, X, (const mzd_t**)t0);

	for(rci_t i=0; i<f->F->nrows; i++) {
		mzd_free(t0[i]);
	}

	m4ri_mm_free(t0);
	m4ri_mm_free(t1);
	m4ri_mm_free(t2);
}



int *crt_init(const deg_t f_len, const deg_t g_len) {
	int *p_best = (int*)m4ri_mm_calloc(M4RIE_CRT_LEN, sizeof(int));
	int  c_best = f_len * g_len;

	int *p = (int*)m4ri_mm_calloc(M4RIE_CRT_LEN, sizeof(int));

	for(deg_t omega=0; omega<8; omega++) {

		deg_t deg_need = f_len+g_len-1-omega;
		deg_t deg_have = 0;
		deg_t deg_poly = 1;

		p[0] = omega;
		for(deg_t d=1; d<M4RIE_CRT_LEN; d++)
			p[d] = 0;

		while (deg_have < deg_need) {
			p[deg_poly] = irreducible_polynomials[deg_poly][0];
			if (deg_have + deg_poly * p[deg_poly] < deg_need) {
				deg_have += deg_poly * p[deg_poly];
			} else {
				deg_t deg_diff = deg_need - deg_have;
				p[deg_poly] = ceil(deg_diff/(double)deg_poly);
				deg_have += deg_poly * p[deg_poly];
			}
			deg_poly ++;
		}

		deg_t deg_diff = deg_have - deg_need;
		if (deg_diff && p[deg_diff] > 0) {
			p[deg_diff]--;
			deg_have -= deg_diff;
		}

		int c = costs[p[0]];
		for(deg_t d=1; d<M4RIE_CRT_LEN; d++)
			c += costs[d] * p[d];

		if (c < c_best) {
			for(deg_t d=0; d<M4RIE_CRT_LEN; d++)
				p_best[d] = p[d];
			c_best = c;
		}
	}  
	m4ri_mm_free(p);
	return p_best;
}

mzd_t *_small_multiplication_map(const deg_t degree) {
	mzd_t *A;
	switch(degree) {
		case 1:
			A = mzd_init(1, 1);
			A->rows[0][0] = 0x1;
			return A;
		case 2:
			A = mzd_init(3, 2);
			A->rows[0][0] = 0x1; A->rows[1][0] = 0x3; A->rows[2][0] = 0x2;
			return A;
		case 3:
			A = mzd_init(6, 3);
			A->rows[0][0] = 0x1; A->rows[1][0] = 0x2; A->rows[2][0] = 0x4;
			A->rows[3][0] = 0x3; A->rows[4][0] = 0x5; A->rows[5][0] = 0x6;
			return A;
		case 4:
			A = mzd_init(9, 4);
			A->rows[0][0] = 0x1; A->rows[1][0] = 0x2; A->rows[2][0] = 0x4; A->rows[3][0] = 0x8; 
			A->rows[4][0] = 0xf; A->rows[5][0] = 0x3; A->rows[6][0] = 0x5; A->rows[7][0] = 0xa; A->rows[8][0] = 0xc;
			return A;
		case 5:
			A = mzd_init(13, 5);
			A->rows[ 0][0] = 0x01; A->rows[ 1][0] = 0x02; A->rows[ 2][0] = 0x08; A->rows[ 3][0] = 0x10; 
			A->rows[ 4][0] = 0x11; A->rows[ 5][0] = 0x03; A->rows[ 6][0] = 0x18; A->rows[ 7][0] = 0x16;
			A->rows[ 8][0] = 0x0d; A->rows[ 9][0] = 0x1b; A->rows[10][0] = 0x17; A->rows[11][0] = 0x1d;
			A->rows[12][0] = 0x1f;
			return A;
		case 6:
			A = mzd_init(17,  6);
			A->rows[ 0][0] = 0x01;     A->rows[ 1][0] = 0x02;     A->rows[ 2][0] = 0x10;      A->rows[ 3][0] = 0x20;
			A->rows[ 4][0] = 0x30;     A->rows[ 5][0] = 0x03;     A->rows[ 6][0] = 0x18;      A->rows[ 7][0] = 0x06;
			A->rows[ 8][0] = 0x12;     A->rows[ 9][0] = 0x0c;     A->rows[10][0] = 0x38;      A->rows[11][0] = 0x07;
			A->rows[12][0] = 0x29;     A->rows[13][0] = 0x25;     A->rows[14][0] = 0x2d;      A->rows[15][0] = 0x1b;
			A->rows[16][0] = 0x3f;
			return A;
		case 7:
			A = mzd_init(22,  7);
			A->rows[ 0][0] = 0x7f;     A->rows[ 1][0] = 0x6e;     A->rows[ 2][0] = 0x3b;     A->rows[ 3][0] = 0x5d;
			A->rows[ 4][0] = 0x6d;     A->rows[ 5][0] = 0x5b;     A->rows[ 6][0] = 0x36;     A->rows[ 7][0] = 0x03;
			A->rows[ 8][0] = 0x05;     A->rows[ 9][0] = 0x11;     A->rows[10][0] = 0x0a;     A->rows[11][0] = 0x44;
			A->rows[12][0] = 0x28;     A->rows[13][0] = 0x50;     A->rows[14][0] = 0x60;     A->rows[15][0] = 0x01;
			A->rows[16][0] = 0x02;     A->rows[17][0] = 0x04;     A->rows[18][0] = 0x08;     A->rows[19][0] = 0x10;
			A->rows[20][0] = 0x20;     A->rows[21][0] = 0x40;     
			return A;
		case 8:
			A = mzd_init(27,  8);
			A->rows[ 0][0] = 0x01;     A->rows[ 1][0] = 0x02;     A->rows[ 2][0] = 0x04;     A->rows[ 3][0] = 0x08;
			A->rows[ 4][0] = 0x10;     A->rows[ 5][0] = 0x20;     A->rows[ 6][0] = 0x40;     A->rows[ 7][0] = 0x80;
			A->rows[ 8][0] = 0x05;     A->rows[ 9][0] = 0x0c;     A->rows[10][0] = 0x44;     A->rows[11][0] = 0x30;
			A->rows[12][0] = 0xc0;     A->rows[13][0] = 0x0f;     A->rows[14][0] = 0xa0;     A->rows[15][0] = 0x11;
			A->rows[16][0] = 0x0a;     A->rows[17][0] = 0xaa;     A->rows[18][0] = 0x03;     A->rows[19][0] = 0x50;
			A->rows[20][0] = 0x22;     A->rows[21][0] = 0xff;     A->rows[22][0] = 0x55;     A->rows[23][0] = 0x33;
			A->rows[24][0] = 0xcc;     A->rows[25][0] = 0x88;     A->rows[26][0] = 0xf0;
			return A;
		case 9:
			A = mzd_init(31,  9);
			A->rows[ 0][0] = 0x100;     A->rows[ 1][0] = 0x001;     A->rows[ 2][0] = 0x002;     A->rows[ 3][0] = 0x003;
			A->rows[ 4][0] = 0x155;     A->rows[ 5][0] = 0x0aa;     A->rows[ 6][0] = 0x1ff;     A->rows[ 7][0] = 0x16d;
			A->rows[ 8][0] = 0x1b6;     A->rows[ 9][0] = 0x0db;     A->rows[10][0] = 0x0e9;     A->rows[11][0] = 0x13a;
			A->rows[12][0] = 0x074;     A->rows[13][0] = 0x1d3;     A->rows[14][0] = 0x09d;     A->rows[15][0] = 0x14e;
			A->rows[16][0] = 0x0b9;     A->rows[17][0] = 0x172;     A->rows[18][0] = 0x05c;     A->rows[19][0] = 0x1cb;
			A->rows[20][0] = 0x0e5;     A->rows[21][0] = 0x12e;     A->rows[22][0] = 0x191;     A->rows[23][0] = 0x0b2;
			A->rows[24][0] = 0x164;     A->rows[25][0] = 0x0c8;     A->rows[26][0] = 0x08f;     A->rows[27][0] = 0x123;
			A->rows[28][0] = 0x0f5;     A->rows[29][0] = 0x07a;     A->rows[30][0] = 0x1ac;
			return A;
		case 10:
			A = mzd_init(36, 10);
			A->rows[ 0][0] = 0x200;     A->rows[ 1][0] = 0x001;     A->rows[ 2][0] = 0x3ff;     A->rows[ 3][0] = 0x36d;
			A->rows[ 4][0] = 0x1b6;     A->rows[ 5][0] = 0x2db;     A->rows[ 6][0] = 0x0e9;     A->rows[ 7][0] = 0x13a;
			A->rows[ 8][0] = 0x274;     A->rows[ 9][0] = 0x1d3;     A->rows[10][0] = 0x29d;     A->rows[11][0] = 0x34e;
			A->rows[12][0] = 0x0b9;     A->rows[13][0] = 0x172;     A->rows[14][0] = 0x25c;     A->rows[15][0] = 0x1cb;
			A->rows[16][0] = 0x2e5;     A->rows[17][0] = 0x32e;     A->rows[18][0] = 0x191;     A->rows[19][0] = 0x2b2;
			A->rows[20][0] = 0x164;     A->rows[21][0] = 0x2c8;     A->rows[22][0] = 0x08f;     A->rows[23][0] = 0x323;
			A->rows[24][0] = 0x0f5;     A->rows[25][0] = 0x07a;     A->rows[26][0] = 0x3ac;     A->rows[27][0] = 0x2f1;
			A->rows[28][0] = 0x1e2;     A->rows[29][0] = 0x3c4;     A->rows[30][0] = 0x178;     A->rows[31][0] = 0x1af;
			A->rows[32][0] = 0x313;     A->rows[33][0] = 0x135;     A->rows[34][0] = 0x09a;     A->rows[35][0] = 0x2bc;
			return A;
		case 11:
			A = mzd_init(40, 11);
			A->rows[ 0][0] = 0x001;     A->rows[ 1][0] = 0x36d;     A->rows[ 2][0] = 0x5b6;     A->rows[ 3][0] = 0x6db;
			A->rows[ 4][0] = 0x555;     A->rows[ 5][0] = 0x2aa;     A->rows[ 6][0] = 0x7ff;     A->rows[ 7][0] = 0x400;
			A->rows[ 8][0] = 0x200;     A->rows[ 9][0] = 0x600;     A->rows[10][0] = 0x4e9;     A->rows[11][0] = 0x53a;
			A->rows[12][0] = 0x274;     A->rows[13][0] = 0x1d3;     A->rows[14][0] = 0x69d;     A->rows[15][0] = 0x74e;
			A->rows[16][0] = 0x4b9;     A->rows[17][0] = 0x172;     A->rows[18][0] = 0x65c;     A->rows[19][0] = 0x5cb;
			A->rows[20][0] = 0x2e5;     A->rows[21][0] = 0x72e;     A->rows[22][0] = 0x591;     A->rows[23][0] = 0x6b2;
			A->rows[24][0] = 0x564;     A->rows[25][0] = 0x2c8;     A->rows[26][0] = 0x48f;     A->rows[27][0] = 0x323;
			A->rows[28][0] = 0x0f5;     A->rows[29][0] = 0x47a;     A->rows[30][0] = 0x7ac;     A->rows[31][0] = 0x2f1;
			A->rows[32][0] = 0x5e2;     A->rows[33][0] = 0x3c4;     A->rows[34][0] = 0x578;     A->rows[35][0] = 0x1af;
			A->rows[36][0] = 0x713;     A->rows[37][0] = 0x135;     A->rows[38][0] = 0x09a;     A->rows[39][0] = 0x6bc;
			return A;
		case 12:
			A = mzd_init(45, 12);
			A->rows[ 0][0] = 0xb6d;     A->rows[ 1][0] = 0xdb6;     A->rows[ 2][0] = 0x6db;     A->rows[ 3][0] = 0x001;
			A->rows[ 4][0] = 0x002;     A->rows[ 5][0] = 0x003;     A->rows[ 6][0] = 0x555;     A->rows[ 7][0] = 0xaaa;
			A->rows[ 8][0] = 0xfff;     A->rows[ 9][0] = 0x800;     A->rows[10][0] = 0x400;     A->rows[11][0] = 0xc00;
			A->rows[12][0] = 0x4e9;     A->rows[13][0] = 0xd3a;     A->rows[14][0] = 0xa74;     A->rows[15][0] = 0x9d3;
			A->rows[16][0] = 0xe9d;     A->rows[17][0] = 0x74e;     A->rows[18][0] = 0x591;     A->rows[19][0] = 0xeb2;
			A->rows[20][0] = 0xd64;     A->rows[21][0] = 0xac8;     A->rows[22][0] = 0xc8f;     A->rows[23][0] = 0xb23;
			A->rows[24][0] = 0x8f5;     A->rows[25][0] = 0x47a;     A->rows[26][0] = 0x7ac;     A->rows[27][0] = 0xaf1;
			A->rows[28][0] = 0x5e2;     A->rows[29][0] = 0xbc4;     A->rows[30][0] = 0xd78;     A->rows[31][0] = 0x9af;
			A->rows[32][0] = 0xf13;     A->rows[33][0] = 0x135;     A->rows[34][0] = 0x89a;     A->rows[35][0] = 0x6bc;
			A->rows[36][0] = 0x631;     A->rows[37][0] = 0xa52;     A->rows[38][0] = 0x294;     A->rows[39][0] = 0x318;
			A->rows[40][0] = 0xdef;     A->rows[41][0] = 0xc63;     A->rows[42][0] = 0x4a5;     A->rows[43][0] = 0x94a;
			A->rows[44][0] = 0x18c;
			return A;
		case 13:
			A = mzd_init(49, 13);
			A->rows[ 0][0] = 0x0001;     A->rows[ 1][0] = 0x1b6d;     A->rows[ 2][0] = 0x0db6;     A->rows[ 3][0] = 0x16db;
			A->rows[ 4][0] = 0x1555;     A->rows[ 5][0] = 0x0aaa;     A->rows[ 6][0] = 0x1fff;     A->rows[ 7][0] = 0x1000;
			A->rows[ 8][0] = 0x0800;     A->rows[ 9][0] = 0x1800;     A->rows[10][0] = 0x14e9;     A->rows[11][0] = 0x1d3a;
			A->rows[12][0] = 0x1a74;     A->rows[13][0] = 0x09d3;     A->rows[14][0] = 0x0e9d;     A->rows[15][0] = 0x074e;
			A->rows[16][0] = 0x1cb9;     A->rows[17][0] = 0x1972;     A->rows[18][0] = 0x0e5c;     A->rows[19][0] = 0x05cb;
			A->rows[20][0] = 0x12e5;     A->rows[21][0] = 0x172e;     A->rows[22][0] = 0x1591;     A->rows[23][0] = 0x1eb2;
			A->rows[24][0] = 0x1d64;     A->rows[25][0] = 0x1ac8;     A->rows[26][0] = 0x0c8f;     A->rows[27][0] = 0x0b23;
			A->rows[28][0] = 0x08f5;     A->rows[29][0] = 0x047a;     A->rows[30][0] = 0x07ac;     A->rows[31][0] = 0x1af1;
			A->rows[32][0] = 0x15e2;     A->rows[33][0] = 0x0bc4;     A->rows[34][0] = 0x0d78;     A->rows[35][0] = 0x09af;
			A->rows[36][0] = 0x0f13;     A->rows[37][0] = 0x1135;     A->rows[38][0] = 0x189a;     A->rows[39][0] = 0x06bc;
			A->rows[40][0] = 0x0631;     A->rows[41][0] = 0x0a52;     A->rows[42][0] = 0x1294;     A->rows[43][0] = 0x0318;
			A->rows[44][0] = 0x1def;     A->rows[45][0] = 0x0c63;     A->rows[46][0] = 0x14a5;     A->rows[47][0] = 0x094a;
			A->rows[48][0] = 0x118c;
			return A;
		case 14:
			A = mzd_init(55, 14);
			A->rows[ 0][0] = 0x0001;     A->rows[ 1][0] = 0x1b6d;     A->rows[ 2][0] = 0x2db6;     A->rows[ 3][0] = 0x36db;
			A->rows[ 4][0] = 0x1555;     A->rows[ 5][0] = 0x2aaa;     A->rows[ 6][0] = 0x3fff;     A->rows[ 7][0] = 0x34e9;
			A->rows[ 8][0] = 0x1d3a;     A->rows[ 9][0] = 0x3a74;     A->rows[10][0] = 0x29d3;     A->rows[11][0] = 0x0e9d;
			A->rows[12][0] = 0x274e;     A->rows[13][0] = 0x1cb9;     A->rows[14][0] = 0x3972;     A->rows[15][0] = 0x2e5c;
			A->rows[16][0] = 0x25cb;     A->rows[17][0] = 0x32e5;     A->rows[18][0] = 0x172e;     A->rows[19][0] = 0x3591;
			A->rows[20][0] = 0x1eb2;     A->rows[21][0] = 0x3d64;     A->rows[22][0] = 0x3ac8;     A->rows[23][0] = 0x2c8f;
			A->rows[24][0] = 0x2b23;     A->rows[25][0] = 0x08f5;     A->rows[26][0] = 0x247a;     A->rows[27][0] = 0x07ac;
			A->rows[28][0] = 0x1af1;     A->rows[29][0] = 0x35e2;     A->rows[30][0] = 0x2bc4;     A->rows[31][0] = 0x0d78;
			A->rows[32][0] = 0x09af;     A->rows[33][0] = 0x2f13;     A->rows[34][0] = 0x3135;     A->rows[35][0] = 0x389a;
			A->rows[36][0] = 0x26bc;     A->rows[37][0] = 0x0631;     A->rows[38][0] = 0x0a52;     A->rows[39][0] = 0x1294;
			A->rows[40][0] = 0x2318;     A->rows[41][0] = 0x3def;     A->rows[42][0] = 0x0c63;     A->rows[43][0] = 0x14a5;
			A->rows[44][0] = 0x294a;     A->rows[45][0] = 0x318c;     A->rows[46][0] = 0x2000;     A->rows[47][0] = 0x1000;
			A->rows[48][0] = 0x0800;     A->rows[49][0] = 0x0400;     A->rows[50][0] = 0x3c00;     A->rows[51][0] = 0x3000;
			A->rows[52][0] = 0x2800;     A->rows[53][0] = 0x1400;     A->rows[54][0] = 0x0c00;
			return A;
		case 15:
			A = mzd_init(60, 15);
			A->rows[ 0][0] = 0x0001;     A->rows[ 1][0] = 0x7fff;     A->rows[ 2][0] = 0x5b6d;     A->rows[ 3][0] = 0x6db6;
			A->rows[ 4][0] = 0x36db;     A->rows[ 5][0] = 0x4000;     A->rows[ 6][0] = 0x2000;     A->rows[ 7][0] = 0x6000;
			A->rows[ 8][0] = 0x74e9;     A->rows[ 9][0] = 0x1d3a;     A->rows[10][0] = 0x3a74;     A->rows[11][0] = 0x69d3;
			A->rows[12][0] = 0x4e9d;     A->rows[13][0] = 0x274e;     A->rows[14][0] = 0x5cb9;     A->rows[15][0] = 0x3972;
			A->rows[16][0] = 0x2e5c;     A->rows[17][0] = 0x65cb;     A->rows[18][0] = 0x72e5;     A->rows[19][0] = 0x172e;
			A->rows[20][0] = 0x7591;     A->rows[21][0] = 0x1eb2;     A->rows[22][0] = 0x3d64;     A->rows[23][0] = 0x7ac8;
			A->rows[24][0] = 0x2c8f;     A->rows[25][0] = 0x6b23;     A->rows[26][0] = 0x48f5;     A->rows[27][0] = 0x647a;
			A->rows[28][0] = 0x47ac;     A->rows[29][0] = 0x1af1;     A->rows[30][0] = 0x35e2;     A->rows[31][0] = 0x6bc4;
			A->rows[32][0] = 0x4d78;     A->rows[33][0] = 0x09af;     A->rows[34][0] = 0x2f13;     A->rows[35][0] = 0x7135;
			A->rows[36][0] = 0x789a;     A->rows[37][0] = 0x26bc;     A->rows[38][0] = 0x4631;     A->rows[39][0] = 0x4a52;
			A->rows[40][0] = 0x5294;     A->rows[41][0] = 0x6318;     A->rows[42][0] = 0x3def;     A->rows[43][0] = 0x0c63;
			A->rows[44][0] = 0x14a5;     A->rows[45][0] = 0x294a;     A->rows[46][0] = 0x318c;     A->rows[47][0] = 0x4d21;
			A->rows[48][0] = 0x1a42;     A->rows[49][0] = 0x7348;     A->rows[50][0] = 0x6690;     A->rows[51][0] = 0x2bb1;
			A->rows[52][0] = 0x5763;     A->rows[53][0] = 0x15d8;     A->rows[54][0] = 0x0576;     A->rows[55][0] = 0x47cd;
			A->rows[56][0] = 0x42bb;     A->rows[57][0] = 0x4857;     A->rows[58][0] = 0x215d;     A->rows[59][0] = 0x3b1f;
			return A;
		case 16:
			A = mzd_init(64, 16);
			A->rows[ 0][0] = 0xdb6d;     A->rows[ 1][0] = 0x6db6;     A->rows[ 2][0] = 0xb6db;     A->rows[ 3][0] = 0x0001;
			A->rows[ 4][0] = 0x0002;     A->rows[ 5][0] = 0x0003;     A->rows[ 6][0] = 0x5555;     A->rows[ 7][0] = 0xaaaa;
			A->rows[ 8][0] = 0xffff;     A->rows[ 9][0] = 0x8000;     A->rows[10][0] = 0x4000;     A->rows[11][0] = 0xc000;
			A->rows[12][0] = 0x74e9;     A->rows[13][0] = 0x9d3a;     A->rows[14][0] = 0x3a74;     A->rows[15][0] = 0xe9d3;
			A->rows[16][0] = 0x4e9d;     A->rows[17][0] = 0xa74e;     A->rows[18][0] = 0x5cb9;     A->rows[19][0] = 0xb972;
			A->rows[20][0] = 0x2e5c;     A->rows[21][0] = 0xe5cb;     A->rows[22][0] = 0x72e5;     A->rows[23][0] = 0x972e;
			A->rows[24][0] = 0xf591;     A->rows[25][0] = 0x1eb2;     A->rows[26][0] = 0x3d64;     A->rows[27][0] = 0x7ac8;
			A->rows[28][0] = 0xac8f;     A->rows[29][0] = 0xeb23;     A->rows[30][0] = 0xc8f5;     A->rows[31][0] = 0x647a;
			A->rows[32][0] = 0x47ac;     A->rows[33][0] = 0x9af1;     A->rows[34][0] = 0x35e2;     A->rows[35][0] = 0x6bc4;
			A->rows[36][0] = 0x4d78;     A->rows[37][0] = 0x89af;     A->rows[38][0] = 0xaf13;     A->rows[39][0] = 0xf135;
			A->rows[40][0] = 0x789a;     A->rows[41][0] = 0x26bc;     A->rows[42][0] = 0xc631;     A->rows[43][0] = 0x4a52;
			A->rows[44][0] = 0x5294;     A->rows[45][0] = 0x6318;     A->rows[46][0] = 0xbdef;     A->rows[47][0] = 0x8c63;
			A->rows[48][0] = 0x94a5;     A->rows[49][0] = 0x294a;     A->rows[50][0] = 0x318c;     A->rows[51][0] = 0xcd21;
			A->rows[52][0] = 0x9a42;     A->rows[53][0] = 0xf348;     A->rows[54][0] = 0xe690;     A->rows[55][0] = 0x2bb1;
			A->rows[56][0] = 0x5763;     A->rows[57][0] = 0x15d8;     A->rows[58][0] = 0x8576;     A->rows[59][0] = 0xc7cd;
			A->rows[60][0] = 0x42bb;     A->rows[61][0] = 0x4857;     A->rows[62][0] = 0x215d;     A->rows[63][0] = 0xbb1f;
			return A;

		default:
			m4ri_die("only degrees up to 16 are implemented but got degree %d\n", degree);
	}
	return NULL;
}

/**
 * \param length The length of the polynomial we want to reduce
 * \param poly A polynomial
 * \param d The degree of poly
 */

mzd_t *_crt_modred_mat(const deg_t length, const word poly, const deg_t d) {
	mzd_t *A = mzd_init(d, length);

	/* (x-infinity)^d */
	if (poly == 0) {
		for(deg_t i=0; i<d; i++) 
			mzd_write_bit(A, i, length-i-1, 1);
		return A;
	}

	mzd_t *f = mzd_init(1, length);
	mzd_t *t = mzd_init(1, length);

	for(deg_t i=0; i<length; i++) {
		mzd_set_ui(f, 0);
		f->rows[0][i/m4ri_radix] = __M4RI_TWOPOW(i%m4ri_radix);
		word ii = i;
		while(ii >= d) {
			/* f ^= gf2x_mul((1ULL<<(ii-d)), poly, length); */
			mzd_set_ui(t, 0);
			mzd_xor_bits(t, 0, ii-d, d+1, poly);

			mzd_add(f, f, t);

			/* ii = gf2x_deg(f); */
			ii = 0;
			for(wi_t j=f->width-1; j>=0; j--) {
				if (f->rows[0][j]) {
					ii = gf2x_deg(f->rows[0][j]) + m4ri_radix*j;
					break;
				}
			}
		}
		for(deg_t j=0; j<= ii; j++) 
			mzd_write_bit(A, j, i, (f->rows[0][j/m4ri_radix]>>(j%m4ri_radix)) & 0x1);
	}
	return A;
}

blm_t *_blm_finish_polymult(const gf2e *ff, blm_t *f) {
	assert( (f != NULL) & (f->F != NULL) & (f->G != NULL) );

	const rci_t m = f->F->nrows;
	const rci_t c_nrows = f->F->ncols + f->G->ncols - 1;

	mzd_t *H = mzd_init(c_nrows, m);

	mzd_t *F_T = mzd_transpose(NULL, f->F);
	mzd_t *G_T = mzd_transpose(NULL, f->G);

	mzd_t *C = mzd_init(m, m);

	/* we find a full rank matrix */

	word v = 0;
	word w = 0;
	rci_t r = 0;
	rci_t rank = 0;

	mzd_t *pivots = mzd_init(m, m4ri_radix*2);

	mzp_t *P = mzp_init(m);
	mzp_t *Q = mzp_init(m);

	while(rank < m) {
		/* x^v = (0, 0, ... , 1, ..., 0, 0) * F_T -> select row v */
		for(wi_t j=0; j< C->width; j++)
			C->rows[r][j] = F_T->rows[v][j] & G_T->rows[w][j];

		pivots->rows[r][0] = v;
		pivots->rows[r][1] = w;

		w++;
		if (w == f->G->ncols) {
			v++;
			if (v == f->F->ncols)
				v = 0;
			w = v;
		}
		r++;
		if (r == C->nrows) {
			mzd_t *D = mzd_copy(NULL, C);
			rank = mzd_ple(D, P, Q, 0);
			mzd_apply_p_left(pivots, P);
			mzd_apply_p_left(C, P);
			mzd_free(D);

			if (rank < m)
				r = rank;
		}
	}
	mzd_free(F_T);
	mzd_free(G_T);
	mzp_free(P);
	mzp_free(Q);

	for(r=0; r<m; r++) {
		/* x^v = (0, 0, ... , 1, ..., 0, 0) * F_T -> select row v */
		v = pivots->rows[r][0];
		w = pivots->rows[r][1];
		for(wi_t j=0; j< C->width; j++)
			C->rows[r][j] = F_T->rows[v][j] & G_T->rows[w][j];
	}

	// This should be replaced by TRSM calls
	mzd_t *D = mzd_inv_m4ri(NULL, C, 0); 
	mzd_free(C);
	mzd_t *DT = mzd_transpose(NULL, D);
	mzd_free(D);

	mzd_t *a = mzd_init(1, m);
	mzd_t *b = mzd_init(1, H->ncols);

	for(rci_t i=0; i<H->nrows; i++) {
		mzd_set_ui(a, 0);
		for(rci_t j=0; j<m; j++) {
			v = pivots->rows[j][0];
			w = pivots->rows[j][1];
			if ((v+w) == i)
				mzd_write_bit(a, 0, j, 1);
		}
		mzd_mul(b, a, DT, 0);

		/* copy result to H */
		for(rci_t j=0; j<H->ncols; j++)
			mzd_write_bit(H, i, j, mzd_read_bit(b, 0, j ));
	}
	mzd_free(a);
	mzd_free(b);
	mzd_free(pivots);

	if (ff == NULL) {
		f->H = H;
	} else { 
		mzd_t *N = _crt_modred_mat(H->nrows, ff->minpoly, ff->degree);
		f->H = mzd_mul(NULL, N, H, 0);
		mzd_free(N);
		mzd_free(H);
	}
	return f;
}

const int costs[17] = {0, 1, 3, 6, 9, 13, 17, 22, 27, 31, 36, 40, 45, 49, 55, 60, 64};
//const int costs[17] = {0, 1, 3, 6, 9, 13, 15, 22, 24, 30, 33, 39, 42, 48, 51, 54, 60}; /* best possible */

blm_t *blm_init_crt(const gf2e *ff, const deg_t f_ncols, const deg_t g_ncols, const int *p, int djb) {
	blm_t *f = m4ri_mm_malloc(sizeof(blm_t));

	// iterator over irreducible polynomials
	int *p_it = (int*)m4ri_mm_calloc(sizeof(int), M4RIE_CRT_LEN); 

	mzd_t *M, *T;

	word poly = 0;

	rci_t m = costs[p[0]];
	for(int d=1; d<M4RIE_CRT_LEN; d++)
		m += costs[d] * p[d];

	f->F = mzd_init(m, f_ncols);
	f->f = NULL;
	f->G = mzd_init(m, g_ncols);
	f->g = NULL;

	rci_t r = 0;


	/**
	 * 1) We construct maps F,G which combine modular reduction to degree d and the linear map
	 *    required for multiplying modulo a polynomial of degree d.
	 */

	/**
	 * 1.1) We deal with (x+infinity)^omega first
	 */

	if(p[0] != 0) {
		deg_t d = p[0];
		mzd_t *N  = _small_multiplication_map(d);
		M = _crt_modred_mat(f_ncols, poly, d);
		T = mzd_init_window(f->F, r, 0, r + costs[d], f_ncols); 
		mzd_mul(T, N, M, 0);
		mzd_free(T);
		mzd_free(M);

		M = _crt_modred_mat(g_ncols, poly, d);
		T = mzd_init_window(f->G, r, 0, r + costs[d], g_ncols); 
		mzd_mul(T, N, M, 0);
		mzd_free(T);
		mzd_free(M);

		mzd_free(N);

		r += costs[d];
	}

	/**
	 * 1.2) We deal with regular polynomial which are co-prime
	 */

	for(deg_t d=1; d<M4RIE_CRT_LEN; d++) {
		if (p[d] == 0)
			continue;

		mzd_t *N  = _small_multiplication_map(d);

		for(int i=0; i<p[d]; i++) {
			if (p_it[d] < irreducible_polynomials[d][0]) {
				poly = irreducible_polynomials[d][ 1 + p_it[d] ];
				p_it[d]++;
			} else if (d/2 && p_it[d/2] < irreducible_polynomials[d/2][0]) {
				/** the minimal polynomial is a square */
				poly = irreducible_polynomials[d/2][ 1 + p_it[d/2]];
				p_it[d/2]++;
				poly = gf2x_mul(poly, poly, d/2+1);
			} else if (d/4 && p_it[d/4] < irreducible_polynomials[d/4][0]) {
				/** the minimal polynomial is a fourth power */
				poly = irreducible_polynomials[d/4][1 + p_it[d/4]];
				p_it[d/4]++;
				poly = gf2x_mul(poly, poly, d/4+1);
				poly = gf2x_mul(poly, poly, d/2+1);
			} else if (d/8 && p_it[d/8] < irreducible_polynomials[d/8][0]) {
				/** the minimal polynomial is an eigth power */
				poly = irreducible_polynomials[d/8][p_it[d/8]+ 1];
				p_it[d/8]++;
				poly = gf2x_mul(poly, poly, d/8+1);
				poly = gf2x_mul(poly, poly, d/4+1);
				poly = gf2x_mul(poly, poly, d/2+1);
			} else {
				m4ri_die("Degree %d is not implemented\n", d);
			}
			M = _crt_modred_mat(f_ncols, poly, d);
			T = mzd_init_window(f->F, r, 0, r + costs[d], f_ncols); 
			mzd_mul(T, N, M, 0);
			mzd_free(T);
			mzd_free(M);

			M = _crt_modred_mat(g_ncols, poly, d);
			T = mzd_init_window(f->G, r, 0, r + costs[d], g_ncols); 
			mzd_mul(T, N, M, 0);
			mzd_free(T);
			mzd_free(M);

			r += costs[d];
		}

		mzd_free(N);
	}

	m4ri_mm_free(p_it);

	/**
	 * 2) We solve for H as we know poly(c) and (F*vec(a) x G*vec(b)). We pick points poly(a) = x^v,
	 * poly(b) = x^w (hence: poly(c) = x^(v+w)).
	 */
	_blm_finish_polymult(ff, f);
	f->h = NULL;

	/**
	 * 3) We compile DJB maps if asked for.
	 */

	if (djb)
		_blm_djb_compile(f);

	return f;
}


void blm_free(blm_t *f) {
	mzd_free(f->F);
	mzd_free(f->G);
	mzd_free(f->H);
	if (f->f != f->g)
		djb_free(f->g);
	djb_free(f->f);
	djb_free(f->h);

	m4ri_mm_free(f);
}
