#ifndef VANE_QUAD_H
#define VANE_QUAD_H

struct quad_mp {
    int ***P;
    int ***Q;
    int *state;
    int *keystream;
    int nvar;
    int klength;
};

typedef struct quad_mp quad_mp_t;

unsigned long long rdtsc(void);

void init_mp(quad_mp_t *, int);
void generate_mpcyclic(quad_mp_t *);
int *evaluate_cyclic(quad_mp_t *);
void generate_keystream_cyclic(quad_mp_t *, int);
void free_quad_mp(quad_mp_t *);
#endif
