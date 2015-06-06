/**
 * Based on John Underhill's implementation
 *
 * Serpent's block cipher
 */


#ifndef	SERPENT_H
#define	SERPENT_h

#include <sys/types.h>
#include <sha2.h>

enum ctx_type {
	SHA512,
	CTX_MAX
};

typedef void (*ctx_init)(void *);
typedef void (*ctx_update)(void *, const void *, size_t);
typedef void (*ctx_final)(uint8_t *, void *);

struct shx {
	SHA2_CTX ctx;
	ctx_init ci;
	ctx_update cu;
	ctx_final cf;
	enum ctx_type ct;
	size_t digestlen;
	size_t blocklen;
	size_t expkeylen;
	int32_t *expkey;
	int rounds;
};

typedef struct shx shx_t;

int shx_init(struct shx *, int, enum ctx_type, const uint8_t *, size_t);
void shx_decrypt(struct shx *, const uint8_t *, size_t, uint8_t *, size_t);
void shx_encrypt(struct shx *, const uint8_t *, size_t, uint8_t *, size_t);
void shx_dispose(struct shx *);

#endif
