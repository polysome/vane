#include <skein.h>
#include <serpent/shx.h>

struct cascade_cipher {
	struct shx sctx;
	Skein1024_Ctxt_t tctx;
	const uint8_t *key;
	size_t keylen;
};

typedef struct cascade_cipher cascade_cipher_t;

int cascade_cipher_init(struct cascade_cipher *, const uint8_t *, size_t);
int cascade_cipher_encrypt(struct cascade_cipher *, const uint8_t *, size_t, uint8_t *, size_t);
int cascade_cipher_decrypt(struct cascade_cipher *, const uint8_t *, size_t, uint8_t *, size_t);
void cascade_cipher_dispose(struct cascade_cipher *);
