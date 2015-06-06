#include <stdlib.h>
#include <string.h>

#include <cascade_cipher/cascade.h>

#define	MAX(a, b)	(a < b ? b : a)
#define	MIN(a, b)	(a < b ? a : b)

int
cascade_cipher_init(struct cascade_cipher *ctx, const uint8_t *key, size_t keylen)
{
	if (shx_init(&ctx->sctx, 128, SHA512, key, keylen) != 0)
		return (-1);
	if (Skein1024_Init(&ctx->tctx, 1024) != 0)
		return (-1);
	if (Skein1024_Update(&ctx->tctx, key, keylen) != 0)
		return (-1);

	ctx->key = key;
	ctx->keylen = keylen;

	return (0);
}

int
cascade_cipher_encrypt(struct cascade_cipher *ctx, const uint8_t *input,
		size_t inputlen, uint8_t *output, size_t outputlen)
{
	uint8_t *skoutput = malloc(SKEIN1024_STATE_WORDS * outputlen);
	size_t i;

	if (Skein1024_Final(&ctx->tctx, skoutput) != 0)
		return (-1);

	shx_encrypt(&ctx->sctx, input, inputlen, output, outputlen);	

	for (i = 0; i < outputlen; i ++)
		output[i] ^= skoutput[i * SKEIN1024_STATE_WORDS];
   free(skoutput);

	return (0);
}

int
cascade_cipher_decrypt(struct cascade_cipher *ctx, const uint8_t *input,
		size_t inputlen, uint8_t *output, size_t outputlen)
{
	uint8_t *skoutput = malloc(SKEIN1024_STATE_WORDS * outputlen);
	uint8_t *inputtmp;
	size_t i;

	if (Skein1024_Final(&ctx->tctx, skoutput) != 0)
		return (-1);

	inputtmp = malloc(outputlen * sizeof(*inputtmp));
	if (inputtmp == 0)
		return (-1);
	memcpy(inputtmp, input, inputlen);

	for (i = 0; i < outputlen; i ++)
		inputtmp[i] ^= skoutput[i * SKEIN1024_STATE_WORDS];

	shx_decrypt(&ctx->sctx, (const uint8_t *)inputtmp, inputlen, output, outputlen);	

	free(inputtmp);
   free(skoutput);
	return (0);
}

void
cascade_cipher_dispose(struct cascade_cipher *ctx)
{
	shx_dispose(&ctx->sctx);
}
