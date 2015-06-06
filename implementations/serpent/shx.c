#include "shx.h"

#include <err.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#define	BTOI32(src, off)		(((src[0 + off]) << 24) |		\
								((src[1 + off]) << 16) |		\
								((src[2 + off]) << 8) |			\
								(src[3 + off])) 

#define	I32TOB(src, dst, off)	dst[off + 3] = src;				\
								dst[off + 2] = src >> 8;		\
								dst[off + 1] = src >> 16;		\
								dst[off] = src >> 24

#define	TORIGHT(x, b)			((int32_t)((uint32_t)x >> b) | (x << (32 - b)))
#define	TOLEFT(x, b)			((x << b) | (int32_t)((uint32_t)x >>  (32 - b)))

static int
xmalloc(void **p, size_t s) {
    assert(p != NULL);
    if (s > 0 && errno == 0) {
        *p = malloc(s);
        if (*p == NULL) {
            errno = ENOMEM;
            return (-1);
        }
        return (0);
    }

    fprintf(stderr, "malloc failed (size %ld / errno %d)\n", s, errno);

    return (-1);
}

struct ctx_data {
	ctx_init ci;
	ctx_update cu;
	ctx_final cf;
	size_t digestlen;
	size_t blocklen;
};

static const struct ctx_data cd[] = {
	{	.ci = (ctx_init)SHA512Init,
		.cu = (ctx_update)SHA512Update,
		.cf = (ctx_final)SHA512Final,
		.digestlen = SHA512_DIGEST_LENGTH,
		.blocklen = SHA512_BLOCK_LENGTH   },
	{	0, 0, 0, 0, 0   }
};

int shx_ctx_init(struct shx *, const uint8_t *, size_t);
int shx_expand_key(struct shx *, const uint8_t *, size_t);

int
shx_init(struct shx *s, int rounds, enum ctx_type ct,
		const uint8_t *key, size_t keylen)
{
    if (ct >= CTX_MAX)
        return (-1);

    errno = 0;

    s->ci = cd[ct].ci;
    s->cu = cd[ct].cu;
    s->cf = cd[ct].cf;
    s->digestlen = cd[ct].digestlen;
    s->blocklen = cd[ct].blocklen;
    s->ct = ct;
    s->rounds = rounds;

    return (shx_ctx_init(s, key, keylen));
}

int
shx_ctx_init(struct shx *s, const uint8_t *key, size_t keylen)
{
    s->ci(&s->ctx);
    return (shx_expand_key(s, key, keylen));
}

static inline void
shx_sb0(int32_t *r)
{
    int32_t t[4] = { 0 };
    t[0] = r[0] ^ r[3];
    t[1] = r[2] ^ t[0];
    t[2] = r[1] ^ t[1];
    r[3] = (r[0] & r[3]) ^ t[2];
    t[3] = r[0] ^ (r[1] & t[0]);
    r[2] = t[2] ^ (r[2] | t[3]);
    r[0] = r[3] & (t[1] ^ t[3]);
    r[1] = ~t[1] & r[0];
    r[0] ^= ~t[3];
}

static inline void
shx_ib0(int32_t *r)
{
    int32_t t[5] = { 0 };
    t[0] = ~r[0];
    t[1] = r[0] ^ r[1];
    t[2] = r[3] ^ (t[0] | t[1]);
    t[3] = r[2] ^ t[2];
    r[2] = t[1] ^ t[3];
    t[4] = t[0] ^ (r[3] & t[1]);
    r[1] = t[2] ^ (r[2] & t[4]);
    r[3] = (r[0] & t[2]) ^ (t[3] | r[1]);
    r[0] = r[3] ^ (t[3] ^ t[4]);
}

static inline void
shx_sb1(int32_t *r)
{
    int32_t t[6] = { 0 };
    t[0] = r[1] ^ ~r[0];
    t[1] = r[2] ^ (r[0] | t[0]);
    r[2] = r[3] ^ t[1];
    t[2] = r[1] ^ (r[3] | t[0]);
    t[3] = t[0] ^ r[2];
    r[3] = t[3] ^ (t[1] & t[2]);
    t[4] = t[1] ^ t[2];
    r[1] = r[3] ^ t[4];
    r[0] = t[1] ^ (t[3] & t[4]);
}

static inline void 
shx_ib1(int32_t *r)
{
    int32_t t[6] = { 0 };
    t[0] = r[0] ^ r[2];
    t[1] = r[0] ^ (r[1] & t[0]);
    t[2] = t[0] ^ t[1];
    r[3] = r[2] ^ t[2];
    t[3] = r[1] ^ (t[0] & t[1]);
    r[1] = t[1] ^ (r[3] | t[3]);
    t[4] = ~r[1];
    t[5] = r[3] ^ t[3];
    r[0] = t[4] ^ t[5];
    r[2] = t[2] ^ (t[4] | t[5]);
}

static inline void
shx_sb2(int32_t *r)
{
    int32_t t[6] = { 0 };
    t[0] = ~r[0];
    t[1] = r[1] ^ r[2];
    t[2] = t[1] ^ (r[2] & t[0]);
    t[3] = r[2] ^ t[0];
    t[4] = r[1] & (r[2] & r[3]);
    t[5] = t[3] ^ t[4];
    r[2] = r[0] ^ ((r[3] | t[4]) & (t[2] | t[3]));
    r[1] = (t[1] ^ t[5]) ^ (r[2] ^ (r[3] | t[0]));
    r[0] = t[2];
    r[3] = t[5];
}

static inline void
shx_ib2(int32_t *r)
{
    int32_t t[7] = { 0 };
    t[0] = r[1] ^ r[3];
    t[1] = r[0] ^ r[2];
    t[2] = r[2] ^ t[1];
    t[3] = r[0] | ~t[0];
    r[0] = t[1] ^ (r[1] & t[2]);
    t[4] = t[0] ^ (t[1] | (r[3] ^ t[3]));
    t[5] = ~t[2];
    t[6] = r[0] | t[4];
    r[1] = t[5] ^ t[6];
    r[2] = (r[3] & t[5]) ^ (t[1] ^ t[6]);
    r[3] = t[4];
}

static inline void
shx_sb3(int32_t *r)
{
    int32_t t[6] = { 0 };
    t[0] = r[0] ^ r[1];
    t[1] = r[0] | r[3];
    t[2] = r[2] ^ r[3];
    t[3] = (r[0] & r[2]) | (t[0] & t[1]);
    r[2] = t[2] ^ t[3];
    t[4] = t[3] ^ (r[1] ^ t[1]);
    r[0] = t[0] ^ (t[2] & t[4]);
    t[5] = r[2] & r[0];
    r[3] = (r[1] | r[3]) ^ (t[2] ^ t[5]);
    r[1] = t[4] ^ t[5];
}

static inline void
shx_ib3(int32_t *r)
{
    int32_t t[5] = { 0 };
    t[0] = r[1] ^ r[2];
    t[1] = r[0] ^ (r[1] & t[0]);
    t[2] = r[3] | t[1];
    t[3] = r[3] ^ (t[0] | t[2]);
    r[2] = (r[2] ^ t[1]) ^ t[3];
    t[4] = (r[0] | r[1]) ^ t[3];
    r[0] = t[0] ^ t[2];
    r[3] = t[1] ^ (r[0] & t[4]);
    r[1] = r[3] ^ (r[0] ^ t[4]);
}

static inline void
shx_sb4(int32_t *r)
{
    int32_t t[7] = { 0 };
    t[0] = r[0] ^ r[3];
    t[1] = r[2] ^ (r[3] & t[0]);
    t[2] = r[1] | t[1];
    r[3] = t[0] ^ t[2];
    t[3] = ~r[1];
    t[4] = t[1] ^ (t[0] | t[3]);
    t[5] = t[0] ^ t[3];
    t[6] = (r[0] & t[4]) ^ (t[3] & t[5]);
    r[1] = (r[0] ^ t[1]) ^ (t[5] & t[6]);
    r[0] = t[4];
    r[2] = t[6];
}

static inline void
shx_ib4(int32_t *r)
{
    int32_t t[6] = { 0 };
    t[0] = r[1] ^ (r[0] & (r[2] | r[3]));
    t[1] = r[2] ^ (r[0] & t[0]);	
    t[2] = r[3] ^ t[1];
    t[3] = ~r[0];
    t[4] = t[0] ^ (t[1] & t[2]);
    t[5] = r[3] ^ (t[2] | t[3]);
    r[1] = t[2];
    r[0] = t[4] ^ t[5];
    r[2] = (t[0] & t[5]) ^ (t[2] ^ t[3]);
    r[3] = t[4];
}

static inline void
shx_sb5(int32_t *r)
{
    int32_t t[7] = { 0 };
    t[0] = ~r[0];
    t[1] = r[0] ^ r[1];
    t[2] = r[0] ^ r[3];
    t[3] = (r[2] ^ t[0]) ^ (t[1] | t[2]);
    t[4] = r[3] & t[3];
    t[5] = t[4] ^ (t[1] ^ t[3]);
    t[6] = t[2] ^ (t[0] | t[3]);
    r[2] = (t[1] | t[4]) ^ t[6];
    r[3] = (r[1] ^ t[4]) ^ (t[5] & t[6]);
    r[0] = t[3];
    r[1] = t[5];
}

static inline void
shx_ib5(int32_t *r)
{
    int32_t t[7] = { 0 };
    t[0] = ~r[2];
    t[1] = r[3] ^ (r[1] & t[0]);
    t[2] = r[0] & t[1];
    t[3] = t[2] ^ (r[1] ^ t[0]);
    t[4] = r[1] | t[3];
    t[5] = t[1] ^ (r[0] & t[4]);
    t[6] = r[0] | r[3];
    r[2] = (r[1] & t[6]) ^ (t[2] | (r[0] ^ r[2]));
    r[0] = t[6] ^ (t[0] ^ t[4]);
    r[1] = t[5];
    r[3] = t[3];
}

static inline void
shx_sb6(int32_t *r)
{
    int32_t t[6] = { 0 };
    t[0] = r[0] ^ r[3];
    t[1] = r[1] ^ t[0];
    t[2] = r[2] ^ (~r[0] | t[0]);
    r[1] ^= t[2];
    t[3] = t[0] | r[1];
    t[4] = r[3] ^ (t[0] | r[1]);
    r[2] = t[1] ^ (t[2] & t[4]);
    t[5] = t[2] ^ t[4];
    r[0] = r[2] ^ t[5];
    r[3] = ~t[2] ^ (t[1] & t[5]);
}

static inline void
shx_ib6(int32_t *r)
{
    int32_t t[8] = { 0 };
    t[0] = ~r[0];
    t[1] = r[0] ^ r[1];
    t[2] = r[2] ^ t[1];
    t[3] = r[3] ^ (r[2] | t[0]);
    t[4] = t[2] ^ t[3];
    t[5] = t[1] ^ (t[2] & t[3]);
    t[6] = t[3] ^ (r[1] | t[5]);
    t[7] = r[1] | t[6];
    r[0] = t[5] ^ t[7];
    r[2] = (r[3] & t[0]) ^ (t[2] ^ t[7]);
    r[1] = t[4];
    r[3] = t[6];
}

static inline void
shx_sb7(int32_t *r)
{
    int32_t t[5] = { 0 };
    t[0] = r[1] ^ r[2];
    t[1] = r[3] & (r[0] | r[1]);
    t[2] = r[0] ^ t[1];
    r[1] ^= (t[2] & (r[3] | t[0]));
    t[3] = t[0] ^ (r[0] & t[2]);
    t[4] = t[2] ^ (t[1] | r[1]);
    r[2] = t[1] ^ (t[3] & t[4]);
    r[0] = ~t[4] ^ (t[3] & r[2]);
    r[3] = t[3];
}

static inline void
shx_ib7(int32_t *r)
{
    int32_t t[5] = { 0 };
    t[0] = r[2] | (r[0] & r[1]);
    t[1] = r[3] & (r[0] | r[1]);
    t[2] = t[0] ^ t[1];
    t[3] = r[1] ^ t[1];
    r[1] = r[0] ^ (t[3] | (t[2] ^ ~r[3]));
    t[4] = (r[2] ^ t[3]) ^ (r[3] | r[1]);
    r[2] = (t[0] ^ r[1]) ^ (t[4] ^ (r[0] & t[2]));
    r[0] = t[4];
    r[3] = t[2];
}

static inline void
shx_linear_transform(int32_t *r)
{
    int32_t x[4] = { 0 };
    x[0] = TOLEFT(r[0], 13);
    x[2] = TOLEFT(r[2], 3);
    x[1] = r[1] ^ x[0] ^ x[2];
    x[3] = r[3] ^ x[2] ^ x[0] << 3;

    r[1] = TOLEFT(x[1], 1);
    r[3] = TOLEFT(x[3], 7);
    r[0] = TOLEFT(x[0] ^ r[1] ^ r[3], 5);
    r[2] = TOLEFT(x[2] ^ r[3] ^ (r[1] << 7), 22);
}

static inline void
shx_invert_transform(int32_t *r)
{
    int32_t x[4] = { 0 };
    x[2] = TORIGHT(r[2], 22) ^ r[3] ^ (r[1] << 7);
    x[0] = TORIGHT(r[0], 5) ^ r[1] ^ r[3];
    x[3] = TORIGHT(r[3], 7);
    x[1] = TORIGHT(r[1], 1);

    r[3] = x[3] ^ x[2] ^ x[0] << 3;
    r[1] = x[1] ^ x[0] ^ x[2];
    r[2] = TORIGHT(x[2], 3);
    r[0] = TORIGHT(x[0], 13);
}


void
shx_decrypt(struct shx *s, const uint8_t *src, size_t srcoff,
        uint8_t *dst, size_t dstoff)
{
    int32_t r[4] = { 0 };
    size_t keyctr = s->expkeylen - 1;

    r[3] = s->expkey[keyctr--] ^ BTOI32(src, srcoff);
    r[2] = s->expkey[keyctr--] ^ BTOI32(src, srcoff + 4);
    r[1] = s->expkey[keyctr--] ^ BTOI32(src, srcoff + 8);
    r[0] = s->expkey[keyctr--] ^ BTOI32(src, srcoff + 12);

    while (keyctr > 4) {
        shx_ib7(r);
        r[3] ^= s->expkey[keyctr--];
        r[2] ^= s->expkey[keyctr--];
        r[1] ^= s->expkey[keyctr--];
        r[0] ^= s->expkey[keyctr--];
        shx_invert_transform(r);

        shx_ib6(r);
        r[3] ^= s->expkey[keyctr--];
        r[2] ^= s->expkey[keyctr--];
        r[1] ^= s->expkey[keyctr--];
        r[0] ^= s->expkey[keyctr--];
        shx_invert_transform(r);

        shx_ib5(r);
        r[3] ^= s->expkey[keyctr--];
        r[2] ^= s->expkey[keyctr--];
        r[1] ^= s->expkey[keyctr--];
        r[0] ^= s->expkey[keyctr--];
        shx_invert_transform(r);

        shx_ib4(r);
        r[3] ^= s->expkey[keyctr--];
        r[2] ^= s->expkey[keyctr--];
        r[1] ^= s->expkey[keyctr--];
        r[0] ^= s->expkey[keyctr--];
        shx_invert_transform(r);

        shx_ib3(r);
        r[3] ^= s->expkey[keyctr--];
        r[2] ^= s->expkey[keyctr--];
        r[1] ^= s->expkey[keyctr--];
        r[0] ^= s->expkey[keyctr--];
        shx_invert_transform(r);

        shx_ib2(r);
        r[3] ^= s->expkey[keyctr--];
        r[2] ^= s->expkey[keyctr--];
        r[1] ^= s->expkey[keyctr--];
        r[0] ^= s->expkey[keyctr--];
        shx_invert_transform(r);

        shx_ib1(r);
        r[3] ^= s->expkey[keyctr--];
        r[2] ^= s->expkey[keyctr--];
        r[1] ^= s->expkey[keyctr--];
        r[0] ^= s->expkey[keyctr--];
        shx_invert_transform(r);

        shx_ib0(r);

        if (keyctr > 4) {
            r[3] ^= s->expkey[keyctr--];
            r[2] ^= s->expkey[keyctr--];
            r[1] ^= s->expkey[keyctr--];
            r[0] ^= s->expkey[keyctr--];
            shx_invert_transform(r);
        }
    }

    keyctr--;
    I32TOB(r[3] ^ s->expkey[keyctr], dst, dstoff);
    keyctr--;
    I32TOB(r[2] ^ s->expkey[keyctr], dst, dstoff + 4);
    keyctr--;
    I32TOB(r[1] ^ s->expkey[keyctr], dst, dstoff + 8);
    I32TOB(r[0] ^ s->expkey[keyctr], dst, dstoff + 12);
}

void
shx_encrypt(struct shx *s, const uint8_t *src, size_t srcoff, 
        uint8_t *dst, size_t dstoff)
{
    int32_t r[4] = { 0 };
    size_t keyctr = 0;

    r[0] = BTOI32(src, srcoff + 12);
    r[1] = BTOI32(src, srcoff + 8);
    r[2] = BTOI32(src, srcoff + 4);
    r[3] = BTOI32(src, srcoff);

    while (keyctr < s->expkeylen - 4) {
        r[0] ^= s->expkey[keyctr++];
        r[1] ^= s->expkey[keyctr++];
        r[2] ^= s->expkey[keyctr++];
        r[3] ^= s->expkey[keyctr++];
        shx_sb0(r);
        shx_linear_transform(r);

        r[0] ^= s->expkey[keyctr++];
        r[1] ^= s->expkey[keyctr++];
        r[2] ^= s->expkey[keyctr++];
        r[3] ^= s->expkey[keyctr++];
        shx_sb1(r);
        shx_linear_transform(r);

        r[0] ^= s->expkey[keyctr++];
        r[1] ^= s->expkey[keyctr++];
        r[2] ^= s->expkey[keyctr++];
        r[3] ^= s->expkey[keyctr++];
        shx_sb2(r);
        shx_linear_transform(r);

        r[0] ^= s->expkey[keyctr++];
        r[1] ^= s->expkey[keyctr++];
        r[2] ^= s->expkey[keyctr++];
        r[3] ^= s->expkey[keyctr++];
        shx_sb3(r);
        shx_linear_transform(r);

        r[0] ^= s->expkey[keyctr++];
        r[1] ^= s->expkey[keyctr++];
        r[2] ^= s->expkey[keyctr++];
        r[3] ^= s->expkey[keyctr++];
        shx_sb4(r);
        shx_linear_transform(r);

        r[0] ^= s->expkey[keyctr++];
        r[1] ^= s->expkey[keyctr++];
        r[2] ^= s->expkey[keyctr++];
        r[3] ^= s->expkey[keyctr++];
        shx_sb5(r);
        shx_linear_transform(r);

        r[0] ^= s->expkey[keyctr++];
        r[1] ^= s->expkey[keyctr++];
        r[2] ^= s->expkey[keyctr++];
        r[3] ^= s->expkey[keyctr++];
        shx_sb6(r);
        shx_linear_transform(r);

        r[0] ^= s->expkey[keyctr++];
        r[1] ^= s->expkey[keyctr++];
        r[2] ^= s->expkey[keyctr++];
        r[3] ^= s->expkey[keyctr++];
        shx_sb7(r);

        if (keyctr < s->expkeylen - 4)
            shx_linear_transform(r);
    }

    keyctr++;
    I32TOB(s->expkey[keyctr] ^ r[0], dst, dstoff + 12);
    keyctr++;
    I32TOB(s->expkey[keyctr] ^ r[1], dst, dstoff + 8);
    keyctr++;
    I32TOB(s->expkey[keyctr] ^ r[2], dst, dstoff + 4);
    I32TOB(s->expkey[keyctr] ^ r[3], dst, dstoff);
}

void
shx_decrypt_block(struct shx *s, const uint8_t *src, uint8_t *dst)
{
    shx_decrypt(s, src, 0, dst, 0);
}

void
shx_encrypt_block(struct shx *s, const uint8_t *src, uint8_t *dst)
{
    shx_encrypt(s, src, 0, dst, 0);
}

int
shx_expand_key(struct shx *s, const uint8_t *key, size_t len)
{
    assert(len > s->digestlen);
    size_t keysize = 4 * (s->rounds + 1);
    size_t keybytes = keysize * 4;
    size_t saltsize = len - s->digestlen;
    uint8_t *rawkey;
    uint8_t *hdfkey;
    uint8_t *hdfsalt;

    s->expkeylen = 0;
    s->expkey = NULL;

    if (saltsize % s->blocklen != 0)
        saltsize = saltsize - (saltsize % s->blocklen);

    if (xmalloc((void **)&hdfsalt, saltsize * sizeof(*hdfsalt)) != 0) {
        fprintf(stderr, "hdfsalt failed (%ld).\n", saltsize);
        return (-1);
    }
    if (xmalloc((void **)&hdfkey, s->digestlen * sizeof(*hdfkey)) != 0) {
        fprintf(stderr, "hdfkey failed.\n");
        goto bad1;
    }
    if (xmalloc((void **)&rawkey, keybytes * sizeof(*rawkey)) != 0) {
        fprintf(stderr, "rawkey failed.\n");
        goto bad2;
    }

    memcpy(hdfkey, key, s->digestlen);	
    memcpy(hdfsalt, key + s->digestlen, saltsize);

    s->cu(&s->ctx, hdfkey, s->digestlen);
    s->cu(&s->ctx, hdfsalt, saltsize);
    s->cf(rawkey, &s->ctx);

    if (xmalloc((void **)&s->expkey, keysize * sizeof(*s->expkey)) != 0) {
        fprintf(stderr, "expkey failed.\n");
        goto bad3;
    }
    s->expkeylen = keysize;
    memcpy(s->expkey, rawkey, keysize);

    free(rawkey);
    free(hdfkey);
    free(hdfsalt);

    return (0);

bad3:
    free(rawkey);
bad2:
    free(hdfkey);
bad1:
    free(hdfsalt);

    return (-1);
}

void
shx_dispose(struct shx *s)
{
    if (s) {
        free(s->expkey);
        s = NULL;
    }
}
