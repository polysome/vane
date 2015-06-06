
#include <string.h>
#include <skein.h>
#include <threefishApi.h>


/*****************************  Skein_256 ******************************/
inline
void Skein_256_Process_Block(Skein_256_Ctxt_t *ctx, const uint8_t *blkPtr,
                             size_t blkCnt, size_t byteCntAdd)
{
    ThreefishKey_t key;
    uint64_t tweak[2];
    int i;
    uint64_t  w[SKEIN_256_STATE_WORDS];           /* local copy of input block */
    uint64_t words[3];

    Skein_assert(blkCnt != 0);                  /* never call with blkCnt == 0! */
    tweak[0] = ctx->h.T[0];
    tweak[1] = ctx->h.T[1];

    do  {
        uint64_t carry = byteCntAdd;

        words[0] = tweak[0] & 0xffffffffL;
        words[1] = ((tweak[0] >> 32) & 0xffffffffL);
        words[2] = (tweak[1] & 0xffffffffL);

        for (i = 0; i < 3; i++) {
            carry += words[i];
            words[i] = carry;
            carry >>= 32;
        }        
        tweak[0] = words[0] & 0xffffffffL;
        tweak[0] |= (words[1] & 0xffffffffL) << 32;
        tweak[1] |= words[2] & 0xffffffffL;

        threefishSetKey(&key, Threefish256, ctx->X, tweak);

        Skein_Get64_LSB_First(w, blkPtr, SKEIN_256_STATE_WORDS);   /* get input block in little-endian format */

        threefishEncryptBlockWords(&key, w, ctx->X);

        blkPtr += SKEIN_256_BLOCK_BYTES;

        /* do the final "feedforward" xor, update context chaining vars */
        ctx->X[0] = ctx->X[0] ^ w[0];
        ctx->X[1] = ctx->X[1] ^ w[1];
        ctx->X[2] = ctx->X[2] ^ w[2];
        ctx->X[3] = ctx->X[3] ^ w[3];

        tweak[1] &= ~SKEIN_T1_FLAG_FIRST;
    } while (--blkCnt);

    ctx->h.T[0] = tweak[0];
    ctx->h.T[1] = tweak[1];
}

inline
void Skein_512_Process_Block(Skein_512_Ctxt_t *ctx, const uint8_t *blkPtr,
                             size_t blkCnt, size_t byteCntAdd)
{
    ThreefishKey_t key;
    uint64_t tweak[2];
    int i;
    uint64_t words[3];
    uint64_t  w[SKEIN_512_STATE_WORDS];           /* local copy of input block */

    Skein_assert(blkCnt != 0);                  /* never call with blkCnt == 0! */
    tweak[0] = ctx->h.T[0];
    tweak[1] = ctx->h.T[1];

    do  {
        uint64_t carry = byteCntAdd;

        words[0] = tweak[0] & 0xffffffffL;
        words[1] = ((tweak[0] >> 32) & 0xffffffffL);
        words[2] = (tweak[1] & 0xffffffffL);

        for (i = 0; i < 3; i++) {
            carry += words[i];
            words[i] = carry;
            carry >>= 32;
        }        
        tweak[0] = words[0] & 0xffffffffL;
        tweak[0] |= (words[1] & 0xffffffffL) << 32;
        tweak[1] |= words[2] & 0xffffffffL;

        threefishSetKey(&key, Threefish512, ctx->X, tweak);

        Skein_Get64_LSB_First(w, blkPtr, SKEIN_512_STATE_WORDS);   /* get input block in little-endian format */

        threefishEncryptBlockWords(&key, w, ctx->X);

        blkPtr += SKEIN_512_BLOCK_BYTES;

        /* do the final "feedforward" xor, update context chaining vars */
        ctx->X[0] = ctx->X[0] ^ w[0];
        ctx->X[1] = ctx->X[1] ^ w[1];
        ctx->X[2] = ctx->X[2] ^ w[2];
        ctx->X[3] = ctx->X[3] ^ w[3];
        ctx->X[4] = ctx->X[4] ^ w[4];
        ctx->X[5] = ctx->X[5] ^ w[5];
        ctx->X[6] = ctx->X[6] ^ w[6];
        ctx->X[7] = ctx->X[7] ^ w[7];

        tweak[1] &= ~SKEIN_T1_FLAG_FIRST;
    } while (--blkCnt);

    ctx->h.T[0] = tweak[0];
    ctx->h.T[1] = tweak[1];
}

inline
void Skein1024_Process_Block(Skein1024_Ctxt_t *ctx, const uint8_t *blkPtr,
                              size_t blkCnt, size_t byteCntAdd)
{
    ThreefishKey_t key;
    uint64_t tweak[2];
    int i;
    uint64_t words[3];
    uint64_t  w[SKEIN1024_STATE_WORDS];           /* local copy of input block */

    Skein_assert(blkCnt != 0);                  /* never call with blkCnt == 0! */
    tweak[0] = ctx->h.T[0];
    tweak[1] = ctx->h.T[1];

    do  {
        uint64_t carry = byteCntAdd;

        words[0] = tweak[0] & 0xffffffffL;
        words[1] = ((tweak[0] >> 32) & 0xffffffffL);
        words[2] = (tweak[1] & 0xffffffffL);

        for (i = 0; i < 3; i++) {
            carry += words[i];
            words[i] = carry;
            carry >>= 32;
        }        
        tweak[0] = words[0] & 0xffffffffL;
        tweak[0] |= (words[1] & 0xffffffffL) << 32;
        tweak[1] |= words[2] & 0xffffffffL;

        threefishSetKey(&key, Threefish1024, ctx->X, tweak);

        Skein_Get64_LSB_First(w, blkPtr, SKEIN1024_STATE_WORDS);   /* get input block in little-endian format */

        threefishEncryptBlockWords(&key, w, ctx->X);

        blkPtr += SKEIN1024_BLOCK_BYTES;

        /* do the final "feedforward" xor, update context chaining vars */
        ctx->X[ 0] = ctx->X[ 0] ^ w[ 0];
        ctx->X[ 1] = ctx->X[ 1] ^ w[ 1];
        ctx->X[ 2] = ctx->X[ 2] ^ w[ 2];
        ctx->X[ 3] = ctx->X[ 3] ^ w[ 3];
        ctx->X[ 4] = ctx->X[ 4] ^ w[ 4];
        ctx->X[ 5] = ctx->X[ 5] ^ w[ 5];
        ctx->X[ 6] = ctx->X[ 6] ^ w[ 6];
        ctx->X[ 7] = ctx->X[ 7] ^ w[ 7];
        ctx->X[ 8] = ctx->X[ 8] ^ w[ 8];
        ctx->X[ 9] = ctx->X[ 9] ^ w[ 9];
        ctx->X[10] = ctx->X[10] ^ w[10];
        ctx->X[11] = ctx->X[11] ^ w[11];
        ctx->X[12] = ctx->X[12] ^ w[12];
        ctx->X[13] = ctx->X[13] ^ w[13];
        ctx->X[14] = ctx->X[14] ^ w[14];
        ctx->X[15] = ctx->X[15] ^ w[15];

        tweak[1] &= ~SKEIN_T1_FLAG_FIRST;
    } while (--blkCnt);

    ctx->h.T[0] = tweak[0];
    ctx->h.T[1] = tweak[1];
}
