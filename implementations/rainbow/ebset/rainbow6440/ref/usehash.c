
#include "rainbow.h"

#if defined(__CRYPTO_EBATS__)
#include "crypto_sign.h"
#endif

#include "crypto_hash_sha256.h"
#include "sizes.h"




int crypto_sign_keypair( unsigned char *pk, unsigned char *sk )
{
	return genkey_pack(pk,sk);
}

int crypto_sign(
  unsigned char *sm,unsigned long long *smlen,
  const unsigned char *m,unsigned long long mlen,
  const unsigned char *sk
)
{
  unsigned char h[32];
  int i;

  crypto_hash_sha256(h,m,mlen);

  /* if (SHORTHASH_BYTES < 32) return -1; */
  i = sign_bin( sm , sk , h );
  *smlen = SIGNATURE_BYTES;
  if (i < 0) return i;

  if (*smlen != SIGNATURE_BYTES) return -1;
  for (i = 0;i < mlen;++i) {
    sm[*smlen] = m[i];
    ++*smlen;
  }
  return 0;
}

int crypto_sign_open(
  unsigned char *m,unsigned long long *mlen,
  const unsigned char *sm,unsigned long long smlen,
  const unsigned char *pk
)
{
  unsigned char h[32];
  int i;

  if (smlen < SIGNATURE_BYTES) return -100;
  crypto_hash_sha256(h,sm+SIGNATURE_BYTES,smlen-SIGNATURE_BYTES);
  /* if (SHORTHASH_BYTES < 32) return -1; */
  int r =  verify_bin( h , pk , sm );
  for (i = SIGNATURE_BYTES;i < smlen;++i) m[i - SIGNATURE_BYTES] = sm[i];
  *mlen = smlen - SIGNATURE_BYTES;
  return r;
}
