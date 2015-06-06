#include <serpent/shx.h>

#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

int
main(int argc, char *argv[])
{
	struct shx s;
	unsigned long x[32];
	const uint8_t *key = "0000000000000000000000000000000000000000000000000000000000000000"
                        "0000000000000000000000000000000000000000000000000000000000000000"                      
                        "0000000000000000000000000000000000000000000000000000000000000000"                      
                        "0000000000000000000000000000000000000000000000000000000000000000";                     
	size_t i = 0;
	for (; i < 32; i ++)
		x[i] = arc4random() % (1 << 31);
	if (shx_init(&s, 128, SHA512, key, strlen(key)) != -1) {
       shx_encrypt(&s, (uint8_t *)x, 32, (uint8_t *)x, 32);
       shx_decrypt(&s, (uint8_t *)x, 32, (uint8_t *)x, 32);
       shx_dispose(&s);
   } else {
       fprintf(stderr, "serpent could not init\n");
   }
	return (0);
}
