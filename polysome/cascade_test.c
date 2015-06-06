#include <cascade_cipher/cascade.h>

int
main(int argc, char *argv[])
{
	struct cascade_cipher cc;
	const uint8_t key[16] = { 0 };
	const uint8_t *input = "1234567890";
	uint8_t output[10] = { 0 };
	uint8_t final[10] = { 0 };

	if (cascade_cipher_init(&cc, key, 16) != - 1) {
		cascade_cipher_encrypt(&cc, input, 10, output, 10);
		cascade_cipher_decrypt(&cc, output, 10, final, 10);
		cascade_cipher_dispose(&cc);
	}

	return (0);
}
