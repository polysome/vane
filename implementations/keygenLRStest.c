#include "keygen.h"

int main(int argc, char *argv[])
{
	generate_key();
	free_all_publickey();
	free_all_privatekey();

	return (0);
}
