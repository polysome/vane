#include <stdio.h>
#include <polysome/polysome.h>

void malloc_dump(int);

int
main(int argc, char *argv[])
{
	if (argc < 2)
		return (-1);

	const char *mod = argv[1];	
	int fd = -1;

	init_modules();
	polysome *p = get_module(mod);
	if (p == NULL) {
		printf("%s invalid\n", mod);
	} else {
		printf("%s valid\n", p->name);
		ovvar_t ov;
		memset(&ov, 0, sizeof(ov));

		init_ov(&ov);
		printf("path %s\n", ov.path);
		
		p->ss(&ov);
		printf("path %s be read\n", ov.fpath);
		printf("ov.a: %d, ov.b: %d\n", ov.a, ov.b);
		p->gs(&ov);

		for (rci_t i = 0; i < (ov.a + ov.b); i ++)
			printf("%ld ", ov.signature[i]);

		printf("\n");

		free_ov(&ov);
	}

	free_modules();

	return (0);
}
