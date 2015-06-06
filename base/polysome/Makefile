LIB=			vane
SHLIB_MAJOR=		0
SHLIB_MINOR=		8

SRCS+=		   abc/abc.c
SRCS+=			serpent/shx.c

SRCS+=			cascade_cipher/cascade.c

SRCS+=			skein/threefishApi.c skein/threefish256Block.c
SRCS+=			skein/threefish512Block.c skein/threefish1024Block.c
SRCS+=			skein/skein_block.c skein/skein.c

SRCS+=			quad/quad256.c

CFLAGS+=		-I${.CURDIR}/. -I${.CURDIR}/abc -I${.CURDIR}/skein/include -I${.CURDIR}/quad
CFLAGS+=		-Wall -std=c99 -Wunused-variable
CLEANFILES+= 		${.CURDIR}/keygenLRStest
CLEANFILES+= 		${.CURDIR}/test
CLEANFILES+=		${.CURDIR}/shx_test
CLEANFILES+=		${.CURDIR}/cascade_test
CLEANFILES+=		${.CURDIR}/quad_test
CLEANFILES+=		${.CURDIR}/*.core

DEBUGLIBS=		yes

keygenLRStest:
			cc -g -O0 -shared -o libkeygenLRS.so basic.c keygenLRS.c -fpie -pie -I. -L. -lm4ri -lm4rie
			cc -g -O0 -o keygenLRStest keygenLRStest.c -std=c99 -fpie -pie -I. -L. -lm4ri -lm4rie

test:
			cc -g -O0 -o test test.c -std=c99 -fpie -pie -I. -L. -lvane -lm
shx_test:
			cc -g -O0 -o shx_test shx_test.c ${CFLAGS} -I. -L. -lvane_g -lm

cascade_test:
			cc -g -O0 -o cascade_test cascade_test.c ${CFLAGS} -I. -L. -lvane_g -lm

quad_test:
			cc -g -O0 -o quad_test quad_test.c ${CFLAGS} -I. -L. -lvane_g -lm

.include <bsd.lib.mk>

.include <bsd.subdir.mk>
