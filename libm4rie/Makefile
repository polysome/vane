LIB=			m4rie
SHLIB_MAJOR=	1
SHLIB_MINOR=	0
SRCS=			blm.c conversion.c conversion_cling16.c conversion_cling8.c conversion_slice16.c		\
				conversion_slice8.c echelonform.c gf2e.c karatsuba.c mzd_poly.c mzd_ptr.c mzd_slice.c	\
				mzed.c newton_john.c ple.c strassen.c trsm.c

HDRS=			m4rie/blm.h m4rie/conversion.h m4rie/echelonform.h m4rie/gf2e.h m4rie/gf2x.h			\
				m4rie/m4ri_functions.h m4rie/m4rie.h m4rie/mzd_poly.h m4rie/mzd_ptr.h m4rie/mzd_slice.h	\
				m4rie/mzed.h m4rie/newton_john.h m4rie/permutation.h m4rie/ple.h m4rie/strassen.h		\
				m4rie/trsm.h

CFLAGS+=		-I${.CURDIR} -I${.CURDIR}/m4rie -I${.CURDIR}/../libm4ri -std=c99
.PATH:			${.CURDIR}/m4rie

DEBUGLIBS=		yes

includes:
	@cd ${.CURDIR}; for i in $(HDRS); do \
	    j="cmp -s $$i ${DESTDIR}/usr/include/$$i || \
		${INSTALL} ${INSTALL_COPY} -o ${BINOWN} -g ${BINGRP} -m 444 $$i \
		${DESTDIR}/usr/include/m4rie"; \
	    echo $$j; \
	    eval "$$j"; \
	done

beforeinstall:
		${INSTALL} -o ${BINOWN} -g ${BINGRP} -m 755 -d /usr/include/m4rie


.include <bsd.lib.mk>
