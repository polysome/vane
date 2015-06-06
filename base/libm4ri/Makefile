#	$OpenBSD: Makefile,v 1.17 2014/11/02 13:00:00 dcarlier Exp $
LIB=			m4ri
SHLIB_MAJOR=	1
SHLIB_MINOR=	0
SRCS=			brilliantrussian.c djb.c echelonform.c graycode.c io.c					\
				misc.c mmc.c mp.c mzd.c mzp.c ple.c ple_russian.c solve.c				\
				strassen.c triangular.c triangular_russian.c

HDRS=			m4ri/brilliantrussian.h m4ri/debug_dump.h m4ri/djb.h m4ri/echelonform.h	\
				m4ri/graycode.h m4ri/io.h m4ri/m4ri.h m4ri/m4ri_config.h m4ri/misc.h	\
				m4ri/mmc.h m4ri/mp.h m4ri/mzd.h m4ri/mzp.h m4ri/parity.h m4ri/ple.h		\
				m4ri/ple_russian.h m4ri/ple_russian_template.h m4ri/solve.h 			\
				m4ri/strassen.h m4ri/triangular.h m4ri/triangular_russian.h m4ri/xor.h	\
				m4ri/xor_template.h

CFLAGS+=		-I${.CURDIR} -I${.CURDIR}/m4ri -std=c99

DEBUGLIBS=		yes

.PATH:			${.CURDIR}/m4ri
includes:
	@cd ${.CURDIR}; for i in $(HDRS); do \
	    j="cmp -s $$i ${DESTDIR}/usr/include/$$i || \
		${INSTALL} ${INSTALL_COPY} -o ${BINOWN} -g ${BINGRP} -m 444 $$i \
		${DESTDIR}/usr/include/m4ri"; \
	    echo $$j; \
	    eval "$$j"; \
	done

beforeinstall:
		${INSTALL} -o ${BINOWN} -g ${BINGRP} -m 755 -d /usr/include/m4ri

.include <bsd.lib.mk>
