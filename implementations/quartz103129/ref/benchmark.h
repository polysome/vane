/*
 * $Id: benchmark.h 1271 2008-06-08 08:06:13Z owenhsin $
 */

#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>

/* Copied from http://en.wikipedia.org/wiki/RDTSC */
#ifdef __cplusplus
extern "C"
{
#endif
static inline uint64_t rdtsc() {
	uint32_t lo, hi;
	/* We cannot use "=A", since this would use %rax on x86_64 */
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	return (uint64_t)hi << 32 | lo;
}
#ifdef __cplusplus
}
#endif

#define RECMAX 6

#define BENCHMARK(bm,call) do { \
		bm_start(&(bm)); \
		call; \
		bm_stop(&(bm)); \
	} while (0)

struct benchmark {
	uint64_t start;
	uint64_t stop;
	double record[RECMAX];
	double acc;
	int currec;
	int count;
};

static inline void
bm_init(struct benchmark *bm)
{
	memset(bm, 0, sizeof(*bm));
}

static inline void
bm_start(struct benchmark *bm)
{
	bm->start = rdtsc();
}

static inline void
bm_stop(struct benchmark *bm)
{
	bm->stop = rdtsc();
	bm->record[bm->currec] = bm->stop - bm->start;
	bm->acc += bm->record[bm->currec];
	bm->currec = (bm->currec + 1) % RECMAX;
	++bm->count;
}

static inline void
bm_dump(char *buf, size_t bufsize, const struct benchmark *bm)
{
	int i;
	size_t len;

	len = snprintf(buf, bufsize, "%.0lf (%d):", bm->acc/bm->count, bm->count);
	buf += len;
	bufsize -= len;
	for (i = 0; i < RECMAX; ++i) {
		len = snprintf(buf, bufsize, " %.0lf", bm->record[i]);
		buf += len;
		bufsize -= len;
	}
}

#endif /* BENCHMARK_H */
