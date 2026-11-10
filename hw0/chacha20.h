#ifndef CHACHA20_H
#define CHACHA20_H

#include <stdint.h>

#define CHACHA20_KEYBYTES 32
#define CHACHA20_NONCEBYTES 8

int crypto_stream_chacha20(unsigned char *c, unsigned long long clen,
                           const unsigned char *n, const unsigned char *k);

/*
 * The ChaCha20 quarter round (RFC 8439, Section 2.1). It updates the four
 * 32-bit words *a, *b, *c, *d in place:
 *
 *   *a += *b;  *d ^= *a;  *d <<<= 16;
 *   *c += *d;  *b ^= *c;  *b <<<= 12;
 *   *a += *b;  *d ^= *a;  *d <<<=  8;
 *   *c += *d;  *b ^= *c;  *b <<<=  7;
 *
 * where <<<= is a rotate left and additions are modulo 2^32.
 *
 * quarterround() is the one you implement, in Armv7-M assembly, in
 * quarterround.S. crypto_stream_chacha20() calls it, so ChaCha20 only produces
 * the correct keystream once your implementation is right.
 *
 * quarterround_ref() is the provided C reference. The tests compare your
 * assembly against it. Do not modify it.
 */
void quarterround(uint32_t *a, uint32_t *b, uint32_t *c, uint32_t *d);
void quarterround_ref(uint32_t *a, uint32_t *b, uint32_t *c, uint32_t *d);

#endif
