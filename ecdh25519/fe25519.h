#ifndef FE25519_H
#define FE25519_H

#include <stdint.h>

/* Five radix-2^51 limbs cover the 255-bit field. */
typedef struct
{
  uint64_t v[5];
}
fe25519;

extern const fe25519 fe25519_zero;
extern const fe25519 fe25519_one;
extern const fe25519 fe25519_two;

extern const fe25519 fe25519_sqrtm1;
extern const fe25519 fe25519_msqrtm1;
extern const fe25519 fe25519_m1;

/* Fixed-size Cortex-M4 helper for the 16-digit field multiplication. */
extern void _fe25519_mul_coeff(uint64_t *out, const uint64_t *a,
                              const uint64_t *b);
extern void _fe25519_add(fe25519 *r, const fe25519 *x, const fe25519 *y);
extern void _fe25519_carry(fe25519 *r);
extern void _fe25519_cmov(fe25519 *r, const fe25519 *x, unsigned b);

void fe25519_freeze(fe25519 *r);

void fe25519_unpack(fe25519 *r, const unsigned char x[32]);

void fe25519_pack(unsigned char r[32], const fe25519 *x);

int fe25519_iszero(const fe25519 *x);

int fe25519_isone(const fe25519 *x);

int fe25519_isnegative(const fe25519 *x);

int fe25519_iseq(const fe25519 *x, const fe25519 *y);

int fe25519_iseq_vartime(const fe25519 *x, const fe25519 *y);

void fe25519_cmov(fe25519 *r, const fe25519 *x, unsigned char b);

void fe25519_setone(fe25519 *r);

void fe25519_setzero(fe25519 *r);

void fe25519_neg(fe25519 *r, const fe25519 *x);

unsigned char fe25519_getparity(const fe25519 *x);

void fe25519_add(fe25519 *r, const fe25519 *x, const fe25519 *y);

void fe25519_double(fe25519 *r, const fe25519 *x);
void fe25519_triple(fe25519 *r, const fe25519 *x);

void fe25519_sub(fe25519 *r, const fe25519 *x, const fe25519 *y);

void fe25519_mul(fe25519 *r, const fe25519 *x, const fe25519 *y);

void fe25519_square(fe25519 *r, const fe25519 *x);

void fe25519_invert(fe25519 *r, const fe25519 *x);

void fe25519_pow2523(fe25519 *r, const fe25519 *x);

void fe25519_invsqrt(fe25519 *r, const fe25519 *x);

void fe25519_print(const fe25519 *x);

#endif
