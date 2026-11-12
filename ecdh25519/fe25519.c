#include <stdio.h>
#include "fe25519.h"

/* p = 2^255 - 19, represented in radix 2^51. */
#define FE_MASK ((UINT64_C(1) << 51) - 1)

const fe25519 fe25519_zero = {{0, 0, 0, 0, 0}};
const fe25519 fe25519_one = {{1, 0, 0, 0, 0}};
const fe25519 fe25519_two = {{2, 0, 0, 0, 0}};
const fe25519 fe25519_sqrtm1 = {{0x61b274a0ea0b0, 0xd5a5fc8f189d, 0x7ef5e9cbd0c60, 0x78595a6804c9e, 0x2b8324804fc1d}};
const fe25519 fe25519_msqrtm1 = {{0x1e4d8b5f15f3d, 0x72a5a0370e762, 0x10a16342f39f, 0x7a6a597fb361, 0x547cdb7fb03e2}};
const fe25519 fe25519_m1 = {{0x7ffffffffffec, 0x7ffffffffffff, 0x7ffffffffffff, 0x7ffffffffffff, 0x7ffffffffffff}};



void fe25519_cmov(fe25519 *r, const fe25519 *x, unsigned char b)
{
  _fe25519_cmov(r, x, b);
}

void fe25519_freeze(fe25519 *r)
{
  uint64_t q, mask;
  unsigned i;
  _fe25519_carry(r);
  /* q is one exactly when the normalized value is at least p. */
  q = (r->v[0] + 19) >> 51;
  for (i = 1; i < 4; ++i)
    q = (r->v[i] + q) >> 51;
  q = (r->v[4] + q) >> 51;
  mask = UINT64_C(0) - q;

  uint64_t sub = FE_MASK-18;
  r->v[0] -= mask & sub;
  sub = FE_MASK;
  for (i = 1; i < 5; ++i) {
    r->v[i] -= mask & sub;
  }
}

void fe25519_unpack(fe25519 *r, const unsigned char x[32])
{
  unsigned i;
  for (i = 0; i < 5; ++i)
    r->v[i] = 0;
  for (i = 0; i < 255; ++i)
    r->v[i / 51] |= (uint64_t)((x[i >> 3] >> (i & 7)) & 1) << (i % 51);
}

void fe25519_pack(unsigned char r[32], const fe25519 *x)
{
  fe25519 y = *x;
  unsigned i;
  fe25519_freeze(&y);
  for (i = 0; i < 32; ++i)
    r[i] = 0;
  for (i = 0; i < 255; ++i)
    r[i >> 3] |= (unsigned char)(((y.v[i / 51] >> (i % 51)) & 1) << (i & 7));
}

int fe25519_iseq(const fe25519 *x, const fe25519 *y)
{
  fe25519 a = *x;
  fe25519 b = *y;
  unsigned int d = 0;
  fe25519_freeze(&a);
  fe25519_freeze(&b);
  for (unsigned i = 0; i < 5; ++i)
    d |= (a.v[i] ^ b.v[i]);
  return d == 0;
}
int fe25519_iseq_vartime(const fe25519 *x, const fe25519 *y)
{
  return fe25519_iseq(x, y);
}

int fe25519_iszero(const fe25519 *x)
{
  return fe25519_iseq(x, &fe25519_zero);
}

int fe25519_isone(const fe25519 *x)
{
  return fe25519_iseq(x, &fe25519_one);
}

int fe25519_isnegative(const fe25519 *x)
{
  fe25519 t = *x;
  fe25519_freeze(&t);
  return (int)(t.v[0] & 1);
}

unsigned char fe25519_getparity(const fe25519 *x)
{
  return (unsigned char)fe25519_isnegative(x);
}

void fe25519_setone(fe25519 *r)
{
  *r = fe25519_one;
}

void fe25519_setzero(fe25519 *r)
{
  *r = fe25519_zero;
}

void fe25519_add(fe25519 *r, const fe25519 *x, const fe25519 *y)
{
  _fe25519_add(r, x, y);
  _fe25519_carry(r);
}

void fe25519_double(fe25519 *r, const fe25519 *x)
{
  fe25519_add(r, x, x);
}

void fe25519_triple(fe25519 *r, const fe25519 *x)
{
  fe25519_double(r, x);
  fe25519_add(r, r, x);
}

void fe25519_sub(fe25519 *r, const fe25519 *x, const fe25519 *y)
{
  unsigned i;
  r->v[0] = x->v[0] + 2 * (FE_MASK - 18) - y->v[0];
  for (i = 1; i < 5; ++i)
    r->v[i] = x->v[i] + 2 * FE_MASK - y->v[i];
  _fe25519_carry(r);
}

void fe25519_neg(fe25519 *r, const fe25519 *x)
{
  fe25519_sub(r, &fe25519_zero, x);
}

void fe25519_mul(fe25519 *r, const fe25519 *x, const fe25519 *y)
{
  /*
   * Cortex-M4 compilers do not provide a native 128-bit integer type. Split
   * the radix-2^51 limbs into alternating 26- and 25-bit limbs so every
   * product and convolution sum fits in uint64_t.
   */
  uint64_t a[10], b[10], h[10] = {0};
  unsigned i, j;
  for (i = 0; i < 5; ++i) {
    a[2 * i] = x->v[i] & ((UINT64_C(1) << 26) - 1);
    a[2 * i + 1] = x->v[i] >> 26;
    b[2 * i] = y->v[i] & ((UINT64_C(1) << 26) - 1);
    b[2 * i + 1] = y->v[i] >> 26;
  }
  for (i = 0; i < 10; ++i)
    for (j = 0; j < 10; ++j) {
      unsigned k = i + j;
      uint64_t odd_pair_mask = UINT64_C(0) - (uint64_t)(i & j & 1U);
      /* Odd limbs start one bit later in the alternating radix. */
      uint64_t k_mask = UINT64_C(0) - (uint64_t)((k - 10) >> 31);
      uint64_t fold = UINT64_C(19) ^ (odd_pair_mask & (UINT64_C(38) ^ UINT64_C(19)));
      uint64_t scale = UINT64_C(1) ^ (odd_pair_mask & (UINT64_C(2) ^ UINT64_C(1)));
      uint64_t mult = fold ^ (k_mask & (fold ^ scale));
      signed index = (k-10) ^ (k_mask & ((k-10) ^ k));
      h[index] += (mult * a[i] * b[j]);
    }

  /* carry in radix 2^26, 2^25, ...; 2^255 folds back as 19. */
  for (unsigned round = 0; round < 3; ++round) {
    for (i = 0; i < 9; ++i) {
      unsigned odd_mask = 0U - (i & 1U);
      unsigned bits = 26U ^ (odd_mask & (26U ^ 25U));
      uint64_t c = h[i] >> bits;
      h[i] &= (UINT64_C(1) << bits) - 1;
      h[i + 1] += c;
    }
    h[0] += (h[9] >> 25) * 19;
    h[9] &= (UINT64_C(1) << 25) - 1;
  }
  r->v[0] = h[0] | (h[1] << 26);
  r->v[1] = h[2] | (h[3] << 26);
  r->v[2] = h[4] | (h[5] << 26);
  r->v[3] = h[6] | (h[7] << 26);
  r->v[4] = h[8] | (h[9] << 26);
}
void fe25519_square(fe25519 *r, const fe25519 *x)
{
  fe25519_mul(r, x, x);
}

void fe25519_pow2523(fe25519 *r, const fe25519 *x)
{
  fe25519 z2, z9, z11, z2_5_0, z2_10_0, z2_20_0, z2_50_0, z2_100_0, t;
  int i;

  fe25519_square(&z2, x);
  fe25519_square(&t, &z2);
  fe25519_square(&t, &t);
  fe25519_mul(&z9, &t, x);
  fe25519_mul(&z11, &z9, &z2);
  fe25519_square(&t, &z11);
  fe25519_mul(&z2_5_0, &t, &z9);
  fe25519_square(&t, &z2_5_0);
  for (i = 1; i < 5; ++i)
    fe25519_square(&t, &t);
  fe25519_mul(&z2_10_0, &t, &z2_5_0);
  fe25519_square(&t, &z2_10_0);
  for (i = 1; i < 10; ++i)
    fe25519_square(&t, &t);
  fe25519_mul(&z2_20_0, &t, &z2_10_0);
  fe25519_square(&t, &z2_20_0);
  for (i = 1; i < 20; ++i)
    fe25519_square(&t, &t);
  fe25519_mul(&t, &t, &z2_20_0);
  fe25519_square(&t, &t);
  for (i = 1; i < 10; ++i)
    fe25519_square(&t, &t);
  fe25519_mul(&z2_50_0, &t, &z2_10_0);
  fe25519_square(&t, &z2_50_0);
  for (i = 1; i < 50; ++i)
    fe25519_square(&t, &t);
  fe25519_mul(&z2_100_0, &t, &z2_50_0);
  fe25519_square(&t, &z2_100_0);
  for (i = 1; i < 100; ++i)
    fe25519_square(&t, &t);
  fe25519_mul(&t, &t, &z2_100_0);
  fe25519_square(&t, &t);
  for (i = 1; i < 50; ++i)
    fe25519_square(&t, &t);
  fe25519_mul(&t, &t, &z2_50_0);
  fe25519_square(&t, &t);
  fe25519_square(&t, &t);
  fe25519_mul(r, &t, x);
}

void fe25519_invert(fe25519 *r, const fe25519 *x)
{
  fe25519 t, x2, x3;
  fe25519_pow2523(&t, x);
  fe25519_square(&t, &t);
  fe25519_square(&t, &t);
  fe25519_square(&t, &t);
  fe25519_square(&x2, x);
  fe25519_mul(&x3, &x2, x);
  fe25519_mul(r, &t, &x3);
}

void fe25519_invsqrt(fe25519 *r, const fe25519 *x)
{
  fe25519 den2, den3, den4, den6, chk, t, t2;
  int b;
  fe25519_square(&den2, x);
  fe25519_mul(&den3, &den2, x);
  fe25519_square(&den4, &den2);
  fe25519_mul(&den6, &den2, &den4);
  fe25519_mul(&t, &den6, x);
  fe25519_pow2523(&t, &t);
  fe25519_mul(&t, &t, &den3);
  fe25519_square(&chk, &t);
  fe25519_mul(&chk, &chk, x);
  fe25519_mul(&t2, &t, &fe25519_sqrtm1);
  b = 1 - fe25519_isone(&chk);
  fe25519_cmov(&t, &t2, (unsigned char)b);
  *r = t;
}

void fe25519_print(const fe25519 *x)
{
  unsigned char b[32];
  int i;
  fe25519_pack(b, x);
  for (i = 0; i < 32; ++i)
    printf("%02x", b[i]);
}
