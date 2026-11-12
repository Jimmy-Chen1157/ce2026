#include "montgomery.h"
#include "fe25519.h"

/* Constant-time swap of field elements. b must be 0 or 1. */
static void fe_cswap(fe25519 *x, fe25519 *y, unsigned char b)
{
  uint64_t mask = UINT64_C(0) - (uint64_t)(b & 1U);
  unsigned i;

  for (i = 0; i < 5; ++i) {
    uint64_t t = mask & (x->v[i] ^ y->v[i]);
    x->v[i] ^= t;
    y->v[i] ^= t;
  }
}

void montgomery_scalarmult(unsigned char out[32],
                           const unsigned char scalar[32],
                           const unsigned char u[32])
{
  static const fe25519 a24 = {{121665, 0, 0, 0, 0}};
  unsigned char k[32];
  fe25519 x1, x2, z2, x3, z3;
  fe25519 a, aa, b, bb, e, c, d, da, cb, t0, t1;
  unsigned char swap = 0;
  int pos;

  /* Copy before writing output so the API also permits in-place use. */
  for (pos = 0; pos < 32; ++pos)
    k[pos] = scalar[pos];
  k[0] &= 248;
  k[31] &= 127;
  k[31] |= 64;

  fe25519_unpack(&x1, u);
  x2 = fe25519_one;
  z2 = fe25519_zero;
  x3 = x1;
  z3 = fe25519_one;

  /* Montgomery ladder from RFC 7748, with a fixed iteration count. */
  for (pos = 254; pos >= 0; --pos) {
    unsigned char bit = (unsigned char)((k[pos >> 3] >> (pos & 7)) & 1U);
    swap ^= bit;
    fe_cswap(&x2, &x3, swap);
    fe_cswap(&z2, &z3, swap);
    swap = bit;

    fe25519_add(&a, &x2, &z2);
    fe25519_square(&aa, &a);
    fe25519_sub(&b, &x2, &z2);
    fe25519_square(&bb, &b);
    fe25519_sub(&e, &aa, &bb);
    fe25519_add(&c, &x3, &z3);
    fe25519_sub(&d, &x3, &z3);
    fe25519_mul(&da, &d, &a);
    fe25519_mul(&cb, &c, &b);

    fe25519_add(&t0, &da, &cb);
    fe25519_square(&x3, &t0);
    fe25519_sub(&t0, &da, &cb);
    fe25519_square(&t0, &t0);
    fe25519_mul(&z3, &x1, &t0);

    fe25519_mul(&x2, &aa, &bb);
    fe25519_mul(&t1, &a24, &e);
    fe25519_add(&t1, &aa, &t1);
    fe25519_mul(&z2, &e, &t1);
  }

  fe_cswap(&x2, &x3, swap);
  fe_cswap(&z2, &z3, swap);
  fe25519_invert(&z2, &z2);
  fe25519_mul(&x2, &x2, &z2);
  fe25519_pack(out, &x2);
}

void montgomery_scalarmult_base(unsigned char out[32],
                                const unsigned char scalar[32])
{
  static const unsigned char base[32] = {9};
  montgomery_scalarmult(out, scalar, base);
}
