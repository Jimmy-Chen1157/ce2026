#include <stdio.h>  //XXX: DEBUG
#include <assert.h> //XXX: DEBUG
#include "group.h"

/*
 * Arithmetic on the twisted Edwards curve -x^2 + y^2 = 1 + dx^2y^2
 * with d = -(121665/121666) = 37095705934669439343138083508754565189542113879843219016388785533085940283555
 * Base point: (15112221349535400772501151409588531511454012693041857206046113283949847762202,46316835694926478169428394003475163141307993866256225615783033603165251855960);
 */

/* d */
static const fe25519 ge25519_ecd = {{UINT64_C(0x34dca135978a3), UINT64_C(0x1a8283b156ebd), UINT64_C(0x5e7a26001c029), UINT64_C(0x739c663a03cbb), UINT64_C(0x52036cee2b6ff)}};
/* 2*d */
static const fe25519 ge25519_ec2d = {{UINT64_C(0x69b9426b2f159), UINT64_C(0x35050762add7a), UINT64_C(0x3cf44c0038052), UINT64_C(0x6738cc7407977), UINT64_C(0x2406d9dc56dff)}};

static const fe25519 ge25519_magic = {{UINT64_C(0x702557fa2bf03), UINT64_C(0x514b7d1a82cc6), UINT64_C(0x7f89efd8b43a7), UINT64_C(0x1aef49ec23700), UINT64_C(0x79376fa30500)}};

const group_ge group_ge_neutral = {{{0}},
                                   {{1}},
                                   {{1}},
                                   {{0}}};

#define ge25519_p3 group_ge

typedef struct
{
  fe25519 x;
  fe25519 z;
  fe25519 y;
  fe25519 t;
} ge25519_p1p1;

typedef struct
{
  fe25519 x;
  fe25519 y;
  fe25519 z;
} ge25519_p2;

typedef struct
{
  fe25519 x;
  fe25519 y;
} ge25519_aff;


static void p1p1_to_p2(ge25519_p2 *r, const ge25519_p1p1 *p)
{
  fe25519_mul(&r->x, &p->x, &p->t);
  fe25519_mul(&r->y, &p->y, &p->z);
  fe25519_mul(&r->z, &p->z, &p->t);
}

static void p1p1_to_p3(ge25519_p3 *r, const ge25519_p1p1 *p)
{
  p1p1_to_p2((ge25519_p2 *)r, p);
  fe25519_mul(&r->t, &p->x, &p->y);
}

static void add_p1p1(ge25519_p1p1 *r, const ge25519_p3 *p, const ge25519_p3 *q)
{
  fe25519 a, b, c, d, t;

  fe25519_sub(&a, &p->y, &p->x); /* A = (Y1-X1)*(Y2-X2) */
  fe25519_sub(&t, &q->y, &q->x);
  fe25519_mul(&a, &a, &t);
  fe25519_add(&b, &p->x, &p->y); /* B = (Y1+X1)*(Y2+X2) */
  fe25519_add(&t, &q->x, &q->y);
  fe25519_mul(&b, &b, &t);
  fe25519_mul(&c, &p->t, &q->t); /* C = T1*k*T2 */
  fe25519_mul(&c, &c, &ge25519_ec2d);
  fe25519_mul(&d, &p->z, &q->z); /* D = Z1*2*Z2 */
  fe25519_add(&d, &d, &d);
  fe25519_sub(&r->x, &b, &a); /* E = B-A */
  fe25519_sub(&r->t, &d, &c); /* F = D-C */
  fe25519_add(&r->z, &d, &c); /* G = D+C */
  fe25519_add(&r->y, &b, &a); /* H = B+A */
}

/* See http://www.hyperelliptic.org/EFD/g1p/auto-twisted-extended-1.html#doubling-dbl-2008-hwcd */
static void dbl_p1p1(ge25519_p1p1 *r, const ge25519_p2 *p)
{
  fe25519 a,b,c,d;
  fe25519_square(&a, &p->x);
  fe25519_square(&b, &p->y);
  fe25519_square(&c, &p->z);
  fe25519_add(&c, &c, &c);
  fe25519_neg(&d, &a);

  fe25519_add(&r->x, &p->x, &p->y);
  fe25519_square(&r->x, &r->x);
  fe25519_sub(&r->x, &r->x, &a);
  fe25519_sub(&r->x, &r->x, &b);
  fe25519_add(&r->z, &d, &b);
  fe25519_sub(&r->t, &r->z, &c);
  fe25519_sub(&r->y, &d, &b);
}

const group_ge group_ge_base = {
  {{UINT64_C(0x62d608f25d51a), UINT64_C(0x412a4b4f6592a), UINT64_C(0x75b7171a4b31d), UINT64_C(0x1ff60527118fe), UINT64_C(0x216936d3cd6e5)}},
  {{UINT64_C(0x6666666666658), UINT64_C(0x4cccccccccccc), UINT64_C(0x1999999999999), UINT64_C(0x3333333333333), UINT64_C(0x6666666666666)}},
  {{1, 0, 0, 0, 0}},
  {{UINT64_C(0x68ab3a5b7dda3), UINT64_C(0xeea2a5eadbb), UINT64_C(0x2af8df483c27e), UINT64_C(0x332b375274732), UINT64_C(0x67875f0fd78b7)}}
};



/* Packing and unpacking is using Hamburg's "ristretto" approach. See
 * https://eprint.iacr.org/2015/673 and
 * https://ristretto.group/
 */
int group_ge_unpack(group_ge *r, const unsigned char x[GROUP_GE_PACKEDBYTES])
{
  fe25519 s, s2, chk, yden, ynum, yden2, xden2, isr, xdeninv, ydeninv, t;
  int ret;
  unsigned char b;

  fe25519_unpack(&s, x);

  /* s = cls.bytesToGf(s,mustBePositive=True) */
  ret = fe25519_isnegative(&s);

  /* yden     = 1-a*s^2    // 1+s^2 */
  /* ynum     = 1+a*s^2    // 1-s^2 */
  fe25519_square(&s2, &s);
  fe25519_add(&yden,&fe25519_one,&s2);
  fe25519_sub(&ynum,&fe25519_one,&s2);

  /* yden_sqr = yden^2 */
  /* xden_sqr = a*d*ynum^2 - yden_sqr */
  fe25519_square(&yden2, &yden);
  fe25519_square(&xden2, &ynum);
  fe25519_mul(&xden2, &xden2, &ge25519_ecd); // d*ynum^2
  fe25519_add(&xden2, &xden2, &yden2); // d*ynum2+yden2
  fe25519_neg(&xden2, &xden2); // -d*ynum2-yden2

  /* isr = isqrt(xden_sqr * yden_sqr) */
  fe25519_mul(&t, &xden2, &yden2);
  fe25519_invsqrt(&isr, &t);

  //Check inverse square root!
  fe25519_square(&chk, &isr);
  fe25519_mul(&chk, &chk, &t);

  ret |= !fe25519_isone(&chk);

  /* xden_inv = isr * yden */
  fe25519_mul(&xdeninv, &isr, &yden);


  /* yden_inv = xden_inv * isr * xden_sqr */
  fe25519_mul(&ydeninv, &xdeninv, &isr);
  fe25519_mul(&ydeninv, &ydeninv, &xden2);

  /* x = 2*s*xden_inv */
  fe25519_mul(&r->x, &s, &xdeninv);
  fe25519_double(&r->x, &r->x);

  /* if negative(x): x = -x */
  b = fe25519_isnegative(&r->x);
  fe25519_neg(&t, &r->x);
  fe25519_cmov(&r->x, &t, b);


  /* y = ynum * yden_inv */
  fe25519_mul(&r->y, &ynum, &ydeninv);

  r->z = fe25519_one;

  /* if cls.cofactor==8 and (negative(x*y) or y==0):
       raise InvalidEncodingException("x*y is invalid: %d, %d" % (x,y)) */
  fe25519_mul(&r->t, &r->x, &r->y);
  ret |= fe25519_isnegative(&r->t);
  ret |= fe25519_iszero(&r->y);


  // Zero all coordinates of point for invalid input; produce invalid point
  fe25519_cmov(&r->x, &fe25519_zero, ret);
  fe25519_cmov(&r->y, &fe25519_zero, ret);
  fe25519_cmov(&r->z, &fe25519_zero, ret);
  fe25519_cmov(&r->t, &fe25519_zero, ret);

  return -ret;
}

/* Packing and unpacking is using Hamburg's "ristretto" approach. See
 * https://eprint.iacr.org/2015/673 and
 * https://ristretto.group/
 */
void group_ge_pack(unsigned char r[GROUP_GE_PACKEDBYTES], const group_ge *x)
{
  fe25519 d, u1, u2, isr, i1, i2, zinv, deninv, nx, ny, s;
  unsigned char b;

  /* u1    = mneg*(z+y)*(z-y) */
  fe25519_add(&d, &x->z, &x->y);
  fe25519_sub(&u1, &x->z, &x->y);
  fe25519_mul(&u1, &u1, &d);

  /* u2    = x*y # = t*z */
  fe25519_mul(&u2, &x->x, &x->y);

  /* isr   = isqrt(u1*u2^2) */
  fe25519_square(&isr, &u2);
  fe25519_mul(&isr, &isr, &u1);
  fe25519_invsqrt(&isr, &isr);

  /* i1    = isr*u1 # sqrt(mneg*(z+y)*(z-y))/(x*y) */
  fe25519_mul(&i1, &isr, &u1);

  /* i2    = isr*u2 # 1/sqrt(a*(y+z)*(y-z)) */
  fe25519_mul(&i2, &isr, &u2);

  /* z_inv = i1*i2*t # 1/z */
  fe25519_mul(&zinv, &i1, &i2);
  fe25519_mul(&zinv, &zinv, &x->t);

  /* if negative(t*z_inv):
       x,y = y*self.i,x*self.i
       den_inv = self.magic * i1 */
  fe25519_mul(&d, &zinv, &x->t);
  b = !fe25519_isnegative(&d);

  fe25519_mul(&nx, &x->y, &fe25519_sqrtm1);
  fe25519_mul(&ny, &x->x, &fe25519_sqrtm1);
  fe25519_mul(&deninv, &ge25519_magic, &i1);

  fe25519_cmov(&nx, &x->x, b);
  fe25519_cmov(&ny, &x->y, b);
  fe25519_cmov(&deninv, &i2, b);

  /* if negative(x*z_inv): y = -y */
  fe25519_mul(&d, &nx, &zinv);
  b = fe25519_isnegative(&d);
  fe25519_neg(&d, &ny);

  fe25519_cmov(&ny, &d, b);

  /* s = (z-y) * den_inv */
  fe25519_sub(&s, &x->z, &ny);
  fe25519_mul(&s, &s, &deninv);

  /* return self.gfToBytes(s,mustBePositive=True) */
  b = fe25519_isnegative(&s);
  fe25519_neg(&d, &s);

  fe25519_cmov(&s, &d, b);

  fe25519_pack(r, &s);
}

void group_ge_add(group_ge *r, const group_ge *x, const group_ge *y)
{
  ge25519_p1p1 t;
  add_p1p1(&t, x, y);
  p1p1_to_p3(r,&t);
}

void group_ge_double(group_ge *r, const group_ge *x)
{
  ge25519_p1p1 t;
  dbl_p1p1(&t, (ge25519_p2 *)x);
  p1p1_to_p3(r,&t);
}

