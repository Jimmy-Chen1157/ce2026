#include "group.h"
#include "smult.h"


int crypto_scalarmult(unsigned char *ss, const unsigned char *sk, const unsigned char *pk)
{
  group_ge p, r0, r1, sum, doubled;
  unsigned char t[32];
  int i;

  for(i=0;i<32;i++) {
    t[i] = sk[i];
  }

  t[0] &= 248;
  t[31] &= 127;
  t[31] |= 64;

  if(group_ge_unpack(&p, pk)) return -1;

  /* Edwards-coordinate ladder: fixed work and no scalar-dependent branch. */
  r0 = group_ge_neutral;
  r1 = p;
  for(i=255;i>=0;i--)
  {
    unsigned char bit = (unsigned char)((t[i >> 3] >> (i & 7)) & 1);
    group_ge_cswap(&r0, &r1, bit);
    group_ge_add(&sum, &r0, &r1);
    group_ge_double(&doubled, &r0);
    r1 = sum;
    r0 = doubled;
    group_ge_cswap(&r0, &r1, bit);
  }

  group_ge_pack(ss, &r0);
  return 0;
}

int crypto_scalarmult_base(unsigned char *pk, const unsigned char *sk)
{
  group_ge table[16], acc, selected;
  unsigned char t[32];
  int i, j;

  for(i = 0; i < 32; i++) t[i] = sk[i];
  t[0] &= 248;
  t[31] &= 127;
  t[31] |= 64;

  /* Build a small fixed-base table, then process four scalar bits at a time.
   * The table scan and all group operations have fixed control flow. */
  table[0] = group_ge_neutral;
  table[1] = group_ge_base;
  for(i = 2; i < 16; i++)
    group_ge_add(&table[i], &table[i - 1], &group_ge_base);

  acc = group_ge_neutral;
  for(i = 63; i >= 0; i--)
  {
    unsigned char digit = (unsigned char)((t[(i * 4) >> 3] >> ((i * 4) & 7)) & 15);
    selected = group_ge_neutral;
    for(j = 0; j < 16; j++) {
      unsigned char match = (unsigned char)(digit == (unsigned char)j);
      fe25519_cmov(&selected.x, &table[j].x, match);
      fe25519_cmov(&selected.y, &table[j].y, match);
      fe25519_cmov(&selected.z, &table[j].z, match);
      fe25519_cmov(&selected.t, &table[j].t, match);
    }
    if(i != 63) {
      group_ge_double(&acc, &acc);
      group_ge_double(&acc, &acc);
      group_ge_double(&acc, &acc);
      group_ge_double(&acc, &acc);
    }
    group_ge_add(&acc, &acc, &selected);
  }
  group_ge_pack(pk, &acc);
  return 0;
}


