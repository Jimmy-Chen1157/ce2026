#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "chacha20.h"
#include "hal.h"

#include "testvectors.inc"

/* RFC 8439, Section 2.1.1: test vector for a single quarter round. */
static const uint32_t qr_in[4]  = {0x11111111, 0x01020304, 0x9b8d6f43, 0x01234567};
static const uint32_t qr_out[4] = {0xea2a92f4, 0xcb1cf8ce, 0x4581472e, 0x5881c4bb};

#define NRANDOM 256

static unsigned char out[OUTLEN];

/* Deterministic xorshift32, so any failure is reproducible. */
static uint32_t rng_state = 0x2545f491;
static uint32_t rng_next(void)
{
  uint32_t x = rng_state;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  rng_state = x;
  return x;
}

static void print_state(const char *label, const uint32_t x[4])
{
  char outstr[128];
  sprintf(outstr, "  %s a=%08lx b=%08lx c=%08lx d=%08lx\n", label,
          (unsigned long)x[0], (unsigned long)x[1],
          (unsigned long)x[2], (unsigned long)x[3]);
  hal_send_str(outstr);
}

/* Test 1: the provided C reference against RFC 8439. This passes out of the
   box. If it ever fails, the code we gave you is broken, not yours. */
static int test_reference(void)
{
  uint32_t x[4];

  hal_send_str("\n=== Test 1: C reference against RFC 8439 ===\n");
  memcpy(x, qr_in, sizeof x);
  quarterround_ref(&x[0], &x[1], &x[2], &x[3]);

  if(memcmp(x, qr_out, sizeof x) != 0)
  {
    hal_send_str("quarterround_ref does not match RFC 8439!\n");
    print_state("expected", qr_out);
    print_state("got     ", x);
    return -1;
  }

  hal_send_str("quarterround_ref matches RFC 8439\n");
  return 0;
}

/* Test 2: your assembly against the same single-round vector. */
static int test_asm_vector(void)
{
  uint32_t x[4];

  hal_send_str("\n=== Test 2: your quarterround against RFC 8439 ===\n");
  memcpy(x, qr_in, sizeof x);
  quarterround(&x[0], &x[1], &x[2], &x[3]);

  if(memcmp(x, qr_out, sizeof x) != 0)
  {
    hal_send_str("quarterround does not match RFC 8439!\n");
    print_state("expected", qr_out);
    print_state("got     ", x);
    return -1;
  }

  hal_send_str("quarterround matches RFC 8439\n");
  return 0;
}

/* Test 3: your assembly must agree with the C reference on many inputs, not
   just one vector. Catches rotate amounts that happen to be right on the
   single RFC vector. */
static int test_asm_equivalence(void)
{
  uint32_t xr[4], xa[4];
  int i, j;

  hal_send_str("\n=== Test 3: your quarterround against C on random inputs ===\n");

  for(i = 0; i < NRANDOM; i++)
  {
    for(j = 0; j < 4; j++)
    {
      xr[j] = xa[j] = rng_next();
    }

    quarterround_ref(&xr[0], &xr[1], &xr[2], &xr[3]);
    quarterround(&xa[0], &xa[1], &xa[2], &xa[3]);

    if(memcmp(xr, xa, sizeof xr) != 0)
    {
      char outstr[64];
      sprintf(outstr, "mismatch on random input %d!\n", i);
      hal_send_str(outstr);
      print_state("C  ", xr);
      print_state("asm", xa);
      return -1;
    }
  }

  hal_send_str("quarterround agrees with C on all random inputs\n");
  return 0;
}

/* Test 4: the whole cipher. ChaCha20 calls your quarterround 8 times per
   double round, so this only passes if it is correct in context. */
static int test_chacha20(void)
{
  int i;

  hal_send_str("\n=== Test 4: ChaCha20 keystream test vector ===\n");
  crypto_stream_chacha20(out, OUTLEN, nonce, key);

  for(i = 0; i < OUTLEN; i++)
  {
    if(cmp[i] != out[i])
    {
      char outstr[96];
      sprintf(outstr, "keystream differs at byte %d: expected %02x, got %02x\n",
              i, cmp[i], out[i]);
      hal_send_str(outstr);
      return -1;
    }
  }

  hal_send_str("ChaCha20 keystream matches the test vector\n");
  return 0;
}

static void run_speed(void)
{
  char outstr[128];
  uint64_t oldcount, newcount;

  hal_send_str("\n=== Benchmarks ===\n");

  oldcount = hal_get_time();
  crypto_stream_chacha20(out, OUTLEN, nonce, key);
  newcount = hal_get_time();
  sprintf(outstr, "cycles for %d bytes of ChaCha20 keystream: ", OUTLEN);
  hal_send_str(outstr);
#ifdef MPS2_AN386
  (void)oldcount; (void)newcount;
  sprintf(outstr, "[cycle counts not meaningful in qemu emulation]\n");
#else
  sprintf(outstr, "%llu\n", newcount-oldcount);
#endif
  hal_send_str(outstr);

  hal_send_str("Benchmarks completed!\n");
}

int main(void)
{
  int test_result = 0;

  hal_setup(CLOCK_FAST);

  hal_send_str("\nHomework 0: ChaCha20 quarter round in Armv7-M assembly\n");

  test_result |= test_reference();
  test_result |= test_asm_vector();
  test_result |= test_asm_equivalence();
  test_result |= test_chacha20();

  run_speed();

  if(test_result != 0) {
    hal_send_str("\n*** TEST FAILED ***\n");
    return -1;
  }

  hal_send_str("\n*** ALL GOOD ***\n");
  return 0;
}
