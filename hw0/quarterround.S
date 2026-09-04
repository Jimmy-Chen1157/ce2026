.syntax unified
.cpu cortex-m4
.thumb

// Homework 0: implement the ChaCha20 quarter round in Armv7-M assembly.
//
// void quarterround(uint32_t *a, uint32_t *b, uint32_t *c, uint32_t *d);
//
// From RFC 8439, Section 2.1:
//
//   1.  a += b; d ^= a; d <<<= 16;
//   2.  c += d; b ^= c; b <<<= 12;
//   3.  a += b; d ^= a; d <<<= 8;
//   4.  c += d; b ^= c; b <<<= 7;
//
// where <<<= denotes a rotate left, and additions are modulo 2^32.
//
// Check your work with:
//   make run-qemu PLATFORM=qemu

.global quarterround
.type quarterround, %function
.align 2
quarterround:
    push {r4-r7}
    ldr r4, [r0]                // a
    ldr r5, [r1]                // b
    ldr r6, [r2]                // c
    ldr r7, [r3]                // d

    // TODO: your implementation goes here.
    //       a = r4, b = r5, c = r6, d = r7

    str r4, [r0]
    str r5, [r1]
    str r6, [r2]
    str r7, [r3]
    pop {r4-r7}
    bx lr
