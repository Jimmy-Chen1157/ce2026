.syntax unified
.thumb    
.text


.macro FE_MUL4_TERM AREG, ABASE, BREG, BBASE, OREG, OBASE
    ldr.w     r0, [\AREG, #(\ABASE + fe_mul_i * 8)]
    ldr.w     r1, [\BREG, #(\BBASE + fe_mul_j * 8)]
    umull     r2, r3, r0, r1
    ldr.w     r1, [\OREG, #(\OBASE + (fe_mul_i + fe_mul_j) * 8)]
    adds.w    r1, r1, r2
    str.w     r1, [\OREG, #(\OBASE + (fe_mul_i + fe_mul_j) * 8)]
    ldr.w     r1, [\OREG, #(\OBASE + (fe_mul_i + fe_mul_j) * 8 + 4)]
    adc.w     r1, r1, r3
    str.w     r1, [\OREG, #(\OBASE + (fe_mul_i + fe_mul_j) * 8 + 4)]
.endm

.macro FE_SCHOOLBOOK4 AREG, ABASE, BREG, BBASE, OREG, OBASE
.set fe_mul_i, 0
.rept 4
.set fe_mul_j, 0
.rept 4
    FE_MUL4_TERM \AREG, \ABASE, \BREG, \BBASE, \OREG, \OBASE
.set fe_mul_j, fe_mul_j + 1
.endr
.set fe_mul_i, fe_mul_i + 1
.endr
.endm

/* 8x8 Karatsuba, using r5/r6 as inputs, r7 as output, and r8 as the
 * 256-byte scratch area: as[4], bs[4], z0[8], z2[8], zs[8]. */
.macro FE_KARATSUBA8
    movs.w      r0, #0
.set fe_k8_zero_off, 0
.rept 64
    str.w     r0, [r8, #fe_k8_zero_off]
.set fe_k8_zero_off, fe_k8_zero_off + 4
.endr

.set fe_k8_i, 0
.rept 4
    ldr.w     r0, [r5, #(fe_k8_i * 8)]
    ldr.w     r1, [r5, #(fe_k8_i * 8 + 32)]
    add.w     r0, r0, r1
    str.w     r0, [r8, #(fe_k8_i * 8)]
    ldr.w     r0, [r6, #(fe_k8_i * 8)]
    ldr.w     r1, [r6, #(fe_k8_i * 8 + 32)]
    add.w     r0, r0, r1
    str.w     r0, [r8, #(fe_k8_i * 8 + 32)]
.set fe_k8_i, fe_k8_i + 1
.endr

    FE_SCHOOLBOOK4 r5, 0,  r6, 0,  r8, 64
    FE_SCHOOLBOOK4 r5, 32, r6, 32, r8, 128
    FE_SCHOOLBOOK4 r8, 0,  r8, 32, r8, 192

    movs.w      r0, #0
.set fe_k8_outzero, 0
.rept 32
    str.w     r0, [r7, #fe_k8_outzero]
.set fe_k8_outzero, fe_k8_outzero + 4
.endr

    /* out = z0 + (zs-z0-z2)*x^4 + z2*x^8. */
.set fe_k8_i, 0
.rept 8
    ldr.w     r0, [r8, #(64 + fe_k8_i * 8)]
    ldr.w     r1, [r8, #(68 + fe_k8_i * 8)]
    str.w     r0, [r7, #(fe_k8_i * 8)]
    str.w     r1, [r7, #(fe_k8_i * 8 + 4)]
    ldr.w     r0, [r8, #(128 + fe_k8_i * 8)]
    ldr.w     r1, [r8, #(132 + fe_k8_i * 8)]
    str.w     r0, [r7, #(64 + fe_k8_i * 8)]
    str.w     r1, [r7, #(68 + fe_k8_i * 8)]
.set fe_k8_i, fe_k8_i + 1
.endr

.set fe_k8_i, 0
.rept 8
    ldr.w     r0, [r8, #(192 + fe_k8_i * 8)]
    ldr.w     r1, [r8, #(64 + fe_k8_i * 8)]
    subs.w    r0, r0, r1
    ldr.w     r2, [r8, #(196 + fe_k8_i * 8)]
    ldr.w     r3, [r8, #(68 + fe_k8_i * 8)]
    sbc.w     r2, r2, r3
    ldr.w     r1, [r8, #(128 + fe_k8_i * 8)]
    subs.w    r0, r0, r1
    ldr.w     r3, [r8, #(132 + fe_k8_i * 8)]
    sbc.w     r2, r2, r3
    ldr.w     r1, [r7, #(32 + fe_k8_i * 8)]
    adds.w    r1, r1, r0
    str.w     r1, [r7, #(32 + fe_k8_i * 8)]
    ldr.w     r1, [r7, #(36 + fe_k8_i * 8)]
    adc.w     r1, r1, r2
    str.w     r1, [r7, #(36 + fe_k8_i * 8)]
.set fe_k8_i, fe_k8_i + 1
.endr
.endm

/* Specialized for n=16: top-level 16x16 Karatsuba, three 8x8 Karatsuba
 * products, and nine 4x4 schoolbook products. Compile-time expansion avoids
 * all run-time loop/conditional branches. Frame offsets: z0=0, z2=128,
 * zs=256, root sums=384/448, reusable 8x8 scratch=512. */
.align  2
.global _fe25519_mul_coeff
.type   _fe25519_mul_coeff, %function
.thumb_func
_fe25519_mul_coeff:
    stmdb.w   sp!, {r4-r11, lr}
    sub.w     sp, sp, #772
    mov.w       r4, r0                  /* output */
    mov.w       r9, r1                  /* a */
    mov.w       r10, r2                 /* b */
    movs.w      r2, #0

.set fe_root_i, 0
.rept 8
    ldr.w     r0, [r9, #(fe_root_i * 8)]
    ldr.w     r1, [r9, #(64 + fe_root_i * 8)]
    add.w     r0, r0, r1
    str.w     r0, [sp, #(384 + fe_root_i * 8)]
    str.w     r2, [sp, #(388 + fe_root_i * 8)]
    ldr.w     r0, [r10, #(fe_root_i * 8)]
    ldr.w     r1, [r10, #(64 + fe_root_i * 8)]
    add.w     r0, r0, r1
    str.w     r0, [sp, #(448 + fe_root_i * 8)]
    str.w     r2, [sp, #(452 + fe_root_i * 8)]
.set fe_root_i, fe_root_i + 1
.endr

    mov.w       r5, r9
    mov.w       r6, r10
    mov.w       r7, sp
    add.w     r8, sp, #512
    FE_KARATSUBA8

    add.w     r5, r9, #64
    add.w     r6, r10, #64
    add.w     r7, sp, #128
    add.w     r8, sp, #512
    FE_KARATSUBA8

    add.w     r5, sp, #384
    add.w     r6, sp, #448
    add.w     r7, sp, #256
    add.w     r8, sp, #512
    FE_KARATSUBA8

    movs.w      r0, #0
.set fe_root_zero, 0
.rept 64
    str.w     r0, [r4, #fe_root_zero]
.set fe_root_zero, fe_root_zero + 4
.endr

.set fe_root_i, 0
.rept 16
    ldr.w     r0, [sp, #(fe_root_i * 8)]
    ldr.w     r1, [sp, #(fe_root_i * 8 + 4)]
    str.w     r0, [r4, #(fe_root_i * 8)]
    str.w     r1, [r4, #(fe_root_i * 8 + 4)]
    ldr.w     r0, [sp, #(128 + fe_root_i * 8)]
    ldr.w     r1, [sp, #(132 + fe_root_i * 8)]
    ldr.w     r2, [r4, #(128 + fe_root_i * 8)]
    adds.w    r2, r2, r0
    str.w     r2, [r4, #(128 + fe_root_i * 8)]
    ldr.w     r2, [r4, #(132 + fe_root_i * 8)]
    adc.w     r2, r2, r1
    str.w     r2, [r4, #(132 + fe_root_i * 8)]
.set fe_root_i, fe_root_i + 1
.endr

.set fe_root_i, 0
.rept 16
    ldr.w     r0, [sp, #(256 + fe_root_i * 8)]
    ldr.w     r1, [sp, #(fe_root_i * 8)]
    subs.w    r0, r0, r1
    ldr.w     r2, [sp, #(260 + fe_root_i * 8)]
    ldr.w     r3, [sp, #(fe_root_i * 8 + 4)]
    sbc.w     r2, r2, r3
    ldr.w     r1, [sp, #(128 + fe_root_i * 8)]
    subs.w    r0, r0, r1
    ldr.w     r3, [sp, #(132 + fe_root_i * 8)]
    sbc.w     r2, r2, r3
    ldr.w     r1, [r4, #(64 + fe_root_i * 8)]
    adds.w    r1, r1, r0
    str.w     r1, [r4, #(64 + fe_root_i * 8)]
    ldr.w     r1, [r4, #(68 + fe_root_i * 8)]
    adc.w     r1, r1, r2
    str.w     r1, [r4, #(68 + fe_root_i * 8)]
.set fe_root_i, fe_root_i + 1
.endr

    add.w     sp, sp, #772
    ldmia.w   sp!, {r4-r11, lr}
    bx        lr
.size _fe25519_mul_coeff, .-_fe25519_mul_coeff
.purgem FE_KARATSUBA8
.purgem FE_SCHOOLBOOK4
.purgem FE_MUL4_TERM

.align  2
.global _fe25519_add
.type   _fe25519_add, %function
_fe25519_add:
    stmdb.w   sp!, {r4-r7}
    mov.w     r4, r0                  /* output */
    mov.w     r5, r1                  /* x */
    mov.w     r6, r2                  /* y */

.rept 5
    ldr     r0, [r5], #4
    ldr     r1, [r6], #4
    adds.w    r0, r0, r1
    ldr     r2, [r5], #4
    ldr     r3, [r6], #4
    adcs.w    r2, r2, r3
    str     r0, [r4], #4
    str     r2, [r4], #4
.endr

    ldmia.w  sp!, {r4-r7}
    bx      lr

.align  2
.global _fe25519_carry
.type   _fe25519_carry, %function
_fe25519_carry:
    stmdb.w   sp!, {r4-r8}
    mov.w     r4, r0
.rept 3
.rept 4
    ldr.w       r5, [r4, #4]          /* c = r->v[i] (upper) */
    ldr.w       r6, [r4, #8]         /* r6 = r->v[i+1] (lower) */
    ldr.w       r7, [r4, #12]          /* r7 = r->v[i+1] (upper) */
    mov         r8, r5, lsl #13
    mov         r8, r8, lsr #13 
    str.w       r8, [r4, #4]
    mov         r5, r5, lsr #19       /* c = c >> 51 */
    adds.w      r6, r6, r5            
    adcs.w      r7, r7, #0             /* r->v[i+1] = r->v[i+1] + c */
    str.w       r6, [r4, #8]
    str.w       r7, [r4, #12]
    add       r4, r4, #8
.endr
    ldr.w       r5, [r4, #4]          /* c = r->v[i] */
    mov         r8, r5, lsl #13
    mov         r8, r8, lsr #13 
    str.w       r8, [r4, #4]
    mov         r5, r5, lsr #19       /* c = c >> 51 */
    ldr.w       r6, [r0]              /* r6 = r->v[0] (lower) */
    ldr.w       r7, [r0, #4]          /* r7 = r->v[0] (upper) */
    adds.w      r6, r6, r5, lsl #4
    adcs.w      r7, r7, #0
    adds.w      r6, r6, r5, lsl #1
    adcs.w      r7, r7, #0
    adds.w      r6, r6, r5
    adcs.w      r7, r7, #0
    str.w       r6, [r0]
    str.w       r7, [r0, #4]
    mov        r4, r0
.endr
    ldmia.w  sp!, {r4-r8}
    bx lr

.align  2
.global _fe25519_cmov 
.type   _fe25519_cmov, %function
_fe25519_cmov:
    stmdb.w   sp!, {r4-r8}
    rsb       r2, r2, #0  /* mask */

.rept 5
    ldr.w       r4, [r0]
    ldr.w       r5, [r0, #4]
    ldr.w       r6, [r1]
    ldr.w       r7, [r1, #4]

    eor       r8, r4, r6
    and       r8, r8, r2
    eor       r4, r4, r8

    eor       r8, r5, r7
    and       r8, r8, r2
    eor       r5, r5, r8

    str.w     r4, [r0], #4
    str.w     r5, [r0], #4
    add       r1, r1, #8
.endr

    ldmia.w   sp!, {r4-r8}
    bx        lr
