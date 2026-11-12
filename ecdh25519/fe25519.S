.syntax unified
.thumb    
.text

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
    stmdb.w   sp!, {r4-r9}
    rsb       r9, r2, #0  /* mask */

.rept 5
    ldr.w       r4, [r0]
    ldr.w       r5, [r0, #4]
    ldr.w       r6, [r1]
    ldr.w       r7, [r1, #4]

    eor       r8, r4, r6
    and       r8, r8, r9
    eor       r4, r4, r8

    eor       r8, r5, r7
    and       r8, r8, r9
    eor       r5, r5, r8

    str.w     r4, [r0], #4
    str.w     r5, [r0], #4
    add       r1, r1, #8
.endr

    ldmia.w   sp!, {r4-r9}
    bx        lr
