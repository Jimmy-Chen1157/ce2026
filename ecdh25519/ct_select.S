/* Constant-time selectors for the radix-2^51 field and extended points.
 * Cortex-M4 is a 32-bit core, so operate on each 64-bit limb as two words.
 */
	.syntax unified
	.thumb
	.text

	.global fe25519_cmov
	.type fe25519_cmov, %function
	.thumb_func
fe25519_cmov:
	/* void fe25519_cmov(fe25519 *r, const fe25519 *x, unsigned char b) */
	push {r4, lr}
	ands r2, r2, #1
	rsb r4, r2, #0             /* mask = 0 - (b & 1) */
	movs r2, #10               /* five 64-bit limbs */
1:
	ldr r3, [r0]
	ldr r12, [r1]
	eors r12, r3
	ands r12, r4
	eors r3, r12
	str r3, [r0]
	adds r0, #4
	adds r1, #4
	subs r2, #1
	bne 1b
	pop {r4, pc}
	.size fe25519_cmov, .-fe25519_cmov

	.global group_ge_cswap
	.type group_ge_cswap, %function
	.thumb_func
group_ge_cswap:
	/* void group_ge_cswap(group_ge *x, group_ge *y, unsigned char b) */
	push {r4, lr}
	ands r2, r2, #1
	rsb r4, r2, #0             /* mask = 0 - (b & 1) */
	movs r2, #40               /* four coordinates, five 64-bit limbs */
2:
	ldr r3, [r0]
	ldr r12, [r1]
	eors r3, r12
	ands r3, r4                /* d = mask & (x ^ y) */
	ldr r12, [r0]
	eors r12, r3
	str r12, [r0]
	ldr r12, [r1]
	eors r12, r3
	str r12, [r1]
	adds r0, #4
	adds r1, #4
	subs r2, #1
	bne 2b
	pop {r4, pc}
	.size group_ge_cswap, .-group_ge_cswap

	.section .note.GNU-stack,"",%progbits
