#ifndef MONTGOMERY_H
#define MONTGOMERY_H

/* RFC 7748 X25519 using 32-byte little-endian inputs and outputs.
 * The scalar is clamped internally. The input u-coordinate's top bit is
 * ignored, as required by RFC 7748.
 */
void montgomery_scalarmult(unsigned char out[32],
                           const unsigned char scalar[32],
                           const unsigned char u[32]);

/* X25519 with the standard base u-coordinate 9. */
void montgomery_scalarmult_base(unsigned char out[32],
                                const unsigned char scalar[32]);

#endif
