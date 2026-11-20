/*
 * BBP hexadecimal digit extractor for pi.
 *
 * Extract up to 16 consecutive hexadecimal digits starting at any
 * zero-based position after the radix point, without computing or
 * storing the preceding digits.
 *
 * Fixed-point arithmetic, with 128 fractional bits, avoids floating-
 * point rounding.  Integer arithmetic gives a conservative error bound;
 * if a requested digit lies too close to a boundary, fail instead of
 * printing digits that cannot be certified.
 *
 * Compile: cc -O3 -std=gnu11 -Wall -Wextra -Wpedantic bbp.c -o bbp
 * Parallel: cc -O3 -march=native -fopenmp bbp.c -o bbp
 */

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>

#if !defined(__SIZEOF_INT128__)
#error "This program needs GCC/Clang unsigned __int128 support"
#endif

typedef __uint128_t u128;

/* For moduli up to 2^32, the product of two residues fits in uint64_t. */
static inline uint64_t pow16_mod_small(uint64_t exp, uint64_t m) {
    uint64_t acc = 1 % m;
    uint64_t base = 16 % m;
    while (exp) {
        if (exp & 1) acc = (acc * base) % m;
        exp >>= 1;
        if (exp) base = (base * base) % m;
    }
    return acc;
}

static inline uint64_t pow16_mod_large(uint64_t exp, uint64_t m) {
    uint64_t acc = 1 % m;
    uint64_t base = 16 % m;
    while (exp) {
        if (exp & 1) acc = (uint64_t)(((u128)acc * base) % m);
        exp >>= 1;
        if (exp) base = (uint64_t)(((u128)base * base) % m);
    }
    return acc;
}

static inline uint64_t pow16_mod(uint64_t exp, uint64_t m) {
    return m <= UINT32_MAX ? pow16_mod_small(exp, m)
                           : pow16_mod_large(exp, m);
}

/* Return floor((r/m)*2^128), where 0 <= r < m.
 * Four base-2^32 division steps avoid an out-of-range 128-bit shift.
 * Small moduli use hardware 64-bit division instead of software 128-bit
 * division. The result is exact, not a floating-point approximation.
 */
static inline u128 fraction128(uint64_t r, uint64_t m) {
    u128 out = 0;
    if (m <= UINT32_MAX) {
        for (unsigned i = 0; i != 4; ++i) {
            uint64_t numerator = r << 32;
            uint64_t digit = numerator / m;
            r = numerator % m;
            out = (out << 32) | digit;
        }
    } else {
        for (unsigned i = 0; i != 4; ++i) {
            u128 numerator = (u128)r << 32;
            uint64_t digit = (uint64_t)(numerator / m);
            r = (uint64_t)(numerator % m);
            out = (out << 32) | digit;
        }
    }
    return out;
}

/* The BBP identity is:
 * pi = sum_{k>=0} 16^-k * (4/(8k+1) - 2/(8k+4)
 *                              - 1/(8k+5) - 1/(8k+6)).
 * Multiplying by 16^offset and discarding integer parts leaves the
 * desired hexadecimal digits. Reduction modulo 2^128 is intentional.
 */
static inline u128 prefix_term(uint64_t power, uint64_t k) {
    uint64_t base = 8 * k;
    const uint64_t m1 = base + 1, m4 = base + 4;
    const uint64_t m5 = base + 5, m6 = base + 6;
    u128 a = fraction128(pow16_mod(power, m1), m1);
    u128 b = fraction128(pow16_mod(power, m4), m4);
    u128 c = fraction128(pow16_mod(power, m5), m5);
    u128 d = fraction128(pow16_mod(power, m6), m6);
    return 4 * a - 2 * b - c - d;
}

static u128 extract_fraction(uint64_t offset) {
    u128 sum = 0;
    /* The reduction is exact (unsigned modulo 2^128), so summation order
     * has no effect; OpenMP can safely parallelize this loop. */
#ifdef _OPENMP
#pragma omp parallel for reduction(+:sum) schedule(static)
#endif
    for (uint64_t k = 0; k <= offset; ++k)
        sum += prefix_term(offset - k, k);

    /* Tail: for r>=1, 16^-r/(8*(offset+r)+j). Each contribution
     * floor(2^(128-4r)/m) is exact. For r>=32, the omitted weighted
     * tail is less than one unit of 2^-128. */
    for (unsigned r = 1; r < 32; ++r) {
        unsigned bits = 128 - 4 * r;
        u128 numerator = (u128)1 << bits;
        uint64_t base = 8 * (offset + r);
        sum += 4 * (numerator / (base + 1));
        sum -= 2 * (numerator / (base + 4));
        sum -= numerator / (base + 5);
        sum -= numerator / (base + 6);
    }
    return sum;
}

static int parse_u64(const char *str, uint64_t *out) {
    if (!str[0] || str[0] == '-' || str[0] == '+') return 0;
    errno = 0;
    char *end;
    unsigned long long value = strtoull(str, &end, 10);
    if (errno || *end) return 0;
    *out = (uint64_t)value;
    return 1;
}

int main(int argc, char **argv) {
    uint64_t offset, requested = 16;
    if ((argc != 2 && argc != 3) || !parse_u64(argv[1], &offset) ||
        (argc == 3 && !parse_u64(argv[2], &requested)) ||
        requested < 1 || requested > 16) {
        fprintf(stderr, "Usage: %s HEX_OFFSET [COUNT=16]  (1 <= COUNT <= 16)\n"
                        "HEX_OFFSET is zero-based after the radix point.\n", argv[0]);
        return 2;
    }
    /* Keep the denominators and error bound within their integer widths. */
    if (offset > (UINT64_MAX - 6) / 8 - 32) {
        fprintf(stderr, "Offset exceeds supported integer range\n");
        return 2;
    }
    u128 value = extract_fraction(offset);
    unsigned shift = 128 - 4 * (unsigned)requested;
    u128 mask = ((u128)1 << shift) - 1;
    u128 low = value & mask;
    u128 unit = (u128)1 << shift;

    /* Per j, each prefix/tail truncated division contributes <1 unit
     * of 2^-128. Sum of absolute BBP weights is 8, and the remaining
     * infinite tail contributes less than one additional unit.
     * This deliberately overestimates the actual error.
     */
    u128 error = 8 * ((u128)offset + 32) + 1;
    if (low <= error || unit - low <= error) {
        fprintf(stderr, "Cannot certify %" PRIu64 " digits at offset %" PRIu64
                        " with 128-bit fixed-point precision\n", requested, offset);
        return 1;
    }
    uint64_t digits = (uint64_t)(value >> shift);
    printf("%0*" PRIX64 "\n", (int)requested, digits);
    return 0;
}
