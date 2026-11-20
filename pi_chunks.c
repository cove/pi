/* Compute pi to arbitrary decimal precision in base-10000 chunks.
 * Every arithmetic operation uses bounded uint64_t intermediates;
 * no floating point or integer wider than a machine register is needed.
 */
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BASE UINT64_C(10000)
#define GUARD_DIGITS 24
#define MAX_DIGITS 10000

typedef uint32_t limb;

/* Most-significant limb first: a[0] is the integer part. */
static void divide(limb *a, size_t count, uint32_t d) {
    uint64_t remainder = 0;
    for (size_t i = 0; i < count; ++i) {
        uint64_t v = remainder * BASE + a[i];
        a[i] = (limb)(v / d);
        remainder = v % d;
    }
}

static int nonzero(const limb *a, size_t count) {
    for (size_t i = 0; i < count; ++i)
        if (a[i]) return 1;
    return 0;
}

static void add(limb *a, const limb *b, size_t count) {
    uint32_t carry = 0;
    for (size_t i = count; i-- > 0;) {
        uint32_t v = a[i] + b[i] + carry;
        a[i] = v % (uint32_t)BASE;
        carry = v / (uint32_t)BASE;
    }
}

/* Requires a >= b. */
static void subtract(limb *a, const limb *b, size_t count) {
    uint32_t borrow = 0;
    for (size_t i = count; i-- > 0;) {
        uint32_t v = b[i] + borrow;
        borrow = a[i] < v;
        a[i] = (a[i] + (uint32_t)BASE - v) % (uint32_t)BASE;
    }
}

static void multiply(limb *a, size_t count, uint32_t m) {
    uint64_t carry = 0;
    for (size_t i = count; i-- > 0;) {
        uint64_t v = (uint64_t)a[i] * m + carry;
        a[i] = (limb)(v % BASE);
        carry = v / BASE;
    }
}

/* arctan(1/x) = 1/x - 1/(3*x^3) + 1/(5*x^5) - ... */
static void arctan_inverse(limb *sum, limb *term, limb *piece,
                           size_t count, uint32_t x) {
    memset(sum, 0, count * sizeof(*sum));
    memset(term, 0, count * sizeof(*term));
    term[0] = 1;
    divide(term, count, x);
    for (uint32_t odd = 1; nonzero(term, count); odd += 2) {
        memcpy(piece, term, count * sizeof(*piece));
        divide(piece, count, odd);
        if ((odd & 2U) == 0) add(sum, piece, count);
        else subtract(sum, piece, count);
        divide(term, count, x * x);
    }
}

int main(int argc, char **argv) {
    char *end;
    unsigned long digits = 100;
    if (argc > 2) {
        fprintf(stderr, "usage: %s [decimal-digits]\n", argv[0]);
        return 1;
    }
    if (argc == 2) {
        digits = strtoul(argv[1], &end, 10);
        if (!argv[1][0] || *end || digits < 1 || digits > MAX_DIGITS) {
            fprintf(stderr, "digits must be between 1 and %d\n", MAX_DIGITS);
            return 1;
        }
    }
    size_t count = 1 + (digits + GUARD_DIGITS + 3) / 4;
    limb *a = calloc(count, sizeof(*a));
    limb *b = calloc(count, sizeof(*b));
    limb *term = calloc(count, sizeof(*term));
    limb *piece = calloc(count, sizeof(*piece));
    if (!a || !b || !term || !piece) {
        fprintf(stderr, "out of memory\n");
        free(a); free(b); free(term); free(piece);
        return 1;
    }
    arctan_inverse(a, term, piece, count, 5);
    arctan_inverse(b, term, piece, count, 239);
    multiply(a, count, 16);
    multiply(b, count, 4);
    subtract(a, b, count);  /* Machin's formula: pi = 16 atan(1/5) - 4 atan(1/239) */

    printf("%" PRIu32 ".", a[0]);
    for (size_t i = 1, printed = 0; printed < digits; ++i) {
        char group[5];
        snprintf(group, sizeof(group), "%04" PRIu32, a[i]);
        size_t take = digits - printed < 4 ? digits - printed : 4;
        fwrite(group, 1, take, stdout);
        printed += take;
    }
    putchar('\n');
    free(a); free(b); free(term); free(piece);
    return 0;
}
