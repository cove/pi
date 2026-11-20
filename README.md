# pi

A small C program that computes decimal digits of π using base-10,000 integer chunks. It implements Machin's formula:

\[
\pi = 16\arctan(1/5) - 4\arctan(1/239).
\]

Each chunk stores four decimal digits in a `uint32_t`; arithmetic uses bounded `uint64_t` intermediates instead of a big-integer library.

## Build

```sh
cc -O2 -std=c11 -Wall -Wextra -Wpedantic pi_chunks.c -o pi_chunks
```

## Run

```sh
./pi_chunks          # 100 digits after the decimal point
./pi_chunks 1000     # 1,000 digits after the decimal point
```

The current implementation accepts 1–10,000 decimal digits and allocates storage dynamically with `calloc`. It computes all digits from the beginning; it does **not** provide random access to a distant position. For hexadecimal digit extraction at a distant position, a BBP-based program would be a different implementation.
