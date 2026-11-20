# pi

Two standalone C programs to compute π, using integer arithmetic instead of a big-number library.

## BBP: jump to any hexadecimal position

`bbp.c` implements the **Bailey–Borwein–Plouffe digit-extraction algorithm**, which calculates up to 16 hexadecimal digits of π starting at a specified position **without computing or storing preceding digits**.

```sh
make                 # optimized serial build
./bbp 0              # 243F6A8885A308D3 (first 16 digits after 3.)
./bbp 1000000        # 6C65E52CB4593500
./bbp 1000000 8      # 6C65E52C (8 digits)
```

The offset is **zero-based**, counting hexadecimal digits after the radix point: `0` means the first digit (`2`); `1000000` means position **1,000,001**. The second argument selects 1–16 consecutive hexadecimal digits and defaults to 16.

For multi-core CPUs with GCC/Clang and OpenMP, use:

```sh
make fast
OMP_NUM_THREADS=5 ./bbp 10000000 16
```

### How it works

The formula is

```text
π = Σ(k=0..∞) 16⁻ᵏ [4/(8k+1) − 2/(8k+4) − 1/(8k+5) − 1/(8k+6)].
```

The program calculates modular powers `16^(offset-k) mod (8k+j)` for terms before the requested position, then evaluates a short, rapidly decreasing tail. It uses a 128-bit fixed-point fraction and a conservative truncation-error bound; if it cannot certify the requested digits, it reports an error instead of printing an unreliable answer. Moduli fitting in 32 bits take a faster 64-bit arithmetic path. Summation is exact modulo `2^128`, so OpenMP can parallelize it without rounding issues.

- **Space:** O(1) per worker, independent of requested position.
- **Time:** O(offset × log(offset)) for a single requested window. Huge offsets still require substantial computation.
- **Portability:** 64-bit GCC or Clang with `unsigned __int128`; OpenMP is optional.
- **Output:** hexadecimal, not decimal. Computing all preceding digits with this algorithm is inefficient; use a full-precision π algorithm for that purpose.

### Tests

```sh
make test
```

The regression tests compare known hexadecimal windows, including the millionth offset from [Bailey's BBP paper](https://www.davidhbailey.com/dhbpapers/bbp-alg.pdf).

## Sequential decimal digits

`pi_chunks.c` uses Machin's formula, with base-10,000 chunks, for computing a *prefix* of π in decimal:

```sh
cc -O2 -std=c11 pi_chunks.c -o pi_chunks
./pi_chunks 1000
```

Its current version is limited to 10,000 decimal digits; it is separate from the random-access hexadecimal BBP implementation.
