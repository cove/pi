CC ?= cc
CFLAGS ?= -O3 -std=gnu11 -Wall -Wextra -Wpedantic

.PHONY: all fast test clean
all: bbp

bbp: bbp.c
	$(CC) $(CFLAGS) $< -o $@

fast:
	$(CC) $(CFLAGS) -march=native -fopenmp bbp.c -o bbp

test:
	python3 tests/test_bbp.py

clean:
	rm -f bbp
