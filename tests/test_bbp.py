#!/usr/bin/env python3
"""Regression tests using known hexadecimal digits of pi (stdlib only)."""
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
CHECKS = {
    0: "243F6A8885A308D3",
    1: "43F6A8885A308D31",
    16: "13198A2E03707344",
    31: "4A4093822299F31D",
    100: "29B7C97C50DD3F84",
    1000: "49F1C09B075372C9",
    10000: "8AC8FCFB8016CBDB",
    1000000: "6C65E52CB4593500",  # Bailey's BBP reference, 1,000,001st digit
}


def run(binary, *args):
    return subprocess.run(
        [str(binary), *(str(x) for x in args)],
        text=True, capture_output=True, check=False,
    )


def verify(binary):
    for offset, expected in CHECKS.items():
        result = run(binary, offset, 16)
        assert result.returncode == 0, (offset, result.stderr)
        assert result.stdout.strip() == expected, (offset, result.stdout, expected)
        for count in (1, 8, 15):
            partial = run(binary, offset, count)
            assert partial.returncode == 0
            assert partial.stdout.strip() == expected[:count]

    for args in ((), ("invalid",), ("-1",), (0, 0), (0, 17),
                 ("18446744073709551615",), (0, "2.5")):
        assert run(binary, *args).returncode != 0, args


with tempfile.TemporaryDirectory() as temp:
    target = Path(temp) / "bbp"
    subprocess.run(
        [os.environ.get("CC", "cc"), "-O3", "-std=gnu11", "-Wall", "-Wextra",
         "-Wpedantic", str(ROOT / "bbp.c"), "-o", str(target)],
        check=True,
    )
    verify(target)
print("BBP tests passed")
