#!/usr/bin/env python3
"""Export the leading Schwarzschild psi0 coefficients from a Wolfram WXF file.

The source WXF uses the mostly-plus convention.  The emitted table uses the
mostly-minus convention used by the C++ code and by BondiGauge_Adjusted.nb.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from wolframclient.deserializers import binary_deserialize
from wolframclient.language.expression import WLFunction


def complex_parts(value):
    if isinstance(value, WLFunction):
        if str(value.head) != "Complex" or len(value.args) != 2:
            raise ValueError(f"unsupported Wolfram expression: {value!r}")
        return value.args[0], value.args[1]
    return value, 0


def negate_decimal_text(value) -> str:
    """Flip a serialized decimal without rounding it through Decimal context."""
    text = str(value)
    if text.startswith("-"):
        return text[1:]
    if value == 0:
        return text
    return "-" + text


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    data = binary_deserialize(args.input.read_bytes())
    if data["FailedModes"]:
        raise RuntimeError(f"WXF contains failed modes: {data['FailedModes']}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        stream.write("# GHZ_PSI0_LM_V1\n")
        stream.write("# signature=mostly-minus\n")
        stream.write("# source_signature=mostly-plus\n")
        stream.write("# source_to_output_factor=-1\n")
        stream.write("# mass=1\n")
        stream.write("# orbital_radius=10\n")
        stream.write(f"# ell_min={data['LMin']}\n")
        stream.write(f"# ell_max={data['ellMax']}\n")
        stream.write(f"# fit_r_min={data['rFitMin']}\n")
        stream.write(f"# fit_r_max={data['rFitMax']}\n")
        stream.write(f"# fit_leading_power={data['p']}\n")
        stream.write(f"# working_precision={data['WorkingPrecision']}\n")
        writer = csv.writer(stream)
        writer.writerow(("ell", "m", "psi0_5o_real", "psi0_5o_imag"))
        for ell, m in data["Modes"]:
            real, imag = complex_parts(data["LeadingCoefficients"][(ell, m)])
            writer.writerow(
                (ell, m, negate_decimal_text(real), negate_decimal_text(imag))
            )


if __name__ == "__main__":
    main()
