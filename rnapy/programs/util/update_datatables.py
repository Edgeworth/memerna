# Copyright 2023 Eliot Courtney.
from pathlib import Path

import cloup

AU_PENALTY = 0.45
GU_PENALTY = 0.45


def update_internal_loop_end_penalties(
    data: str,
    remove_au: bool = True,
    remove_gu: bool = True,
) -> str:
    out = ""
    for line in data.splitlines():
        s, estr = line.split()
        e = float(estr)
        mid = 0
        if len(s) == 6:
            mid = 3
        elif len(s) == 7:
            mid = 3
        elif len(s) == 8:
            mid = 4
        else:
            raise ValueError(f"invalid internal loop: {s}")
        if remove_au and (s[0] == "A" and s[-1] == "U" or s[0] == "U" and s[-1] == "A"):
            e -= AU_PENALTY
        if remove_gu and (s[0] == "G" and s[-1] == "U" or s[0] == "U" and s[-1] == "G"):
            e -= GU_PENALTY
        if remove_au and (
            s[mid - 1] == "A" and s[mid] == "U" or s[mid - 1] == "U" and s[mid] == "A"
        ):
            e -= AU_PENALTY
        if remove_gu and (
            s[mid - 1] == "G" and s[mid] == "U" or s[mid - 1] == "U" and s[mid] == "G"
        ):
            e -= GU_PENALTY
        out += f"{s} {e:.2f}\n"
    return out


@cloup.command()
@cloup.option(
    "-i",
    "--input",
    "inp",
    type=cloup.Path(file_okay=False, exists=True, resolve_path=True, path_type=Path),
    required=True,
)
@cloup.option(
    "-o",
    "--output",
    "out",
    type=cloup.Path(file_okay=False, exists=True, writable=True, resolve_path=True, path_type=Path),
    required=True,
)
def update_datatables(inp: Path, out: Path) -> None:
    (out / "internal_1x1.data").write_text(
        update_internal_loop_end_penalties((inp / "internal_1x1.data").read_text()),
    )
    (out / "internal_1x2.data").write_text(
        update_internal_loop_end_penalties((inp / "internal_1x2.data").read_text()),
    )
    (out / "internal_2x2.data").write_text(
        update_internal_loop_end_penalties((inp / "internal_2x2.data").read_text()),
    )
