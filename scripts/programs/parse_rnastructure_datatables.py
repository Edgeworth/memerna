# Copyright 2016 Eliot Courtney.
import itertools
from pathlib import Path
import re
from typing import Tuple

import cloup

MAX = 0x0F0F0F0F
ORDER = "ACGU"


def parse_number(val: str, default: int = MAX) -> int:
    if val == ".":
        return default
    res = int(val.replace(".", ""))
    if float(val) * 10 != res:
        raise ValueError(f"invalid number: {val}")
    return res


# Converts ordering of dangles.
def parse_dangle_file(data: str) -> Tuple[str, str]:
    lines = [i.strip() for i in re.sub(r" +", " ", data).split("\n")]
    outputs = ["", ""]
    idx = 0
    for i in range(len(lines)):  # pylint: disable=consider-using-enumerate
        if re.match(r"(\s*5\' --> 3\'\s*){4}", lines[i]):
            output_idx = 1
            if lines[i + 1].count("X") > 0:
                output_idx = 0
            values = [parse_number(i.strip()) for i in lines[i + 4].split()]
            for m in range(4):
                for c in range(4):
                    outputs[
                        output_idx
                    ] += f"{ORDER[idx % 4]}{ORDER[c]}{ORDER[m]} {values[m * 4 + c]}\n"
            idx += 1
    return (outputs[0], outputs[1])


def parse_2x2_file(data: str) -> str:
    lines = [i.strip() for i in re.sub(r" +", " ", data).split("\n")]
    output = ""
    idx = 0
    for i in range(len(lines)):  # pylint: disable=consider-using-enumerate
        if re.match(r"(\s*3\' <-- 5\'\s*){4}", lines[i]):
            matrix_lines = [[j.strip() for j in i.split()] for i in lines[i + 1 : i + 5]]
            for m in range(4):
                for r in range(4):
                    for c in range(4):
                        val = parse_number(matrix_lines[r][m * 4 + c])
                        output += f"{ORDER[idx]}{ORDER[r]}{ORDER[c]}{ORDER[m]} {val}\n"
            idx += 1
    return output


def parse_1x1_internal_loop(data: str) -> str:
    lines = [i.strip() for i in re.sub(r" +", " ", data).split("\n")]
    output = ""
    for i in range(len(lines)):  # pylint: disable=consider-using-enumerate
        if re.match(r"(\s*5\' --> 3\'\s*){6}", lines[i]):
            t3prime = re.findall(r"[GUAC]", lines[i + 2])
            t5prime = re.findall(r"[GUAC]", lines[i + 3])
            if len(t3prime) != 12:
                raise ValueError(f"invalid t3prime: {t3prime}")
            if len(t5prime) != 12:
                raise ValueError(f"invalid t5prime: {t5prime}")
            matrix_lines = [[j.strip() for j in i.split()] for i in lines[i + 6 : i + 10]]
            for m in range(6):
                for r in range(4):
                    for c in range(4):
                        val = parse_number(matrix_lines[r][m * 4 + c])
                        output += f"{t3prime[2 * m]}{ORDER[r]}{t3prime[2 * m + 1]}{t5prime[2 * m + 1]}{ORDER[c]}{t5prime[2 * m]} {val}\n"
    return output


def parse_1x2_internal_loop(data: str) -> str:
    lines = [i.strip() for i in re.sub(r" +", " ", data).split("\n")]
    output = ""
    for i in range(len(lines)):  # pylint: disable=consider-using-enumerate
        if re.match(r"(\s*5\' --> 3\'\s*){6}", lines[i]):
            t3prime = re.findall(r"[GUAC]", lines[i + 2])
            t5prime = re.findall(r"[GUAC]", lines[i + 3])
            if len(t3prime) != 12:
                raise ValueError(f"invalid t3prime: {t3prime}")
            if len(t5prime) != 12:
                raise ValueError(f"invalid t5prime: {t5prime}")

            extra = re.findall(r"[GUAC]", lines[i + 4])
            if len(extra) != 6:
                raise ValueError(f"invalid extra: {extra}")
            matrix_lines = [[j.strip() for j in i.split()] for i in lines[i + 6 : i + 10]]
            for m in range(6):
                for r in range(4):
                    for c in range(4):
                        val = parse_number(matrix_lines[r][m * 4 + c])
                        output += f"{t3prime[2 * m]}{ORDER[r]}{t3prime[2 * m + 1]}{t5prime[2 * m + 1]}{extra[m]}{ORDER[c]}{t5prime[2 * m]} {val}\n"
    return output


def parse_2x2_internal_loop(data: str) -> str:
    lines = [i.strip() for i in re.sub(r" +", " ", data).split("\n")]
    output = ""
    # Skip first example.
    for i in range(15, len(lines)):
        if re.match(r"\s*5\' ------> 3\'\s*", lines[i]):
            t3prime = re.findall(r"[GUAC]", lines[i + 1])
            t5prime = re.findall(r"[GUAC]", lines[i + 2])
            if len(t3prime) != 2:
                raise ValueError(f"invalid t3prime: {t3prime}")
            if len(t5prime) != 2:
                raise ValueError(f"invalid t5prime: {t5prime}")

            matrix_lines = [[j.strip() for j in i.split()] for i in lines[i + 4 : i + 20]]
            for x1, x2, y1, y2 in itertools.product(range(4), range(4), range(4), range(4)):
                val = parse_number(matrix_lines[4 * x1 + x2][4 * y1 + y2])
                output += f"{t3prime[0]}{ORDER[x1]}{ORDER[y1]}{t3prime[1]}{t5prime[1]}{ORDER[y2]}{ORDER[x2]}{t5prime[0]} {val}\n"
    return output


def parse_map_file(data: str) -> str:
    m = re.findall(r"([GUAC]+)\s*(\S+)", data)
    return "".join(f"{i[0]} {i[1].replace('.', '')}\n" for i in m)


def parse_loop_file(data: str) -> Tuple[str, str, str]:
    m = re.findall(r"(\d+)\s+([0-9.\-+]+)\s+([0-9.\-+]+)\s+([0-9.\-+]+)", data)
    internal, bulge, hairpin = "", "", ""
    for i in m:
        internal += f"{i[0]} {parse_number(i[1], 0)}\n"
        bulge += f"{i[0]} {parse_number(i[2], 0)}\n"
        hairpin += f"{i[0]} {parse_number(i[3], 0)}\n"
    return (internal, bulge, hairpin)


# Outputs AXYA number
def parse_stack_txt(data: str) -> str:
    return parse_2x2_file(data)


def parse_terminal_txt(data: str) -> str:
    return parse_2x2_file(data)


@cloup.command()
@cloup.option(
    "-i",
    "--input",
    type=cloup.Path(file_okay=False, exists=True, path_type=Path),
    required=True,
)
@cloup.option(
    "-o",
    "--output",
    type=cloup.Path(file_okay=False, exists=True, writable=True, path_type=Path),
    required=True,
)
def parse_rnastructure_datatables(inp: Path, out: Path) -> None:
    (out / "hairpin.data").write_text(
        parse_map_file((inp / "triloop.txt").read_text())
        + parse_map_file((inp / "tloop.txt").read_text())
        + parse_map_file((inp / "hexaloop.txt").read_text()),
    )
    (out / "stacking.data").write_text(parse_stack_txt((inp / "stack.txt").read_text()))
    (out / "terminal.data").write_text(parse_terminal_txt((inp / "tstack.txt").read_text()))

    internal, bulge, hairpin = parse_loop_file((inp / "loop.txt").read_text())
    (out / "internal_initiation.data").write_text(internal)
    (out / "bulge_initiation.data").write_text(bulge)
    (out / "hairpin_initiation.data").write_text(hairpin)
    (out / "internal_1x1.data").write_text(
        parse_1x1_internal_loop((inp / "int11.txt").read_text()),
    )
    (out / "internal_1x2.data").write_text(
        parse_1x2_internal_loop((inp / "int21.txt").read_text()),
    )
    (out / "internal_2x2.data").write_text(
        parse_2x2_internal_loop((inp / "int22.txt").read_text()),
    )

    dangle3, dangle5 = parse_dangle_file((inp / "dangle.txt").read_text())
    (out / "dangle3.data").write_text(dangle3)
    (out / "dangle5.data").write_text(dangle5)
