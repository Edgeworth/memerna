# Copyright 2016 Eliot Courtney.
from pathlib import Path
import re
import cloup

from scripts.util.file import read_file, write_file


MAX = 0x0F0F0F0F
ORDER = "ACGU"


def parse_number(val, default=MAX):
    if val == ".":
        return default
    else:
        assert float(val) * 10 == int(val.replace(".", ""))
        return int(val.replace(".", ""))


# Converts ordering of dangles.
def parse_dangle_file(data):
    lines = [i.strip() for i in re.sub(r" +", " ", data).split("\n")]
    outputs = ["", ""]
    idx = 0
    for i in range(len(lines)):
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
    return outputs


def parse_2x2_file(data):
    lines = [i.strip() for i in re.sub(r" +", " ", data).split("\n")]
    output = ""
    idx = 0
    for i in range(len(lines)):
        if re.match(r"(\s*3\' <-- 5\'\s*){4}", lines[i]):
            matrix_lines = [[j.strip() for j in i.split()] for i in lines[i + 1 : i + 5]]
            for m in range(4):
                for r in range(4):
                    for c in range(4):
                        val = parse_number(matrix_lines[r][m * 4 + c])
                        output += f"{ORDER[idx]}{ORDER[r]}{ORDER[c]}{ORDER[m]} {val}\n"
            idx += 1
    return output


def parse_1x1_internal_loop(data):
    lines = [i.strip() for i in re.sub(r" +", " ", data).split("\n")]
    output = ""
    for i in range(len(lines)):
        if re.match(r"(\s*5\' --> 3\'\s*){6}", lines[i]):
            t3prime = re.findall(r"[GUAC]", lines[i + 2])
            t5prime = re.findall(r"[GUAC]", lines[i + 3])
            assert len(t3prime) == 12
            assert len(t5prime) == 12
            matrix_lines = [[j.strip() for j in i.split()] for i in lines[i + 6 : i + 10]]
            for m in range(6):
                for r in range(4):
                    for c in range(4):
                        val = parse_number(matrix_lines[r][m * 4 + c])
                        output += f"{t3prime[2 * m]}{ORDER[r]}{t3prime[2 * m + 1]}{t5prime[2 * m + 1]}{ORDER[c]}{t5prime[2 * m]} {val}\n"
    return output


def parse_1x2_internal_loop(data):
    lines = [i.strip() for i in re.sub(r" +", " ", data).split("\n")]
    output = ""
    for i in range(len(lines)):
        if re.match(r"(\s*5\' --> 3\'\s*){6}", lines[i]):
            t3prime = re.findall(r"[GUAC]", lines[i + 2])
            t5prime = re.findall(r"[GUAC]", lines[i + 3])
            assert len(t3prime) == 12
            assert len(t5prime) == 12
            extra = re.findall(r"[GUAC]", lines[i + 4])
            assert len(extra) == 6
            matrix_lines = [[j.strip() for j in i.split()] for i in lines[i + 6 : i + 10]]
            for m in range(6):
                for r in range(4):
                    for c in range(4):
                        val = parse_number(matrix_lines[r][m * 4 + c])
                        output += "%s%s%s%s%s%s%s %d\n" % (
                            t3prime[2 * m],
                            ORDER[r],
                            t3prime[2 * m + 1],
                            t5prime[2 * m + 1],
                            extra[m],
                            ORDER[c],
                            t5prime[2 * m],
                            val,
                        )
    return output


def parse_2x2_internal_loop(data):
    lines = [i.strip() for i in re.sub(r" +", " ", data).split("\n")]
    output = ""
    # Skip first example.
    for i in range(15, len(lines)):
        if re.match(r"\s*5\' ------> 3\'\s*", lines[i]):
            t3prime = re.findall(r"[GUAC]", lines[i + 1])
            t5prime = re.findall(r"[GUAC]", lines[i + 2])
            assert len(t3prime) == 2
            assert len(t5prime) == 2
            matrix_lines = [[j.strip() for j in i.split()] for i in lines[i + 4 : i + 20]]
            for x1 in range(4):
                for x2 in range(4):
                    for y1 in range(4):
                        for y2 in range(4):
                            val = parse_number(matrix_lines[4 * x1 + x2][4 * y1 + y2])
                            output += f"{t3prime[0]}{ORDER[x1]}{ORDER[y1]}{t3prime[1]}{t5prime[1]}{ORDER[y2]}{ORDER[x2]}{t5prime[0]} {val}\n"
    return output


def parse_map_file(data):
    m = re.findall(r"([GUAC]+)\s*(\S+)", data)
    return "".join(f"{i[0]} {i[1].replace('.', '')}\n" for i in m)


def parse_loop_file(data):
    m = re.findall(r"(\d+)\s+([0-9.\-+]+)\s+([0-9.\-+]+)\s+([0-9.\-+]+)", data)
    internal, bulge, hairpin = "", "", ""
    for i in m:
        internal += f"{i[0]} {parse_number(i[1], 0)}\n"
        bulge += f"{i[0]} {parse_number(i[2], 0)}\n"
        hairpin += f"{i[0]} {parse_number(i[3], 0)}\n"
    return (internal, bulge, hairpin)


# Outputs AXYA number
def parse_stack_txt(data):
    return parse_2x2_file(data)


def parse_terminal_txt(data):
    return parse_2x2_file(data)


@cloup.command()
@cloup.option(
    "-i", "--input", type=cloup.Path(file_okay=False, exists=True, path_type=Path), required=True
)
@cloup.option(
    "-o",
    "--output",
    type=cloup.Path(file_okay=False, exists=True, writable=True, path_type=Path),
    required=True,
)
def parse_rnastructure_datatables(input, output):
    write_file(
        output / "hairpin.data",
        parse_map_file(read_file(input / "triloop.txt"))
        + parse_map_file(read_file(input / "tloop.txt"))
        + parse_map_file(read_file(input / "hexaloop.txt")),
    )
    write_file(output / "stacking.data", parse_stack_txt(read_file(input / "stack.txt")))
    write_file(output / "terminal.data", parse_terminal_txt(read_file(input / "tstack.txt")))

    internal, bulge, hairpin = parse_loop_file(read_file(input / "loop.txt"))
    write_file(output / "internal_initiation.data", internal)
    write_file(output / "bulge_initiation.data", bulge)
    write_file(output / "hairpin_initiation.data", hairpin)
    write_file(
        output / "internal_1x1.data",
        parse_1x1_internal_loop(read_file(input / "int11.txt")),
    )
    write_file(
        output / "internal_1x2.data",
        parse_1x2_internal_loop(read_file(input / "int21.txt")),
    )
    write_file(
        output / "internal_2x2.data",
        parse_2x2_internal_loop(read_file(input / "int22.txt")),
    )

    dangle3, dangle5 = parse_dangle_file(read_file(input / "dangle.txt"))
    write_file(output / "dangle3.data", dangle3)
    write_file(output / "dangle5.data", dangle5)
