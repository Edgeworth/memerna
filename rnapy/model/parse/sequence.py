# Copyright 2022 Eliot Courtney.
import string


def default_brackets() -> list[str]:
    return ["()", "[]", "{}", "<>"] + [
        f"{a}{b}" for a, b in zip(string.ascii_uppercase, string.ascii_lowercase, strict=True)
    ]


def seq_to_primary(seq: str) -> str:
    # Only consider fully determined sequences.
    seq = seq.upper().replace("T", "U")
    if not all(i in "GUAC" for i in seq):
        raise ValueError(f"Invalid character in sequence: {seq}")
    return seq


def primary_to_seq(r: str) -> str:
    return r


def validate_secondary(s: list[int]) -> None:
    for i, pair in enumerate(s):
        if pair < -1 or pair >= len(s):
            raise ValueError(f"Invalid secondary structure: pair index {pair} at position {i}.")
        if pair == i:
            raise ValueError("Invalid secondary structure: self-pairing.")
        if pair != -1 and s[pair] != i:
            raise ValueError(f"Invalid secondary structure: pair {i}-{pair} is not reciprocal.")


def db_to_secondary(db: str) -> list[int]:
    brackets = default_brackets()
    opening = {v[0] for v in brackets}
    closing = {v[1]: v[0] for v in brackets}
    stack: dict[str, list[int]] = {}
    s = [-1 for i in range(len(db))]
    for i, v in enumerate(db):
        if v in opening:
            stack.setdefault(v, [])
            stack[v].append(i)
        elif v in closing:
            o = closing[v]
            if o not in stack or not stack[o]:
                raise ValueError(f"Invalid character in DB: {v}")
            pair = stack[o].pop()
            s[pair] = i
            s[i] = pair
        elif v != ".":
            raise ValueError(f"Invalid character in DB: {v}")
    if any(stack.values()):
        raise ValueError("Invalid character in DB: unmatched opening bracket")
    return s


def _crossing_graph(s: list[int]) -> tuple[list[tuple[int, int]], list[set[int]]]:
    pairs = [(i, pair) for i, pair in enumerate(s) if pair > i]
    adjacent = [set[int]() for _ in pairs]
    for node, (start, end) in enumerate(pairs):
        for other in range(node):
            other_end = pairs[other][1]
            if start < other_end < end:
                adjacent[node].add(other)
                adjacent[other].add(node)
    return pairs, adjacent


def _dsatur_colors(s: list[int]) -> tuple[list[int], int]:
    pairs, adjacent = _crossing_graph(s)

    pair_colors = [-1] * len(pairs)
    neighbor_colors = [set[int]() for _ in pairs]
    uncolored = set(range(len(pairs)))
    while uncolored:
        node = max(uncolored, key=lambda v: (len(neighbor_colors[v]), len(adjacent[v]), -v))
        color = 0
        while color in neighbor_colors[node]:
            color += 1
        pair_colors[node] = color
        uncolored.remove(node)
        for other in adjacent[node] & uncolored:
            neighbor_colors[other].add(color)

    colors = [-1] * len(s)
    for (start, end), color in zip(pairs, pair_colors, strict=True):
        colors[start] = color
        colors[end] = color
    return colors, max(pair_colors, default=-1) + 1


def secondary_to_db(s: list[int]) -> str:
    brackets = default_brackets()
    validate_secondary(s)

    colors, color_count = _dsatur_colors(s)
    if color_count > len(brackets):
        raise ValueError("Invalid secondary structure: not enough bracket types.")

    popularity: dict[int, int] = {}
    for i, pair in enumerate(s):
        if pair != -1 and pair < i:
            popularity.setdefault(colors[i], 0)
            popularity[colors[i]] += 1
    stacks_by_pop = sorted(popularity.keys(), key=lambda x: popularity[x], reverse=True)
    stack_map = {v: i for i, v in enumerate(stacks_by_pop)}
    return "".join(
        "." if pair == -1 else brackets[stack_map[colors[i]]][int(pair < i)]
        for i, pair in enumerate(s)
    )
