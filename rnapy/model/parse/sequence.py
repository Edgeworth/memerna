# Copyright 2022 Eliot Courtney.
import string


def default_brackets() -> list[str]:
    return ["()", "[]", "{}", "<>"] + [
        f"{a}{b}" for a, b in zip(string.ascii_lowercase, string.ascii_uppercase, strict=True)
    ]


def seq_to_primary(seq: str) -> str:
    # Only consider fully determined sequences.
    seq = seq.upper()
    if not all(i in "GUAC" for i in seq):
        raise ValueError(f"Invalid character in sequence: {seq}")
    return seq


def primary_to_seq(r: str) -> str:
    return r


def db_to_secondary(db: str, brackets: list[str] | None = None) -> list[int]:
    brackets = brackets or default_brackets()
    opening = {v[0] for v in brackets}
    closing = {v[1]: v[0] for v in brackets}
    stack: dict[str, list[int]] = {}
    s = [-1 for i in range(len(db))]
    for i, v in enumerate(db):
        if v in opening:
            stack.setdefault(v, [])
            stack[v].append(i)
        elif v in closing:
            pair = stack[closing[v]].pop()
            s[pair] = i
            s[i] = pair
        elif v != ".":
            raise ValueError(f"Invalid character in DB: {v}")
    return s


def secondary_to_db(s: list[int], brackets: list[str] | None = None) -> str:
    """Handles pseudoknots."""
    brackets = brackets or default_brackets()
    inter = []
    popularity: dict[int, int] = {}
    stacks: list[list[int]] = []
    for i, p in enumerate(s):
        if p == i:
            raise ValueError("Invalid secondary structure: self-pairing.")
        if p == -1:
            inter.append((-1, 0))
        elif p > i:
            best = len(stacks)
            for sid, stack in enumerate(stacks):
                if not stack or p < stack[-1]:
                    best = sid
            if best == len(stacks):
                stacks.append([])
            stacks[best].append(p)
            inter.append((best, 0))
        else:
            best = -1
            for sid, stack in enumerate(stacks):
                if stack and stack[-1] == i:
                    best = sid
            if best == -1:
                raise ValueError("Invalid secondary structure: unpaired closing bracket.")
            stacks[best].pop()
            inter.append((best, 1))
            popularity.setdefault(best, 0)
            popularity[best] += 1
    db = ""
    stacks_by_pop = sorted(popularity.keys(), key=lambda x: popularity[x], reverse=True)
    stack_map = {v: i for i, v in enumerate(stacks_by_pop)}
    for sid, bid in inter:
        if sid == -1:
            db += "."
        else:
            db += brackets[stack_map[sid]][bid]
    return db
