# Copyright 2022 Eliot Courtney.
from collections import deque
import re
import string

from scripts.model.rna import Rna

BRACKETS = ["()", "[]", "{}", "<>"]
BRACKETS += ["%s%s" % i for i in zip(string.ascii_lowercase, string.ascii_uppercase)]


def seq_to_primary(seq: str):
    # Only consider fully determined sequences.
    seq = seq.upper()
    assert all(i in "GUAC" for i in seq)
    return seq


def primary_to_seq(r):
    return r


def db_to_secondary(db: str):
    opening = {v[0] for v in BRACKETS}
    closing = {v[1]: v[0] for v in BRACKETS}
    stack = {}
    s = [-1 for i in range(len(db))]
    for i, v in enumerate(db):
        if v in opening:
            stack.setdefault(v, [])
            stack[v].append(i)
        elif v in closing:
            pair = stack[closing[v]].pop()
            s[pair] = i
            s[i] = pair
        else:
            assert v == "."
    return s


# Handles pseudoknots.
def secondary_to_db(s):
    inter = []
    popularity = {}
    stacks = []
    for i, p in enumerate(s):
        assert p != i
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
            assert best != -1
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
            db += BRACKETS[stack_map[sid]][bid]
    return db


def rnas_from_ct_file(data: str):
    q = deque(data.strip().split("\n"))
    rnas = []
    while len(q) > 0:
        length = int(re.search(r"(\d+)", q[0].strip()).group(0))
        subdata = f"{q[0]}\n"
        q.popleft()
        for i in range(length):
            subdata += f"{q[0]}\n"
            q.popleft()
        rnas.append(Rna.from_ct_file(subdata))
    return rnas
