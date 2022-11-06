# Copyright 2022 Eliot Courtney.
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Iterator, Sequence
import pandas as pd

from rnapy.model.rna import Rna


@dataclass
class Dataset:
    name: str
    dfs: dict[str, pd.DataFrame] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.dfs)

    def __iter__(self) -> Iterator[str]:
        return iter(self.dfs)

    def __getitem__(self, key: str) -> pd.DataFrame:
        return self.dfs[key]

    def concat(self, other: Dataset) -> Dataset:
        # Concat shared dataframes:
        shared = self.keys() & other.keys()
        dfs = {k: pd.concat([self[k], other[k]]) for k in shared}
        # Add unique dataframes:
        dfs.update({k: v for k, v in self.dfs.items() if k not in shared})
        dfs.update({k: v for k, v in other.dfs.items() if k not in shared})
        return Dataset(self.name, dfs)

    def keys(self) -> set[str]:
        return set(self.dfs)

    def exclude(self, names: str | Sequence[str]) -> Dataset:
        if isinstance(names, str):
            names = [names]
        return Dataset(self.name, {k: v for k, v in self.dfs.items() if k not in names})


@dataclass
class RnaAccuracy:
    ppv: float
    sensitivity: float

    def fscore(self) -> float:
        fscore = 2.0 * self.sensitivity * self.ppv
        if fscore != 0:
            fscore /= self.ppv + self.sensitivity
        return fscore

    @staticmethod
    def from_rna(true: Rna, pred: Rna) -> "RnaAccuracy":
        assert pred.s is not None and true.s is not None
        num_predicted = sum(1 for i in pred.s if i != -1) / 2
        num_true = sum(1 for i in true.s if i != -1) / 2
        num_correct = 0.0
        for i, pair in enumerate(pred.s):
            if pair != -1 and true.s[i] == pair:
                num_correct += 1
        num_correct /= 2
        assert num_correct <= num_predicted and num_correct <= num_true
        ppv, sensitivity = num_correct, num_correct
        if num_correct != 0:
            ppv /= num_predicted
            sensitivity /= num_true
        return RnaAccuracy(ppv=ppv, sensitivity=sensitivity)

    def __str__(self) -> str:
        return (
            f"F-Score: {self.fscore:.2f} - PPV: {self.ppv:.2f}"
            f" - Sensitivity: {self.sensitivity:.2f}"
        )
