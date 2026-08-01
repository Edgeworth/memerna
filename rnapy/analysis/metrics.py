# Copyright 2022 Eliot Courtney.
from __future__ import annotations

from dataclasses import dataclass
from typing import override

from rnapy.model.rna import Rna


@dataclass
class RnaAccuracy:
    ppv: float
    sensitivity: float

    def f1_score(self) -> float:
        fscore = 2.0 * self.sensitivity * self.ppv
        if fscore != 0:
            fscore /= self.ppv + self.sensitivity
        return fscore

    @staticmethod
    def from_rna(true: Rna, pred: Rna) -> RnaAccuracy:
        if pred.s is None:
            raise ValueError("Predicted structure is None.")
        if true.s is None:
            raise ValueError("True structure is None.")
        if len(true.s) != len(pred.s):
            raise ValueError("True and predicted structures must have the same length.")

        if true.r is not None and pred.r is not None and true.r != pred.r:
            raise ValueError("True and predicted primary sequences must match.")

        num_predicted = sum(1 for i in pred.s if i != -1) / 2
        num_true = sum(1 for i in true.s if i != -1) / 2
        num_correct = 0.0
        for i, pair in enumerate(pred.s):
            if pair != -1 and true.s[i] == pair:
                num_correct += 1
        num_correct /= 2
        if num_correct > num_predicted:
            raise ValueError("More correct pairs than predicted pairs.")
        if num_correct > num_true:
            raise ValueError("More correct pairs than true pairs.")
        if num_true == 0 and num_predicted == 0:
            return RnaAccuracy(ppv=1.0, sensitivity=1.0)
        ppv = num_correct / num_predicted if num_predicted else 0.0
        sensitivity = num_correct / num_true if num_true else 0.0
        return RnaAccuracy(ppv=ppv, sensitivity=sensitivity)

    @override
    def __str__(self) -> str:
        return (
            f"F1-Score: {self.f1_score():.2f} - PPV: {self.ppv:.2f}"
            f" - Sensitivity: {self.sensitivity:.2f}"
        )
