# Copyright 2022 Eliot Courtney.
from __future__ import annotations

from dataclasses import dataclass

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
        ppv, sensitivity = num_correct, num_correct
        if num_correct != 0:
            ppv /= num_predicted
            sensitivity /= num_true
        return RnaAccuracy(ppv=ppv, sensitivity=sensitivity)

    def __str__(self) -> str:
        return (
            f"F1-Score: {self.f1_score:.2f} - PPV: {self.ppv:.2f}"
            f" - Sensitivity: {self.sensitivity:.2f}"
        )
