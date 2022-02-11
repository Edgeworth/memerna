# Copyright 2022 Eliot Courtney.
from dataclasses import dataclass
from scripts.model.rna import Rna


@dataclass
class RnaAccuracy:
    ppv: float
    sensitivity: float

    def __init__(self, ppv: float, sensitivity: float):
        self.ppv = ppv
        self.sensitivity = sensitivity

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
        return f"F-Score: {self.fscore:.2f} - PPV: {self.ppv:.2f} - Sensitivity: {self.sensitivity:.2f}"
