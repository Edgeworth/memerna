# Copyright 2022 Eliot Courtney.
from re import L


class RnaAccuracy:
    def __init__(self, ppv, sensitivity):
        self.ppv = ppv
        self.sensitivity = sensitivity
        self.fscore = 2.0 * sensitivity * ppv
        if self.fscore != 0:
            self.fscore /= ppv + sensitivity

    @staticmethod
    def from_rna(true, predicted):
        num_predicted = sum(1 for i in predicted.pairs if i != -1) / 2
        num_true = sum(1 for i in true.pairs if i != -1) / 2
        num_correct = 0
        for i, pair in enumerate(predicted.pairs):
            if pair != -1 and true.pairs[i] == pair:
                num_correct += 1
        num_correct /= 2
        assert num_correct <= num_predicted and num_correct <= num_true
        ppv, sensitivity = num_correct, num_correct
        if num_correct != 0:
            ppv /= num_predicted
            sensitivity /= num_true
        return RnaAccuracy(ppv=ppv, sensitivity=sensitivity)

    def __str__(self):
        return f"F-Score: {self.fscore:.2f} - PPV: {self.ppv:.2f} - Sensitivity: {self.sensitivity:.2f}"
