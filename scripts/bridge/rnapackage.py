# Copyright 2022 Eliot Courtney.

from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from scripts.model.config import EnergyCfg, SuboptCfg
from scripts.model.rna import Rna
from scripts.util.command import CmdLimits, run_cmd


@dataclass
class RnaPackage:
    path: Path
    limits: CmdLimits = CmdLimits()

    def efn(self, rna: Rna, cfg: EnergyCfg):
        pass

    def fold(self, rna: Rna, cfg: EnergyCfg):
        pass

    def partition(self, rna: Rna, cfg: EnergyCfg):
        pass

    def subopt(self, rna: Rna, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg):
        pass

    def _run_cmd(
        self,
        *cmd,
        input: Optional[str] = None,
        return_stdout: bool = True,
        stdout_path: Optional[Path] = None,
    ):
        return run_cmd(
            *cmd,
            input=input,
            return_stdout=return_stdout,
            stdout_path=stdout_path,
            limits=self.limits,
            cwd=self.path,
        )
