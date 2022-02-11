# Copyright 2022 Eliot Courtney.
from dataclasses import dataclass
from dataclasses import field
from pathlib import Path

from scripts.model.config import EnergyCfg
from scripts.model.config import SuboptCfg
from scripts.model.rna import Rna
from scripts.util.command import CmdLimits
from scripts.util.command import run_cmd
from scripts.util.command import CmdResult
from typing import Tuple


@dataclass
class RnaPackage:
    path: Path
    limits: CmdLimits = field(default_factory=CmdLimits)
    env: dict[str, str] = field(default_factory=dict)

    def efn(self, rna: Rna, cfg: EnergyCfg) -> Tuple[float, CmdResult]:
        pass

    def fold(self, rna: Rna, cfg: EnergyCfg) -> Tuple[Rna, CmdResult]:
        pass

    def partition(self, rna: Rna, cfg: EnergyCfg) -> None:
        pass

    def subopt(
        self, rna: Rna, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg
    ) -> Tuple[list[Rna], CmdResult]:
        pass

    def _run_cmd(
        self,
        *cmd: str,
        inp: str | None = None,
        return_stdout: bool = True,
        stdout_path: Path | None = None,
    ) -> CmdResult:
        return run_cmd(
            *cmd,
            inp=inp,
            return_stdout=return_stdout,
            stdout_path=stdout_path,
            limits=self.limits,
            cwd=self.path,
            extra_env=self.env,
        )
