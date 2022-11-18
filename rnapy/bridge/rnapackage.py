# Copyright 2022 Eliot Courtney.
from dataclasses import dataclass
from dataclasses import field
from pathlib import Path

from rnapy.model.model_cfg import EnergyCfg
from rnapy.model.model_cfg import SuboptCfg
from rnapy.model.rna import Rna
from rnapy.util.command import CmdLimits
from rnapy.util.command import CmdResult
from rnapy.util.command import run_cmd


@dataclass
class RnaPackage:
    path: Path
    limits: CmdLimits = field(default_factory=CmdLimits)
    env: dict[str, str] = field(default_factory=dict)

    def name(self) -> str:
        raise NotImplementedError

    def efn(self, rna: Rna, cfg: EnergyCfg) -> tuple[float, CmdResult]:
        raise NotImplementedError

    def fold(self, rna: Rna, cfg: EnergyCfg) -> tuple[Rna, CmdResult]:
        raise NotImplementedError

    def partition(self, rna: Rna, cfg: EnergyCfg) -> None:
        raise NotImplementedError

    def subopt(
        self,
        rna: Rna,
        energy_cfg: EnergyCfg,
        subopt_cfg: SuboptCfg,
    ) -> tuple[list[Rna], CmdResult]:
        raise NotImplementedError

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
