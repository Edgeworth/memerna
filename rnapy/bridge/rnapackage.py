# Copyright 2022 Eliot Courtney.
from dataclasses import dataclass, field
from decimal import Decimal
from pathlib import Path

from rnapy.model.model_cfg import EnergyCfg, SuboptCfg
from rnapy.model.rna import Rna
from rnapy.util.command import CmdLimits, CmdResult, run_cmd
from rnapy.util.util import strict_merge


@dataclass
class RnaPackage:
    path: Path
    limits: CmdLimits = field(default_factory=CmdLimits)
    env: dict[str, str] = field(default_factory=dict)

    def desc(self, *, energy_cfg: EnergyCfg | None, subopt_cfg: SuboptCfg | None) -> dict[str, str]:
        desc = {"package_name": self.package_name()}
        if energy_cfg is not None:
            desc = strict_merge(desc, energy_cfg.desc())
        if subopt_cfg is not None:
            desc = strict_merge(desc, subopt_cfg.desc())
        return desc

    def package_name(self) -> str:
        raise NotImplementedError

    def efn(self, rna: Rna, cfg: EnergyCfg) -> tuple[Decimal, CmdResult]:
        raise NotImplementedError

    def fold(self, rna: Rna, cfg: EnergyCfg) -> tuple[Rna, CmdResult]:
        raise NotImplementedError

    def partition(self, rna: Rna, cfg: EnergyCfg) -> None:
        raise NotImplementedError

    def subopt(
        self, rna: Rna, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg
    ) -> tuple[list[Rna] | int, CmdResult]:
        raise NotImplementedError

    def _run_cmd(
        self,
        *cmd: str,
        stdout_to_str: bool,
        stdin_inp: str | None = None,
        stdout_path: Path | None = None,
    ) -> CmdResult:
        return run_cmd(
            *cmd,
            stdin_inp=stdin_inp,
            stdout_to_str=stdout_to_str,
            stdout_path=stdout_path,
            limits=self.limits,
            cwd=self.path,
            extra_env=self.env,
        )
