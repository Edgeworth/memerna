# Copyright 2022 Eliot Courtney.
import re
import tempfile
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import override

from rnapy.bridge.rnapackage import RnaPackage
from rnapy.model.model_cfg import CtdCfg, EnergyCfg, LonelyPairs, SuboptCfg
from rnapy.model.parse.rna_parser import RnaParser
from rnapy.model.rna import Rna
from rnapy.util.command import CmdResult
from rnapy.util.util import fast_linecount


@dataclass
class RNAstructure(RnaPackage):
    def __post_init__(self) -> None:
        self.env["DATAPATH"] = str(self.path / "data_tables")

    def check_energy_cfg(self, cfg: EnergyCfg) -> None:
        if cfg.lonely_pairs != LonelyPairs.HEURISTIC:
            raise NotImplementedError(
                "RNAstructure does not support modifying lonely pairs behavior"
            )
        if cfg.ctd != CtdCfg.ALL:
            raise NotImplementedError("RNAstructure does not support turning off CTDs")
        if cfg.model is not None:
            raise NotImplementedError("RNAstructure energy model configuration not supported")

    def check_subopt_cfg(self, cfg: SuboptCfg) -> None:
        if cfg.delta is None:
            raise NotImplementedError("RNAstructure does not support subopt without delta")
        if cfg.strucs is not None:
            raise NotImplementedError(
                "RNAstructure does not support subopt with max number of structures"
            )
        if cfg.time_secs is not None:
            raise NotImplementedError("RNAstructure does not support subopt with max time")
        if cfg.algorithm is not None:
            raise NotImplementedError("RNAstructure does not support custom subopt algorithm")

    @override
    def name(self) -> str:
        return "RNAstructure"

    @override
    def efn(self, rna: Rna, cfg: EnergyCfg) -> tuple[Decimal, CmdResult]:
        self.check_energy_cfg(cfg)
        with tempfile.NamedTemporaryFile("w") as fin, tempfile.NamedTemporaryFile("r") as fout:
            fin.write(RnaParser.to_ct_file(rna))
            fin.flush()
            # Note that not giving the -s flag doesn't make it logarithmic.
            # RNAstructure 5.8 adds the logarithmic and asymmetry models together in this case.
            # RNAstructure also uses a coefficient of -6 for the number of branches, rather than
            # the fitted -9.
            res = self._run_cmd("./exe/efn2", "-s", fin.name, fout.name, stdout_to_str=True)
            match = re.search(r"[eE]nergy = (.+)", fout.read())
            if match is None:
                raise RuntimeError(f"Could not find energy in {fout.read()}")
            energy = Decimal(match.group(1))
        return energy, res

    @override
    def fold(self, rna: Rna, cfg: EnergyCfg) -> tuple[Rna, CmdResult]:
        self.check_energy_cfg(cfg)
        with tempfile.NamedTemporaryFile("w") as fin, tempfile.NamedTemporaryFile("r") as fout:
            fin.write(RnaParser.to_seq_file(rna))
            fin.flush()
            res = self._run_cmd("./exe/Fold", "-mfe", fin.name, fout.name, stdout_to_str=True)
            predicted = RnaParser.from_any_file(fout.read())
        return predicted, res

    @override
    def partition(self, _: Rna, cfg: EnergyCfg) -> None:
        self.check_energy_cfg(cfg)
        raise NotImplementedError

    @override
    def subopt(
        self, rna: Rna, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg
    ) -> tuple[list[Rna] | int, CmdResult]:
        self.check_energy_cfg(energy_cfg)
        self.check_subopt_cfg(subopt_cfg)
        with tempfile.NamedTemporaryFile("w") as fin, tempfile.NamedTemporaryFile("r") as fout:
            fin.write(RnaParser.to_seq_file(rna))
            fin.flush()

            if subopt_cfg.delta is None:
                raise ValueError("SuboptCfg.delta must be set")
            res = self._run_cmd(
                "./exe/AllSub",
                "-a",
                f"{subopt_cfg.delta}",
                fin.name,
                fout.name,
                stdout_to_str=False,
            )

            if subopt_cfg.count_only:
                count = fast_linecount(Path(fout.name))
                divisor = (
                    len(rna) + 1
                )  # ct file has one line per nucleotide, plus one for the header
                assert (
                    count % divisor == 0
                ), f"Expected count to be a multiple of {divisor}, got {count}"
                return count // divisor, res

            output = fout.read()
            # TODO(3): does not extract energy yet
            subopts = RnaParser.multi_from_ct_file(output)
        return subopts, res

    def __str__(self) -> str:
        return self.name()
