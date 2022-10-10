# Copyright 2022 Eliot Courtney.
from dataclasses import dataclass
import re
import tempfile

from rnapy.bridge.rnapackage import RnaPackage
from rnapy.model.model_cfg import CtdCfg
from rnapy.model.model_cfg import EnergyCfg
from rnapy.model.model_cfg import SuboptCfg
from rnapy.model.parse.rna_parser import RnaParser
from rnapy.model.rna import Rna
from rnapy.util.command import CmdResult


@dataclass
class RNAstructure(RnaPackage):
    def __post_init__(self) -> None:
        self.env["DATAPATH"] = str(self.path / "data_tables")

    def check_energy_cfg(self, cfg: EnergyCfg) -> None:
        if cfg.lonely_pairs:
            raise NotImplementedError("RNAstructure does not support turning on lonely pairs")
        if cfg.ctd != CtdCfg.ALL:
            raise NotImplementedError("RNAstructure does not support turning off CTDs")

    def check_subopt_cfg(self, cfg: SuboptCfg) -> None:
        if cfg.delta is None:
            raise NotImplementedError("RNAstructure does not support subopt without delta")
        if cfg.strucs is not None:
            raise NotImplementedError(
                "RNAstructure does not support subopt with max number of structures",
            )

    def name(self) -> str:
        return "RNAstructure"

    def efn(self, rna: Rna, cfg: EnergyCfg) -> tuple[float, CmdResult]:
        self.check_energy_cfg(cfg)
        with tempfile.NamedTemporaryFile("w") as fin, tempfile.NamedTemporaryFile("r") as fout:
            fin.write(RnaParser.to_ct_file(rna))
            fin.flush()
            # Note that not giving the -s flag doesn't make it logarithmic.
            # RNAstructure 5.8 adds the logarithmic and asymmetry models together in this case.
            # RNAstructure also uses a coefficient of -6 for the number of branches, rather than
            # the fitted -9.
            res = self._run_cmd("./exe/efn2", "-s", fin.name, fout.name)
            match = re.search(r"[eE]nergy = (.+)", fout.read())
            assert match is not None
            energy = float(match.group(1))
        return energy, res

    def fold(self, rna: Rna, cfg: EnergyCfg) -> tuple[Rna, CmdResult]:
        self.check_energy_cfg(cfg)
        with tempfile.NamedTemporaryFile("w") as fin, tempfile.NamedTemporaryFile("r") as fout:
            fin.write(RnaParser.to_seq_file(rna))
            fin.flush()
            res = self._run_cmd("./exe/Fold", "-mfe", fin.name, fout.name)
            predicted = RnaParser.from_any_file(fout.read())
        return predicted, res

    def partition(self, rna: Rna, cfg: EnergyCfg) -> None:
        self.check_energy_cfg(cfg)
        raise NotImplementedError

    def subopt(
        self,
        rna: Rna,
        energy_cfg: EnergyCfg,
        subopt_cfg: SuboptCfg,
    ) -> tuple[list[Rna], CmdResult]:
        self.check_energy_cfg(energy_cfg)
        self.check_subopt_cfg(subopt_cfg)
        with tempfile.NamedTemporaryFile("w") as fin, tempfile.NamedTemporaryFile("r") as fout:
            fin.write(RnaParser.to_seq_file(rna))
            fin.flush()

            assert subopt_cfg.delta is not None
            res = self._run_cmd(
                "./exe/AllSub",
                "-a",
                f"{subopt_cfg.delta / 10.0:.2f}",
                fin.name,
                fout.name,
            )
            output = fout.read()
            # TODO(3): does not extract energy yet
            subopts = RnaParser.multi_from_ct_file(output)
        return subopts, res

    def __str__(self) -> str:
        return "RNAstructure"
