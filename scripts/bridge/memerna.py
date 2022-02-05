from pathlib import Path

from scripts.bridge.rnapackage import RnaPackage
from scripts.model.config import EnergyCfg, SuboptCfg
from scripts.model.parse import db_to_secondary
from scripts.model.rna import Rna
from scripts.util.command import run_cmd


class MemeRna(RnaPackage):
    def __init__(self, path: Path):
        self.path = path

    def energy_cfg_args(self, cfg: EnergyCfg):
        args = []
        if cfg.lonely_pairs:
            args.append("--lonely-pairs")
        else:
            args.append("--no-lonely-pairs")
        args += ["--ctd", f"{cfg.ctd}"]
        return args

    def efn(self, rna: Rna, cfg: EnergyCfg):
        args = self.energy_cfg_args(cfg)
        res = run_cmd("efn", *args, rna.r, rna.db(), cwd=self.path)
        energy = float(res.stdout.splitlines()[0].strip()) / 10.0
        return energy, res

    def fold(self, rna: Rna, cfg: EnergyCfg):
        args = self.energy_cfg_args(energy_cfg)
        res = run_cmd("fold", *args, rna.r, cwd=self.path)
        _, db, _ = res.stdout.strip().split("\n")
        return Rna(rna.name, rna.r, db_to_secondary(db)), res

    def partition(self, rna: Rna, cfg: EnergyCfg):
        args = self.energy_cfg_args(energy_cfg)
        raise NotImplementedError

    def subopt(self, rna: Rna, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg):
        args = self.energy_cfg_args(energy_cfg)
        if delta >= 0:
            args += ["-delta", str(delta)]
        if max_num >= 0:
            args += ["-num", str(max_num)]
        assert len(args) > 0
        with tempfile.NamedTemporaryFile("r") as out:
            res = run_cmd(
                "subopt",
                *args,
                rna.seq,
                record_stdout=out.name,
                limits=limits,
                cwd=self.path,
            )
            if num_only:
                res2 = run_command("wc", "-l", out.name, record_stdout=True)
                num_lines = int(res2.stdout.strip().split(" ")[0])
                retval = num_lines - 1
            else:
                assert not self.ctd_output
                retval = []
                output = out.read()
                for i in output.splitlines()[:-1]:
                    energy, db = re.split(r"\s+", i.strip())
                    retval.append(
                        (float(energy) / 10.0, Rna.from_name_seq_db(rna.name, rna.seq, db)),
                    )
        return retval, res
