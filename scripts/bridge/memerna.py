from pathlib import Path

from scripts.bridge.rnapackage import RnaPackage
from scripts.model.config import EnergyCfg, SuboptCfg
from scripts.model.rna import Rna


class MemeRna(RnaPackage):
    def __init__(self, loc: Path):
        self.loc = loc

    def efn(self, rna: Rna):
        prev_dir = os.getcwd()
        os.chdir(self.loc)
        res = run_command(os.path.join("efn"), rna.seq, rna.db(), record_stdout=True)
        os.chdir(prev_dir)
        energy = float(res.stdout.strip().split(" ")[1]) / 10.0
        return energy, res

    def fold(self, rna: Rna, cfg: EnergyCfg):
        prev_dir = os.getcwd()
        os.chdir(self.loc)
        res = run_command(os.path.join("fold"), rna.seq, record_stdout=True)
        os.chdir(prev_dir)
        _, db, _ = res.stdout.strip().split("\n")
        predicted = Rna.from_name_seq_db(rna.name, rna.seq, db.strip())
        return predicted, res

    def subopt(self, rna: Rna, energy_cfg: EnergyCfg, subopt_cfg: SuboptCfg):
        maxdelta_args = []
        if delta >= 0:
            maxdelta_args += ["-delta", str(delta)]
        if max_num >= 0:
            maxdelta_args += ["-num", str(max_num)]
        assert len(maxdelta_args) > 0
        with tempfile.NamedTemporaryFile("r") as out:
            prev_dir = os.getcwd()
            os.chdir(self.loc)
            res = try_command(
                os.path.join("subopt"),
                *self.subopt_args,
                *maxdelta_args,
                rna.seq,
                record_stdout=out.name,
                limits=limits,
            )
            os.chdir(prev_dir)
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

    def __str__(self):
        return self.name
