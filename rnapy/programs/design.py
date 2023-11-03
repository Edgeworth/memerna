# Copyright 2022 Eliot Courtney.
from pathlib import Path

import cloup

from rnapy.design.models.ff.simple import SimpleFF
from rnapy.design.models.transformer.mlm import MLMTransformer
from rnapy.design.rna.pipeline import RnaPipeline
from rnapy.design.rna.pipeline_cfg import RnaPipelineCfg


@cloup.command()
@cloup.option(
    "--output-path",
    default=Path.home() / "bin" / "design",
    type=cloup.Path(
        dir_okay=True,
        file_okay=False,
        exists=False,
        writable=True,
        resolve_path=True,
        path_type=Path,
    ),
)
@cloup.option(
    "--checkpoint-path",
    required=False,
    type=cloup.Path(
        dir_okay=False,
        file_okay=True,
        exists=True,
        writable=False,
        resolve_path=True,
        path_type=Path,
    ),
)
@cloup.option("--model", default="simple", type=cloup.Choice(["simple", "transformer"]))
def design(output_path: Path, checkpoint_path: Path | None, model: str) -> None:
    cfg = RnaPipelineCfg(activation="gelu")

    match model:
        case "simple":
            cfg.model_class = SimpleFF
            cfg.batch_size = 1024 * 16
        case "transformer":
            cfg.model_class = MLMTransformer
            cfg.batch_size = 128
    pl = RnaPipeline(cfg=cfg, output_path=output_path, checkpoint_path=checkpoint_path)

    if checkpoint_path:
        # print(pipeline.predict_oneshot("....((..((....))....))..."))

        # e.g. part of rd0260: GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUC
        print(pl.predict_db_only("(((((((((((.((...((((....))))..)).)))..((((..((((....))))...))"))

    else:
        pl.train(5000)
