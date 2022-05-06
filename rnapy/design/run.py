from pathlib import Path

from rnapy.design.rna.pipeline import RnaPipeline


def run(output_path: Path, checkpoint_path: Path | None) -> None:
    pipeline = RnaPipeline(output_path=output_path, checkpoint_path=checkpoint_path)

    if checkpoint_path:
        pipeline.predict("....((..((....))....))...")
    else:
        pipeline.train(5000)
