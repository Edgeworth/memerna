from pathlib import Path

from rnapy.design.rna.pipeline import RnaPipeline


def run(output_path: Path, checkpoint_path: Path | None) -> None:
    pipeline = RnaPipeline(output_path=output_path, checkpoint_path=checkpoint_path)

    if checkpoint_path:
        pipeline.predict(
            """The 2021 World Snooker Championship was a professional snooker
            tournament that took place from 17 April to 3 May at the Crucible
            Theatre in Sheffield, England. It was the 45th consecutive year the
            World Snooker Championship was held at the Crucible Theatre and was
            the 15th and final ranking event of the 2020-21 snooker season.
            It was organised by the World Snooker Tour, a subsidiary of the
            World Professional Billiards and Snooker Association.""",
        )
    else:
        pipeline.train(5000)
