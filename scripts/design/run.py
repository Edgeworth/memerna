from pathlib import Path

import click
from scripts.design.fashion_model import FashionModel
from scripts.design.harness.config import TrainConfig
from scripts.design.harness.trainer import Trainer
from scripts.design.rna.pipeline import RnaPipeline
import torchvision
from torchvision.transforms import ToTensor


def run_fashion(output_path: Path) -> None:
    click.echo(f"Running fashion at {output_path}")
    # Download training data from open datasets.
    train_data = torchvision.datasets.FashionMNIST(
        root=output_path,
        train=True,
        download=True,
        transform=ToTensor(),
    )

    # Download validaiton data from open datasets.
    valid_data = torchvision.datasets.FashionMNIST(
        root=output_path,
        train=False,
        download=True,
        transform=ToTensor(),
    )

    cfg = TrainConfig(
        model_name="FashionModel",
        output_path=output_path,
        profile=False,
    )
    trainer = Trainer(
        model=FashionModel(),
        train_data=train_data,
        valid_data=valid_data,
        cfg=cfg,
    )
    trainer.run(5)


def run(output_path: Path, checkpoint_path: Path | None) -> None:
    # run_fashion(output_path)
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
