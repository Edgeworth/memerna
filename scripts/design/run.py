from pathlib import Path

import click
from scripts.design.fashion_model import FashionModel
from scripts.design.language import LanguagePipeline, run_transformer
from scripts.design.trainer import Trainer
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

    trainer = Trainer(
        model=FashionModel(), train_data=train_data, valid_data=valid_data, path=output_path
    )
    trainer.run(5)


def run(output_path: Path) -> None:
    # run_fashion(output_path)
    pipeline = LanguagePipeline(output_path=output_path)
    pipeline.run(5)
