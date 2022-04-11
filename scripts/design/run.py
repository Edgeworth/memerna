from pathlib import Path
from torchvision import datasets
from torchvision.transforms import ToTensor
import click
from scripts.design.trainer import Trainer


def run_design(output_path: Path) -> None:
    click.echo(f"Running design at {output_path}")
    # Download training data from open datasets.
    train_data = datasets.FashionMNIST(
        root=output_path,
        train=True,
        download=True,
        transform=ToTensor(),
    )

    # Download validaiton data from open datasets.
    valid_data = datasets.FashionMNIST(
        root=output_path,
        train=False,
        download=True,
        transform=ToTensor(),
    )

    trainer = Trainer(train_data, valid_data, output_path)
    trainer.run(5)
