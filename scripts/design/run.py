from pathlib import Path

import click
from scripts.design.trainer import Trainer
from scripts.design.transformer.positional_encoder import PositionalEncoder
import torch
from torchvision import datasets
from torchvision.transforms import ToTensor


def run_fashion(output_path: Path) -> None:
    click.echo(f"Running fashion at {output_path}")
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


def run(output_path: Path) -> None:
    # run_fashion(output_path)
    # run_transformer(output_path)
    enc = PositionalEncoder(d_emb=4, max_seq_len=100)
    a = enc(torch.zeros(5, 3, 4))
    print(a)
