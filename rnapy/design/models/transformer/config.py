from dataclasses import dataclass


@dataclass
class MLMTransformerConfig:
    mask_prop: float = 0.15
    """Proportion of masked bases to generate in the training data"""

    mask_correct_prop: float = 0.15
    """Proportion of masked bases to not mask and just pass through"""

    mask_incorrect_prop: float = 0.15
    """Proportion of masked bases to randomly change in the input"""

    d_emb: int = 512
    """Dimension of the embedding used for the transformer"""

    nhead: int = 8
    """Number of heads used for the transformer"""

    num_encoder_layers: int = 12
    """Number of layers in the encoder"""

    num_decoder_layers: int = 12
    """Number of layers in the decoder"""
