from dataclasses import dataclass


@dataclass
class RnaPipelineConfig:
    struc_len: int = 128
    """Length of the primary structures used to generate training data"""

    max_seq_len: int = 32
    """Maximum length of the sequences used to train the transformer"""

    batch_size: int = 64
    """Number of sequences to train on at a time"""

    unk_proportion: float = 0.2
    """Proportion of unknown bases to generate in the training data"""

    d_emb: int = 16
    """Dimension of the embedding used for the transformer"""

    nhead: int = 8
    """Number of heads used for the transformer"""

    num_encoder_layers: int = 6
    """Number of layers in the encoder"""

    num_decoder_layers: int = 6
    """Number of layers in the decoder"""

    dim_feedforward: int = 2048
    """Dimension of the feedforward layers in the transformer"""
