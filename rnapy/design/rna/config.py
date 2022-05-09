from dataclasses import dataclass

from rnapy.design.rna.tensor import BasicRnaTensor
from rnapy.design.rna.tensor import RnaTensor


@dataclass
class RnaPipelineConfig:
    train_num_struc: int = 100000
    """Number of structures to generate training data from"""

    valid_num_struc: int = 5000
    """Number of structures to generate validation data from"""

    struc_len: int = 30
    """Length of the primary structures used to generate training data"""

    max_seq_len: int = 32
    """Maximum length of the sequences used to train the transformer"""

    batch_size: int = 32
    """Number of sequences to train on at a time"""

    mask_prop: float = 0.15
    """Proportion of masked bases to generate in the training data"""

    mask_correct_prop: float = 0.15
    """Proportion of masked bases to not mask and just pass through"""

    mask_incorrect_prop: float = 0.15
    """Proportion of masked bases to randomly change in the input"""

    activation: str = "gelu"
    """Activation function to use in the pipeline"""

    d_emb: int = 256
    """Dimension of the embedding used for the transformer"""

    nhead: int = 8
    """Number of heads used for the transformer"""

    num_encoder_layers: int = 12
    """Number of layers in the encoder"""

    num_decoder_layers: int = 12
    """Number of layers in the decoder"""

    dim_feedforward: int = 2048
    """Dimension of the feedforward layers in the transformer"""

    tensor: RnaTensor = BasicRnaTensor()