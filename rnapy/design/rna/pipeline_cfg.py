from dataclasses import dataclass

from rnapy.design.harness.model import Model
from rnapy.design.models.transformer.transformer_cfg import MLMTransformerCfg
from rnapy.design.rna.tensor import ChunkedRnaTensor, RnaTensor


@dataclass
class RnaPipelineCfg:
    train_num_struc: int = 10000
    """Number of structures to generate training data from"""

    valid_num_struc: int = 5000
    """Number of structures to generate validation data from"""

    struc_len: int = 128
    """Length of the primary structures used to generate training data"""

    max_seq_len: int = 32
    """Maximum length of the sequences used to train the model"""

    batch_size: int = 128
    """Number of sequences to train on at a time"""

    activation: str = "gelu"
    """Activation function to use in the pipeline"""

    dim_feedforward: int = 2048
    """Dimension of the feedforward layers used (e.g. controls transformer feedforward)"""

    tensor: RnaTensor = ChunkedRnaTensor(1)

    mlm: MLMTransformerCfg = MLMTransformerCfg()

    model_class: type[Model] = Model  # Needs to be set.
