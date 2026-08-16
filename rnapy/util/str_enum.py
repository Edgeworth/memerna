# Copyright 2026 Eliot Courtney.
from enum import Enum


class StrEnum(str, Enum):
    def __str__(self) -> str:
        return self.value
