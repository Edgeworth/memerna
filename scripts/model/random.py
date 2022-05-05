import random


class RandomRna:
    @staticmethod
    def primary(length: int) -> str:
        """Uniform random sampling of primary structures."""
        return "".join(random.choice("ACGU") for _ in range(length))
