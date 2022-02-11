def human_size(b: int, binary: bool = True) -> str:
    def fmt(f: float) -> str:
        return (f"{f:.2f}").rstrip("0").rstrip(".")

    units = ["B", "KiB", "MiB", "GiB"]
    base = 1024
    if not binary:
        units = ["B", "KB", "MB", "GB"]
        base = 1000

    cur = float(b)
    for unit in units[:-1]:
        if abs(cur) < base:
            return f"{fmt(cur)} {unit}"
        cur /= base
    return f"{fmt(cur)} {units[-1]}"
