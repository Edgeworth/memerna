# TODO

- Add support for different kinds of partition function rescaling.
  - log space
    - different logsumexp implementations
  - per nucleotide scaling
  - per nuc scaling with log-space fallback.

- need a type to represent the boltzenergy but whether it's in log space or not + the operations.
- BoltzEnergy<ScalingKind>

- BoltzEnergyBuffer<ScalingKind, Size> - buffers up to Size BoltzEnergy values, then adds them together.
  - try for 2,4,8,16 etc.

LogSumExp kinds:
- BoostLSE
- OnlineLSE - rescale when max value increases
- MaxLSE - just return the maximum value
- ComponentLSE<ShiftKind, LogKind, ExpKind, FilterKind, SummationKind>

ShiftKind:
- None
- MaxShift - compute max value and subtract from all values, compute m + log(sum(exp(x_i - m)))

LogKind:
- standard library
- TaylorApprox<N> - decompose into m + log1p(exp(x*))
- BithackTaylorApprox<N> - decompose S = 2^k * u, log(S) = k*log(2) + log1p(u-1)
- PadeApprox<N>
- BithackPadeApprox<N>
- MinimaxApprox<N>
- Lookup table idea?

ExpKind:
- standard library
- TaylorApprox<N>
- BithackTaylorApprox<N>
- PadeApprox<N>
- BithackPadeApprox<N>
- MinimaxApprox<N>

FilterKind:
- None
- Cutoff(x) - only include terms where x_i > x (should be combined with max shifting)

SummationKind:
- LinearScan
- TreeReduction

ScalingKind
- None
- Log(LSEKind)
- PerNucleotide(Factor)
- PerNucleotideWithFallback(Factor, Log)
