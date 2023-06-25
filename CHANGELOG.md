
# 2023/06

Added suboptimal folding for t22:
- Uses same Expansion struct as for randomised traceback
- Pairs are handled explicitly in the State rather than implicitly like in t04
  by looking at the DP tables (for DP_P)

Added SHAPE support for t22:
- add pf_paired and pf_unpaired vectors to the energy model - these are
  pseudofree energy penalties for each nucleotide
- Each substructure handles paired pseudofree energy for its opening pair.
- Closing pairs are the opening pair of the next substructure, and are handled
  by it.
- Unlike RNAstructure:
  - unpaired pseudofree energy is applied for special hairpins
    hairpins.
  - paired pseudofree energy is not double counted inside helices or missed for
    lonely pairs
  - paired pseudofree energy is counted across a size 1 bulge loop
