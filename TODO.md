### Zuker optimisation ideas
- list of actual base pairs (skip things)? memory
- magic data structures to reduce multiloop computation to log time
- Python code to generate C code to do energy for optimisation. - probs not useful
- Bake AUGU penalties etc into things
- Create second set of energy functions  and zuker
-- optimised versions and easy to understand versions
- Vectorisation (SSE)
- Sliding window, reduce memory usage - needs st en not sz st?

### Todo
- Lyngso optimisation?
- Speed up internal loops somehow
- American fuzzy lop?
- Download and test that sparse folding program from that paper
- Grep codebase for todos
- Fuzz tables!!
- Add choices to argparse.

### Fuzzing
- Randomised energy model fuzzing
-- Plus comparison to rnastructure / rnark
- Run on digitalocean?
- Run on max computer
