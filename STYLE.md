Some notes on the style used in memerna:

# C++
- Format according to the .clang-format
- Prefer composition and data oriented APIs. Use inheritance carefully, if at
   all.
- Prefer pass by value and move to pass by const-ref and copy.
- Prefer non-const pointer for output parameters, const-ref for params,
   const-pointer for optional references.

# Python
- Format according to black -l 100

# Memerna
There are two conceptual APIs:
- A high level API, using Context, which makes it easy to do regular RNA
  operations, such as MFE folding with the Turner 2004 model. These
  bundle the outputs into *Result structs.
- A low level API, where each thing does one thing. These take their inputs as
  basic building block data structures and shouldn't take high level API types.
