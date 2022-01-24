Some notes on the coding style used in memerna:

# C++
- Format according to the .clang-format
- Prefer composition and data oriented APIs. Use inheritance carefully, if at
   all.
- Prefer pass by value and move to pass by const-ref and copy.
- Prefer non-const pointer for output parameters, const-ref for params,
   const-pointer for optional references.

# Python
- Format according to black -l 100
