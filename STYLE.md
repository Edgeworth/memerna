Some notes on the coding style used in memerna:

1. Prefer composition and data oriented APIs. Use inheritance carefully, if at
   all.
2. Prefer pass by value and move to pass by const-ref and copy.
3. Prefer non-const pointer for output parameters, const-ref for params,
   const-pointer for optional references.
