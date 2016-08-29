### Zuker optimisation ideas
- Sliding window, reduce memory usage - needs st en not sz st?

### Todo
- American fuzzy lop?
- Wutchys
- Grep codebase for todos

### Wuchty testing plan
- Implement reporting stacking for wuchty
- List of CTDs (enum), put into energy/ etc. allow specifying custom stacking for branches
- Write bruteforce version - needs to enumerate all stackings as well
-- current bruteforce runs only up to N = 25, which is at most ~11 branches 
-- will need way to run through every possible
- Compare to brute force
- Check MFE is always top result
- Check for no duplicated structures (must include stacking)
- Check efn gives same energy value for each thing
- All of above for random energy models
