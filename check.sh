#!/bin/sh

set -xe  # exit on failure and print commands

# Check various build configurations and run tests. Good to run before push.
# Uses :::+ to match mpfr-rnastructure and no-mpfr-no-rnastructure
parallel --progress --jobs 8 ./build.py --test {} ">" /dev/null ::: \
  --type=debug --type=relwithdebinfo ::: --sanitizer=asan \
  --sanitizer=tsan --sanitizer=ubsan ::: \
  --rnastructure --no-rnastructure :::+ --mpfr --no-mpfr ::: \
  --compiler=clang
