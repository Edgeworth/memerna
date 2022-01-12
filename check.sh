#!/bin/sh

set -xe  # exit on failure and print commands

# Check various build configurations and run tests. Good to run before push.
parallel --progress --jobs 1 ./build.py --test --type {} ">" /dev/null ::: \
  debug relwithdebinfo ::: '' --use-rnastructure ::: '' --use-mpfr ::: '' --use-clang
