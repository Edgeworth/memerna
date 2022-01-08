#!/bin/sh

set -xe  # exit on failure and print commands

# Check various build configurations and run tests. Good to run before push.
parallel --progress --load 100% ./build.py --type {} ">" /dev/null ::: \
  debug relwithdebinfo ::: '' --use-rnastructure ::: '' --use-mpfr ::: '' --use-clang
