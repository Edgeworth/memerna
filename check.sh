#!/bin/sh

set -xe  # exit on first failure and print commands

# Check various build configurations and run tests. Good to run before push.
for type in debug relwithdebinfo
do
  ./build.py --type $type
  ./build.py --type $type --use-rnastructure
  ./build.py --type $type --use-mpfr
  ./build.py --type $type --use-rnastructure --use-mpfr
  ./build.py --type $type --use-clang
  ./build.py --type $type --use-clang --use-rnastructure
  ./build.py --type $type --use-clang --use-mpfr
  ./build.py --type $type --use-clang --use-rnastructure --use-mpfr
done
