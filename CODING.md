# Style notes

Set up pre-commit.

## C++

- C++20
- Format according to the .clang-format
- Prefer composition and data oriented APIs. Use inheritance carefully, if at
  all.
- Prefer non-const pointer for output parameters, const-ref for params,
  const-pointer for optional references.

## Python

- Python 3.10
- Format according to black -l 100

## Directories

- cmake: CMake scripts
- extern: external projects and data
- data: energy model data for memerna
- examples: various dot-bracket example folded RNAs
- scripts: scripts for various things
- src: source
- tests: tests

### memerna API notes

There are two conceptual APIs:

- A high level API, using Context, which makes it easy to do regular RNA
  operations, such as MFE folding with the Turner 2004 model. These
  bundle the outputs into \*Result structs.
- A low level API, where each thing does one thing. These take their inputs as
  basic building block data structures and shouldn't take high level API types.

If a submodule needs some configuration, it should declare a \*Cfg struct, e.g.
CtxCfg, EnergyCfg, etc. If it wants to be used from the command line it should
also declare Opt options and a FromArgParse static method on the config struct,
in a config.h and config.cpp file.

In general, Cfg structs should not include other Cfg structs. If code needs
multiple Cfg structs, take them separately. This prevents
structs from being included multiple times and potentially conflicting. Maybe
an exception if it's a leaf Cfg struct that won't be included in anything else.

config.h should declare a RegisterOpts method which registers options for
that Cfg struct. It should also call RegisterOpts on any other modules it
would need, if it wanted to be invoked on the command line. e.g. Subopt
should call RegisterOpts for energy and model.

Some notes on the code structure:

## rnapy directory

Contains Python code. Organised as follows:

- analysis: Library for analysing data generated from RNA programs.
- bridge: Library for running various RNA programs.
- build: Library for building and running memerna programs.
- data: Library for getting and manipulating RNA data.
- programs: Programs invoking various scripts functionality.
- util: Misc common utilities.

## tests directory

Structure should mirror the source directory.

## Releasing

Update version numbers in CMakeLists.txt, pyproject.toml, and util/version.h.

Create a new branch from main - branch has major.minor version number.

```bash
git checkout -b release/0.2
```

Add a tag for the release, using the format `release-tag/0.2.0` - it contains
the full version number, including patch.

```bash
git tag release-tag/0.2.0
git push origin release/0.2
git push origin release-tag/0.2.0
```

Then create a release on GitHub, using the tag `release-tag/0.2.0` and the
title "v0.2.0" and description "Release 0.2.0".
