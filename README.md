# Molecular Dynamics Microframework

[![Release][release-badge]][release-url]
[![Build Status][build-badge]][build-url]
[![Boost License][license-badge]][license-url]
![C++14,17,20][cxx-badge]

Header-only molecular dynamics microframework for C++.

[cxx-badge]: https://img.shields.io/badge/C%2B%2B-14%2F17%2F20-orange.svg
[license-badge]: https://img.shields.io/badge/license-Boost-blue.svg
[license-url]: https://github.com/snsinfu/micromd/blob/master/LICENSE.txt
[build-badge]: https://github.com/snsinfu/micromd/workflows/test/badge.svg
[build-url]: https://github.com/snsinfu/micromd/actions?query=workflow%3Atest
[release-badge]: https://img.shields.io/github/release/snsinfu/micromd.svg
[release-url]: https://github.com/snsinfu/micromd/releases

## Install

### Option 1: Single-header build

Go to the [release page][release-url] and download `md.hpp` into your include
directory. You can then use micromd by `#include <md.hpp>`.

### Option 2: Git submodule

Clone the repository as a git submodule (change the submodule path as you like):

```
git submodule add https://github.com/snsinfu/micromd submodules/github.com/snsinfu/micromd
```

Then, add this to your g++/clang++ flags (`CXXFLAGS` if you use Makefile):

```
-isystem submodules/github.com/snsinfu/micromd/include
```

You can then use micromd by `#include <md.hpp>`.

## Test

```console
git clone https://github.com/snsinfu/micromd
cd micromd/tests
make
```

## License

Boost v1.
