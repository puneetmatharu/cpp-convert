# cpp-fpdiff

## Building

To build the code you will need [CMake](https://cmake.org/download/)

```bash
cmake -B build && cmake --build build
```

If you have [Ninja](https://ninja-build.org), you can specify it as the build generator by running the following command instead:

```bash
cmake -G Ninja -B build && cmake --build build
```

If you are on an Apple MacBook with an M-series chip, you should run the command below instead to ensure the build is optimised for the Arm architecture:

```bash
CMAKE_APPLE_SILICON_PROCESSOR="arm64" cmake -G Ninja -B build && cmake --build build
```

In the build directory, you should now see the executable `oomph_convert`.

## Benchmark

TODO: Fill this in...

## License

The files `oomph-convert.py` and `oomph_utilities.h` have been borrowed from the [`oomph-lib`](https://github.com/oomph-lib/oomph-lib) repository and thus falls under the licensing conditions of that repository (see [LICENSE](./LICENCE)). The remaining code also falls under the purview of this license.
