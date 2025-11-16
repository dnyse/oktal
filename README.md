# AdvPT Programming Project - Oktal

## Building the Project

### Prerequisites
- CMake 3.31 or higher
- C++23 compatible compiler.

### Compilation Steps

1. **Navigate to the project root:**
```bash
   cd oktal
```

2. **Generate the build system:**
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_CXX_COMPILER=clang++-19
```

3. **Build the project:**
```bash
   cmake --build build -j8
```

### Running Tests

After building, run the test suite:
```bash
cd build
ctest
```

To run specific tests by name:
```bash
ctest -R <test-name-pattern>
```

### Build Options

- To disable the test suite: `cmake -S . -B build -DOKTAL_BUILD_TESTSUITE=OFF`

### Cleaning

- **Clean build artifacts:** `cmake --build build --target clean`
- **Clean and rebuild:** `cmake --build build --clean-first`
- **Fresh configuration:** `cmake --fresh -S . -B build`

### Linter
Run clang-tid:
```bash
cd build
run-clang-tidy-19 -p .
```

## Install dependecies via Spack 
Environments can be created from the manifest file (the user input), the lockfile, or the entire environment at once. Create an environment from a manifest using:
```bash
spack env create oktal ./.spack-env/spack.yaml
# or
spack env create oktal ./.spack-env/spack.lock
```
To activate an environment, use the following command:
```bash
spack env activate oktal
```
You can add packages to the environment with
```bash
spack add <package>
spack concretize
```
In addition to adding individual specs to an environment, one can install the entire environment at once using the command
```bash
spack install
```
### Independent Environments
Independent environments can be located in any directory outside of Spack.
```bash
spack env activate --create ./spack-env
```
### Environment Views
Spack Environments can have an associated filesystem view, which is a directory with a more traditional structure <view>/bin, <view>/lib, <view>/include in which all files of the installed packages are linked.
```bash
.spack-env
├── spack.lock
├── spack.yaml
└── view -> /home/$USER/University/EUMaster4HPC/ws2025-group-45/.spack-env/._view/ttdiwszqvfbboietfm7h7tjuqqlfkxka
```
The Spack docs can be found [here](https://spack.readthedocs.io/en/latest/environments.html#independent-environments).
