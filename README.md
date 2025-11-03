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
   cmake -S . -B build
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
