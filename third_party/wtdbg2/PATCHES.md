# wtdbg2 — Patches for macOS ARM (Apple Silicon)

This directory contains wtdbg2 v2.5 (https://github.com/ruanjue/wtdbg2)
with the following modifications to enable building on macOS ARM (aarch64).

The original code is licensed under GPL-3.0 (see LICENSE.txt).

## Changes

### 1. SSE to ARM NEON translation
- **ksw.c, poacns.h**: Wrapped `#include <emmintrin.h>` / `<tmmintrin.h>` with
  `#ifdef __aarch64__` guards to use `sse2neon.h` on ARM platforms.
- **sse2neon.h**: Added from https://github.com/DLTcollab/sse2neon (MIT licensed).

### 2. Bug fix: undefined identifier on macOS
- **kbm.c:326**: Changed `_proc_deamon->ncpu` to `_sig_proc_deamon->ncpu`.
  The variable `_proc_deamon` is only defined inside the thread function scope;
  the file-scope pointer is `_sig_proc_deamon`.

### 3. Makefile: cross-platform support
- Removed x86-only flags (`-mpopcnt`, `-msse4.2`) on non-x86 platforms.
- Removed `-lrt` on macOS (not available on Darwin).
- Changed `-O4` to `-O3` and `-g3` to `-g` for compiler compatibility.
- Added `-std=gnu11` for gcc 15+ compatibility.
- Changed `CC :=` to `CC ?=` to allow override from parent Makefile.
