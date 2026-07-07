# C-style parsing code refactor (libmbutil, libmbc)

Analysis and modernization of the parsing code in `libraries/libmbutil` and
`libraries/libmbc`, targeting error sources linked to C-style char scanning
and copying functions and fixed-length `char*` buffers, converted to modern
C++ (`std::string`) where applicable.

## Sources of error identified

### libmbutil (C++ parsing code)

1. **`LowParser::sCurrWordBuf`** (`parser.cc`) — a manually doubled `char*`
   buffer grown with `SAFENEWARR`/`memcpy`. On EOF in mid-word, `PackWords()`
   returned *without null-terminating*, so `sGetWord()` could return stale
   content from a previous, longer word.
2. **`HighParser::sStringBuf[8192]` / `sStringBufWithSpaces[8192]`**
   (`parser.h`) — fixed buffers behind `GetString()`, `GetStringWithDelims()`
   and `ParseWord()`:
   - `GetString()` and `ParseWord()` *silently truncated* input beyond
     8 KB (`BUFSIZ`);
   - `GetStringWithDelims()` threw a hard "End-of-buffer" error for long
     strings;
   - the `,`/`;` "null string" path of `GetStringWithDelims()` returned the
     *unmodified previous buffer contents* (stale garbage) instead of the
     documented null result.
3. **`IncludeParser`** (`parsinc.cc`) — `sCurrPath`/`sInitialPath`/`sCurrFile`
   were raw `char*` managed with `SAFESTRDUP`/`SAFEDELETEARR`, with
   `sInitialPath` *aliasing* the first `sCurrPath` allocation: after
   `Close()`, `sInitialPath` dangled; in non-`USE_INCLUDE_PARSER` builds
   `GetLineData()` called `strcmp(NULL, NULL)` (undefined behavior).
   `GetFileName()` failed outright for resolved paths ≥ 8192 characters.
4. **`expand_environment()` / `resolve_filename()`** (`parsinc.cc`) — fixed
   `subst[10]` substitution array (hard failure on more than 8 `$VAR`
   substitutions in a file name) plus hand-computed `memcpy` length
   arithmetic; `resolve_filename()` used `strcpy` concatenation and a fragile
   pointer-range test to decide buffer ownership.
5. **`MathParser`** (`mathp.cc`):
   - the number-token scanner used `char s[BUFSIZ]` with an artificial
     "value too long" failure;
   - `readplugin()` (both the expression-evaluator and the `DO_NOT_USE_EE`
     copies) accumulated arguments in `char buf[BUFSIZ]` (hard
     "buffer overflow" error) and leaked the strdup'ed `argv` strings on
     error paths;
   - the `DO_NOT_USE_EE` copy of `readplugin()` referenced an undeclared
     variable `c` instead of `cIn` — that build variant **did not compile**
     (verified against the pristine tree);
   - `trim_arg()` read past index 0 on an all-whitespace string.

### libmbc (C coupling API)

`mbc.c`/`sock.c` are part of the public **C** API, so `std::string` is not
applicable; issues were fixed in C idiom instead:

6. **`winsock_err_string(int, char*)`** (`sock.c`) — `strcpy` into an unsized
   caller buffer with **no `default` case**: an unknown error code left the
   caller's `char msg[100]` uninitialized, which was then printed.
7. **`sock_err_string()`** (`sock.c`) — guarded by `#ifdef _win32`
   (lowercase, never defined), so Windows got `strerror()` for winsock error
   codes; the intended branch returned a pointer to a **stack-local buffer**
   (dangling pointer).

## Changes

- `parser.h` / `parser.cc`: `LowParser` word buffer and both `HighParser`
  string buffers are now `std::string` — no length limits, no silent
  truncation, always terminated. Public `const char*` return types are
  preserved (pointers into the member strings), so none of the external call
  sites change.
- `parsinc.h` / `parsinc.cc`: path/file members and `MyInput` use
  `std::string` (fixes the `sInitialPath` aliasing dangle);
  `expand_environment()`/`resolve_filename()` rewritten around `std::string`
  with no substitution limit; `GetFileName()` no longer has a size cap.
- `mathp.h` / `mathp.cc`: number scanner and both `readplugin()` copies use
  `std::string`/`std::vector<std::string>` (a small `std::vector<char*>` shim
  preserves the legacy `PlugIn::Read(argc, argv)` interface); `trim_arg()`
  takes `std::string&`; the `c`/`cIn` bug and the error-path leaks are gone.
- `sock.h` / `sock.c`, `mbc.c` and call sites (`usesock.cc`, `dataman6.cc`,
  `socketstreamdrive.cc`, `strext.cc`): `sock_err_string()` and
  `winsock_err_string()` now return `const char*` string literals (with an
  "Unknown winsock error" default) — no copying, no caller buffers.

### Deliberate behavior change

`GetStringWithDelims()` on an empty argument (`,`/`;`) now returns an empty
string instead of stale buffer contents.

## Verification

- Every modified file syntax-checked with `g++ -std=c++20` / `gcc` against a
  stub `mbconfig.h` (unconfigured tree), including the `-DDEBUG` and
  `-DDO_NOT_USE_EE` variants; the pristine tree compiles identically except
  for the pre-existing `DO_NOT_USE_EE` breakage, which is now fixed.
- An end-to-end test linking the real modified sources exercised: C-style
  comments, hex escapes (`\41`), keywords, ints/reals, a 20 000-character
  quoted string (fatal before), and an `include: "$TESTDIR/sub.mbd"` whose
  `set:` statement is visible in the symbol table afterwards — all pass,
  clean under AddressSanitizer + UBSan.

## Possible follow-ups (out of scope)

`fn_UNIX.cc`, `crypt.cc`, `mbstrbuf.cc` and the `SAFESTRDUP` macros in
`mynewmem.h` still use C-style copying, though they are not parsing code.
