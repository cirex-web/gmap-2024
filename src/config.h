/* src/config.h.  Generated from config.h.in by configure.  */
/* src/config.h.in.  Generated from configure.ac by autoheader.  */

/* Define if building universal (internal helper macro) */
/* #undef AC_APPLE_UNIVERSAL_BUILD */

/* Define to 1 if host is Power8 */
/* #undef AX_HOST_POWER8 */

/* Define to 1 if using 'alloca.c'. */
/* #undef C_ALLOCA */

/* Define to 1 if you have 'alloca', as a function or macro. */
/* #undef HAVE_ALLOCA */

/* Define to 1 if <alloca.h> works. */
/* #undef HAVE_ALLOCA_H */

/* Define to 1 if bsr command is available in assembly */
/* #undef HAVE_ASM_BSR */

/* Define to 1 if __builtin_clz works. */
/* #undef HAVE_BUILTIN_CLZ */

/* Define to 1 if __builtin_clzll works. */
/* #undef HAVE_BUILTIN_CLZLL */

/* Define to 1 if __builtin_ctz works. */
/* #undef HAVE_BUILTIN_CTZ */

/* Define to 1 if __builtin_ctzll works. */
/* #undef HAVE_BUILTIN_CTZLL */

/* Define to 1 if __builtin_popcount works. */
/* #undef HAVE_BUILTIN_POPCOUNT */

/* Define to 1 if you have a working bzlib library. */
#define HAVE_BZLIB 1

/* Define to 1 if the system has the type `caddr_t'. */
#define HAVE_CADDR_T 1

/* Define to 1 if you have the `ceil' function. */
#define HAVE_CEIL 1

/* Define to 1 if you have the <dirent.h> header file, and it defines `DIR'.
   */
#define HAVE_DIRENT_H 1

/* Define to 1 if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H 1

/* Define to 1 if you have the `floor' function. */
#define HAVE_FLOOR 1

/* Define to 1 if fseeko (and presumably ftello) exists and is declared. */
#define HAVE_FSEEKO 1

/* Define to 1 if you have the `index' function. */
#define HAVE_INDEX 1

/* Define to 1 if compiler supports extern inline */
#define HAVE_INLINE 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the <limits.h> header file. */
#define HAVE_LIMITS_H 1

/* Define to 1 if you have the `log' function. */
#define HAVE_LOG 1

/* Define to 1 if you support Intel intrinsic _lzcnt instruction */
/* #undef HAVE_LZCNT */

/* Define to 1 if you have the `madvise' function. */
#define HAVE_MADVISE 1

/* Define to 1 if MADV_DONTNEED available for madvise. */
#define HAVE_MADVISE_MADV_DONTNEED 1

/* Define to 1 if MADV_RANDOM available for madvise */
#define HAVE_MADVISE_MADV_RANDOM 1

/* Define to 1 if MADV_SEQUENTIAL available for madvise */
#define HAVE_MADVISE_MADV_SEQUENTIAL 1

/* Define to 1 if MADV_WILLNEED available for madvise. */
#define HAVE_MADVISE_MADV_WILLNEED 1

/* Define to 1 if you have the `memcpy' function. */
#define HAVE_MEMCPY 1

/* Define to 1 if you have the `memmove' function. */
#define HAVE_MEMMOVE 1

/* Define to 1 if you have the `memset' function. */
#define HAVE_MEMSET 1

/* Define to 1 if you have a working 'mmap' system call. */
#define HAVE_MMAP 1

/* Define to 1 if MAP_FILE available for mmap. */
#define HAVE_MMAP_MAP_FILE 1

/* Define to 1 if MAP_PRIVATE available for mmap */
#define HAVE_MMAP_MAP_PRIVATE 1

/* Define to 1 if MAP_SHARED available for mmap */
#define HAVE_MMAP_MAP_SHARED 1

/* Define to 1 if MAP_VARIABLE available for mmap. */
/* #undef HAVE_MMAP_MAP_VARIABLE */

/* Define to 1 if _mm_extract_epi64 intrinsic is available. */
/* #undef HAVE_MM_EXTRACT_EPI64 */

/* Define to 1 if you support Intel intrinsic _mm_popcnt_u32/64 instructions
   */
/* #undef HAVE_MM_POPCNT */

/* Define to 1 if _mm_popcnt_u64 intrinsic is available. */
/* #undef HAVE_MM_POPCNT_U64 */

/* Define to 1 if you have the `munmap' function. */
#define HAVE_MUNMAP 1

/* Define to 1 if you have the <ndir.h> header file, and it defines `DIR'. */
/* #undef HAVE_NDIR_H */

/* Define to 1 if you support Intel intrinsic _pext instruction */
/* #undef HAVE_PEXT */

/* Define to 1 if you support Intel intrinsic _popcnt instruction */
/* #undef HAVE_POPCNT */

/* Define to 1 if you have the `pow' function. */
#define HAVE_POW 1

/* Define if you have POSIX threads libraries and header files. */
#define HAVE_PTHREAD 1

/* Define to 1 if you have the `rint' function. */
#define HAVE_RINT 1

/* Define to 1 if you have the `semctl' function. */
#define HAVE_SEMCTL 1

/* Define to 1 if you have the `semget' function. */
#define HAVE_SEMGET 1

/* Define to 1 if you have the `semop' function. */
#define HAVE_SEMOP 1

/* Define to 1 if you have the `shmat' function. */
#define HAVE_SHMAT 1

/* Define to 1 if you have the `shmctl' function. */
#define HAVE_SHMCTL 1

/* Define to 1 if you have the `shmdt' function. */
#define HAVE_SHMDT 1

/* Define to 1 if you have the `shmget' function. */
#define HAVE_SHMGET 1

/* Define to 1 if SHM_NORESERVE available for shmget. */
/* #undef HAVE_SHM_NORESERVE */

/* Define to 1 if you have the `sigaction' function. */
#define HAVE_SIGACTION 1

/* Define to 1 if you have the `stat64' function. */
#define HAVE_STAT64 1

/* Define to 1 if you have the <stddef.h> header file. */
#define HAVE_STDDEF_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strtoul' function. */
#define HAVE_STRTOUL 1

/* have struct stat64 */
/* #undef HAVE_STRUCT_STAT64 */

/* Define to 1 if your compiler has STTNI commands. */
/* #undef HAVE_STTNI */

/* Define to 1 if you have the `sysconf' function. */
#define HAVE_SYSCONF 1

/* Define to 1 if you have the `sysctl' function. */
#define HAVE_SYSCTL 1

/* Define to 1 if you have the <sys/dir.h> header file, and it defines `DIR'.
   */
/* #undef HAVE_SYS_DIR_H */

/* Define to 1 if you have the <sys/ndir.h> header file, and it defines `DIR'.
   */
/* #undef HAVE_SYS_NDIR_H */

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you support Intel intrinsic _tzcnt instruction */
/* #undef HAVE_TZCNT */

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have a working zlib library. */
#define HAVE_ZLIB 1

/* Define to 1 if your zlib library has a gzbuffer function. */
#define HAVE_ZLIB_GZBUFFER 1

/* Define MAP_FAILED here if not available otherwise. */
/* #undef MAP_FAILED */

/* Define to 1 if your C compiler doesn't accept -c and -o together. */
/* #undef NO_MINUS_C_MINUS_O */

/* Name of package */
#define PACKAGE "gmap"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "Thomas Wu <twu@gene.com>"

/* Define to the full name of this package. */
#define PACKAGE_NAME "gmap"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "gmap 2024-02-22"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "gmap"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "2024-02-22"

/* pagesize is available via sysconf */
#define PAGESIZE_VIA_SYSCONF 1

/* pagesize is available via sysctl */
#define PAGESIZE_VIA_SYSCTL 1

/* Define to necessary symbol if this constant uses a non-standard name on
   your system. */
/* #undef PTHREAD_CREATE_JOINABLE */

/* The size of `off_t', as computed by sizeof. */
#define SIZEOF_OFF_T 8

/* The size of `unsigned long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG 8

/* The size of `unsigned long long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG_LONG 8

/* Define to 1 if your compiler requires an immediate in sse2 _mm_slli_epi32
   command. */
#define SSE2_SLLI_CONST_IMM8 1

/* If using the C implementation of alloca, define if you know the
   direction of stack growth for your system; otherwise it will be
   automatically deduced at runtime.
	STACK_DIRECTION > 0 => grows toward higher addresses
	STACK_DIRECTION < 0 => grows toward lower addresses
	STACK_DIRECTION = 0 => direction of growth unknown */
/* #undef STACK_DIRECTION */

/* Define to 1 if all of the C90 standard headers exist (not just the ones
   required in a freestanding environment). This macro is provided for
   backward compatibility; new code need not use it. */
#define STDC_HEADERS 1

/* Define this if we can use the "b" mode for fopen safely. */
#define USE_FOPEN_BINARY 1

/* Define this if we can use the "t" mode for fopen safely. */
#define USE_FOPEN_TEXT 1

/* Define to 1 to use Intel SIMD intrinsics rather than SIMDe */
/* #undef USE_INTEL_INTRINSICS */

/* Version number of package */
#define VERSION "2024-02-22"

/* Define WORDS_BIGENDIAN to 1 if your processor stores words with the most
   significant byte first (like Motorola and SPARC, unlike Intel). */
#if defined AC_APPLE_UNIVERSAL_BUILD
# if defined __BIG_ENDIAN__
#  define WORDS_BIGENDIAN 1
# endif
#else
# ifndef WORDS_BIGENDIAN
/* #  undef WORDS_BIGENDIAN */
# endif
#endif

/* Number of bits in a file offset, on hosts where this is settable. */
/* #undef _FILE_OFFSET_BITS */

/* Define to 1 to make fseeko visible on some hosts (e.g. glibc 2.2). */
/* #undef _LARGEFILE_SOURCE */

/* Define for large files, on AIX-style hosts. */
/* #undef _LARGE_FILES */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to `long int' if <sys/types.h> does not define. */
/* #undef off_t */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to empty if the keyword `volatile' does not work. Warning: valid
   code using `volatile' can become incorrect without. Disable with care. */
/* #undef volatile */
