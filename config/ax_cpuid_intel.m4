AC_DEFUN([AX_CPUID_INTEL],
[AC_REQUIRE([AC_PROG_CC])
AC_LANG_PUSH([C])
AC_CHECK_HEADER([immintrin.h],[
# Test for SSE2 support
  AC_MSG_CHECKING(for sse2 support)
  AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[#include <immintrin.h>]],
                         [[return _may_i_use_cpu_feature(_FEATURE_SSE2) ? 0 : 9;]])],
        [AC_MSG_RESULT(yes)
         ax_cv_cpu_has_sse2_ext=yes],
	[AC_MSG_RESULT(no)])

# Test for SSSE3 support
  AC_MSG_CHECKING(for ssse3 support)
  AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[#include <immintrin.h>]],
                         [[return _may_i_use_cpu_feature(_FEATURE_SSSE3) ? 0 : 9;]])],
        [AC_MSG_RESULT(yes)
         ax_cv_cpu_has_ssse3_ext=yes],
	[AC_MSG_RESULT(no)])

# Test for SSE4.1 support
  AC_MSG_CHECKING(for sse4.1 support)
  AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[#include <immintrin.h>]],
                         [[return _may_i_use_cpu_feature(_FEATURE_SSE4_1) ? 0 : 9;]])],
        [AC_MSG_RESULT(yes)
         ax_cv_cpu_has_sse41_ext=yes],
	[AC_MSG_RESULT(no)])

# Test for SSE4.2 support
  AC_MSG_CHECKING(for sse4.2 support)
  AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[#include <immintrin.h>]],
                         [[return _may_i_use_cpu_feature(_FEATURE_SSE4_2) ? 0 : 9;]])],
        [AC_MSG_RESULT(yes)
         ax_cv_cpu_has_sse42_ext=yes],
	[AC_MSG_RESULT(no)])

# Test for AVX2 support
  AC_MSG_CHECKING(for avx2 support)
  AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[#include <immintrin.h>]],
                         [[return _may_i_use_cpu_feature(_FEATURE_AVX2 | _FEATURE_FMA | _FEATURE_BMI | _FEATURE_LZCNT | _FEATURE_MOVBE) ? 0 : 9;]])],
        [AC_MSG_RESULT(yes)
         ax_cv_cpu_has_avx2_ext=yes],
	[AC_MSG_RESULT(no)])

# Test for AVX512 (F and CD) support
  AC_MSG_CHECKING(for avx512 support)
  AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[#include <immintrin.h>]],
                         [[return _may_i_use_cpu_feature(_FEATURE_AVX512F | _FEATURE_AVX512CD) ? 0 : 9;]])],
        [AC_MSG_RESULT(yes)
         ax_cv_cpu_has_avx512_ext=yes],
	[AC_MSG_RESULT(no)])

# Test for AVX512BW and VL support (not in Knights Landing or Knights Mill)
  AC_MSG_CHECKING(for avx512bw and avx512vl support)
  AC_RUN_IFELSE(
	[AC_LANG_PROGRAM([[#include <immintrin.h>]],
                         [[return _may_i_use_cpu_feature(_FEATURE_AVX512BW | _FEATURE_AVX512VL) ? 0 : 9;]])],
        [AC_MSG_RESULT(yes)
         ax_cv_cpu_has_avx512bw_ext=yes],
	[AC_MSG_RESULT(no)])
],[])
AC_LANG_POP([C])
])

