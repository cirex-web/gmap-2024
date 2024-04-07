static char rcsid[] = "$Id: cpuid.c 226315 2023-02-28 18:22:20Z twu $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "cpuid.h"
#include <stdio.h>
#include <string.h>
#include <sys/utsname.h>


#if !defined(__x86_64__) && !defined(__i386__)
/* Non-Intel machine; checking to see if it is an Arm machine */

void
CPUID_support (bool *arm_support_p,
	       bool *sse2_support_p, bool *ssse3_support_p, bool *sse4_1_support_p, bool *sse4_2_support_p,
	       bool *avx2_support_p, bool *avx512_support_p, bool *avx512bw_support_p) {
  struct utsname utsname;
  int rc;

  if ((rc = uname(&utsname)) != 0) {
    fprintf(stderr,"System call to uname fails with rc %d\n",rc);

  } else if (!strncmp(utsname.machine,"arm64",5)) {
    *arm_support_p = true;

  } else {
    *arm_support_p = false;
  }

  *sse2_support_p = false;
  *ssse3_support_p = false;
  *sse4_1_support_p = false;
  *sse4_2_support_p = false;
  *avx2_support_p = false;
  *avx512_support_p = false;
  *avx512bw_support_p = false;

  return;
}


#elif defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 1300)

#include <immintrin.h>

void
CPUID_support (bool *arm_support_p,
	       bool *sse2_support_p, bool *ssse3_support_p, bool *sse4_1_support_p, bool *sse4_2_support_p,
	       bool *avx2_support_p, bool *avx512_support_p, bool *avx512bw_support_p) {

  *arm_support_p = false;

  *sse2_support_p = _may_i_use_cpu_feature(_FEATURE_SSE2);
  *ssse3_support_p = _may_i_use_cpu_feature(_FEATURE_SSSE3);
  *sse4_1_support_p = _may_i_use_cpu_feature(_FEATURE_SSE4_1);
  *sse4_2_support_p = _may_i_use_cpu_feature(_FEATURE_SSE4_2);
  *avx2_support_p = _may_i_use_cpu_feature(_FEATURE_AVX2 | _FEATURE_FMA | _FEATURE_BMI | _FEATURE_LZCNT | _FEATURE_MOVBE);
  *avx512_support_p = _may_i_use_cpu_feature(_FEATURE_AVX512F | _FEATURE_AVX512CD);
  *avx512bw_support_p = _may_i_use_cpu_feature(_FEATURE_AVX512BW | _FEATURE_AVX512VL);

  return;
}


#else  /* non-Intel compiler */

#define EAX 0
#define EBX 1
#define ECX 2
#define EDX 3

#include <stdint.h>
#if defined(_MSC_VER)
#include <intrin.h>
#endif

static void
run_cpuid (uint32_t eax, uint32_t ecx, uint32_t *abcd) {
#if defined(_MSC_VER)
  __cpuidex(abcd, eax, ecx);
#else
  uint32_t ebx, edx;
#if defined(__i386__) && defined(__PIC__)
  /* in case of PIC under 32-bit, EBX cannot be clobbered */
  __asm__ ("movl %%ebx, %%edi \n\t cpuid \n\t xchgl %%ebx, %%edi"
	   : "=D" (ebx), "+a" (eax), "+c" (ecx), "=d" (edx));
#else
  __asm__ ("cpuid" : "+b" (ebx), "+a" (eax), "+c" (ecx), "=d" (edx));
#endif
  abcd[EAX] = eax, abcd[EBX] = ebx, abcd[ECX] = ecx, abcd[EDX] = edx;
#endif
}

static int
check_xcr0_ymm () {
  uint32_t xcr0;
  uint32_t ymm_xmm = ((1 << 2) | (1 << 1));
#if defined(_MSC_VER)
  xcr0 = (uint32_t)_xgetbv(0);	/* min VS2010 SP1 compiler is required */
#else
  __asm__ ("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx");
#endif
  return ((xcr0 & ymm_xmm) == ymm_xmm);	/* checking if xmm and ymm state are enabled in XCR0 */
}

static int
check_xcr0_zmm () {
  uint32_t xcr0;
  uint32_t zmm_ymm_xmm = ((7 << 5) | (1 << 2) | (1 << 1));
#if defined(_MSC_VER)
  xcr0 = (uint32_t)_xgetbv(0);	/* min VS2010 SP1 compiler is required */
#else
  __asm__ ("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx");
#endif
  return ((xcr0 & zmm_ymm_xmm) == zmm_ymm_xmm);
}


void
CPUID_support (bool *arm_support_p,
	       bool *sse2_support_p, bool *ssse3_support_p, bool *sse4_1_support_p, bool *sse4_2_support_p,
	       bool *avx2_support_p, bool *avx512_support_p, bool *avx512bw_support_p) {
  uint32_t abcd[4];
  uint32_t sse2_mask = (1 << 26); /* edx */
  uint32_t ssse3_mask = (1 << 9); /* ecx */
  uint32_t sse4_1_mask = (1 << 19); /* ecx */
  uint32_t sse4_2_mask = (1 << 20); /* ecx */
#if 0
  uint32_t popcnt_mask = (1 << 23); /* ecx */
#endif
  uint32_t fma_movbe_osxsave_mask = ((1 << 12) | (1 << 22) | (1 << 27)); /* ecx */
  uint32_t avx2_bmi12_mask = ((1 << 5) | (1 << 3) | (1 << 8)); /* ebx */
  uint32_t lzcnt_mask = (1 << 5); /* ecx */

  uint32_t osxsave_mask = (1 << 27); /* ecx */
  uint32_t avx512_mask = ((1 << 16) | (1 << 28));  /* ebx */
  uint32_t avx512bw_vl_mask = ((1 << 30) | (1 << 31));  /* ebx */


  *arm_support_p = false;

  run_cpuid(1, 0, abcd);
#ifdef MAIN
  printf("CPUID          1, 0 returns EAX %08X EBX %08X ECX %08X EDX %08X\n",abcd[EAX],abcd[EBX],abcd[ECX],abcd[EDX]);
#endif

  *sse2_support_p = ((abcd[EDX] & sse2_mask) == sse2_mask) ? true : false;
  *ssse3_support_p = ((abcd[ECX] & ssse3_mask) == ssse3_mask) ? true : false;
  *sse4_1_support_p = ((abcd[ECX] & sse4_1_mask) == sse4_1_mask) ? true : false;
  *sse4_2_support_p = ((abcd[ECX] & sse4_2_mask) == sse4_2_mask) ? true : false;
  /* *popcnt_support_p = ((abcd[ECX] & popcnt_mask) == popcnt_mask); */

  if ((abcd[ECX] & fma_movbe_osxsave_mask) != fma_movbe_osxsave_mask) {
    *avx2_support_p = false;
  } else if (!check_xcr0_ymm()) {
    *avx2_support_p = false;
  } else {
    run_cpuid(7, 0, abcd);
#ifdef MAIN
    printf("CPUID          7, 0 returns EAX %08X EBX %08X ECX %08X EDX %08X\n",abcd[EAX],abcd[EBX],abcd[ECX],abcd[EDX]);
#endif

    if ((abcd[EBX] & avx2_bmi12_mask) != avx2_bmi12_mask) {
      *avx2_support_p = false;
    } else {
      run_cpuid(0x80000001, 0, abcd);
#ifdef MAIN
      printf("CPUID 0x80000001, 0 returns EAX %08X EBX %08X ECX %08X EDX %08X\n",abcd[EAX],abcd[EBX],abcd[ECX],abcd[EDX]);
#endif

      if ((abcd[ECX] & lzcnt_mask) != lzcnt_mask) {
	*avx2_support_p = false;
      } else {
	*avx2_support_p = true;
      }
    }
  }

  run_cpuid(1, 0, abcd);
  if ((abcd[ECX] & osxsave_mask) != osxsave_mask) {
    *avx512_support_p = false;
    *avx512bw_support_p = false;
  } else if (!check_xcr0_zmm()) {
    *avx512_support_p = false;
    *avx512bw_support_p = false;
  } else {
    if ((abcd[EBX] & avx512_mask) != avx512_mask) {
      *avx512_support_p = true;	/* Should be false, but book/Web examples skip this check */
    } else {
      *avx512_support_p = true;
    }
    *avx512bw_support_p = ((abcd[EBX] & avx512bw_vl_mask) == avx512bw_vl_mask) ? true : false;
  }

  return;
}

#endif	/* non-Intel compiler */


#ifdef MAIN
int
main (int argc, char *argv[]) {
  (void)(argc);
  (void)(argv);
  
  bool arm_support_p;
  bool sse2_support_p;
  bool ssse3_support_p;
  bool sse4_1_support_p;
  bool sse4_2_support_p;
  bool avx2_support_p;
  bool avx512_support_p;
  bool avx512bw_support_p;


  CPUID_support(&arm_support_p,
		&sse2_support_p,&ssse3_support_p,&sse4_1_support_p,&sse4_2_support_p,
		&avx2_support_p,&avx512_support_p,&avx512bw_support_p);

  printf("arm support: %d\n",arm_support_p);
  printf("sse2 support: %d\n",sse2_support_p);
  printf("ssse3 support: %d\n",ssse3_support_p);
  printf("sse4.1 support: %d\n",sse4_1_support_p);
  printf("sse4.2 support: %d\n",sse4_2_support_p);
  printf("avx2 support: %d\n",avx2_support_p);
  printf("avx512 support: %d\n",avx512_support_p);
  printf("avx512bw support: %d\n",avx512bw_support_p);

  return 0;
}
#endif
