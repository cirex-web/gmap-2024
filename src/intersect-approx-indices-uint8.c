static char rcsid[] = "$Id: 7ec5b6f1007e74013548f020555c35441c4842ad $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* Code is modified from the SIMDCompressionAndIntersection package,
   by Daniel Lemire, Leonid Boytsov, and Nathan Kurz, file
   intersection.cpp.  */

/**
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 *
 */

/* Minimum requirement for SIMD is SSE4_2 for _mm_cmpgt_epi64 */

#include "assert.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include "mem.h"
#include "bool.h"
#include "intersect-approx-indices-uint8.h"

#include "simd.h"

#if 0
#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#ifdef HAVE_SSE4_1
#include <smmintrin.h>
#endif
#ifdef HAVE_AVX2
#include <immintrin.h>
#endif
#endif


#define INCLUDE_EXACT 1

#define EPI16_MAX 0x8000	     /* 4 nibbles * 4 bits/nibble */
#define EPI32_MAX 0x80000000	     /* 8 nibbles * 4 bits/nibble*/
#define EPI64_MAX 0x8000000000000000 /* 16 nibbles * 4 bits/nibble */
/*                  0123456701234567 */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Intersect_approx_old */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

/* Check SIMD against non-SIMD */
#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif

/* binary search */
#ifdef DEBUG20
#define debug20(x) x
#else
#define debug20(x)
#endif



#if defined(HAVE_AVX2) && defined(DEBUG)
static void
print_vector_signed (__m256i x) {
  printf("%10d %10d %10d %10d %10d %10d %10d %10d\n",
	 _mm256_extract_epi32(x,0) + EPI32_MAX,
	 _mm256_extract_epi32(x,1) + EPI32_MAX,
	 _mm256_extract_epi32(x,2) + EPI32_MAX,
	 _mm256_extract_epi32(x,3) + EPI32_MAX,
	 _mm256_extract_epi32(x,4) + EPI32_MAX,
	 _mm256_extract_epi32(x,5) + EPI32_MAX,
	 _mm256_extract_epi32(x,6) + EPI32_MAX,
	 _mm256_extract_epi32(x,7) + EPI32_MAX);
  return;
}

static void
print_vector_unsigned (__m256i x) {
  printf("%10u %10u %10u %10u %10u %10u %10u %10u\n",
	 _mm256_extract_epi32(x,0),
	 _mm256_extract_epi32(x,1),
	 _mm256_extract_epi32(x,2),
	 _mm256_extract_epi32(x,3),
	 _mm256_extract_epi32(x,4),
	 _mm256_extract_epi32(x,5),
	 _mm256_extract_epi32(x,6),
	 _mm256_extract_epi32(x,7));
  return;
}

static void
print_vector_hex (__m256i x) {
  printf("%08X %08X %08X %08X %08X %08X %08X %08X\n",
	 _mm256_extract_epi32(x,0),
	 _mm256_extract_epi32(x,1),
	 _mm256_extract_epi32(x,2),
	 _mm256_extract_epi32(x,3),
	 _mm256_extract_epi32(x,4),
	 _mm256_extract_epi32(x,5),
	 _mm256_extract_epi32(x,6),
	 _mm256_extract_epi32(x,7));
  return;
}
#endif


/* Fast scalar scheme designed by N. Kurz. */
/* Delta is required to be added to A to become equivalent with B */
/* Desired result should be A + diagtermA (== B + diagtermB) */
static int *
scalar (int *indices, int **init_indices, int **end_indices, const int A_index, const int B_index,
	const Univcoord_T *A, const Univcoord_T *endA,
	const Univcoord_T *B, const Univcoord_T *endB, int deltaA,
	const Univcoord_T below_slop, const Univcoord_T above_slop,
	bool A_first_p) {

  const Univcoord_T *initA = A;
  const Univcoord_T *initB = B;

  const Univcoord_T *B_save;
  int allocation, nindices;
  int *more_indices;

  debug(printf("Entered scalar with %lu and %lu entries\n",endA - A,endB - B));

  while (A < endA) {
    /* Advance B until (A - below_slop) */
    while (B < endB && (*B) + below_slop < (*A) + deltaA) {
      B++;
    }

    B_save = B;
    debug(printf("Advancing B from %u up to %u\n",*B,(*A) + deltaA + above_slop));
    while (B < endB && (*B) <= (*A) + deltaA + above_slop) {
      debug(printf("(1) Generating output %ld %ld\n",
		   A_index + (A - initA),B_index + (B - initB)));
      if (indices >= (*end_indices) /* same as indices + 2 > (*end_indices) */) {
	allocation = (*end_indices) - (*init_indices);
	more_indices = (int *) MALLOC(2 * allocation * sizeof(int));

	nindices = indices - (*init_indices);
	memcpy(more_indices,*init_indices,nindices * sizeof(int));

	FREE(*init_indices);
	*init_indices = more_indices;
	*end_indices = &(more_indices[2 * allocation]);
	indices = &(more_indices[nindices]);
      }
	
#ifdef INCLUDE_EXACT
      if (A_first_p == true) {
	*indices++ = A_index + (A - initA);
	*indices++ = B_index + (B - initB);
      } else {
	*indices++ = B_index + (B - initB);
	*indices++ = A_index + (A - initA);
      }
#else
      if ((*A) + deltaA != (*B)) {
	if (A_first_p == true) {
	  *indices++ = A_index + (A - initA);
	  *indices++ = B_index + (B - initB);
	} else {
	  *indices++ = B_index + (B - initB);
	  *indices++ = A_index + (A - initA);
	}
      }
#endif
      assert(indices <= *end_indices);

      B++;
    }
    B = B_save;

    A++;
  }

  return indices;
}


/* Performs scalar against a single value of A */
static inline int *
scalar_one_set (int *indices, int **init_indices, int **end_indices, const int A_index, const int B_index,
		const Univcoord_T *A, const Univcoord_T *B, const Univcoord_T *endB,
		int deltaA, const Univcoord_T below_slop, const Univcoord_T above_slop,
		bool A_first_p) {

  const Univcoord_T *initA = A;
  const Univcoord_T *initB = B;

  int allocation, nindices;
  int *more_indices;

  debug(printf("Entered scalar_one_set with %lu entries\n",endB - B));

  debug(printf("Advancing B from %u up to %u\n",*B,(*A) + deltaA - below_slop));
  while (B < endB && (*B) + below_slop < (*A) + deltaA) {
    B++;
  }

  while (B < endB && (*B) <= (*A) + deltaA + above_slop) {
    debug(printf("(2) Generating output %ld %ld\n",
		 A_index + (A - initA),B_index + (B - initB)));

    if (indices >= (*end_indices) /* same as indices + 2 > (*end_indices) */) {
      allocation = (*end_indices) - (*init_indices);
      more_indices = (int *) MALLOC(2 * allocation * sizeof(int));
      
      nindices = indices - (*init_indices);
      memcpy(more_indices,*init_indices,nindices * sizeof(int));

      FREE(*init_indices);
      *init_indices = more_indices;
      *end_indices = &(more_indices[2 * allocation]);
      indices = &(more_indices[nindices]);
    }

#ifdef INCLUDE_EXACT
    if (A_first_p == true) {
      *indices++ = A_index + (A - initA);
      *indices++ = B_index + (B - initB);
    } else {
      *indices++ = B_index + (B - initB);
      *indices++ = A_index + (A - initA);
    }
#else
    if ((*A) + deltaA != (*B)) {
      if (A_first_p == true) {
	*indices++ = A_index + (A - initA);
	*indices++ = B_index + (B - initB);
      } else {
	*indices++ = B_index + (B - initB);
	*indices++ = A_index + (A - initA);
      }
    }
#endif
    assert(indices <= *end_indices);

    B++;
  }

  return indices;
}


#ifdef __GNUC__
#define COMPILER_LIKELY(x) __builtin_expect((x), 1)
#define COMPILER_RARELY(x) __builtin_expect((x), 0)
#else
#define COMPILER_LIKELY(x) x
#define COMPILER_RARELY(x) x
#endif

/**
 * Intersections scheme designed by N. Kurz that works very
 * well when intersecting an array with another where the density
 * differential is small (between 2 to 10).
 *
 * It assumes that lenRare <= lenFreq.
 *
 * Note that this is not symmetric: flipping the rare and freq pointers
 * as well as lenRare and lenFreq could lead to significant performance
 * differences.
 *
 * The matchOut pointer can safely be equal to the rare pointer.
 *
 */

#ifdef HAVE_AVX512
#define NELTS 8		/* 512 / 64 bits per UINT8 */
#define MAX_IDX 65536		/* 2^16 */
#else
#define NELTS 4			/* 256 / 64 per UINT 8*/
#define MAX_IDX 256		/* 2^8 */
#endif

static int match_start[MAX_IDX];
static int match_n[MAX_IDX];

static void
initialize_match_values () {
  int high_bit, low_bit, bit;
  int idx;

  for (idx = 0; idx < MAX_IDX; idx++) {
    match_start[idx] = -1;
    match_n[idx] = 0;
  }

  for (high_bit = 0; high_bit < NELTS; high_bit++) {
    for (low_bit = 0; low_bit <= high_bit; low_bit++) {
      idx = 0;
      for (bit = low_bit; bit <= high_bit; bit++) {
	idx |= (1 << bit);
      }
      match_start[idx] = low_bit;
      match_n[idx] = high_bit - low_bit + 1;
      debug(printf("%d %d %d\n",idx,match_start[idx],match_n[idx]));
      assert(match_n[idx] <= NELTS);
    }
  }
      
  return;
}


#ifdef HAVE_SSE4_2
static int
binary_search_univcoord (int lowi, int highi, const Univcoord_T *positions, Univcoord_T goal) {
  int middlei;

  debug20(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug20(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,positions[lowi],middlei,positions[middlei],
		   highi,positions[highi],goal));
    if (goal < positions[middlei]) {
      highi = middlei;
    } else if (goal > positions[middlei]) {
      lowi = middlei + 1;
    } else {
      debug20(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }

  debug20(printf("binary search returns %d\n",highi));
  return highi;
}
#endif


#ifdef HAVE_SSE4_2
static int *
v1 (int *indices, int **init_indices, int **end_indices, int rare_index, int freq_index,
    const Univcoord_T *rare, int lenRare,
    const Univcoord_T *freq, int lenFreq, int delta_rare,
    const Univcoord_T below_slop, const Univcoord_T above_slop, bool rare_first_p) {

  const Univcoord_T *initRare = rare;
  const Univcoord_T *initFreq = freq;
  const Univcoord_T *endRare = &(rare[lenRare]);
  const Univcoord_T *endFreq = &(freq[lenFreq]);

  const Univcoord_T kRareSpace = 0;	/* 0. type was uint64_t */

#ifndef INCLUDE_EXACT
  Univcoord_T valFreq;
#endif
  Univcoord_T valRare;
  Univcoord_T goal;
  int npositions_rare, npositions_freq, j;

  int i;
  const Univcoord_T below_slop_plus_1 = below_slop + 1; /* Because only cmpgt is available, not cmpge */
  const Univcoord_T above_slop_plus_1 = above_slop + 1; /* Because only cmpgt is available, not cmpge */
  int allocation, nindices;
  int *more_indices;

#ifdef HAVE_AVX512
  const Univcoord_T kProbe = (0 + 1) * 2 * 4; /* 8 64-mers in 512 bits */
  const Univcoord_T kFreqSpace = 2 * 4 * (0 + 1) - 1;

  __m512i Rare_low, Rare_high;
  __m512i F, F_save;
  __m512i _epi64_offset;
  __mmask16 idx, M_above, M_below;
#elif defined(HAVE_AVX2)
  const Univcoord_T kProbe = (0 + 1) * 2 * 2; /* 4 64-mers in 256 bits */
  const Univcoord_T kFreqSpace = 2 * 2 * (0 + 1) - 1;

  __m256i Rare_low, Rare_high;
  __m256i F, F_save;
  __m256i _epi64_offset;
  __m256i M, M_above, M_below;
  int idx;			/* equivalent to __mmask8 */
#else
  const Univcoord_T kProbe = (0 + 1) * 2 * 2; /* 4 64-mers in 128+128 bits */
  const Univcoord_T kFreqSpace = 2 * 2 * (0 + 1) - 1; /* 7 */

  __m128i Rare_low, Rare_high;
  __m128i F0, F1, F0_save, F1_save;
  __m128i _epi64_offset;
  __m128i M, M_above, M_below;
  int idx;			/* equivalent to __mmask8 */
#endif

  const Univcoord_T *stopRare = endRare - kRareSpace;
  const Univcoord_T *stopFreq = endFreq - kFreqSpace;

#ifdef HAVE_AVX512
  _epi64_offset = _mm512_set1_epi64(EPI64_MAX);
#elif defined(HAVE_AVX2)
  _epi64_offset = _mm256_set1_epi64x(EPI64_MAX);
#else
  _epi64_offset = _mm_set1_epi64((__m64) EPI64_MAX);
#endif


  assert(/*lenRare*/(endRare - initRare) <= /*lenFreq*/(endFreq - initFreq));

  /* valRare = (*rare) + delta_rare; -- For comparison, now computed inside condition and inside main loop */
  /* outRare = (*rare) + diagterm_rare; -- For output, not needed */

  /* Avoid issues where valRare - below_slop_plus_1 < 0 */
  while (COMPILER_RARELY(rare < stopRare &&
			 (valRare = (*rare) + delta_rare) < below_slop_plus_1)) {
    indices = scalar_one_set(indices, &(*init_indices), &(*end_indices),
			     rare_index + (rare - initRare), freq_index /*+ (freq - initFreq)*/,
			     rare, freq, endFreq, delta_rare, 
			     below_slop, above_slop, rare_first_p);
    rare++;
    /* valRare = (*rare) + delta_rare; -- For comparison, now computed inside condition and inside main loop */ 
    /* outRare = (*rare) + diagterm_rare; -- For output, not needed */
  }

  if (COMPILER_LIKELY(rare < stopRare)) {
    valRare = (*rare) + delta_rare;

    /* Advance freq */
    npositions_rare = endRare - rare;
    npositions_freq = endFreq - freq;

    if (50 * npositions_rare > npositions_freq) {
      /* Use linear */
      while (freq < stopFreq && /*maxFreq*/freq[kProbe - 1] + below_slop < valRare) {
	debug(printf("Advancing freq by one vector because %u + below_slop < %u\n",
		     freq[kProbe - 1],valRare));
	freq += kProbe;
      }

    } else {
      /* Use galloping search */
      goal = valRare - below_slop;
      j = 1;
      while (j < npositions_freq && freq[j] < goal) {
	j <<= 1;		/* gallop by 2 */
      }
      if (j >= npositions_freq) {
	j = binary_search_univcoord(j >> 1,npositions_freq,freq,goal);
      } else {
	j = binary_search_univcoord(j >> 1,j,freq,goal);
      }
      freq += j;
    }
  }

  if (freq < stopFreq) {
#ifdef HAVE_AVX512
    Rare_low = _mm512_set1_epi64(valRare - below_slop_plus_1 - EPI64_MAX);
    Rare_high = _mm512_set1_epi64(valRare + above_slop_plus_1 - EPI64_MAX);

    F = _mm512_loadu_si512((const __m512i *)(freq));
    F = _mm512_sub_epi64(F, _epi64_offset);
#elif defined(HAVE_AVX2)
    Rare_low = _mm256_set1_epi64x(valRare - below_slop_plus_1 - EPI64_MAX);
    Rare_high = _mm256_set1_epi64x(valRare + above_slop_plus_1 - EPI64_MAX);

    F = _mm256_loadu_si256((const __m256i *)(freq));
    F = _mm256_sub_epi64(F, _epi64_offset);
#else
    Rare_low = _mm_set1_epi64((__m64) (valRare - below_slop_plus_1 - EPI64_MAX));
    Rare_high = _mm_set1_epi64((__m64) (valRare + above_slop_plus_1 - EPI64_MAX));

    F0 = _mm_lddqu_si128((const __m128i *)(freq));
    F0 = _mm_sub_epi64(F0, _epi64_offset);
    F1 = _mm_lddqu_si128((const __m128i *)(freq + 2));
    F1 = _mm_sub_epi64(F1, _epi64_offset);
#endif
  }

  while (rare < stopRare && freq < stopFreq) {
    debug(printf("rare is at item %lu with value %u\n",rare - initRare,valRare));
    debug(printf("freq is at item %lu with values %u..%u\n",freq - initFreq,freq[0],freq[kProbe - 1]));

    const Univcoord_T *freq_save = freq;

    /* First iteration on freq */
#ifdef HAVE_AVX512
    F_save = F;

    M_above = _mm512_cmpgt_epi64_mask(F, Rare_low);
    M_below = _mm512_cmplt_epi64_mask(F, Rare_high);
    idx = M_above & M_below;
#elif defined(HAVE_AVX2)
    F_save = F;

    M_above = _mm256_cmpgt_epi64(F, Rare_low);
    M_below = _mm256_cmpgt_epi64(Rare_high, F); /* cmplt not available */
    M = _mm256_and_si256(M_above,M_below);
    idx = _mm256_movemask_pd((__m256d) M);
#else
    F0_save = F0;
    F1_save = F1;

    M_above = _mm_cmpgt_epi64(F1, Rare_low);
    M_below = _mm_cmpgt_epi64(Rare_high, F1); /* cmplt not available */
    M = _mm_and_si128(M_above,M_below);
    idx = _mm_movemask_pd((__m128d) M);

    M_above = _mm_cmpgt_epi64(F0, Rare_low);
    M_below = _mm_cmpgt_epi64(Rare_high, F0); /* cmplt not available */
    M = _mm_and_si128(M_above,M_below);
    idx = (idx << 2) + _mm_movemask_pd((__m128d) M);
#endif

    debug(printf("Got idx %d => %d x %d\n",idx,match_start[idx],match_n[idx]));

    if (indices + 2*match_n[idx] > (*end_indices)) {
      allocation = (*end_indices) - (*init_indices);
      more_indices = (int *) MALLOC(2 * allocation * sizeof(int));
      
      nindices = indices - (*init_indices);
      memcpy(more_indices,*init_indices,nindices * sizeof(int));
      
      FREE(*init_indices);
      *init_indices = more_indices;
      *end_indices = &(more_indices[2 * allocation]);
      indices = &(more_indices[nindices]);
      assert(indices + 2*match_n[idx] <= (*end_indices));
    }

    const Univcoord_T *freq_ptr = freq + match_start[idx];
    for (i = 0; i < match_n[idx]; i++, freq_ptr++) {
      debug(printf("(3) Generating output %ld %ld\n",
		   rare_index + (rare - initRare), freq_index + (freq_ptr - initFreq)));

#ifdef INCLUDE_EXACT
      if (rare_first_p == true) {
	*indices++ = rare_index + (rare - initRare);
	*indices++ = freq_index + (freq_ptr - initFreq);
      } else {
	*indices++ = freq_index + (freq_ptr - initFreq);
	*indices++ = rare_index + (rare - initRare);
      }
#else
      if (valRare != (valFreq = *freq_ptr++)) {
	if (rare_first_p == true) {
	  *indices++ = rare_index + (rare - initRare);
	  *indices++ = freq_index + (freq_ptr - initFreq);
	} else {
	  *indices++ = freq_index + (freq_ptr - initFreq);
	  *indices++ = rare_index + (rare - initRare);
	}
      }
#endif
    }

    /* Subsequent iterations on freq */
    while (freq + kProbe < stopFreq && /*minFreq*/freq[kProbe] <= valRare + above_slop) {
      freq += kProbe;
      
#ifdef HAVE_AVX512
      F = _mm512_loadu_si512((const __m512i *)(freq));
      F = _mm512_sub_epi64(F, _epi64_offset);
    
      M_above = _mm512_cmpgt_epi64_mask(F, Rare_low);
      M_below = _mm512_cmplt_epi64_mask(F, Rare_high);
      idx = M_above & M_below;
#elif defined(HAVE_AVX2)
      F = _mm256_loadu_si256((const __m256i *)(freq));
      F = _mm256_sub_epi64(F, _epi64_offset);
    
      M_above = _mm256_cmpgt_epi64(F, Rare_low);
      M_below = _mm256_cmpgt_epi64(Rare_high, F); /* cmplt not available */
      M = _mm256_and_si256(M_above,M_below);

      idx = _mm256_movemask_pd((__m256d) M);
#else
      F0 = _mm_lddqu_si128((const __m128i *)(freq));
      F0 = _mm_sub_epi64(F0, _epi64_offset);
      F1 = _mm_lddqu_si128((const __m128i *)(freq + 2));
      F1 = _mm_sub_epi64(F1, _epi64_offset);

      M_above = _mm_cmpgt_epi64(F1, Rare_low);
      M_below = _mm_cmpgt_epi64(Rare_high, F1); /* cmplt not available */
      M = _mm_and_si128(M_above,M_below);
      idx = _mm_movemask_pd((__m128d) M);

      M_above = _mm_cmpgt_epi64(F0, Rare_low);
      M_below = _mm_cmpgt_epi64(Rare_high, F0); /* cmplt not available */
      M = _mm_and_si128(M_above,M_below);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M);
#endif

      debug(printf("Got idx %d => %d x %d\n",idx,match_start[idx],match_n[idx]));

      if (indices + 2*match_n[idx] > (*end_indices)) {
	allocation = (*end_indices) - (*init_indices);
	more_indices = (int *) MALLOC(2 * allocation * sizeof(int));

	nindices = indices - (*init_indices);
	memcpy(more_indices,*init_indices,nindices * sizeof(int));

	FREE(*init_indices);
	*init_indices = more_indices;
	*end_indices = &(more_indices[2 * allocation]);
	indices = &(more_indices[nindices]);
	assert(indices + 2*match_n[idx] <= (*end_indices));
      }

      const Univcoord_T *freq_ptr = freq + match_start[idx];
      for (i = 0; i < match_n[idx]; i++, freq_ptr++) {
	debug(printf("(4) Generating output %ld %ld\n",
		     rare_index + (rare - initRare), freq_index + (freq_ptr - initFreq)));
#ifdef INCLUDE_EXACT
	if (rare_first_p == true) {
	  *indices++ = rare_index + (rare - initRare);
	  *indices++ = freq_index + (freq_ptr - initFreq);
	} else {
	  *indices++ = freq_index + (freq_ptr - initFreq);
	  *indices++ = rare_index + (rare - initRare);
	}
#else
	if (valRare != (valFreq = *freq_ptr++)) {
	  if (rare_first_p == true) {
	    *indices++ = rare_index + (rare - initRare);
	    *indices++ = freq_index + (freq_ptr - initFreq);
	  } else {
	    *indices++ = freq_index + (freq_ptr - initFreq);
	    *indices++ = rare_index + (rare - initRare);
	  }
	}
#endif
      }
    }
    /* End of iterations on freq */

    if (freq + kProbe >= stopFreq) {
      freq += kProbe;
      indices = scalar_one_set(indices, &(*init_indices), &(*end_indices),
			       rare_index + (rare - initRare), freq_index + (freq - initFreq),
			       rare, freq, endFreq, delta_rare, 
			       below_slop, above_slop, rare_first_p);
    }

    if (++rare < stopRare) {
      valRare = (*rare) + delta_rare; /* For comparison */
      /* outRare = (*rare) + diagterm_rare; -- For output, not needed */

#ifdef HAVE_AVX512
      Rare_low = _mm512_set1_epi64(valRare - below_slop_plus_1 - EPI64_MAX);
      Rare_high = _mm512_set1_epi64(valRare + above_slop_plus_1 - EPI64_MAX);
      F = F_save;
#elif defined(HAVE_AVX2)
      Rare_low = _mm256_set1_epi64x(valRare - below_slop_plus_1 - EPI64_MAX);
      Rare_high = _mm256_set1_epi64x(valRare + above_slop_plus_1 - EPI64_MAX);
      F = F_save;
#else
      Rare_low = _mm_set1_epi64((__m64) (valRare - below_slop_plus_1 - EPI64_MAX));
      Rare_high = _mm_set1_epi64((__m64) (valRare + above_slop_plus_1 - EPI64_MAX));
      F0 = F0_save;
      F1 = F1_save;
#endif

      /* Advance freq */
      freq = freq_save;

      npositions_rare = endRare - rare;
      npositions_freq = endFreq - freq;

      if (50 * npositions_rare > npositions_freq) {
	/* Use linear */
	while (freq < stopFreq && /*maxFreq*/freq[kProbe - 1] + below_slop < valRare) {
	  debug(printf("Advancing freq by %d because %u + below_slop < %u\n",
		       kProbe,freq[kProbe - 1],valRare));
	  freq += kProbe;
	}

      } else {
	/* Use galloping search */
	goal = valRare - below_slop;
	j = 1;
	while (j < npositions_freq && freq[j] < goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= npositions_freq) {
	  j = binary_search_univcoord(j >> 1,npositions_freq,freq,goal);
	} else {
	  j = binary_search_univcoord(j >> 1,j,freq,goal);
	}
	freq += j;
      }

      if (freq < stopFreq && freq != freq_save) {
#ifdef HAVE_AVX512
	F = _mm512_loadu_si512((const __m512i *)(freq));
	F = _mm512_sub_epi64(F, _epi64_offset);
#elif defined(HAVE_AVX2)
	F = _mm256_loadu_si256((const __m256i *)(freq));
	F = _mm256_sub_epi64(F, _epi64_offset);
#else
	F0 = _mm_lddqu_si128((const __m128i *)(freq));
	F0 = _mm_sub_epi64(F0, _epi64_offset);
	F1 = _mm_lddqu_si128((const __m128i *)(freq + 2));
	F1 = _mm_sub_epi64(F1, _epi64_offset);
#endif
      }

      debug(printf("Advancing rare to %lu to compare with stop %lu\n",
		   rare - initRare,stopRare - initRare));
      debug(printf("Restoring freq to %lu\n",freq - initFreq));
    }
  }

  return scalar(indices, &(*init_indices), &(*end_indices),
		rare_index + (rare - initRare), freq_index + (freq - initFreq),
		rare, endRare, freq, endFreq, delta_rare, 
		below_slop, above_slop, rare_first_p);
}
#endif


/**
 * Our main heuristic.
 *
 * The out pointer can be set1 if length1<=length2,
 * or else it can be set2 if length1>length2.
 */
/* set1 + diagterm1 == set2 + diagterm2 */
#ifdef HAVE_SSE4_2
int *
Intersect_approx_indices_uint8 (int *ndiagpairs,
				const Univcoord_T *set1, const int length1, int diagterm1,
				const Univcoord_T *set2, const int length2, int diagterm2,
				const Univcoord_T below_slop, const Univcoord_T above_slop) {
  int *indices, *init_indices, *end_indices;
  int allocation;

  debug(printf("Entered Intersect_approx_indices with %d and %d items\n",length1,length2));

#ifdef DEBUG
  int i;
  for (i = 0; i < length1; i++) {
    printf("%u\n",set1[i]);
  }
  printf("\n");
  for (i = 0; i < length2; i++) {
    printf("%u\n",set2[i]);
  }
  printf("\n");
#endif    

  if ((length1 == 0) || (length2 == 0)) {
    *ndiagpairs = 0;
    debug(printf("Exiting Intersect_approx_indices with 0 items\n"));
    return (int *) NULL;

  } else {
    debug(printf("Using v1 method\n"));
    if (length1 <= length2) {
      if ((allocation = 2 * length1) < 2*NELTS) {
	/* Needs to be greater than 2*match_n */
	allocation = 2*NELTS;
      }
      indices = init_indices = (int *) MALLOC(allocation * sizeof(int));
      end_indices = &(indices[allocation]);
      indices = v1(indices, &init_indices, &end_indices, /*rare_index*/0, /*freq_index*/0,
		     set1, length1, set2, length2,
		     /*delta_rare*/(diagterm1 - diagterm2),
		     below_slop, above_slop, /*rare_first_p*/true);
    } else {
      if ((allocation = 2 * length2) < 2*NELTS) {
	allocation = 2*NELTS;
      }
      indices = init_indices = (int *) MALLOC(allocation * sizeof(int));
      end_indices = &(indices[allocation]);
      indices = v1(indices, &init_indices, &end_indices, /*rare_index*/0, /*freq_index*/0,
		     set2, length2, set1, length1,
		     /*delta_rare*/(diagterm2 - diagterm1),
		     above_slop, below_slop, /*rare_first_p*/false);
    }
    assert(indices <= end_indices);

    *ndiagpairs = (indices - init_indices) / 2;

#ifdef DEBUG
    int k;
    printf("Exiting Intersect_approx_indices with %d items\n",*ndiagpairs);
    for (i = 0, k = 0; i < *ndiagpairs; i++, k += 2) {
      printf("%u %u\n",init_indices[k],init_indices[k+1]);
    }
#endif

    return init_indices;
  }
}

#else

int *
Intersect_approx_indices_uint8 (int *ndiagpairs,
				const Univcoord_T *set1, const int length1, int diagterm1,
				const Univcoord_T *set2, const int length2, int diagterm2,
				const Univcoord_T below_slop, const Univcoord_T above_slop) {
  int *indices, *init_indices, *end_indices;
  int allocation;

  if ((length1 == 0) || (length2 == 0)) {
    *ndiagpairs = 0;
    return (int *) NULL;

  } else {
    if (length1 <= length2) {
      allocation = 2 * length1;
      indices = init_indices = (int *) MALLOC(allocation * sizeof(int));
      end_indices = &(indices[allocation]);
      indices = scalar(indices, &init_indices, &end_indices, /*A_index*/0, /*B_index*/0,
			 set1, /*endA*/&(set1[length1]),
			 set2, /*endB*/&(set2[length2]),
			 /*deltaA*/(diagterm1 - diagterm2),
			 below_slop, above_slop, /*A_first_p*/true);
    } else {
      allocation = 2 * length2;
      indices = init_indices = (int *) MALLOC(allocation * sizeof(int));
      end_indices = &(indices[allocation]);
      indices = scalar(indices, &init_indices, &end_indices, /*A_index*/0, /*B_index*/0,
			 set2, /*endA*/&(set2[length2]),
			 set1, /*endB*/&(set1[length1]),
			 /*deltaA*/(diagterm2 - diagterm1),
			 above_slop, below_slop, /*A_first_p*/false);
    }

    *ndiagpairs = (indices - init_indices) / 2;
    return init_indices;
  }
}
#endif


void
Intersect_approx_indices_uint8_setup () {
  initialize_match_values();
  return;
}  