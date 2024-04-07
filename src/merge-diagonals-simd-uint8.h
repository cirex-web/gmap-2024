/* $Id: merge-diagonals-simd-uint8.h 224755 2021-12-13 00:42:15Z twu $ */
#ifndef MERGE_DIAGONALS_SIMD_UINT8_INCLUDED
#define MERGE_DIAGONALS_SIMD_UINT8_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"		/* For HAVE_64_BIT */
#endif

#include "types.h"
#include "univcoord.h"
#include "mergeinfo.h"

/* kmer-search.c calls Merge_diagonals_large with stream_high_list and
   stream_low_list.  localdb.c calls Merge_diagonals_uint4 and
   Merge_diagonals_uint8 */

#if defined(HAVE_AVX512) || defined(HAVE_AVX2)
extern Univcoord_T *
Merge_diagonals_large (int *nelts1, unsigned char **stream_high_array, UINT4 **stream_low_array,
		       int *streamsize_array, int *diagterm_array, int nstreams,
		       Mergeinfo_uint8_T mergeinfo);
extern Univcoord_T *
Merge_diagonals_uint8 (int *nelts1, Univcoord_T **stream_array, int *streamsize_array,
		       int nstreams, Mergeinfo_uint8_T mergeinfo);
#endif

#endif


