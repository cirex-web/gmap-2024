static char rcsid[] = "$Id: dc1507a14904fd014196e6e5fb60c7a1444f70b8 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "localdb-read.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>		/* For munmap */

#include "assert.h"
#include "mem.h"
#include "genomicpos.h"
#include "genomebits_consec.h"
#include "genomebits_count.h"
#include "genomebits_trim.h"
/* #include "merge-uint4.h" -- for get_by_bisection */
#include "intersect-uint2.h"
#include "mergeinfo.h"

#ifndef LARGE_GENOMES
#include "merge-diagonals-simd-uint4.h"
#elif defined(HAVE_AVX512) || defined(HAVE_AVX2)
#include "merge-diagonals-simd-uint8.h"
#else
#include "merge-diagonals-heap.h"
#endif


/* Print univdiagonals */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Method used */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Details of suffix array search */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* get_exhaustive */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Localdb_get_one */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif

/* Dump each suffix array contents */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

/* Dump saindex */
#ifdef DEBUG9A
#define debug9a(x) x
#else
#define debug9a(x)
#endif

/* Compare specialized Localdb_get against a general method */
#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif


/* It is critical to have the "U" at the end of each of these values.
   Otherwise, a compiler sometimes optimizes to give strange results */
#define SAINDEX16_LENGTH 256U	/* 2^8 */
#define SARRAY16_LENGTH 65536U	/* 2^16 */
#define SARRAY8_LENGTH 256U	/* 2^8 */


#define OLIGOSIZE 4 /* nucleotides for saindex to specify 256 pointers */


static Genome_T genomecomp;
static Univcoord_T genomelength;
static int index1part;
static int index1interval;

void
Localdb_setup (Genome_T genomecomp_in, Univcoord_T genomelength_in,
	       int index1part_in, int index1interval_in) {
  genomecomp = genomecomp_in;
  genomelength = genomelength_in;
  index1part = index1part_in;
  index1interval = index1interval_in;

  return;
}


#define T Localdb_T
struct T {
  uint16_t *saindex16;
  uint16_t *sarray16;
  uint8_t *sarray8;
  uint16_t *sasort16;

  Access_T localdb_access;

  int saindex16_shmid;
  key_t saindex16_key;
  int saindex16_fd;
  size_t saindex16_len;

  int sarray16_shmid;
  key_t sarray16_key;
  int sarray16_fd;
  size_t sarray16_len;

  int sarray8_shmid;
  key_t sarray8_key;
  int sarray8_fd;
  size_t sarray8_len;

  int sasort16_shmid;
  key_t sasort16_key;
  int sasort16_fd;
  size_t sasort16_len;
};


void
Localdb_free (T *old) {

  if ((*old)->localdb_access == ALLOCATED_PRIVATE) {
    FREE_KEEP((*old)->saindex16);
    FREE_KEEP((*old)->sarray16);
    FREE_KEEP((*old)->sarray8);
    FREE_KEEP((*old)->sasort16);

  } else if ((*old)->localdb_access == ALLOCATED_SHARED) {
    Access_deallocate((*old)->saindex16,(*old)->saindex16_shmid,(*old)->saindex16_key);
    Access_deallocate((*old)->sarray16,(*old)->sarray16_shmid,(*old)->sarray16_key);
    Access_deallocate((*old)->sarray8,(*old)->sarray8_shmid,(*old)->sarray8_key);
    Access_deallocate((*old)->sasort16,(*old)->sasort16_shmid,(*old)->sasort16_key);

#ifdef HAVE_MMAP
  } else if ((*old)->localdb_access == MMAPPED) {
    munmap((void *) (*old)->saindex16,(*old)->saindex16_len);
    munmap((void *) (*old)->sarray16,(*old)->sarray16_len);
    munmap((void *) (*old)->sarray8,(*old)->sarray8_len);
    munmap((void *) (*old)->sasort16,(*old)->sasort16_len);
    close((*old)->saindex16_fd);
    close((*old)->sarray16_fd);
    close((*old)->sarray8_fd);
    close((*old)->sasort16_fd);
#endif
  }
  
  FREE(*old);
}


/* Position 10470: CGGTACCC */
/* TODO: Write get char that works on genomebits */

/* Taken from Johannes Fischer, Advanced Text Indexing Techniques, Algorithm 1 */
/* Does not use LCP, so time is O(m * log(n)) */
static void
sarray_search_uint1 (int *matchlength, int *initptr, int *finalptr, uint8_t *SA, int n,
		     char *query, int querylength, int pos5, int pos3, int min_matchlength,
		     Univcoord_T region_term, Compress_T query_compress, bool plusp,
		     Genomebits_T genomebits) {
  int low, high, mid;
  Univcoord_T pos;
  int nmatches = 0;
  int seglength = pos3 - pos5;
  char mismatch_char;

#ifdef DEBUG1
  char *Buffer;
  Buffer = (char *) MALLOC((seglength+2)*sizeof(char));
#endif

  debug1(printf("Entered sarray_search_uint1\n"));

  low = 0;
  high = n;

  while (low < high) {
    /* Compute mid for unsigned ints.  Want floor((low+high)/2). */
    /* low 0, high 2, mid = 1.  low 0, high 3, mid = 1.  low 1, high 2, mid = 1.  low 1, high 3, mid = 2 */
    mid = low + (high - low + 0)/2;

    pos = region_term + SA[mid];

#ifdef DEBUG1
    printf("low %d, high %d, mid %d => position %u\n",low,high,mid,pos);
    Genome_fill_buffer_simple(genomecomp,/*left*/pos,seglength,Buffer);
    if (SA[mid] + seglength >= n) {
      Buffer[n - SA[mid]] = '$';
      Buffer[n - SA[mid] + 1] = '\0';
    }
    printf("%d\t%d\t%s\n",mid,SA[mid],Buffer);
    printf("query\t\t%.*s\n",seglength,&(query[pos5]));
#endif

    nmatches = Genomebits_consecutive_matches_wmm(&mismatch_char,genomebits,query_compress,
						  /*univdiagonal*/pos-pos5+querylength,querylength,
						  pos5,pos3,plusp,/*genestrand*/0);
    if (SA[mid] + nmatches >= n) {
      /* Terminating char is always less than query */
      debug1(printf("nmatches %d, genomic char $\n",nmatches));
      low = mid + 1;
    } else if (nmatches >= seglength) {
      debug1(printf("nmatches %d >= seglength\n",nmatches));
      high = mid;
    } else if (mismatch_char > query[pos5+nmatches]) {
      /* assert(mismatch_char == Genome_bits_get_char(pos+nmatches)); */
      debug1(printf("nmatches %d, genomic char %c > query %c\n",nmatches,
		   mismatch_char,query[pos5+nmatches]));
      high = mid;
    } else {
      /* assert(mismatch_char == Genome_bits_get_char(pos+nmatches)); */
      debug1(printf("nmatches %d, genomic char %c <= query %c\n",nmatches,
		   mismatch_char,query[pos5+nmatches]));
      low = mid + 1;
    }
  }

  if ((*initptr = low) == n) {
    /* Query is higher than any suffix */
    *matchlength = 0;
    *finalptr = *initptr;
  } else {
    pos = region_term + SA[*initptr];
    nmatches = Genomebits_consecutive_matches_wmm(&mismatch_char,genomebits,query_compress,
						  /*univdiagonal*/pos-pos5+querylength,querylength,
						  pos5,pos3,plusp,/*genestrand*/0);
    debug1(printf("Setting initptr to be %d with %d matches\n",low,nmatches));

    if (nmatches < seglength) {
      *matchlength = nmatches;
    } else {
      *matchlength = seglength;
    }

    if (*matchlength < min_matchlength) {
      *finalptr = *initptr;
    } else {
      *finalptr = (*initptr) + 1;
      while (*finalptr < n &&
	     Genomebits_consecutive_matches_wmm(&mismatch_char,genomebits,query_compress,
						/*univdiagonal*/region_term + SA[*finalptr] - pos5 + querylength,
						querylength,pos5,pos3,plusp,/*genestrand*/0) >= *matchlength) {
	*finalptr += 1;
      }
    }
  }


#ifdef DEBUG1
  printf("sarray_search_uint1 returning %d..%d with matchlength %d\n",*initptr,*finalptr,*matchlength);
  for (mid = *initptr; mid < *finalptr; mid++) {
    printf("%u\n",region_term + SA[mid]);
  }
  printf("\n");
#endif

#ifdef DEBUG1
  FREE(Buffer);
#endif

  return;
}


#ifdef DEBUG2
static void
oligo_nt (char nt[], int oligo, int oligosize) {
  int i, j;

  j = oligosize-1;
  for (i = 0; i < oligosize; i++) {
    switch (oligo & 0x03) {
    case 0: nt[j] = 'A'; break;
    case 1: nt[j] = 'C'; break;
    case 2: nt[j] = 'G'; break;
    case 3: nt[j] = 'T'; break;
    }
    oligo >>= 2;
    j--;
  }

  nt[oligosize] = '\0';
  return;
}
#endif


static int
nt_oligo (char *query, int oligosize) {
  int oligo = 0;
  int i;

  for (i = 0; i < oligosize; i++) {
    oligo *= 4;
    
    switch (query[i]) {
    case 'A': break;
    case 'C': oligo += 1; break;
    case 'G': oligo += 2; break;
    case 'T': oligo += 3; break;
    default: oligo += 3; break;	/* Treat X as T */
    }
  }

  return oligo;
}

static unsigned char
nt_oligo_byte (char *query, int oligosize) {
  unsigned char oligo = 0;
  int i;

  for (i = 0; i < oligosize; i++) {
    oligo <<= 2;
    
    switch (query[i]) {
    case 'A': break;
    case 'C': oligo += 1; break;
    case 'G': oligo += 2; break;
    case 'T': oligo += 3; break;
    default: oligo += 3; break;	/* Treat X as T */
    }
  }

  return oligo;
}

static inline unsigned char
revise_oligo_byte (unsigned char oligo, char c) {
  switch (c) {
  case 'A': return (oligo << 2);
  case 'C': return ((oligo << 2) + 1);
  case 'G': return ((oligo << 2) + 2);
  case 'T': return ((oligo << 2) + 3);
  default: return ((oligo << 2) + 3);	/* Treat X as T */
  }
}


static void
sarray_search_uint2 (int *matchlength, int *initptr, int *finalptr, uint16_t *saindex, uint16_t *SA, int n,
		     char *query, int querylength, int pos5, int pos3, int min_matchlength, Univcoord_T region_term,
		     Compress_T query_compress, bool plusp, Genomebits_T genomebits) {
  int low, high, mid;
  int oligo;
  Univcoord_T pos;
  int nmatches = 0;
  int seglength = pos3 - pos5;
  char mismatch_char;

#if defined(DEBUG1) || defined(DEBUG9)
  char *suffix;
  suffix = (char *) MALLOC((seglength+2)*sizeof(char)); /* Add 2 so we can print $ at position n */
#endif

  debug1(printf("Entered sarray_search_uint2\n"));

  
#ifdef DEBUG9
  int i;
  unsigned int k = 0;

#ifdef DEBUG9A
  /* Dump saindex pointers */
  for (k = 0; k < SAINDEX16_LENGTH; k++) {
    printf("%d\t%d\n",k,saindex[k]);
  }
  printf("\n");
#endif

  for (i = 0; i < n; i++) {
    Genome_fill_buffer_simple(genomecomp,/*left*/region_term+SA[i],seglength,suffix);
    if (SA[i] + seglength >= n) {
      suffix[n - SA[i]] = '$';
      suffix[n - SA[i] + 1] = '\0';
    }
    printf("%d\t%d\t%s\n",i,SA[i],suffix);
  }
#endif

#if 0
  for (int i = 0; i < n; i++) {
    if (SA[i] >= n) {
      printf("At i %d, SA[i] is %u, but n is %u\n",i,SA[i],n);
      abort();
    }
  }
#endif

  if (seglength < OLIGOSIZE) {
    low = 0;
    high = n;
  } else if ((oligo = nt_oligo(&(query[pos5]),OLIGOSIZE)) == SAINDEX16_LENGTH - 1) {
    low = saindex[oligo];
    high = n;
  } else {
    low = saindex[oligo];
    high = saindex[oligo+1];
  }
  debug1(printf("query: %.*s => oligo %d, initial low %d, initial high %d\n",
	       OLIGOSIZE,&(query[pos5]),oligo,low,high));

  if (low == high) {
    /* matchlength must be less than OLIGOSIZE */
    *matchlength = 0;
    *initptr = *finalptr = low;
    debug1(printf("sarray_search_uint2 returning %d..%d\n\n",*initptr,*finalptr));
#if defined(DEBUG1) || defined(DEBUG9)
    FREE(suffix);
#endif
    return;
  }

  while (low < high) {
    /* Compute mid for unsigned ints.  Want floor((low+high)/2). */
    /* low 0, high 2, mid = 1.  low 0, high 3, mid = 1.  low 1, high 2, mid = 1.  low 1, high 3, mid = 2 */
    mid = low + (high - low + 0)/2;

    pos = region_term + SA[mid];

#ifdef DEBUG1
    printf("low %d, high %d, mid %d => position %u\n",low,high,mid,pos);
    Genome_fill_buffer_simple(genomecomp,/*left*/pos,seglength,suffix);
    if (SA[mid] + seglength >= n) {
      suffix[n - SA[mid]] = '$';
      suffix[n - SA[mid] + 1] = '\0';
    }
    printf("%d\t%d\t%s\n",mid,SA[mid],suffix);
    printf("query\t\t%.*s\n",seglength,&(query[pos5]));
#endif

    nmatches = Genomebits_consecutive_matches_wmm(&mismatch_char,genomebits,query_compress,
						  /*univdiagonal*/pos-pos5+querylength,querylength,
						  pos5,pos3,plusp,/*genestrand*/0);


    if (SA[mid] + nmatches >= n) {
      /* Terminating char is always less than query */
      debug1(printf("nmatches %d, genomic char $\n",nmatches));
      low = mid + 1;
    } else if (nmatches >= seglength) {
      debug1(printf("nmatches %d >= seglength\n",nmatches));
      high = mid;
    } else if (mismatch_char > query[pos5+nmatches]) {
      /* assert(mismatch_char == Genome_bits_get_char(pos+nmatches)); */
      debug1(printf("nmatches %d, genomic char %c > query %c\n",nmatches,
		   mismatch_char,query[pos5+nmatches]));
      high = mid;
    } else {
      /* assert(mismatch_char == Genome_bits_get_char(pos+nmatches)); */
      debug1(printf("nmatches %d, genomic char %c <= query %c\n",nmatches,
		   mismatch_char,query[pos5+nmatches]));
      low = mid + 1;
    }
  }

  if ((*initptr = low) == n) {
    /* Query is higher than any suffix */
    *matchlength = 0;
    *finalptr = *initptr;
  } else {
    pos = region_term + SA[*initptr];
    nmatches = Genomebits_consecutive_matches_wmm(&mismatch_char,genomebits,query_compress,
						  /*univdiagonal*/pos-pos5+querylength,querylength,
						  pos5,pos3,plusp,/*genestrand*/0);
    debug1(printf("Setting initptr to be %d with %d matches\n",low,nmatches));


    if (nmatches < seglength) {
      *matchlength = nmatches;
    } else {
      *matchlength = seglength;
    }

    if (*matchlength < min_matchlength) {
      *finalptr = *initptr;
    } else {
      *finalptr = (*initptr) + 1;
      while (*finalptr < n &&
	     Genomebits_consecutive_matches_wmm(&mismatch_char,genomebits,query_compress,
						/*univdiagonal*/(region_term + SA[*finalptr])-pos5+querylength,
						querylength,pos5,pos3,
						plusp,/*genestrand*/0) >= *matchlength) {
	*finalptr += 1;
      }
    }
  }


#ifdef DEBUG1
  printf("sarray_search_uint2 returning %d..%d with matchlength %d\n",*initptr,*finalptr,*matchlength);
  for (mid = *initptr; mid < *finalptr; mid++) {
    printf("%u\n",region_term + SA[mid]);
  }
  printf("\n");
#endif
  
#if defined(DEBUG1) || defined(DEBUG9)
  FREE(suffix);
#endif

  return;
}


#if 0
    /* Binary search for finalptr.  ? More accurate than marching because of the terminating character */
    high = n;

    while (low < high) {
      /* Compute mid for unsigned ints.  Want ceil((low+high)/2). */
      /* low 0, high 2, mid = 1.  low 0, high 3, mid = 2.  low 1, high 2, mid = 2.  low 1, high 3, mid = 2 */
      mid = low + (high - low + 1)/2;

      pos = region_term + SA[mid];

#ifdef DEBUG1
      printf("low %d, high %d, mid %d => position %u\n",low,high,mid,pos);
      Genome_fill_buffer_simple(genomecomp,/*left*/pos,seglength,Buffer);
      if (SA[mid] + seglength >= n) {
	Buffer[n - SA[mid]] = '$';
	Buffer[n - SA[mid] + 1] = '\0';
      }
      printf("%d\t%d\t%s\n",mid,SA[mid],Buffer);
      printf("query\t\t%.*s\n",seglength,&(query[pos5]));
#endif

      nmatches = Genomebits_consecutive_matches_wmm(&mismatch_char,genomebits,query_compress,
						    /*univdiagonal*/pos-pos5+querylength,querylength,
						    pos5,pos3,plusp,/*genestrand*/0);
      if (SA[mid] + nmatches >= n) {
	/* Terminating char is always less than query */
	debug1(printf("nmatches %d, genomic char $\n",nmatches));
	low = mid;
      } else if (nmatches >= seglength) {
	debug1(printf("nmatches %d >= seglength\n",nmatches));
	low = mid;
      } else if (Genome_bits_get_char(pos+nmatches) < query[pos5+nmatches]) {
	debug1(printf("nmatches %d, genomic char %c < query %c\n",nmatches,
		     Genome_bits_get_char(pos+nmatches),query[pos5+nmatches]));
	low = mid;
      } else {
	debug1(printf("nmatches %d, genomic char %c >= query %c\n",nmatches,
		     Genome_bits_get_char(pos+nmatches),query[pos5+nmatches]));
	high = mid - 1;
      }
    }
  }
  *finalptr = high;


  low--;
  high = n;

  while (low < high) {
    /* Compute mid for unsigned ints.  Want ceil((low+high)/2). */
    /* low 0, high 2, mid = 1.  low 0, high 3, mid = 2.  low 1, high 2, mid = 2.  low 1, high 3, mid = 2 */
    mid = low + (high - low + 1)/2;

    pos = region_term + SA[mid];
    nmatches = Genomebits_consecutive_matches_wmm(&mismatch_char,genomebits,query_compress,
						  /*univdiagonal*/pos-pos5+querylength,querylength,
						  pos5,pos3,plusp,/*genestrand*/0);
    if (nmatches == seglength || (c = Genome_bits_get_char(pos+nmatches)) < query[pos5+nmatches]) {
      low = mid;
    } else {
      high = mid - 1;
    }
  }
  *finalptr = high;
#endif


static int
univcoord_compare (const void *a, const void *b) {
  Univcoord_T x = * (Univcoord_T *) a;
  Univcoord_T y = * (Univcoord_T *) b;

  if (x < y) {
    return -1;
  } else if (y < x) {
    return 1;
  } else {
    return 0;
  }
}


/* Hard to find maximum matchlength without scanning first */
static int
get_general (int *matchlength, Univcoord_T *diagonals, T this,
	     char *queryptr, int pos5, int pos3, int querylength, int min_matchlength,
	     Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
	     Compress_T query_compress, bool plusp, Genomebits_T genomebits) {
  int region_matchlength;
  int nentries;
  Univcoord_T *out = diagonals;
  Univcoord_T low_position_start, high_position_start, diagonal, last_diagonal, *p, *q;
  int initptr, finalptr, ptr, n;
  uint16_t *saindex16;
  uint16_t *sarray16;
  uint8_t *sarray8;
  Univcoord_T low_regioni, high_regioni, regioni;
  Univcoord_T large_region_term, small_region_term;

  *matchlength = 0;

  if (low_univdiagonal < (Univcoord_T) (querylength - pos5)) {
    low_position_start = 0U;
  } else {
    low_position_start = low_univdiagonal - querylength + pos5;
  }
  if (high_univdiagonal < (Univcoord_T) (querylength - pos5)) {
    high_position_start = 0U;
  } else {
    high_position_start = high_univdiagonal - querylength + pos5;
  }

  low_regioni = low_position_start/SARRAY16_LENGTH;
  high_regioni = high_position_start/SARRAY16_LENGTH;

  debug1(printf("Entering get_general\n"));

  /* low region */
  regioni = low_regioni;
  large_region_term = regioni * SARRAY16_LENGTH;
  if (low_position_start < (small_region_term = large_region_term + SARRAY16_LENGTH - SARRAY8_LENGTH/2)) {
    /* Start with L */
    saindex16 = &(this->saindex16[regioni * SAINDEX16_LENGTH]);
    sarray16 = &(this->sarray16[large_region_term]);
    n = SARRAY16_LENGTH;
    debug1(printf("L* (region %u): %u..%u\n",regioni,large_region_term,large_region_term + n));
    sarray_search_uint2(&region_matchlength,&initptr,&finalptr,saindex16,sarray16,
			n,queryptr,querylength,pos5,pos3,
			min_matchlength,large_region_term,
			query_compress,plusp,genomebits);
    if (region_matchlength < *matchlength) {
      /* Skip */
    } else {
      if (region_matchlength > *matchlength) {
	*matchlength = region_matchlength;
	out = diagonals;		/* Reset output */
      }
      for (ptr = initptr; ptr < finalptr; ptr++) {
	if ((diagonal = large_region_term + sarray16[ptr] - pos5 + querylength) >= low_univdiagonal &&
	    diagonal <= high_univdiagonal) {
	  *out++ = diagonal;
	}
      }
    }
  }

  /* S */
  if (high_regioni > low_regioni) {
    sarray8 = &(this->sarray8[regioni * SARRAY8_LENGTH]);
    n = SARRAY8_LENGTH;
    debug1(printf("S* (region %u): %u..%u\n",regioni,small_region_term,small_region_term + n));
    sarray_search_uint1(&region_matchlength,&initptr,&finalptr,sarray8,
			n,queryptr,querylength,pos5,pos3,
			min_matchlength,small_region_term,
			query_compress,plusp,genomebits);
    if (region_matchlength < *matchlength) {
      /* Skip */
    } else {
      if (region_matchlength > *matchlength) {
	*matchlength = region_matchlength;
	out = diagonals;	/* Reset output */
      }
      for (ptr = initptr; ptr < finalptr; ptr++) {
	if ((diagonal = small_region_term + sarray8[ptr] - pos5 + querylength) >= low_univdiagonal &&
	    diagonal <= high_univdiagonal) {
	  *out++ = diagonal;
	}
      }
    }
  }
  
  /* Middle regions */
  for (regioni = low_regioni + 1; regioni < high_regioni; regioni++) {
    large_region_term += SARRAY16_LENGTH;
    small_region_term += SARRAY16_LENGTH;

    /* L */
    saindex16 = &(this->saindex16[regioni * SAINDEX16_LENGTH]);
    sarray16 = &(this->sarray16[large_region_term]);
    n = SARRAY16_LENGTH;
    debug1(printf("L* (region %u): %u..%u\n",regioni,large_region_term,large_region_term + n));
    sarray_search_uint2(&region_matchlength,&initptr,&finalptr,saindex16,sarray16,
			n,queryptr,querylength,pos5,pos3,
			min_matchlength,large_region_term,
			query_compress,plusp,genomebits);
    if (region_matchlength < *matchlength) {
      /* Skip */
    } else {
      if (region_matchlength > *matchlength) {
	*matchlength = region_matchlength;
	out = diagonals;	/* Reset output */
      }
      for (ptr = initptr; ptr < finalptr; ptr++) {
	diagonal = large_region_term + sarray16[ptr] - pos5 + querylength;
	assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	*out++ = diagonal;
      }
    }
    
    /* S */
    sarray8 = &(this->sarray8[regioni * SARRAY8_LENGTH]);
    n = (small_region_term + SARRAY8_LENGTH < genomelength) ? SARRAY8_LENGTH : genomelength - small_region_term;
    debug1(printf("S* (region %u): %u..%u\n",regioni,small_region_term,small_region_term + n));
    sarray_search_uint1(&region_matchlength,&initptr,&finalptr,sarray8,
			n,queryptr,querylength,pos5,pos3,
			min_matchlength,small_region_term,
			query_compress,plusp,genomebits);
    if (region_matchlength < *matchlength) {
      /* Skip */
    } else {
      if (region_matchlength > *matchlength) {
	*matchlength = region_matchlength;
	out = diagonals;	/* Reset output */
      }
      for (ptr = initptr; ptr < finalptr; ptr++) {
	if ((diagonal = small_region_term + sarray8[ptr] - pos5 + querylength) <= high_univdiagonal) {
	  *out++ = diagonal;
	}
      }
    }
  }

  /* High region */
  large_region_term += SARRAY16_LENGTH;
  if (high_univdiagonal < large_region_term + SARRAY8_LENGTH/2) {
    /* No need for L */
  } else {      
    /* L */
    saindex16 = &(this->saindex16[regioni * SAINDEX16_LENGTH]);
    sarray16 = &(this->sarray16[large_region_term]);
    n = (large_region_term + SARRAY16_LENGTH < genomelength) ? SARRAY16_LENGTH: genomelength - large_region_term;
    debug1(printf("L* (region %u): %u..%u\n",regioni,large_region_term,large_region_term + n));
    sarray_search_uint2(&region_matchlength,&initptr,&finalptr,saindex16,sarray16,
			n,queryptr,querylength,pos5,pos3,
			min_matchlength,large_region_term,
			query_compress,plusp,genomebits);
    if (region_matchlength < *matchlength) {
      /* Skip */
    } else {
      if (region_matchlength > *matchlength) {
	*matchlength = region_matchlength;
	out = diagonals;	/* Reset output */
      }
      for (ptr = initptr; ptr < finalptr; ptr++) {
	if ((diagonal = large_region_term + sarray16[ptr] - pos5 + querylength) <= high_univdiagonal) {
	  *out++ = diagonal;
	}
      }
    }
  }

  if ((nentries = out - diagonals) == 0) {
    return 0;

  } else {
    qsort(diagonals,nentries,sizeof(Univcoord_T),univcoord_compare);

    /* Remove duplicates */
    q = diagonals;
    last_diagonal = *q;
    for (p = &(diagonals[1]); p < &(diagonals[nentries]); p++) {
      if (*p == last_diagonal) {
	/* Skip: repetitive */
      } else {
	*++q = *p;
      }
      last_diagonal = *p;
    }
    nentries = (q - diagonals) + 1;
    return nentries;
  }
}


#ifdef DEBUG0
static void
print_univdiagonals (Univcoord_T *diagonals, int nentries) {
  int i;

  printf("%d entries starting at address %p\n",nentries,diagonals);
  for (i = 0; i < nentries; i++) {
    printf("%u\n",diagonals[i]);
  }
  printf("\n");
    
  return;
}
#endif


#ifdef DEBUG15
static void
check_univdiagonals_general (Univcoord_T *diagonals, int nentries, T this,
			     char *queryptr, int pos5, int pos3, int querylength,
			     Univcoord_T low_position, Univcoord_T high_position,
			     Compress_T query_compress, Genomebits_T genomebits) {
  Univcoord_T *diagonals_alt;
  int nentries_alt, i;
  bool errorp = false;


  diagonals_alt = get_general(&nentries_alt,this,queryptr,
					    pos5,pos3,querylength,
					    low_position,high_position,
					    query_compress,plusp,genomebits);

  if (nentries != nentries_alt) {
    errorp = true;
  } else {
    for (i = 0; i < nentries; i++) {
      if (diagonals[i] != diagonals_alt[i]) {
	errorp = true;
      }
    }
  }

  if (errorp == true) {
    printf("Standard method: %d entries\n",nentries);
    for (i = 0; i < nentries; i++) {
      printf("%u\n",diagonals[i]);
    }
    printf("\n");
    
    printf("General method: %d entries\n",nentries_alt);
    for (i = 0; i < nentries_alt; i++) {
      printf("%u\n",diagonals_alt[i]);
    }
    printf("\n");
    abort();
  }

  FREE(diagonals_alt);
  return;
}
#endif


/* ? Do results need to be sorted */
static int
get_by_sarray (int *matchlength, Univcoord_T *diagonals, T this,
	       char *queryptr, int pos5, int pos3, int querylength, int min_matchlength,
	       Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
	       Compress_T query_compress, bool plusp, Genomebits_T genomebits) {
  int nentries;
  Univcoord_T *out = diagonals;
  Univcoord_T low_position_start, high_position_start, diagonal, last_diagonal, *p, *q;
  int low_initptr, low_finalptr, mid_initptr, mid_finalptr, high_initptr, high_finalptr,
    small0_initptr, small0_finalptr, small1_initptr, small1_finalptr, ptr, n;
  uint16_t *low_saindex16, *mid_saindex16, *high_saindex16;
  uint16_t *low_sarray16, *mid_sarray16, *high_sarray16;
  uint8_t *small0_sarray8, *small1_sarray8;
  Univcoord_T low_regioni, mid_regioni, high_regioni;
  Univcoord_T low_region_term, mid_region_term, high_region_term,
    small0_region_term, small1_region_term;
  int low_matchlength, mid_matchlength, high_matchlength, small0_matchlength, small1_matchlength;

  debug1(printf("Entered get_by_sarray with pos5 %d, pos3 %d, %u..%u\n",
	       pos5,pos3,low_univdiagonal,high_univdiagonal));

  assert(pos5 < pos3);

  if (low_univdiagonal < (Univcoord_T) (querylength - pos5)) {
    low_position_start = 0U;
  } else {
    low_position_start = low_univdiagonal - querylength + pos5;
  }
  if (high_univdiagonal < (Univcoord_T) (querylength - pos5)) {
    high_position_start = 0U;
  } else {
    high_position_start = high_univdiagonal - querylength + pos5;
  }

  low_regioni = low_position_start/SARRAY16_LENGTH;
  high_regioni = high_position_start/SARRAY16_LENGTH;

  debug1(printf("low %u => region %u.  high %u => region %u\n",
	       low_position_start,low_regioni,high_position_start,high_regioni));

  if (high_regioni == low_regioni) {
    debug1(printf("Case L\n"));
    low_region_term = low_regioni * SARRAY16_LENGTH;
    low_saindex16 = &(this->saindex16[low_regioni * SAINDEX16_LENGTH]);
    low_sarray16 = &(this->sarray16[low_region_term]);
    n = (low_region_term + SARRAY16_LENGTH < genomelength) ? SARRAY16_LENGTH : genomelength - low_region_term;
    debug1(printf("L (region %u): %u..%u\n",low_regioni,low_region_term,low_region_term + n));
    sarray_search_uint2(&(*matchlength),&low_initptr,&low_finalptr,low_saindex16,low_sarray16,
			n,queryptr,querylength,pos5,pos3,
			min_matchlength,low_region_term,
			query_compress,plusp,genomebits);

    for (ptr = low_initptr; ptr < low_finalptr; ptr++) {
      if ((diagonal = low_region_term + low_sarray16[ptr] - pos5 + querylength) >= low_univdiagonal &&
	  diagonal <= high_univdiagonal) {
	*out++ = diagonal;
      }
    }

    if ((nentries = out - diagonals) == 0) {
      debug0(print_univdiagonals(diagonals,nentries));
      return 0;

    } else {
      /* One region, so no need to remove duplicates */
#ifdef ALWAYS_SORT
      qsort(diagonals,nentries,sizeof(Univcoord_T),univcoord_compare);
#endif
      debug0(print_univdiagonals(diagonals,nentries));
      return nentries;
    }

  } else if (high_regioni == low_regioni + 1) {
    low_region_term = low_regioni * SARRAY16_LENGTH;
    high_region_term = low_region_term + SARRAY16_LENGTH;

    if (low_position_start >= high_region_term - SARRAY8_LENGTH/2 &&
	high_position_start < high_region_term + SARRAY8_LENGTH/2) {
      debug1(printf("Case S\n"));
      small0_region_term = high_region_term - SARRAY8_LENGTH/2;
      small0_sarray8 = &(this->sarray8[low_regioni * SARRAY8_LENGTH]);
      n = (small0_region_term + SARRAY8_LENGTH < genomelength) ? SARRAY8_LENGTH : genomelength - small0_region_term;
      debug1(printf("S (region %u): %u..%u\n",low_regioni,small0_region_term,small0_region_term + n));
      sarray_search_uint1(&(*matchlength),&small0_initptr,&small0_finalptr,small0_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small0_region_term,
			  query_compress,plusp,genomebits);

      for (ptr = small0_initptr; ptr < small0_finalptr; ptr++) {
	if ((diagonal = small0_region_term + small0_sarray8[ptr] - pos5 + querylength) >= low_univdiagonal &&
	    diagonal <= high_univdiagonal) {
	  *out++ = diagonal;
	}
      }

      if ((nentries = out - diagonals) == 0) {
	debug0(print_univdiagonals(/*diagonals*/NULL,nentries));
	return 0;
      } else {
	/* One region, so no need to remove duplicates */
#ifdef ALWAYS_SORT
	qsort(diagonals,nentries,sizeof(Univcoord_T),univcoord_compare);
#endif
	debug0(print_univdiagonals(diagonals,nentries));
	return nentries;
      }

    } else if (low_position_start >= high_region_term - SARRAY8_LENGTH/2) {
      debug1(printf("Case S,L\n"));
      small0_region_term = high_region_term - SARRAY8_LENGTH/2;
      small0_sarray8 = &(this->sarray8[low_regioni * SARRAY8_LENGTH]);
      n = SARRAY8_LENGTH;
      debug1(printf("S (region %u): %u..%u\n",low_regioni,small0_region_term,small0_region_term + n));
      sarray_search_uint1(&small0_matchlength,&small0_initptr,&small0_finalptr,small0_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small0_region_term,
			  query_compress,plusp,genomebits);

      high_saindex16 = &(this->saindex16[high_regioni * SAINDEX16_LENGTH]);
      high_sarray16 = &(this->sarray16[high_region_term]);
      n = (high_region_term + SARRAY16_LENGTH < genomelength) ? SARRAY16_LENGTH : genomelength - high_region_term;
      debug1(printf("L (region %u): %u..%u\n",high_regioni,high_region_term,high_region_term + n));
      sarray_search_uint2(&high_matchlength,&high_initptr,&high_finalptr,high_saindex16,high_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,high_region_term,
			  query_compress,plusp,genomebits);

      *matchlength = small0_matchlength;
      if (high_matchlength > *matchlength) *matchlength = high_matchlength;

      if (small0_matchlength == *matchlength) {
	for (ptr = small0_initptr; ptr < small0_finalptr; ptr++) {
	  if ((diagonal = small0_region_term + small0_sarray8[ptr] - pos5 + querylength) >= low_univdiagonal) {
	    *out++ = diagonal;
	  }
	}
      }
      if (high_matchlength == *matchlength) {
	for (ptr = high_initptr; ptr < high_finalptr; ptr++) {
	  if ((diagonal = high_region_term + high_sarray16[ptr] - pos5 + querylength) <= high_univdiagonal) {
	    *out++ = diagonal;
	  }
	}
      }

      if ((nentries = out - diagonals) == 0) {
	debug0(print_univdiagonals(/*diagonals*/NULL,nentries));
	return 0;

      } else {
	qsort(diagonals,nentries,sizeof(Univcoord_T),univcoord_compare);
	/* Remove duplicates */
	q = diagonals;
	last_diagonal = *q;
	for (p = &(diagonals[1]); p < &(diagonals[nentries]); p++) {
	  if (*p == last_diagonal) {
	    /* Skip: repetitive */
	  } else {
	    *++q = *p;
	  }
	  last_diagonal = *p;
	}
	nentries = (q - diagonals) + 1;
	debug0(print_univdiagonals(diagonals,nentries));
	return nentries;
      }

    } else if (high_position_start < high_region_term + SARRAY8_LENGTH/2) {
      debug1(printf("Case L,S\n"));
      low_saindex16 = &(this->saindex16[low_regioni * SAINDEX16_LENGTH]);
      low_sarray16 = &(this->sarray16[low_region_term]);
      n = SARRAY16_LENGTH;
      debug1(printf("L (region %u): %u..%u\n",low_regioni,low_region_term,low_region_term + n));
      sarray_search_uint2(&low_matchlength,&low_initptr,&low_finalptr,low_saindex16,low_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,low_region_term,
			  query_compress,plusp,genomebits);

      small0_region_term = high_region_term - SARRAY8_LENGTH/2;
      small0_sarray8 = &(this->sarray8[low_regioni * SARRAY8_LENGTH]);
      n = (small0_region_term + SARRAY8_LENGTH < genomelength) ? SARRAY8_LENGTH : genomelength - small0_region_term;
      debug1(printf("S (region %u): %u..%u\n",low_regioni,small0_region_term,small0_region_term + n));
      sarray_search_uint1(&small0_matchlength,&small0_initptr,&small0_finalptr,small0_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small0_region_term,
			  query_compress,plusp,genomebits);

      *matchlength = low_matchlength;
      if (small0_matchlength > *matchlength) *matchlength = small0_matchlength;

      if (low_matchlength == *matchlength) {
	for (ptr = low_initptr; ptr < low_finalptr; ptr++) {
	  if ((diagonal = low_region_term + low_sarray16[ptr] - pos5 + querylength) >= low_univdiagonal) {
	    *out++ = diagonal;
	  }
	}
      }
      if (small0_matchlength == *matchlength) {
	for (ptr = small0_initptr; ptr < small0_finalptr; ptr++) {
	  if ((diagonal = small0_region_term + small0_sarray8[ptr] - pos5 + querylength) <= high_univdiagonal) {
	    *out++ = diagonal;
	  }
	}
      }

      if ((nentries = out - diagonals) == 0) {
	debug0(print_univdiagonals(/*diagonals*/NULL,nentries));
	return 0;

      } else {
	qsort(diagonals,nentries,sizeof(Univcoord_T),univcoord_compare);
	/* Remove duplicates */
	q = diagonals;
	last_diagonal = *q;
	for (p = &(diagonals[1]); p < &(diagonals[nentries]); p++) {
	  if (*p == last_diagonal) {
	    /* Skip: repetitive */
	  } else {
	    *++q = *p;
	  }
	  last_diagonal = *p;
	}
	nentries = (q - diagonals) + 1;
	debug0(print_univdiagonals(diagonals,nentries));
	return nentries;
      }

    } else {
      debug1(printf("Case L,S,L\n"));
      low_saindex16 = &(this->saindex16[low_regioni * SAINDEX16_LENGTH]);
      low_sarray16 = &(this->sarray16[low_region_term]);
      n = SARRAY16_LENGTH;
      debug1(printf("L (region %u): %u..%u\n",low_regioni,low_region_term,low_region_term + n));
      sarray_search_uint2(&low_matchlength,&low_initptr,&low_finalptr,low_saindex16,low_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,low_region_term,
			  query_compress,plusp,genomebits);

      small0_region_term = high_region_term - SARRAY8_LENGTH/2;
      small0_sarray8 = &(this->sarray8[low_regioni * SARRAY8_LENGTH]);
      n = SARRAY8_LENGTH;
      debug1(printf("S (region %u): %u..%u\n",low_regioni,small0_region_term,small0_region_term + n));
      sarray_search_uint1(&small0_matchlength,&small0_initptr,&small0_finalptr,small0_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small0_region_term,
			  query_compress,plusp,genomebits);

      high_saindex16 = &(this->saindex16[high_regioni * SAINDEX16_LENGTH]);
      high_sarray16 = &(this->sarray16[high_region_term]);
      n = (high_region_term + SARRAY16_LENGTH < genomelength) ? SARRAY16_LENGTH: genomelength - high_region_term;
      debug1(printf("L (region %u): %u..%u\n",high_regioni,high_region_term,high_region_term + n));
      sarray_search_uint2(&high_matchlength,&high_initptr,&high_finalptr,high_saindex16,high_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,high_region_term,
			  query_compress,plusp,genomebits);

      *matchlength = low_matchlength;
      if (small0_matchlength > *matchlength) *matchlength = small0_matchlength;
      if (high_matchlength > *matchlength) *matchlength = high_matchlength;

      if (low_matchlength == *matchlength) {
	for (ptr = low_initptr; ptr < low_finalptr; ptr++) {
	  if ((diagonal = low_region_term + low_sarray16[ptr] - pos5 + querylength) >= low_univdiagonal) {
	    *out++ = diagonal;
	  }
	}
      }
      if (small0_matchlength == *matchlength) {
	for (ptr = small0_initptr; ptr < small0_finalptr; ptr++) {
	  diagonal = small0_region_term + small0_sarray8[ptr] - pos5 + querylength;
	  assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	  *out++ = diagonal;
	}
      }
      if (high_matchlength == *matchlength) {
	for (ptr = high_initptr; ptr < high_finalptr; ptr++) {
	  if ((diagonal = high_region_term + high_sarray16[ptr] - pos5 + querylength) <= high_univdiagonal) {
	    *out++ = diagonal;
	  }
	}
      }

      if ((nentries = out - diagonals) == 0) {
	debug0(print_univdiagonals(/*diagonals*/NULL,nentries));
	return 0;

      } else {
	qsort(diagonals,nentries,sizeof(Univcoord_T),univcoord_compare);
	/* Remove duplicates */
	q = diagonals;
	last_diagonal = *q;
	for (p = &(diagonals[1]); p < &(diagonals[nentries]); p++) {
	  if (*p == last_diagonal) {
	    /* Skip: repetitive */
	  } else {
	    *++q = *p;
	  }
	  last_diagonal = *p;
	}
	nentries = (q - diagonals) + 1;
	debug0(print_univdiagonals(diagonals,nentries));
	return nentries;
      }
    }

  } else if (high_regioni == low_regioni + 2) {
    mid_regioni = low_regioni + 1;
    low_region_term = low_regioni * SARRAY16_LENGTH;
    mid_region_term = low_region_term + SARRAY16_LENGTH;
    high_region_term = mid_region_term + SARRAY16_LENGTH;

    if (low_position_start >= mid_region_term - SARRAY8_LENGTH/2 &&
	high_position_start < high_region_term + SARRAY8_LENGTH/2) {
      debug1(printf("Case S,L,S\n"));
      small0_region_term = mid_region_term - SARRAY8_LENGTH/2;
      small0_sarray8 = &(this->sarray8[low_regioni * SARRAY8_LENGTH]);
      n = SARRAY8_LENGTH;
      debug1(printf("S (region %u): %u..%u\n",low_regioni,small0_region_term,small0_region_term + n));
      sarray_search_uint1(&small0_matchlength,&small0_initptr,&small0_finalptr,small0_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small0_region_term,
			  query_compress,plusp,genomebits);

      mid_saindex16 = &(this->saindex16[mid_regioni * SAINDEX16_LENGTH]);
      mid_sarray16 = &(this->sarray16[mid_region_term]);
      n = SARRAY16_LENGTH;
      debug1(printf("L (region %u): %u..%u\n",mid_regioni,mid_region_term,mid_region_term + n));
      sarray_search_uint2(&mid_matchlength,&mid_initptr,&mid_finalptr,mid_saindex16,mid_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,mid_region_term,
			  query_compress,plusp,genomebits);

      small1_region_term = high_region_term - SARRAY8_LENGTH/2;
      small1_sarray8 = &(this->sarray8[mid_regioni * SARRAY8_LENGTH]);
      n = (small1_region_term + SARRAY8_LENGTH < genomelength) ? SARRAY8_LENGTH : genomelength - small1_region_term;
      debug1(printf("S (region %u): %u..%u\n",mid_regioni,small1_region_term,small1_region_term + n));
      sarray_search_uint1(&small1_matchlength,&small1_initptr,&small1_finalptr,small1_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small0_region_term,
			  query_compress,plusp,genomebits);

      *matchlength = small0_matchlength;
      if (mid_matchlength > *matchlength) *matchlength = mid_matchlength;
      if (small1_matchlength > *matchlength) *matchlength = small1_matchlength;

      if (small0_matchlength == *matchlength) {
	for (ptr = small0_initptr; ptr < small0_finalptr; ptr++) {
	  if ((diagonal = small0_region_term + small0_sarray8[ptr] - pos5 + querylength) >= low_univdiagonal) {
	    *out++ = diagonal;
	  }
	}
      }
      if (mid_matchlength == *matchlength) {
	for (ptr = mid_initptr; ptr < mid_finalptr; ptr++) {
	  diagonal = mid_region_term + mid_sarray16[ptr] - pos5 + querylength;
	  assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	  *out++ = diagonal;
	}
      }
      if (small1_matchlength == *matchlength) {
	for (ptr = small1_initptr; ptr < small1_finalptr; ptr++) {
	  if ((diagonal = small1_region_term + small1_sarray8[ptr] - pos5 + querylength) <= high_univdiagonal) {
	    *out++ = diagonal;
	  }
	}
      }

      if ((nentries = out - diagonals) == 0) {
	debug0(print_univdiagonals(/*diagonals*/NULL,nentries));
	return 0;

      } else {
	qsort(diagonals,nentries,sizeof(Univcoord_T),univcoord_compare);
	/* Remove duplicates */
	q = diagonals;
	last_diagonal = *q;
	for (p = &(diagonals[1]); p < &(diagonals[nentries]); p++) {
	  if (*p == last_diagonal) {
	    /* Skip: repetitive */
	  } else {
	    *++q = *p;
	  }
	  last_diagonal = *p;
	}
	nentries = (q - diagonals) + 1;
	debug0(print_univdiagonals(diagonals,nentries));
	return nentries;
      }

    } else if (low_position_start >= mid_region_term - SARRAY8_LENGTH/2) {
      debug1(printf("Case S,L,S,L\n"));
      small0_region_term = mid_region_term - SARRAY8_LENGTH/2;
      small0_sarray8 = &(this->sarray8[low_regioni * SARRAY8_LENGTH]);
      n = SARRAY8_LENGTH;
      debug1(printf("S (region %u): %u..%u\n",low_regioni,small0_region_term,small0_region_term + n));
      sarray_search_uint1(&small0_matchlength,&small0_initptr,&small0_finalptr,small0_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small0_region_term,
			  query_compress,plusp,genomebits);

      mid_saindex16 = &(this->saindex16[mid_regioni * SAINDEX16_LENGTH]);
      mid_sarray16 = &(this->sarray16[mid_region_term]);
      n = SARRAY16_LENGTH;
      debug1(printf("L (region %u): %u..%u\n",mid_regioni,mid_region_term,mid_region_term + n));
      sarray_search_uint2(&mid_matchlength,&mid_initptr,&mid_finalptr,mid_saindex16,mid_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,mid_region_term,
			  query_compress,plusp,genomebits);

      small1_region_term = high_region_term - SARRAY8_LENGTH/2;
      small1_sarray8 = &(this->sarray8[mid_regioni * SARRAY8_LENGTH]);
      n = SARRAY8_LENGTH;
      debug1(printf("S (region %u): %u..%u\n",mid_regioni,small1_region_term,small1_region_term + n));
      sarray_search_uint1(&small1_matchlength,&small1_initptr,&small1_finalptr,small1_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small1_region_term,
			  query_compress,plusp,genomebits);


      high_saindex16 = &(this->saindex16[high_regioni * SAINDEX16_LENGTH]);
      high_sarray16 = &(this->sarray16[high_region_term]);
      n = (high_region_term + SARRAY16_LENGTH < genomelength) ? SARRAY16_LENGTH: genomelength - high_region_term;
      debug1(printf("L (region %u): %u..%u\n",high_regioni,high_region_term,high_region_term + n));
      sarray_search_uint2(&high_matchlength,&high_initptr,&high_finalptr,high_saindex16,high_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,high_region_term,
			  query_compress,plusp,genomebits);

      *matchlength = small0_matchlength;
      if (mid_matchlength > *matchlength) *matchlength = mid_matchlength;
      if (small1_matchlength > *matchlength) *matchlength = small1_matchlength;
      if (high_matchlength > *matchlength) *matchlength = high_matchlength;

      if (small0_matchlength == *matchlength) {
	for (ptr = small0_initptr; ptr < small0_finalptr; ptr++) {
	  if ((diagonal = small0_region_term + small0_sarray8[ptr] - pos5 + querylength) >= low_univdiagonal) {
	    *out++ = diagonal;
	  }
	}
      }
      if (mid_matchlength == *matchlength) {
	for (ptr = mid_initptr; ptr < mid_finalptr; ptr++) {
	  diagonal = mid_region_term + mid_sarray16[ptr] - pos5 + querylength;
	  assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	  *out++ = diagonal;
	}
      }
      if (small1_matchlength == *matchlength) {
	for (ptr = small1_initptr; ptr < small1_finalptr; ptr++) {
	  diagonal = small1_region_term + small1_sarray8[ptr] - pos5 + querylength;
	  assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	  *out++ = diagonal;
	}
      }
      if (high_matchlength == *matchlength) {
	for (ptr = high_initptr; ptr < high_finalptr; ptr++) {
	  if ((diagonal = high_region_term + high_sarray16[ptr] - pos5 + querylength) <= high_univdiagonal) {
	    *out++ = diagonal;
	  }
	}
      }

      if ((nentries = out - diagonals) == 0) {
	debug0(print_univdiagonals(/*diagonals*/NULL,nentries));
	return 0;

      } else {
	qsort(diagonals,nentries,sizeof(Univcoord_T),univcoord_compare);
	/* Remove duplicates */
	q = diagonals;
	last_diagonal = *q;
	for (p = &(diagonals[1]); p < &(diagonals[nentries]); p++) {
	  if (*p == last_diagonal) {
	    /* Skip: repetitive */
	  } else {
	    *++q = *p;
	  }
	  last_diagonal = *p;
	}
	nentries = (q - diagonals) + 1;
	debug0(print_univdiagonals(diagonals,nentries));
	return nentries;
      }

    } else if (high_position_start < high_region_term + SARRAY8_LENGTH/2) {
      debug1(printf("Case L,S,L,S\n"));
      low_saindex16 = &(this->saindex16[low_regioni * SAINDEX16_LENGTH]);
      low_sarray16 = &(this->sarray16[low_region_term]);
      n = SARRAY16_LENGTH;
      debug1(printf("L (region %u): %u..%u\n",low_regioni,low_region_term,low_region_term + n));
      sarray_search_uint2(&low_matchlength,&low_initptr,&low_finalptr,low_saindex16,low_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,low_region_term,
			  query_compress,plusp,genomebits);

      small0_region_term = mid_region_term - SARRAY8_LENGTH/2;
      small0_sarray8 = &(this->sarray8[low_regioni * SARRAY8_LENGTH]);
      n = SARRAY8_LENGTH;
      debug1(printf("S (region %u): %u..%u\n",low_regioni,small0_region_term,small0_region_term + n));
      sarray_search_uint1(&small0_matchlength,&small0_initptr,&small0_finalptr,small0_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small0_region_term,
			  query_compress,plusp,genomebits);

      mid_saindex16 = &(this->saindex16[mid_regioni * SAINDEX16_LENGTH]);
      mid_sarray16 = &(this->sarray16[mid_region_term]);
      n = SARRAY16_LENGTH;
      debug1(printf("L (region %u): %u..%u\n",mid_regioni,mid_region_term,mid_region_term + n));
      sarray_search_uint2(&mid_matchlength,&mid_initptr,&mid_finalptr,mid_saindex16,mid_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,mid_region_term,
			  query_compress,plusp,genomebits);
      
      small1_region_term = high_region_term - SARRAY8_LENGTH/2;
      small1_sarray8 = &(this->sarray8[mid_regioni * SARRAY8_LENGTH]);
      n = (small1_region_term + SARRAY8_LENGTH < genomelength) ? SARRAY8_LENGTH : genomelength - small1_region_term;
      debug1(printf("S (region %u): %u..%u\n",mid_regioni,small1_region_term,small1_region_term + n));
      sarray_search_uint1(&small1_matchlength,&small1_initptr,&small1_finalptr,small1_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small1_region_term,
			  query_compress,plusp,genomebits);

      *matchlength = low_matchlength;
      if (small0_matchlength > *matchlength) *matchlength = small0_matchlength;
      if (mid_matchlength > *matchlength) *matchlength = mid_matchlength;
      if (small1_matchlength > *matchlength) *matchlength = small1_matchlength;
      
      if (low_matchlength == *matchlength) {
	for (ptr = low_initptr; ptr < low_finalptr; ptr++) {
	  if ((diagonal = low_region_term + low_sarray16[ptr] - pos5 + querylength) >= low_univdiagonal) {
	    *out++ = diagonal;
	  }
	}
      }
      if (small0_matchlength == *matchlength) {
	for (ptr = small0_initptr; ptr < small0_finalptr; ptr++) {
	  diagonal = small0_region_term + small0_sarray8[ptr] - pos5 + querylength;
	  assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	  *out++ = diagonal;
	}
      }
      if (mid_matchlength == *matchlength) {
	for (ptr = mid_initptr; ptr < mid_finalptr; ptr++) {
	  diagonal = mid_region_term + mid_sarray16[ptr] - pos5 + querylength;
	  assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	  *out++ = diagonal;
	}
      }
      if (small1_matchlength == *matchlength) {
	for (ptr = small1_initptr; ptr < small1_finalptr; ptr++) {
	  if ((diagonal = small1_region_term + small1_sarray8[ptr] - pos5 + querylength) <= high_univdiagonal) {
	    *out++ = diagonal;
	  }
	}
      }

      if ((nentries = out - diagonals) == 0) {
	debug0(print_univdiagonals(/*diagonals*/NULL,nentries));
	return nentries;

      } else {
	qsort(diagonals,nentries,sizeof(Univcoord_T),univcoord_compare);
	/* Remove duplicates */
	q = diagonals;
	last_diagonal = *q;
	for (p = &(diagonals[1]); p < &(diagonals[nentries]); p++) {
	  if (*p == last_diagonal) {
	    /* Skip: repetitive */
	  } else {
	    *++q = *p;
	  }
	  last_diagonal = *p;
	}
	nentries = (q - diagonals) + 1;
	debug0(print_univdiagonals(diagonals,nentries));
	return nentries;
      }

    } else {
      debug1(printf("Case L,S,L,S,L\n"));
      low_saindex16 = &(this->saindex16[low_regioni * SAINDEX16_LENGTH]);
      low_sarray16 = &(this->sarray16[low_region_term]);
      n = SARRAY16_LENGTH;
      debug1(printf("L (region %u): %u..%u\n",low_regioni,low_region_term,low_region_term + n));
      sarray_search_uint2(&low_matchlength,&low_initptr,&low_finalptr,low_saindex16,low_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,low_region_term,
			  query_compress,plusp,genomebits);

      small0_region_term = mid_region_term - SARRAY8_LENGTH/2;
      small0_sarray8 = &(this->sarray8[low_regioni * SARRAY8_LENGTH]);
      n = SARRAY8_LENGTH;
      debug1(printf("S (region %u): %u..%u\n",low_regioni,small0_region_term,small0_region_term + n));
      sarray_search_uint1(&small0_matchlength,&small0_initptr,&small0_finalptr,small0_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small0_region_term,
			  query_compress,plusp,genomebits);

      mid_saindex16 = &(this->saindex16[mid_regioni * SAINDEX16_LENGTH]);
      mid_sarray16 = &(this->sarray16[mid_region_term]);
      n = SARRAY16_LENGTH;
      debug1(printf("L (region %u): %u..%u\n",mid_regioni,mid_region_term,mid_region_term + n));
      sarray_search_uint2(&mid_matchlength,&mid_initptr,&mid_finalptr,mid_saindex16,mid_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,mid_region_term,
			  query_compress,plusp,genomebits);
      
      small1_region_term = high_region_term - SARRAY8_LENGTH/2;
      small1_sarray8 = &(this->sarray8[mid_regioni * SARRAY8_LENGTH]);
      n = SARRAY8_LENGTH;
      debug1(printf("S (region %u): %u..%u\n",mid_regioni,small1_region_term,small1_region_term + n));
      sarray_search_uint1(&small1_matchlength,&small1_initptr,&small1_finalptr,small1_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small1_region_term,
			  query_compress,plusp,genomebits);

      high_saindex16 = &(this->saindex16[high_regioni * SAINDEX16_LENGTH]);
      high_sarray16 = &(this->sarray16[high_region_term]);
      n = (high_region_term + SARRAY16_LENGTH < genomelength) ? SARRAY16_LENGTH: genomelength - high_region_term;
      debug1(printf("L (region %u): %u..%u\n",high_regioni,high_region_term,high_region_term + n));
      sarray_search_uint2(&high_matchlength,&high_initptr,&high_finalptr,high_saindex16,high_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,high_region_term,
			  query_compress,plusp,genomebits);

      *matchlength = low_matchlength;
      if (small0_matchlength > *matchlength) *matchlength = small0_matchlength;
      if (mid_matchlength > *matchlength) *matchlength = mid_matchlength;
      if (small1_matchlength > *matchlength) *matchlength = small1_matchlength;
      if (high_matchlength > *matchlength) *matchlength = high_matchlength;

      if (low_matchlength == *matchlength) {
	for (ptr = low_initptr; ptr < low_finalptr; ptr++) {
	  if ((diagonal = low_region_term + low_sarray16[ptr] - pos5 + querylength) >= low_univdiagonal) {
	    *out++ = diagonal;
	  }
	}
      }
      if (small0_matchlength == *matchlength) {
	for (ptr = small0_initptr; ptr < small0_finalptr; ptr++) {
	  diagonal = small0_region_term + small0_sarray8[ptr] - pos5 + querylength;
	  assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	  *out++ = diagonal;
	}
      }
      if (mid_matchlength == *matchlength) {
	for (ptr = mid_initptr; ptr < mid_finalptr; ptr++) {
	  diagonal = mid_region_term + mid_sarray16[ptr] - pos5 + querylength;
	  assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	  *out++ = diagonal;
	}
      }
      if (small1_matchlength == *matchlength) {
	for (ptr = small1_initptr; ptr < small1_finalptr; ptr++) {
	  diagonal = small1_region_term + small1_sarray8[ptr] - pos5 + querylength;
	  assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	  *out++ = diagonal;
	}
      }
      if (high_matchlength == *matchlength) {
	for (ptr = high_initptr; ptr < high_finalptr; ptr++) {
	  if ((diagonal = high_region_term + high_sarray16[ptr] - pos5 + querylength) <= high_univdiagonal) {
	    *out++ = diagonal;
	  }
	}
      }

      if ((nentries = out - diagonals) == 0) {
	debug0(print_univdiagonals(/*diagonals*/NULL,nentries));
	return 0;

      } else {
	qsort(diagonals,nentries,sizeof(Univcoord_T),univcoord_compare);
	/* Remove duplicates */
	q = diagonals;
	last_diagonal = *q;
	for (p = &(diagonals[1]); p < &(diagonals[nentries]); p++) {
	  if (*p == last_diagonal) {
	    /* Skip: repetitive */
	  } else {
	    *++q = *p;
	  }
	  last_diagonal = *p;
	}
	nentries = (q - diagonals) + 1;
	debug0(print_univdiagonals(diagonals,nentries));
	return nentries;
      }
    }

  } else {
#if 0
    /* If the region is larger than 65536, then it is likely coming from end < start */
    fprintf(stderr,"get_by_sarray cannot handle a region this large (%llu) from %llu to %llu => regioni %d..%d\n",
	    (unsigned long long) (high_position_start - low_position_start),
	    (unsigned long long) low_position_start,(unsigned long long) high_position_start,
	    low_regioni,high_regioni);
    abort();
    return 0;
#else
    return get_general(&(*matchlength),&(*diagonals),this,
		       queryptr,pos5,pos3,querylength,min_matchlength,
		       low_univdiagonal,high_univdiagonal,query_compress,
		       plusp,genomebits);
#endif
  }
}


#if 0
/* Overwrites existing array */
static int
remove_duplicate_diagonals (Univcoord_T *diagonals, int ndiagonals) {
  Univcoord_T *out = diagonals;
  int i, j;
  
  i = 0;
  while (i < ndiagonals) {
    *out++ = diagonals[i];
    j = i + 1;
    while (j < ndiagonals && diagonals[j] == diagonals[i]) {
      j++;
    }
    i = j;
  }
  
  return (out - diagonals);
}
#endif


static int
get_by_ends (int *nmismatches, Univcoord_T *diagonals, T this, unsigned short *localdb_alloc,
	     char *queryptr, int pos5, int pos3, int querylength,
	     Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
	     Compress_T query_compress, bool plusp, int genestrand,
	     Genomebits_T genomebits, int nmismatches_allowed) {
  int nentries;
  Univcoord_T *out = diagonals;
  int lowi5, highi5, lowi3, highi3, initptr, finalptr, nptrs5, nptrs3, n;
  uint16_t *saindex16;
  uint16_t *sasort16;
  Univcoord_T low_position_start, high_position_start, low_regioni, high_regioni, regioni;
  Univcoord_T region_term;

  int i;
  unsigned char oligo5, oligo3;
  int ref_nmismatches, best_mismatches, mm;
#ifndef HAVE_SSE2
  unsigned short *ignore_ptr;
#endif
  
#ifdef DEBUG2
  char nt[5];
#endif


  debug2(printf("Entering get_by_ends with %.*s\n",pos3-pos5,&(queryptr[pos5])));  
  debug2(printf("low_univdiagonal %u, high_univdiagonal %u\n",low_univdiagonal,high_univdiagonal));
  assert(pos3 - pos5 >= OLIGOSIZE);

  if (low_univdiagonal < (Univcoord_T) (querylength - pos5)) {
    low_position_start = 0U;
  } else {
    low_position_start = low_univdiagonal - querylength + pos5;
  }
  if (high_univdiagonal < (Univcoord_T) (querylength - pos5)) {
    high_position_start = 0U;
  } else {
    high_position_start = high_univdiagonal - querylength + pos5;
  }

  low_regioni = low_position_start/SARRAY16_LENGTH;
  high_regioni = high_position_start/SARRAY16_LENGTH;

  oligo5 = nt_oligo_byte(&(queryptr[pos5]),OLIGOSIZE);
  oligo3 = nt_oligo_byte(&(queryptr[pos3 - OLIGOSIZE]),OLIGOSIZE);

  /* Gather 4-mer positions for the two end oligos */
  region_term = low_regioni * SARRAY16_LENGTH;
  for (regioni = low_regioni; regioni <= high_regioni; regioni++) {
    saindex16 = &(this->saindex16[regioni * SAINDEX16_LENGTH]);
    sasort16 = &(this->sasort16[region_term]);
    n = (region_term + SARRAY16_LENGTH < genomelength) ? SARRAY16_LENGTH: genomelength - region_term;

    initptr = saindex16[oligo5];
    finalptr = (oligo5 == SAINDEX16_LENGTH - 1) ? n : saindex16[oligo5+1];
    debug2(oligo_nt(nt,oligo5,4));
    debug2(printf("oligo5 %d (%s) => initptr %d, finalptr %d\n",oligo5,nt,initptr,finalptr));

    /* Filter diagonals */
    /* Bound by low_position and high_position */
    lowi5 = initptr;
    while (lowi5 < finalptr && (region_term + sasort16[lowi5]) - pos5 + querylength < low_univdiagonal) {
      lowi5++;
    }
    highi5 = finalptr - 1;
    while (highi5 >= lowi5 && (region_term + sasort16[highi5]) - pos5 + querylength > high_univdiagonal) {
      highi5--;
    }
    nptrs5 = highi5 - lowi5 + 1;


    initptr = saindex16[oligo3];
    finalptr = (oligo3 == SAINDEX16_LENGTH - 1) ? n : saindex16[oligo3+1];
    debug2(oligo_nt(nt,oligo3,4));
    debug2(printf("oligo3 %d (%s) => initptr %d, finalptr %d\n",oligo3,nt,initptr,finalptr));

    /* Filter diagonals */
    /* Bound by low_position and high_position */
    lowi3 = initptr;
    while (lowi3 < finalptr && (region_term + sasort16[lowi3]) - pos5 + querylength < low_univdiagonal) {
      lowi3++;
    }
    highi3 = finalptr - 1;
    while (highi3 >= lowi3 && (region_term + sasort16[highi3]) - pos5 + querylength > high_univdiagonal) {
      highi3--;
    }
    nptrs3 = highi3 - lowi3 + 1;

    if (nptrs5 > 0 && nptrs3 > 0) {
#ifdef HAVE_SSE2
      out += Intersect_uint2(out,localdb_alloc,
			     &(sasort16[lowi5]),nptrs5,/*diagterm1*/querylength - pos5,
			     &(sasort16[lowi3]),nptrs3,/*diagterm2*/querylength - (pos3 - OLIGOSIZE),
			     region_term);
#else
      out += Intersect_uint2_scalar(out,&ignore_ptr,
			     &(sasort16[lowi5]),nptrs5,/*diagterm1*/querylength - pos5,
			     &(sasort16[lowi3]),nptrs3,/*diagterm2*/querylength - (pos3 - OLIGOSIZE),
			     region_term);
#endif
    }

    region_term += SARRAY16_LENGTH;
  }

  if ((nentries = out - diagonals) == 0) {
    return 0;

  } else {
    best_mismatches = querylength;
    for (i = 0; i < out - diagonals; i++) {
      if ((mm = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,/*alt*/NULL,query_compress,
						      /*univdiagonal*/diagonals[i],querylength,
						      pos5,pos3,plusp,genestrand)) < best_mismatches) {
	nentries = 0;
	diagonals[nentries++] = diagonals[i];
	best_mismatches = mm;

      } else if (mm == best_mismatches) {
	diagonals[nentries++] = diagonals[i];
      }
    }

    if (best_mismatches > nmismatches_allowed) {
      return 0;

    } else {
#ifdef DEBUG2 
      printf("\n");
      printf("%d entries:",nentries);
      for (i = 0; i < nentries; i++) {
	printf(" %u",diagonals[i]);
      }
      printf("\n");
#endif

      *nmismatches = best_mismatches;
      return nentries;
    }
  }
}



/* Copied from kmer-search.c */
static int
most_prevalent_count (Univcoord_T *values, int nvalues) {
  int max_count, count;
  Univcoord_T *ptr, *end, *first;

  assert(nvalues > 0);

  ptr = values;
  end = &(values[nvalues]);

  max_count = 1;
  while (ptr < end) {
    first = ptr;
    debug0(printf("DIAGONAL %u\n",*first));
    if (ptr + max_count - 1 >= end) {
      /* End of list fails */
      debug0(printf(" => Goes beyond end of list\n"));
      ptr = end;

    } else if (ptr[max_count - 1] != *first) {
      /* Fails */
      debug0(printf(" => Fails\n"));
      if (ptr[max_count - 1] != ptr[max_count - 2]) {
	/* Can jump immediately */
	ptr += max_count - 1;
      } else {
	/* Advance forward until we see a new value */
	while (ptr < end && *ptr == *first) {
	  ptr++;
	}
      }

    } else {
      /* Contender */
      ptr += max_count;
      while (ptr < end && *ptr == *first) {
	ptr++;
      }
      debug0(printf(" => Count %ld\n",ptr - first));
      if ((count = ptr - first) > max_count) {
	max_count = count;
      }
    }
  }

  return max_count;
}


/* Modified from kmer-search.c */
/* Overwrites values */
static void
find_prevalent (int *nloci, Univcoord_T *values, int nvalues, int subopt_count) {
  Univcoord_T *out, *ptr, *end, *first;
  /* int querystart, queryend; */
  /* int nmismatches5, nmismatches3; */

  out = &(values[0]);		/* overwrite */

  /* Get all values that have subopt_count or more */
  ptr = values;
  end = &(values[nvalues]);
  while (ptr < end) {
    first = ptr;
    debug0(printf("DIAGONAL %u\n",*first));
    if (ptr + subopt_count - 1 >= end) {
      /* End of list fails */
      debug0(printf(" => Goes beyond end of list\n"));
      ptr = end;

    } else if (ptr[subopt_count - 1] != *first) {
      /* Fails */
      debug0(printf(" => Fails\n"));
      if (ptr[subopt_count - 1] != ptr[subopt_count - 2]) {
	/* Can jump immediately */
	ptr += subopt_count - 1;
      } else {
	/* Advance forward until we see a new value */
	while (ptr < end && *ptr == *first) {
	  ptr++;
	}
      }

    } else {
      /* Contender */
      *out++ = *first;
      
      ptr += subopt_count;
      while (ptr < end && *ptr == *first) {
	ptr++;
      }
    }
  }

  *nloci = out - &(values[0]);
  return;
}


#define SUBOPT 2


/* Gets positions of 4-mers from low_position to high_position.  We do
   not worry about the small regions */
static int
get_exhaustive (bool *trimmedp, int *best_matchlength, int *nmismatches, Univcoord_T *diagonals,
		T this, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
		Univcoord_T **streamptr_alloc, int *streamsize_alloc, Mergeinfo_T mergeinfo,
		char *queryptr, int pos5, int pos3, int querylength,
		Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
		Compress_T query_compress, bool plusp, int genestrand, Genomebits_T genomebits,
		bool extend5p, bool trim5p, bool trim3p) {
  int nentries;
  Univcoord_T low_position_start, high_position_start, last_diagonal, *p, *q;
  int initptr, finalptr, ptr, n;
  uint16_t *saindex16;
  uint16_t *sasort16;
  Univcoord_T low_regioni, high_regioni, regioni;
  Univcoord_T region_term;

  List_T allocated = NULL, a;
  int streami = 0;

  /* _stream and _merged are aligned */
  Univcoord_T *merged_start, *_merged, *_stream, *last_space, *out;
  int streamsize, spacesize, last_spacesize;
  int nmerged, merged_low, merged_high;

  int max_count, subopt_count, lowi, highi;
  int pos;
  unsigned char oligo;
  
#ifdef MEASURE_WHOLE_END
  int best_mismatches, ref_nmismatches, mm;
#else
  int matchlength, trimpos5, trimpos3;
  int nmismatches_to_trimpos;
  int best_nmatches, nmatches;
#endif

#ifdef DEBUG2
  char nt[OLIGOSIZE + 1];
#endif


  *trimmedp = false;

  debug2(printf("Entering get_exhaustive with %.*s\n",pos3-pos5,&(queryptr[pos5])));
  last_space = streamspace_alloc;
  last_spacesize = 0;

  if (low_univdiagonal < (Univcoord_T) (querylength - pos5)) {
    low_position_start = 0U;
  } else {
    low_position_start = low_univdiagonal - querylength + pos5;
  }
  if (high_univdiagonal < (Univcoord_T) (querylength - pos5)) {
    high_position_start = 0U;
  } else {
    high_position_start = high_univdiagonal - querylength + pos5;
  }

  low_regioni = low_position_start/SARRAY16_LENGTH;
  high_regioni = high_position_start/SARRAY16_LENGTH;

  region_term = low_regioni * SARRAY16_LENGTH;
  for (regioni = low_regioni; regioni <= high_regioni; regioni++) {
    saindex16 = &(this->saindex16[regioni * SAINDEX16_LENGTH]);
    sasort16 = &(this->sasort16[region_term]);
    n = (region_term + SARRAY16_LENGTH < genomelength) ? SARRAY16_LENGTH: genomelength - region_term;

    if (pos5 + OLIGOSIZE <= pos3) {
      oligo = nt_oligo_byte(&(queryptr[pos5]),OLIGOSIZE - 1);
    }
    for (pos = pos5; pos + OLIGOSIZE <= pos3; pos++) {
      /* printf("Revising oligo with %c\n",queryptr[pos + OLIGOSIZE - 1]); */
      if ((oligo = revise_oligo_byte(oligo,queryptr[pos + OLIGOSIZE - 1])) == SAINDEX16_LENGTH - 1) {
	initptr = saindex16[oligo];
	finalptr = n;
      } else {
	initptr = saindex16[oligo];
	finalptr = saindex16[oligo+1];
      }

      lowi = initptr;
      if (regioni == low_regioni || regioni == low_regioni + 1) { /* Check over two regions, in case we are near a boundary */
	/* Filter diagonals */
	while (lowi < finalptr && (region_term + sasort16[lowi]) - pos + querylength < low_univdiagonal) {
	  lowi++;
	}
      }
      highi = finalptr - 1;
      if (regioni == high_regioni || regioni == high_regioni - 1) { /* Check over two regions, in case we are near a boundary */
	/* Filter diagonals */
	while (highi >= initptr && (region_term + sasort16[highi]) - pos + querylength > high_univdiagonal) {
	  highi--;
	}
      }

      debug2(oligo_nt(nt,oligo,4));
      debug2(printf("oligo %d (%s) => initptr %d, finalptr %d => lowi %d, highi %d\n",
		    oligo,nt,initptr,finalptr,lowi,highi));
      /* lowi and highi are inclusive */
      
      if (highi >= lowi) {
	streamsize = highi - lowi + 1;
	spacesize = streamsize + 1; /* Add extra value for Sedgesort */

	_stream = &(last_space[(last_spacesize + 31)/32 * 32]); /* Find an aligned position */
	if (_stream + spacesize >= streamspace_alloc + streamspace_max_alloc) {
	  /* About to overflow streamspace_alloc, so align instead */
	  debug2(printf("Space of %d would overflow, so allocating instead\n",spacesize));
	  MALLOC_ALIGN(_stream,spacesize*sizeof(Univcoord_T));
	  allocated = List_push(allocated,(void *) _stream);
	} else {
	  last_space = _stream;
	  last_spacesize = spacesize;
	}

	/* assert(streami < ((shortsplicedist + 65536)/65536 + 1) * querylength); */
	out = _stream;
	streamptr_alloc[streami] = _stream;
	streamsize_alloc[streami] = streamsize;

	for (ptr = lowi; ptr <= highi; ptr++) {
	  *out++ = region_term + sasort16[ptr] - pos + querylength; /* diagonal */
	}
	debug2(printf("At querypos %d, %d entries\n",pos,streamsize));

	streami++;
      }
    }

    region_term += SARRAY16_LENGTH;
  }

  /* Merge streams */
  if (streami == 0) {
    return 0;

  } else {
#ifdef LARGE_GENOMES
    _merged = Merge_diagonals_uint8(&nmerged,streamptr_alloc,streamsize_alloc,/*nstreams*/streami,mergeinfo);
#else
    _merged = Merge_diagonals_uint4(&nmerged,streamptr_alloc,streamsize_alloc,/*nstreams*/streami,mergeinfo);
#endif

    for (a = allocated; a != NULL; a = List_next(a)) {
      _stream = (Univcoord_T *) List_head(a);
      FREE_ALIGN(_stream);
    }
    List_free(&allocated);


    /* Remove univdiagonals outside the bounds */
    merged_low = 0;
    while (merged_low < nmerged && _merged[merged_low] < low_univdiagonal) {
      merged_low++;
    }

    merged_high = nmerged - 1;
    while (merged_high >= 0 && _merged[merged_high] > high_univdiagonal) {
      merged_high--;
    }

    merged_start = &(_merged[merged_low]);
    nmerged = merged_high - merged_low + 1;

    max_count = most_prevalent_count(merged_start,nmerged);
    debug2(printf("max_count is %d\n",max_count));
    if ((subopt_count = max_count - SUBOPT) < 1) {
      subopt_count = 1;
    }
    find_prevalent(&nentries,merged_start,nmerged,subopt_count);


  /* Remove repetitive diagonals */
#ifdef DEBUG2
    printf("\n");
    printf("%d entries before removing repetitive:",nentries);
    for (int i = 0; i < nentries; i++) {
      printf(" %u",merged_start[i]);
      assert(merged_start[i] >= low_univdiagonal && merged_start[i] <= high_univdiagonal);
    }
    printf("\n");
#endif

    /* Skip repetitive diagonals and store those with the fewest mismatches */

    q = diagonals;
    last_diagonal = 0;
    best_nmatches = 0;
    *best_matchlength = 0;

    for (p = &(merged_start[0]); p < &(merged_start[nentries]); p++) {
      if (*p == last_diagonal) {
	/* Skip: repetitive */
      } else {
	if (extend5p == true) {
	  /* Perform medial trim first, where we expect mismatches for an insertion */
	  if (trim3p == false) {
	    trimpos3 = pos3;
	  } else {
	    trimpos3 = Genomebits_trim_qend(&nmismatches_to_trimpos,query_compress,genomebits,
					    /*univdiagonal*/(*p),querylength,pos5,pos3,plusp,genestrand);
	  }
      
	  if (trimpos3 <= pos5) {
	    trimpos5 = pos5;
	  } else if (trim5p == false) {
	    trimpos5 = pos5;
	  } else {
	    trimpos5 = Genomebits_trim_qstart(&nmismatches_to_trimpos,query_compress,genomebits,
					      /*univdiagonal*/(*p),querylength,pos5,trimpos3,plusp,genestrand);
	  }

	} else {
	  /* Perform medial trim first, where we expect mismatches for an insertion */
	  if (trim5p == false) {
	    trimpos5 = pos5;
	  } else {
	    /* printf("Calling Genomebits_trim_qstart with %u, pos5 %d, pos3 %d\n",(*p),pos5,pos3); */
	    trimpos5 = Genomebits_trim_qstart(&nmismatches_to_trimpos,query_compress,genomebits,
					      /*univdiagonal*/(*p),querylength,pos5,pos3,plusp,genestrand);
	  }

	  if (trimpos5 >= pos3) {
	    trimpos3 = pos3;
	  } else if (trim3p == false) {
	    trimpos3 = pos3;
	  } else {
	    /* printf("Calling Genomebits_trim_qend with %u, trimpos5 %d, pos3 %d\n",(*p),trimpos5,pos3); */
	    trimpos3 = Genomebits_trim_qend(&nmismatches_to_trimpos,query_compress,genomebits,
					    /*univdiagonal*/(*p),querylength,trimpos5,pos3,plusp,genestrand);
	  }
	}

	if ((matchlength = trimpos3 - trimpos5) <= 0) {
	  /* Skip */
	  debug2(printf("%u has trimpos %d..%d\n",(*p),trimpos5,trimpos3));

	} else if ((nmatches = matchlength - nmismatches_to_trimpos) > best_nmatches) {
	  debug2(printf("%u has %d mismatches to trimpos %d..%d => %d matches, matchlength %d => more so resetting q\n",
			(*p),nmismatches_to_trimpos,trimpos5,trimpos3,nmatches,matchlength));
	  /* New winner */
	  q = diagonals;
	  *q++ = *p;
	  *best_matchlength = matchlength;
	  best_nmatches = nmatches;
	  *trimmedp = false;
	  if (trimpos5 != pos5 || trimpos3 != pos3) {
	    *trimmedp = true;
	  }

	} else if (nmatches == best_nmatches && matchlength < (*best_matchlength)) {
	  debug2(printf("%u has %d mismatches to trimpos %d..%d => %d matches, matchlength %d => shorter so resetting q\n",
			(*p),nmismatches_to_trimpos,trimpos5,trimpos3,nmatches,matchlength));
	  /* New winner */
	  q = diagonals;
	  *q++ = *p;
	  *best_matchlength = matchlength;
	  *trimmedp = false;
	  if (trimpos5 != pos5 || trimpos3 != pos3) {
	    *trimmedp = true;
	  }

	} else if (nmatches == best_nmatches && matchlength == (*best_matchlength)) {
	  debug2(printf("%u has %d mismatches to trimpos %d..%d => %d matches, matchlength %d => same so saving and advancing q\n",
			(*p),nmismatches_to_trimpos,trimpos5,trimpos3,nmatches,matchlength));
	  *q++ = *p;
	  /* Tie */
	  if (trimpos5 != pos5 || trimpos3 != pos3) {
	    *trimmedp = true;
	  }
	  
	} else {
	  debug2(printf("%u has %d mismatches to trimpos %d..%d => %d matches, matchlength %d => fewer so skipping\n",
			(*p),nmismatches_to_trimpos,trimpos5,trimpos3,nmatches,matchlength));
	  /* Worse */
	}
      }
      last_diagonal = *p;
    }

    FREE_ALIGN(_merged);

    assert((*best_matchlength) >= best_nmatches);
    *nmismatches = (*best_matchlength) - best_nmatches;
    nentries = q - diagonals;

#ifdef DEBUG2 
    printf("\n");
    printf("%d entries:",nentries);
    for (int i = 0; i < nentries; i++) {
      printf(" %u",diagonals[i]);
    }
    printf("\n");
#endif
    
    return nentries;
  }
}


#ifdef CHECK_ASSERTIONS
static void
check_ascending (Univcoord_T *coords, int n) {
  Univcoord_T prev_coord;
  int i;

  prev_coord = coords[0];
  for (i = 1; i < n; i++) {
    if (coords[i] <= prev_coord) {
      printf("Expecting ascending, but at %d, got %u <= %u\n",
	     i,coords[i],prev_coord);
      abort();
    }
    prev_coord = coords[i];
  }
 
  return;
}
#endif



/* localdb_alloc is needed for the results of the first and last sasort entries in get_by_ends */
int
Localdb_get (bool *sortedp, bool *trimmedp,
	     int *matchlength, int *nmismatches, Univcoord_T *diagonals_alloc, 
	     T this, unsigned short *localdb_alloc, Stage1_T stage1,
	     char *queryptr, int pos5, int pos3, int querylength,
	     Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
	     Compress_T query_compress, bool plusp, int genestrand,
	     Genomebits_T genomebits, int nmismatches_allowed,
	     bool extend5p, bool trim5p, bool trim3p) {
  int nhits;

  int streamspace_max_alloc = stage1->streamspace_max_alloc;
  Univcoord_T *streamspace_alloc = stage1->streamspace_alloc;
  Univcoord_T **streamptr_alloc = stage1->streamptr_alloc;
  int *streamsize_alloc = stage1->streamsize_alloc;
  Mergeinfo_T mergeinfo = stage1->mergeinfo;

  debug(printf("\nEntered Localdb_get with low_univdiagonal %u, high_univdiagonal %u, plusp %d, pos5 %d, pos3 %d, query %s, compress_fwdp %d\n",
	       low_univdiagonal,high_univdiagonal,plusp,pos5,pos3,queryptr,Compress_fwdp(query_compress)));
  assert(Compress_fwdp(query_compress) == plusp);

  *matchlength = pos3 - pos5;	/* Needed for get_by_ends */

  if (pos5 >= pos3) {
    return 0;

  } else if (high_univdiagonal < low_univdiagonal) {
    return 0;

  } else if ((nhits = get_by_sarray(&(*matchlength),&(*diagonals_alloc),this,
				    queryptr,pos5,pos3,querylength,
				    /*min_matchlength*/(pos3 - pos5),low_univdiagonal,high_univdiagonal,
				    query_compress,plusp,genomebits)) > 0) {
    /* Previously required pos3 - pos5 < index1part + index1interval because otherwise we would have found it from indexdb */
    /* assert(pos3 - pos5 < index1part + index1interval); */

    debug(printf("Localdb_get using get_by_sarray (and no mismatches allowed).  Returning %d hits\n",nhits));
    /* Not guaranteed to be ascending */
    *nmismatches = 0;
    *sortedp = false;
    *trimmedp = false;		/* Entire string matches */
    return nhits;

  } else if (pos3 - pos5 > 2*OLIGOSIZE &&
	     (nhits = get_by_ends(&(*nmismatches),&(*diagonals_alloc),this,localdb_alloc,
				  queryptr,pos5,pos3,querylength,
				  low_univdiagonal,high_univdiagonal,query_compress,
				  plusp,genestrand,genomebits,nmismatches_allowed)) > 0) {
    /* If equal to OLIGOSIZE, we would have found it using get_by_sarray */
    /* If <= 2*OLIGOSIZE, and there is a mismatch we won't find it, since at least one end will have a mismatch */
    debug(printf("Localdb_get using get_by_ends and %d mismatches allowed.  Returning %d hits\n",
		 nmismatches_allowed,nhits));
    *sortedp = true;
    *trimmedp = false;		/* Both ends match */
    return nhits;

  } else {
    debug(printf("Localdb_get using get_exhaustive and %d mismatches allowed.  Returning %d hits\n",
		 (pos3 - pos5)/2,nhits));
    nhits = get_exhaustive(&(*trimmedp),&(*matchlength),&(*nmismatches),&(*diagonals_alloc),this,
			   streamspace_max_alloc,streamspace_alloc,streamptr_alloc,streamsize_alloc,
			   mergeinfo,queryptr,pos5,pos3,querylength,low_univdiagonal,high_univdiagonal,
			   query_compress,plusp,genestrand,genomebits,extend5p,trim5p,trim3p);
    *sortedp = true;
    return nhits;
  }
}


Univcoord_T
Localdb_get_one_low (Univcoord_T *diagonals_alloc, 
		     T this, unsigned short *localdb_alloc, Stage1_T stage1,
		     char *queryptr, int querylength,
		     Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
		     Compress_T query_compress, bool plusp, int genestrand,
		     Genomebits_T genomebits, int nmismatches_allowed) {
  Univcoord_T lowest;
  bool sortedp, trimmedp;
  int matchlength, nmismatches;
  int ndiagonals, i;

  if ((ndiagonals = Localdb_get(&sortedp,&trimmedp,&matchlength,&nmismatches,diagonals_alloc,
				this,localdb_alloc,stage1,queryptr,
				/*pos5*/0,/*pos3*/querylength,querylength,
				low_univdiagonal,high_univdiagonal,
				query_compress,plusp,genestrand,genomebits,nmismatches_allowed,
				/*extend5p*/true,/*trim5p*/true,/*trim3p*/true)) == 0) {
    return (Univcoord_T) 0;

  } else if (ndiagonals == 1) {
    debug5(printf("Returning %u with matchlength %d, nmismatches %d\n",
		  diagonals_alloc[0],matchlength,nmismatches));
    return diagonals_alloc[0];

  } else if (sortedp == true) {
    debug5(printf("Returning %u with matchlength %d, nmismatches %d\n",
		  diagonals_alloc[0],matchlength,nmismatches));
    return diagonals_alloc[0];

  } else {
    lowest = diagonals_alloc[0];
    for (i = 1; i < ndiagonals; i++) {
      if (diagonals_alloc[i] < lowest) {
	lowest = diagonals_alloc[i];
      }
    }

    debug5(printf("Returning %u with matchlength %d, nmismatches %d\n",
		  lowest,matchlength,nmismatches));
    return lowest;
  }
}


Univcoord_T
Localdb_get_one_high (Univcoord_T *diagonals_alloc, 
		     T this, unsigned short *localdb_alloc, Stage1_T stage1,
		     char *queryptr, int querylength,
		     Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
		     Compress_T query_compress, bool plusp, int genestrand,
		     Genomebits_T genomebits, int nmismatches_allowed) {
  Univcoord_T highest;
  bool sortedp, trimmedp;
  int matchlength, nmismatches;
  int ndiagonals, i;

  if ((ndiagonals = Localdb_get(&sortedp,&trimmedp,&matchlength,&nmismatches,diagonals_alloc,
				this,localdb_alloc,stage1,queryptr,
				/*pos5*/0,/*pos3*/querylength,querylength,
				low_univdiagonal,high_univdiagonal,
				query_compress,plusp,genestrand,genomebits,nmismatches_allowed,
				/*extend5p*/true,/*trim5p*/true,/*trim3p*/true)) == 0) {
    
    return (Univcoord_T) 0;

  } else if (ndiagonals == 1) {
    debug5(printf("Returning %u with matchlength %d, nmismatches %d\n",
		  diagonals_alloc[0],matchlength,nmismatches));
    return diagonals_alloc[0];

  } else if (sortedp == true) {
    debug5(printf("Returning %u with matchlength %d, nmismatches %d\n",
		  diagonals_alloc[ndiagonals - 1],matchlength,nmismatches));
    return diagonals_alloc[ndiagonals - 1];

  } else {
    highest = diagonals_alloc[0];
    for (i = 1; i < ndiagonals; i++) {
      if (diagonals_alloc[i] > highest) {
	highest = diagonals_alloc[i];
      }
    }

    debug5(printf("Returning %u with matchlength %d, nmismatches %d\n",
		  highest,matchlength,nmismatches));
    return highest;
  }
}



T
Localdb_new (char *genomesubdir, char *fileroot,
	     Access_mode_T localdb_access, bool sharedp,
	     bool multiple_sequences_p, bool preload_shared_memory_p, bool unload_shared_memory_p) {
  T new;

  char *saindex16_filename, *sarray16_filename, *sarray8_filename, *sasort16_filename;
  char *comma;
  double seconds0, seconds1, seconds2;

  saindex16_filename = (char *) MALLOC((strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				       strlen(".saindex16")+1)*sizeof(char));
  sarray16_filename = (char *) MALLOC((strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				       strlen(".sarray16")+1)*sizeof(char));
  sarray8_filename = (char *) MALLOC((strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				      strlen(".sarray8")+1)*sizeof(char));
  sasort16_filename = (char *) MALLOC((strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				       strlen(".sasort16")+1)*sizeof(char));
  sprintf(saindex16_filename,"%s/%s.saindex16",genomesubdir,fileroot);
  sprintf(sarray16_filename,"%s/%s.sarray16",genomesubdir,fileroot);
  sprintf(sarray8_filename,"%s/%s.sarray8",genomesubdir,fileroot);
  sprintf(sasort16_filename,"%s/%s.sasort16",genomesubdir,fileroot);

  if (Access_file_exists_p(saindex16_filename) == false ||
      Access_file_exists_p(sarray16_filename) == false ||
      Access_file_exists_p(sarray8_filename) == false ||
      Access_file_exists_p(sasort16_filename) == false) {
    fprintf(stderr,"Unable to find localdb files %s, %s, %s, or %s.  Will align without them.\n",
	    saindex16_filename,sarray16_filename,sarray8_filename,sasort16_filename);
    FREE(sasort16_filename);
    FREE(sarray8_filename);
    FREE(sarray16_filename);
    FREE(saindex16_filename);
    return (T) NULL;

  } else if (localdb_access == USE_ALLOCATE) {
    new = (T) MALLOC(sizeof(*new));
    fprintf(stderr,"Allocating memory for localdb...");

    if (
#ifndef HAVE_MMAP
	0 &&
#endif
	multiple_sequences_p == false && preload_shared_memory_p == false && unload_shared_memory_p == false) {
#ifdef HAVE_MMAP
      new->saindex16 = (uint16_t *) Access_mmap(&new->saindex16_fd,&new->saindex16_len,&seconds0,
					       saindex16_filename,/*randomp*/false);
      new->sarray16 = (uint16_t *) Access_mmap(&new->sarray16_fd,&new->sarray16_len,&seconds1,
					       sarray16_filename,/*randomp*/false);
      new->sarray8 = (uint8_t *) Access_mmap(&new->sarray8_fd,&new->sarray8_len,&seconds2,
					     sarray8_filename,/*randomp*/false);
      new->sasort16 = (uint16_t *) Access_mmap(&new->sasort16_fd,&new->sasort16_len,&seconds1,
					       sasort16_filename,/*randomp*/false);
      new->localdb_access = MMAPPED;
#endif
    } else if (sharedp == true) {
      new->saindex16 = (uint16_t *) Access_allocate_shared(&new->localdb_access,&new->saindex16_shmid,&new->saindex16_key,
							  &new->saindex16_fd,&new->saindex16_len,&seconds0,
							  saindex16_filename,sizeof(uint16_t));
      new->sarray16 = (uint16_t *) Access_allocate_shared(&new->localdb_access,&new->sarray16_shmid,&new->sarray16_key,
							  &new->sarray16_fd,&new->sarray16_len,&seconds1,
							  sarray16_filename,sizeof(uint16_t));
      new->sarray8 = (uint8_t *) Access_allocate_shared(&new->localdb_access,&new->sarray8_shmid,&new->sarray8_key,
							  &new->sarray8_fd,&new->sarray8_len,&seconds2,
							  sarray8_filename,sizeof(uint8_t));
      new->sasort16 = (uint16_t *) Access_allocate_shared(&new->localdb_access,&new->sasort16_shmid,&new->sasort16_key,
							  &new->sasort16_fd,&new->sasort16_len,&seconds1,
							  sasort16_filename,sizeof(uint16_t));
    } else {
      new->saindex16 = (uint16_t *) Access_allocate_private(&new->localdb_access,&new->saindex16_len,&seconds0,
							   saindex16_filename,sizeof(uint16_t));
      new->sarray16 = (uint16_t *) Access_allocate_private(&new->localdb_access,&new->sarray16_len,&seconds1,
							   sarray16_filename,sizeof(uint16_t));
      new->sarray8 = (uint8_t *) Access_allocate_private(&new->localdb_access,&new->sarray8_len,&seconds2,
							 sarray8_filename,sizeof(uint8_t));
      new->sasort16 = (uint16_t *) Access_allocate_private(&new->localdb_access,&new->sasort16_len,&seconds1,
							   sasort16_filename,sizeof(uint16_t));
    }

    /* Don't check for sarray8, since it can be empty when genomelength < 65536 */
    if (new->saindex16 == NULL || new->sarray16 == NULL ||
	/*new->sarray8 == NULL ||*/ new->sasort16 == NULL) {
      fprintf(stderr,"insufficient memory\n");
      exit(9);
    } else {
      comma = Genomicpos_commafmt(new->saindex16_len + new->sarray16_len + new->sarray8_len + new->sasort16_len);
      fprintf(stderr,"done (%s bytes",comma);
      FREE(comma);

      if (multiple_sequences_p == true) {
	fprintf(stderr,", %.2f sec",seconds0 + seconds1 + seconds2);
      }
      fprintf(stderr,")\n");
    }

#ifdef HAVE_MMAP
  } else if (localdb_access == USE_MMAP_PRELOAD || localdb_access == USE_MMAP_ONLY) {
    fprintf(stderr,"Memory mapping localdb files...");
    new = (T) MALLOC(sizeof(*new));
    new->saindex16 = (uint16_t *) Access_mmap(&new->saindex16_fd,&new->saindex16_len,&seconds0,
					     saindex16_filename,/*randomp*/true);
    new->sarray16 = (uint16_t *) Access_mmap(&new->sarray16_fd,&new->sarray16_len,&seconds1,
					     sarray16_filename,/*randomp*/true);
    new->sarray8 = (uint8_t *) Access_mmap(&new->sarray8_fd,&new->sarray8_len,&seconds2,
					     sarray8_filename,/*randomp*/true);
    new->sasort16 = (uint16_t *) Access_mmap(&new->sasort16_fd,&new->sasort16_len,&seconds1,
					     sasort16_filename,/*randomp*/true);

    if (new->saindex16 == NULL || new->sarray16 == NULL || new->sarray8 == NULL || new->sasort16 == NULL) {
      fprintf(stderr,"Insufficient memory for mmap (will use disk file instead, but program will be slow)\n");
      new->localdb_access = FILEIO;
    } else {
      comma = Genomicpos_commafmt(new->saindex16_len + new->sarray16_len + new->sarray8_len + new->sasort16_len);
      fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds0 + seconds1 + seconds2);
      FREE(comma);
      new->localdb_access = MMAPPED;
    }
#endif

  } else {
    fprintf(stderr,"Don't recognize localdb_access %d\n",localdb_access);
    abort();
  }

  FREE(sasort16_filename);
  FREE(sarray8_filename);
  FREE(sarray16_filename);
  FREE(saindex16_filename);

  return new;
}
