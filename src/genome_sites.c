static char rcsid[] = "$Id: genome_sites.c 225220 2022-11-01 22:46:01Z twu $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "genome_sites.h"

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>		/* For tolower() */

#include "assert.h"
#include "except.h"
#include "maxent_hr.h"
#include "popcount.h"
#include "dinucl_bits.h"


#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif


/* Splice sites */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

#ifdef DEBUG2
static void
write_chars_comp (UINT4 high, UINT4 low, UINT4 flags) {
  char Buffer[33];
  int i;

  Buffer[32] = '\0';
  /* printf("%08X %08X %08X => ",high,low,flags); */
  for (i = 0; i < 16; i++) {
    switch (low & 3U) {
    case 0U: Buffer[i] = 'A'; break;
    case 1U: Buffer[i] = 'C'; break;
    case 2U: Buffer[i] = 'G'; break;
    case 3U: Buffer[i] = 'T'; break;
    default: abort();
    }
    low >>= 2;
  }
  for ( ; i < 32; i++) {
    switch (high & 3U) {
    case 0U: Buffer[i] = 'A'; break;
    case 1U: Buffer[i] = 'C'; break;
    case 2U: Buffer[i] = 'G'; break;
    case 3U: Buffer[i] = 'T'; break;
    default: abort();
    }
    high >>= 2;
  }
  for (i = 0; i < 32; i++) {
    if ((flags & 1U) == 1U) {
      Buffer[i] = 'N';
    }
    flags >>= 1;
  }

  printf("%s",Buffer);
  return;
}
#endif


#ifdef DEBUG2
static void
Genome_print_blocks (Genomecomp_T *blocks, Univcoord_T startpos, Univcoord_T endpos) {
  /* Chrpos_T length = endpos - startpos; */
  Univcoord_T startblock, endblock, ptr;
  int startdiscard, enddiscard;
  Genomecomp_T high, low, flags;
  int i;

  /* sequence = (char *) CALLOC(length+1,sizeof(char)); */

  ptr = startblock = startpos/32U*3;
  endblock = endpos/32U*3;
  startdiscard = startpos % 32;
  enddiscard = endpos % 32;
  
  /*      high: 9F61B62A  low: 6D68A157  flags: 00000000 */
  printf("                                              \t");
  printf("%llu\t",(unsigned long long) startblock/3*32U);
  for (i = 0; i < startdiscard; i++) {
    printf("*");
  }
  printf("\n");

  for (ptr = startblock ; ptr <= endblock; ptr += 3) {
#ifdef WORDS_BIGENDIAN
    high = Bigendian_convert_uint(blocks[ptr]);
    low = Bigendian_convert_uint(blocks[ptr+1]);
    flags = Bigendian_convert_uint(blocks[ptr+2]);
#else
    high = blocks[ptr]; low = blocks[ptr+1]; flags = blocks[ptr+2];
#endif
    printf("high: %08X  low: %08X  flags: %08X\t",high,low,flags);
    printf("%llu\t",(unsigned long long) ptr/3*32U);
    write_chars_comp(high,low,flags);
    printf("\n");
  }

  /*      high: 9F61B62A  low: 6D68A157  flags: 00000000 */
  printf("                                              \t");
  printf("%llu\t",(unsigned long long) (endblock+3)/3*32U);
  for (i = 0; i < enddiscard; i++) {
    printf(" ");
  }
  for ( ; i < 32; i++) {
    printf("*");
  }
  printf("\n");


  return;
}
#endif


#if 0
static void
Genome_print_blocks_chrpos (Genomecomp_T *blocks, Univcoord_T startpos, Univcoord_T endpos, Univcoord_T chroffset) {
  /* Chrpos_T length = endpos - startpos; */
  Univcoord_T startblock, endblock, ptr;
  int startdiscard, enddiscard;
  Genomecomp_T high, low, flags;
  int i;

  /* sequence = (char *) CALLOC(length+1,sizeof(char)); */

  ptr = startblock = startpos/32U*3;
  endblock = endpos/32U*3;
  startdiscard = startpos % 32;
  enddiscard = endpos % 32;
  
  /*      high: 9F61B62A  low: 6D68A157  flags: 00000000 */
  printf("                                              \t");
  printf("%llu\t%llu\t",(unsigned long long) startblock,(unsigned long long) startblock/3*32U - chroffset);
  for (i = 0; i < startdiscard; i++) {
    printf("*");
  }
  printf("\n");

  for (ptr = startblock ; ptr <= endblock; ptr += 3) {
#ifdef WORDS_BIGENDIAN
    high = Bigendian_convert_uint(blocks[ptr]);
    low = Bigendian_convert_uint(blocks[ptr+1]);
    flags = Bigendian_convert_uint(blocks[ptr+2]);
#else
    high = blocks[ptr]; low = blocks[ptr+1]; flags = blocks[ptr+2];
#endif
    printf("high: %08X  low: %08X  flags: %08X\t",high,low,flags);
    printf("%llu\t%llu\t",(unsigned long long) ptr,(unsigned long long) ptr/3*32U - chroffset);
    write_chars_comp(high,low,flags);
    printf("\n");
  }

  /*      high: 9F61B62A  low: 6D68A157  flags: 00000000 */
  printf("                                              \t");
  printf("%llu\t%llu\t",(unsigned long long) endblock+3,(unsigned long long) (endblock+3)/3*32U - chroffset);
  for (i = 0; i < enddiscard; i++) {
    printf(" ");
  }
  for ( ; i < 32; i++) {
    printf("*");
  }
  printf("\n");


  return;
}
#endif


#define T Genome_T

static T genome;
static T genomealt;


void
Genome_sites_setup (T genome_in, T genomealt_in) {
  genome = genome_in;
  genomealt = genomealt_in;

  return;
}


/*                 76543210 */
#define HIGH_BIT 0x80000000

#define clear_start(diff,startdiscard) (diff & (~0U << (startdiscard)))
#define clear_end(diff,enddiscard) (diff & ~(~0U << (enddiscard)))

#define clear_start_mask(startdiscard) (~0U << (startdiscard))
#define clear_end_mask(enddiscard) (~(~0U << (enddiscard)))

/* Same speed: clear_highbit(diff,relpos) diff -= (HIGH_BIT >> relpos) */
/* Note: xor assumes that bit at relpos was on */
#define clear_highbit(diff,relpos) diff ^= (HIGH_BIT >> relpos)

/* Slower: clear_lowbit(diff,relpos) diff -= (1 << relpos) */
#define clear_lowbit(diff,relpos) diff &= (diff - 1);


/************************************************************************
 *   Splice sites
 ************************************************************************/


static Genomecomp_T
block_find (Genomecomp_T *high_halfsite, Genomecomp_T *low_halfsite, Genomecomp_T *ptr, const Genomecomp_T *splicesite_bits) {
  Genomecomp_T found, compare, flags;

  /* High */
  debug2(printf("Evaluating high %08X and low %08X\n",*ptr,ptr[1]));
#ifdef WORDS_BIGENDIAN
  compare = Bigendian_convert_uint(*ptr++);
#else
  compare = *ptr++;
#endif
  /* Get high_halfsite bit */
  found = splicesite_bits[compare >> 16];
  *high_halfsite = (found & 0x100) >> 8;
  found = (found << 24) | 0x00FFFFFF;
  debug2(printf("  splicesite_bits_3: %08X, high_halfsite %d\n",
		splicesite_bits[compare >> 16] << 24,*high_halfsite));

  found &= (splicesite_bits[compare & 0x0000FFFF] << 16) | 0xFE00FFFF; /* Use FE to allow for high bit */
  debug2(printf("  splicesite_bits_2: %08X\n",splicesite_bits[compare & 0x0000FFFF] << 16));

  /* Low */
#ifdef WORDS_BIGENDIAN
  compare = Bigendian_convert_uint(*ptr++);
#else
  compare = *ptr++;
#endif
  found &= (splicesite_bits[compare >> 16] << 8) | 0xFFFE00FF;
  debug2(printf("  splicesite_bits_1: %08X\n",splicesite_bits[compare >> 16] << 8));

  found &= splicesite_bits[compare & 0x0000FFFF] | 0xFFFFFE00;
  debug2(printf("  splicesite_bits_0: %08X, low_halfsite %d\n",
		splicesite_bits[compare & 0x0000FFFF],found & 0x1));


  /* Flags: N is considered a mismatch */
  debug2(printf("Marking flags: genome %08X ",*ptr));
#ifdef WORDS_BIGENDIAN
  flags = Bigendian_convert_uint(*ptr);
#else
  flags = (*ptr);
#endif

  found &= ~flags;
  found &= ~(flags << 1);

  *low_halfsite = found & 0x00000001;

  *high_halfsite &= ~(flags >> 31);

  /* splicesite_bits are 1-based, so need to reduce values by 1 */
  /* No need to clear top bit */
  found >>= 1;
  debug2(printf(" => found %08X, high_halfsite %d, low_halfsite %d\n",found,*high_halfsite,*low_halfsite));

  return found;
}



static Genomecomp_T
block_find_snp (Genomecomp_T *high_halfsite, Genomecomp_T *low_halfsite, Genomecomp_T *ref_ptr, Genomecomp_T *alt_ptr,
		const Genomecomp_T *splicesite_bits) {
  Genomecomp_T found, ref, alt, flags;

  /* High */
  debug2(printf("Evaluating ref high %08X and low %08X, and alt high %08X and low %08X\n",
		*ref_ptr,ref_ptr[1],*alt_ptr,alt_ptr[1]));
#ifdef WORDS_BIGENDIAN
  ref = Bigendian_convert_uint(*ref_ptr++);
  alt = Bigendian_convert_uint(*alt_ptr++);
#else
  ref = *ref_ptr++;
  alt = *alt_ptr++;
#endif
  /* Get high_halfsite bit */
  found = splicesite_bits[ref >> 16] | splicesite_bits[alt >> 16];
  *high_halfsite = (found & 0x100) >> 8;
  found = (found << 24) | 0x00FFFFFF;
  debug2(printf("  splicesite_bits_3: %08X | %08X, high_halfsite %d\n",
		splicesite_bits[ref >> 16] << 24,splicesite_bits[alt >> 16] << 24,*high_halfsite));

 /* Use FE to allow for high bit */
  found &= ((splicesite_bits[ref & 0x0000FFFF] | splicesite_bits[alt & 0x0000FFFF]) << 16) | 0xFE00FFFF;
  debug2(printf("  splicesite_bits_2: %08X | %08X\n",
		splicesite_bits[ref & 0x0000FFFF] << 16,splicesite_bits[alt & 0x0000FFFF] << 16));

  /* Low */
#ifdef WORDS_BIGENDIAN
  ref = Bigendian_convert_uint(*ref_ptr++);
  alt = Bigendian_convert_uint(*alt_ptr++);
#else
  ref = *ref_ptr++;
  alt = *alt_ptr++;
#endif
  found &= ((splicesite_bits[ref >> 16] | splicesite_bits[alt >> 16]) << 8) | 0xFFFE00FF;
  debug2(printf("  splicesite_bits_1: %08X | %08X\n",
		splicesite_bits[ref >> 16] << 8,splicesite_bits[alt >> 16] << 8));

  found &= (splicesite_bits[ref & 0x0000FFFF] | splicesite_bits[alt & 0x0000FFFF]) | 0xFFFFFE00;
  debug2(printf("  splicesite_bits_0: %08X | %08X, low_halfsite %d\n",
		splicesite_bits[ref & 0x0000FFFF],splicesite_bits[ref & 0x0000FFFF],found & 0x1));


  /* Handle only reference flags, which indicate N */
  /* Reference flags: N is considered a mismatch */
  debug2(printf("Marking flags: genome %08X ",*ref_ptr));
#ifdef WORDS_BIGENDIAN
  flags = Bigendian_convert_uint(*ref_ptr);
#else
  flags = (*ref_ptr);
#endif

  found &= ~flags;
  found &= ~(flags << 1);

  *low_halfsite = found & 0x00000001;

  *high_halfsite &= ~(flags >> 31);

  /* splicesite_bits are 1-based, so need to reduce values by 1 */
  /* No need to clear top bit */
  found >>= 1;
  debug2(printf(" => found %08X, high_halfsite %d, low_halfsite %d\n",found,*high_halfsite,*low_halfsite));

  return found;
}


/* Fills sites and site_knowni with nfound entries.
   Integrates old_knownpos and old_knowni from knownsplicing. */
static int
splicesites (int *sites, int *knowni, int *old_knownpos, int *old_knowni,
	     Univcoord_T left, int pos5, int pos3,
	     const Genomecomp_T *splicesite_bits, int splicepos_offset) {
  int nfound = 0, offset;
  int startdiscard, enddiscard;
  Univcoord_T startblocki, endblocki;
  Genomecomp_T *ptr, *altptr, *end;
  Genomecomp_T found;
  Genomecomp_T high_halfsite, low_halfsite, prev_high_halfsite;
  int pos;
#ifdef HAVE_BUILTIN_CTZ
  int relpos;
#else
  Genomecomp_T lowbit;
#endif

  debug2(
	printf("\n\n");
	printf("Genome (in splicesites):\n");
	Genome_print_blocks(genome->blocks,left+pos5,left+pos3);
	printf("\n");
	);

  /* nshift = left % 32; */

  startblocki = (left+pos5)/32U*3;
  endblocki = (left+pos3)/32U*3;
  startdiscard = (left+pos5) % 32;
  enddiscard = (left+pos3) % 32;

  offset = -startdiscard + pos5 + splicepos_offset;
  
  debug2(printf("left = %llu, pos5 = %d, pos3 = %d, startblocki = %llu, endblocki = %llu\n",
		(unsigned long long) left,pos5,pos3,(unsigned long long) startblocki,(unsigned long long) endblocki));
  debug2(printf("startdiscard = %d, enddiscard = %d\n",startdiscard,enddiscard));

  if (endblocki == startblocki) {
    /* Advance old_knownpos past pos5 */
    while (*old_knownpos < pos5) {
      debug2(printf("Skipping old_knownpos #%d at %d < pos5 %d\n",*old_knowni,*old_knownpos,pos5));
      old_knowni++;
      old_knownpos++;
    }

    if (genomealt != genome) {
      found = block_find_snp(&high_halfsite,&low_halfsite,
			     &(genome->blocks[startblocki]),&(genomealt->blocks[startblocki]),splicesite_bits);
    } else {
      found = block_find(&high_halfsite,&low_halfsite,&(genome->blocks[startblocki]),splicesite_bits);
    }
    found = clear_start(found,startdiscard);
    found = clear_end(found,enddiscard);
    debug2(printf("adding masks %08X and %08x\n",clear_start_mask(startdiscard),clear_end_mask(enddiscard)));

    while (found != 0U) {
#ifdef HAVE_BUILTIN_CTZ
      pos = offset + (relpos = __builtin_ctz(found));
      while (*old_knownpos < pos) {
	debug2(printf("Adding old_knownpos #%d at %d < pos %d\n",*old_knowni,*old_knownpos,pos));
	knowni[nfound] = *old_knowni++;
	sites[nfound++] = *old_knownpos++;
      }
      if (*old_knownpos == pos) {
	debug2(printf("Adding old_knownpos #%d at %d == pos %d\n",*old_knowni,*old_knownpos,pos));
	knowni[nfound] = *old_knowni++;
	sites[nfound++] = *old_knownpos++;
      } else {
	knowni[nfound] = -1;
	sites[nfound++] = pos;
      }
      clear_lowbit(found,relpos);
#else      
      debug2(printf("found is %08X, -found & found is %08X\n",found,-found & found));
      pos = offset + mod_37_bit_position[(lowbit = -found & found) % 37];
      while (*old_knownpos < pos) {
	debug2(printf("Adding old_knownpos #%d at %d < pos %d\n",*old_knowni,*old_knownpos,pos));
	knowni[nfound] = *old_knowni++;
	sites[nfound++] = *old_knownpos++;
      }
      if (*old_knownpos == pos) {
	debug2(printf("Adding old_knownpos #%d at %d == pos %d\n",*old_knowni,*old_knownpos,pos));
	knowni[nfound] = *old_knowni++;
	sites[nfound++] = *old_knownpos++;
      } else {
	knowni[nfound] = -1;
	sites[nfound++] = pos;
      }
      found -= lowbit;
#endif
      debug2(printf("found is %08X => offset %llu + relpos %d\n",found,(unsigned long long) offset,relpos));
    }

    /* Add old_knownpos to pos3 */
    while (*old_knownpos < pos3) {
      debug2(printf("Adding old_knownpos #%d at %d < pos3 %d\n",*old_knowni,*old_knownpos,pos3));
      knowni[nfound] = *old_knowni++;
      sites[nfound++] = *old_knownpos++;
    }

    return nfound;

  } else {
    /* Advance old_knownpos past pos5 */
    while (*old_knownpos < pos5) {
      debug2(printf("Skipping old_knownpos #%d at %d < pos5 %d\n",*old_knowni,*old_knownpos,pos5));
      old_knowni++;
      old_knownpos++;
    }

    /* Startblock */
    if (genomealt != genome) {
      found = block_find_snp(&high_halfsite,&low_halfsite,
			     &(genome->blocks[startblocki]),&(genomealt->blocks[startblocki]),splicesite_bits);
    } else {
      found = block_find(&high_halfsite,&low_halfsite,&(genome->blocks[startblocki]),splicesite_bits);
    }
    found = clear_start(found,startdiscard);
    debug2(printf("adding start mask %08x\n",clear_start_mask(startdiscard)));

    while (found != 0U) {
#ifdef HAVE_BUILTIN_CTZ
      pos = offset + (relpos = __builtin_ctz(found));
      while (*old_knownpos < pos) {
	debug2(printf("Adding old_knownpos #%d at %d < pos %d\n",*old_knowni,*old_knownpos,pos));
	knowni[nfound] = *old_knowni++;
	sites[nfound++] = *old_knownpos++;
      }
      if (*old_knownpos == pos) {
	debug2(printf("Adding old_knownpos #%d at %d == pos %d\n",*old_knowni,*old_knownpos,pos));
	knowni[nfound] = *old_knowni++;
	sites[nfound++] = *old_knownpos++;
      } else {
	knowni[nfound] = -1;
	sites[nfound++] = pos;
      }
      clear_lowbit(found,relpos);
#else
      debug2(printf("found is %08X, -found & found is %08X\n",found,-found & found));
      pos = offset + mod_37_bit_position[(lowbit = -found & found) % 37];
      while (*old_knownpos < pos) {
	debug2(printf("Adding old_knownpos #%d at %d < pos %d\n",*old_knowni,*old_knownpos,pos));
	knowni[nfound] = *old_knowni++;
	sites[nfound++] = *old_knownpos++;
      }
      if (*old_knownpos == pos) {
	debug2(printf("Adding old_knownpos #%d at %d == pos %d\n",*old_knowni,*old_knownpos,pos));
	knowni[nfound] = *old_knowni++;
	sites[nfound++] = *old_knownpos++;
      } else {
	knowni[nfound] = -1;
	sites[nfound++] = pos;
      }
      found -= lowbit;
#endif
      debug2(printf("found is %08X => offset %llu + relpos %d\n",found,(unsigned long long) offset,relpos));
    }

    ptr = &(genome->blocks[startblocki+3]);
    altptr = &(genomealt->blocks[startblocki+3]);
    end = &(genome->blocks[endblocki]);
    offset += 32;
    while (ptr < end) {
      prev_high_halfsite = high_halfsite;
      if (genomealt != genome) {
	found = block_find_snp(&high_halfsite,&low_halfsite,ptr,altptr,splicesite_bits);
      } else {
	found = block_find(&high_halfsite,&low_halfsite,ptr,splicesite_bits);
      }

      if (low_halfsite & prev_high_halfsite) {
	debug2(printf("low_halfsite & prev_high_halfsite => offset %llu - 1\n",(unsigned long long) offset));
	pos = offset - 1;	/* verified that this should be offset - 1 */
	while (*old_knownpos < pos) {
	debug2(printf("Adding old_knownpos #%d at %d < pos %d\n",*old_knowni,*old_knownpos,pos));
	  knowni[nfound] = *old_knowni++;
	  sites[nfound++] = *old_knownpos++;
	}
	if (*old_knownpos == pos) {
	  debug2(printf("Adding old_knownpos #%d at %d == pos %d\n",*old_knowni,*old_knownpos,pos));
	  knowni[nfound] = *old_knowni++;
	  sites[nfound++] = *old_knownpos++;
	} else {
	  knowni[nfound] = -1;
	  sites[nfound++] = pos;
	}
      }

      while (found != 0U) {
#ifdef HAVE_BUILTIN_CTZ
	pos = offset + (relpos = __builtin_ctz(found));
	while (*old_knownpos < pos) {
	  debug2(printf("Adding old_knownpos #%d at %d < pos %d\n",*old_knowni,*old_knownpos,pos));
	  knowni[nfound] = *old_knowni++;
	  sites[nfound++] = *old_knownpos++;
	}
	if (*old_knownpos == pos) {
	  debug2(printf("Adding old_knownpos #%d at %d == pos %d\n",*old_knowni,*old_knownpos,pos));
	  knowni[nfound] = *old_knowni++;
	  sites[nfound++] = *old_knownpos++;
	} else {
	  knowni[nfound] = -1;
	  sites[nfound++] = pos;
	}
	clear_lowbit(found,relpos);
#else
	debug2(printf("found is %08X, -found & found is %08X\n",found,-found & found));
	pos = offset + mod_37_bit_position[(lowbit = -found & found) % 37];
	while (*old_knownpos < pos) {
	  debug2(printf("Adding old_knownpos #%d at %d < pos %d\n",*old_knowni,*old_knownpos,pos));
	  knowni[nfound] = *old_knowni++;
	  sites[nfound++] = *old_knownpos++;
	}
	if (*old_knownpos == pos) {
	  debug2(printf("Adding old_knownpos #%d at %d == pos %d\n",*old_knowni,*old_knownpos,pos));
	  knowni[nfound] = *old_knowni++;
	  sites[nfound++] = *old_knownpos++;
	} else {
	  knowni[nfound] = -1;
	  sites[nfound++] = pos;
	}
	found -= lowbit;
#endif
	debug2(printf("found is %08X => offset %llu + relpos %d\n",found,(unsigned long long) offset,relpos));
      }

      ptr += 3;
      altptr += 3;
      offset += 32;
    }

    /* Endblock */
    prev_high_halfsite = high_halfsite;
    if (genomealt != genome) {
      found = block_find_snp(&high_halfsite,&low_halfsite,ptr,altptr,splicesite_bits);
    } else {
      found = block_find(&high_halfsite,&low_halfsite,ptr,splicesite_bits);
    }
    found = clear_end(found,enddiscard);
    debug2(printf("adding end mask %08x\n",clear_end_mask(enddiscard)));

    if (low_halfsite & prev_high_halfsite) {
      debug2(printf("low_halfsite & prev_high_halfsite => offset %llu - 1\n",(unsigned long long) offset));
      pos = offset - 1;		/* verified that this should be offset - 1 */
      while (*old_knownpos < pos) {
	debug2(printf("Adding old_knownpos #%d at %d < pos %d\n",*old_knowni,*old_knownpos,pos));
	knowni[nfound] = *old_knowni++;
	sites[nfound++] = *old_knownpos++;
      }
      if (*old_knownpos == pos) {
	debug2(printf("Adding old_knownpos #%d at %d == pos %d\n",*old_knowni,*old_knownpos,pos));
	knowni[nfound] = *old_knowni++;
	sites[nfound++] = *old_knownpos++;
      } else {
	knowni[nfound] = -1;
	sites[nfound++] = pos;
      }
    }

    while (found != 0U) {
#ifdef HAVE_BUILTIN_CTZ
      pos = offset + (relpos = __builtin_ctz(found));
      while (*old_knownpos < pos) {
	debug2(printf("Adding old_knownpos #%d at %d < pos %d\n",*old_knowni,*old_knownpos,pos));
	knowni[nfound] = *old_knowni++;
	sites[nfound++] = *old_knownpos++;
      }
      if (*old_knownpos == pos) {
	debug2(printf("Adding old_knownpos #%d at %d == pos %d\n",*old_knowni,*old_knownpos,pos));
	knowni[nfound] = *old_knowni++;
	sites[nfound++] = *old_knownpos++;
      } else {
	knowni[nfound] = -1;
	sites[nfound++] = pos;
      }
      clear_lowbit(found,relpos);
#else
      debug2(printf("found is %08X, -found & found is %08X\n",found,-found & found));
      pos = offset + mod_37_bit_position[(lowbit = -found & found) % 37];
      while (*old_knownpos < pos) {
	debug2(printf("Adding old_knownpos #%d at %d < pos %d\n",*old_knowni,*old_knownpos,pos));
	knowni[nfound] = *old_knowni++;
	sites[nfound++] = *old_knownpos++;
      }
      if (*old_knownpos == pos) {
	debug2(printf("Adding old_knownpos #%d at %d == pos %d\n",*old_knowni,*old_knownpos,pos));
	knowni[nfound] = *old_knowni++;
	sites[nfound++] = *old_knownpos++;
      } else {
	knowni[nfound] = -1;
	sites[nfound++] = pos;
      }
      found -= lowbit;
#endif
      debug2(printf("found is %08X => offset %llu + relpos %d\n",found,(unsigned long long) offset,relpos));
    }

    /* Add old_knownpos to pos3 */
    while (*old_knownpos < pos3) {
      debug2(printf("Adding old_knownpos #%d at %d < pos3 %d\n",*old_knowni,*old_knownpos,pos3));
      knowni[nfound] = *old_knowni++;
      sites[nfound++] = *old_knownpos++;
    }

    return nfound;
  }
}


int
Genome_donor_sites (int *sites, int *knowni, int *old_knownpos, int *old_knowni,
		    Univcoord_T left, int pos5, int pos3) {
  return splicesites(sites,knowni,old_knownpos,old_knowni,left,pos5,pos3,
		     donor_gtgc_bits,/*splicepos_offset*/0);
}

int
Genome_acceptor_sites (int *sites, int *knowni, int *old_knownpos, int *old_knowni,
		       Univcoord_T left, int pos5, int pos3) {
  return splicesites(sites,knowni,old_knownpos,old_knowni,left,
		     pos5 - /*splicepos_offset*/2,pos3 - /*splicepos_offset*/2,
		     acceptor_bits,/*splicepos_offset*/2);
}

int
Genome_antidonor_sites (int *sites, int *knowni, int *old_knownpos, int *old_knowni,
			Univcoord_T left, int pos5, int pos3) {
  return splicesites(sites,knowni,old_knownpos,old_knowni,left,
		     pos5 - /*splicepos_offset*/2,pos3 - /*splicepos_offset*/2,
		     antidonor_acgc_bits,/*splicepos_offset*/2);
}

int
Genome_antiacceptor_sites (int *sites, int *knowni, int *old_knownpos, int *old_knowni,
			   Univcoord_T left, int pos5, int pos3) {
  return splicesites(sites,knowni,old_knownpos,old_knowni,left,pos5,pos3,
		     antiacceptor_bits,/*splicepos_offset*/0);
}


#if 0
static int
splicesites_novel (int *sites, Univcoord_T left, int pos5, int pos3,
		   const Genomecomp_T *splicesite_bits, int splicepos_offset) {
  int nfound = 0, offset;
  int startdiscard, enddiscard;
  Univcoord_T startblocki, endblocki;
  Genomecomp_T *ptr, *altptr, *end;
  Genomecomp_T found;
  Genomecomp_T high_halfsite, low_halfsite, prev_high_halfsite;
  int pos;
#ifdef HAVE_BUILTIN_CTZ
  int relpos;
#else
  Genomecomp_T lowbit;
#endif

  debug2(
	printf("\n\n");
	printf("Genome (in splicesites):\n");
	Genome_print_blocks(genome->blocks,left+pos5,left+pos3);
	printf("\n");
	);

  /* nshift = left % 32; */

  startblocki = (left+pos5)/32U*3;
  endblocki = (left+pos3)/32U*3;
  startdiscard = (left+pos5) % 32;
  enddiscard = (left+pos3) % 32;

  offset = -startdiscard + pos5 + splicepos_offset;
  
  debug2(printf("left = %llu, pos5 = %d, pos3 = %d, startblocki = %llu, endblocki = %llu\n",
		(unsigned long long) left,pos5,pos3,(unsigned long long) startblocki,(unsigned long long) endblocki));
  debug2(printf("startdiscard = %d, enddiscard = %d\n",startdiscard,enddiscard));

  if (endblocki == startblocki) {
    if (genomealt != genome) {
      found = block_find_snp(&high_halfsite,&low_halfsite,
			     &(genome->blocks[startblocki]),&(genomealt->blocks[startblocki]),splicesite_bits);
    } else {
      found = block_find(&high_halfsite,&low_halfsite,&(genome->blocks[startblocki]),splicesite_bits);
    }
    found = clear_start(found,startdiscard);
    found = clear_end(found,enddiscard);
    debug2(printf("adding masks %08X and %08x\n",clear_start_mask(startdiscard),clear_end_mask(enddiscard)));

    while (found != 0U) {
#ifdef HAVE_BUILTIN_CTZ
      pos = offset + (relpos = __builtin_ctz(found));
      sites[nfound++] = pos;
      clear_lowbit(found,relpos);
#else      
      debug2(printf("found is %08X, -found & found is %08X\n",found,-found & found));
      pos = offset + mod_37_bit_position[(lowbit = -found & found) % 37];
      sites[nfound++] = pos;
      found -= lowbit;
#endif
      debug2(printf("found is %08X => offset %llu + relpos %d\n",found,(unsigned long long) offset,relpos));
    }

    return nfound;

  } else {
    /* Startblock */
    if (genomealt != genome) {
      found = block_find_snp(&high_halfsite,&low_halfsite,
			     &(genome->blocks[startblocki]),&(genomealt->blocks[startblocki]),splicesite_bits);
    } else {
      found = block_find(&high_halfsite,&low_halfsite,&(genome->blocks[startblocki]),splicesite_bits);
    }
    found = clear_start(found,startdiscard);
    debug2(printf("adding start mask %08x\n",clear_start_mask(startdiscard)));

    while (found != 0U) {
#ifdef HAVE_BUILTIN_CTZ
      pos = offset + (relpos = __builtin_ctz(found));
      sites[nfound++] = pos;
      clear_lowbit(found,relpos);
#else
      debug2(printf("found is %08X, -found & found is %08X\n",found,-found & found));
      pos = offset + mod_37_bit_position[(lowbit = -found & found) % 37];
      sites[nfound++] = pos;
      found -= lowbit;
#endif
      debug2(printf("found is %08X => offset %llu + relpos %d\n",found,(unsigned long long) offset,relpos));
    }

    ptr = &(genome->blocks[startblocki+3]);
    altptr = &(genomealt->blocks[startblocki+3]);
    end = &(genome->blocks[endblocki]);
    offset += 32;
    while (ptr < end) {
      prev_high_halfsite = high_halfsite;
      if (genomealt != genome) {
	found = block_find_snp(&high_halfsite,&low_halfsite,ptr,altptr,splicesite_bits);
      } else {
	found = block_find(&high_halfsite,&low_halfsite,ptr,splicesite_bits);
      }

      if (low_halfsite & prev_high_halfsite) {
	debug2(printf("low_halfsite & prev_high_halfsite => offset %llu - 1\n",(unsigned long long) offset));
	pos = offset - 1;	/* verified that this should be offset - 1 */
	sites[nfound++] = pos;
      }

      while (found != 0U) {
#ifdef HAVE_BUILTIN_CTZ
	pos = offset + (relpos = __builtin_ctz(found));
	sites[nfound++] = pos;
	clear_lowbit(found,relpos);
#else
	debug2(printf("found is %08X, -found & found is %08X\n",found,-found & found));
	pos = offset + mod_37_bit_position[(lowbit = -found & found) % 37];
	sites[nfound++] = pos;
	found -= lowbit;
#endif
	debug2(printf("found is %08X => offset %llu + relpos %d\n",found,(unsigned long long) offset,relpos));
      }

      ptr += 3;
      altptr += 3;
      offset += 32;
    }

    /* Endblock */
    prev_high_halfsite = high_halfsite;
    if (genomealt != genome) {
      found = block_find_snp(&high_halfsite,&low_halfsite,ptr,altptr,splicesite_bits);
    } else {
      found = block_find(&high_halfsite,&low_halfsite,ptr,splicesite_bits);
    }
    found = clear_end(found,enddiscard);
    debug2(printf("adding end mask %08x\n",clear_end_mask(enddiscard)));

    if (low_halfsite & prev_high_halfsite) {
      debug2(printf("low_halfsite & prev_high_halfsite => offset %llu - 1\n",(unsigned long long) offset));
      pos = offset - 1;		/* verified that this should be offset - 1 */
      sites[nfound++] = pos;
    }

    while (found != 0U) {
#ifdef HAVE_BUILTIN_CTZ
      pos = offset + (relpos = __builtin_ctz(found));
      sites[nfound++] = pos;
      clear_lowbit(found,relpos);
#else
      debug2(printf("found is %08X, -found & found is %08X\n",found,-found & found));
      pos = offset + mod_37_bit_position[(lowbit = -found & found) % 37];
      sites[nfound++] = pos;
      found -= lowbit;
#endif
      debug2(printf("found is %08X => offset %llu + relpos %d\n",found,(unsigned long long) offset,relpos));
    }

    return nfound;
  }
}
#endif


#if 0
int
Genome_donor_sites_novel (int *sites, Univcoord_T left, int pos5, int pos3) {
  return splicesites_novel(sites,left,pos5,pos3,
			   donor_gtgc_bits,/*splicepos_offset*/0);
}

int
Genome_acceptor_sites_novel (int *sites, Univcoord_T left, int pos5, int pos3) {
  return splicesites_novel(sites,left,pos5 - /*splicepos_offset*/2,
			   pos3 - /*splicepos_offset*/2,acceptor_bits,/*splicepos_offset*/2);
}

int
Genome_antidonor_sites_novel (int *sites, Univcoord_T left, int pos5, int pos3) {
  return splicesites_novel(sites,left,pos5 - /*splicepos_offset*/2,
			   pos3 - /*splicepos_offset*/2,antidonor_acgc_bits,/*splicepos_offset*/2);
}

int
Genome_antiacceptor_sites_novel (int *sites, Univcoord_T left, int pos5, int pos3) {
  return splicesites_novel(sites,left,pos5,pos3,
			   antiacceptor_bits,/*splicepos_offset*/0);
}
#endif

