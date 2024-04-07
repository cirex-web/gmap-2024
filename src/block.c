static char rcsid[] = "$Id: block.c 223349 2020-10-28 02:49:25Z twu $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "block.h"

#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#ifdef PMAP
#include "oligop.h"
#else
#include "oligo.h"
#endif
#include "reader.h"		/* also has cDNAEnd_T */

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* In reading inward from both ends, the block5 and block3 share the
   same reader.  But in reading outward from the middle, the block5
   and block3 each have their own reader. */

#define T Block_T
struct T {
  Reader_T reader;

  cDNAEnd_T cdnaend;

  int querylength;
  int last_querypos;
  Oligostate_T last_state;
  Width_T oligosize;

#ifdef PMAP
  Oligospace_T aaindex;
#else
  int leftreadshift;
  Oligospace_T forward;
  Oligospace_T revcomp;
  Oligospace_T oligomask;
#endif

  /* For saving state */
  int last_querypos_save;
  Oligostate_T last_state_save;
#ifdef PMAP
  Oligospace_T aaindex_save;
#else
  Oligospace_T forward_save;
  Oligospace_T revcomp_save;
#endif
};


Reader_T
Block_reader (T this) {
  return this->reader;
}

int
Block_querypos (T this) {
  return this->last_querypos;
}

#ifdef PMAP
Oligospace_T
Block_aaindex (T this) {
 return this->aaindex;
}
#else
Oligospace_T
Block_forward (T this) {
  if (this->cdnaend == FIVE) {
    return this->forward & this->oligomask;
  } else {
    return (this->forward >> this->leftreadshift) & this->oligomask;
  }
}

Oligospace_T
Block_revcomp (T this) {
  if (this->cdnaend == FIVE) {
    return (this->revcomp >> this->leftreadshift) & this->oligomask;
  } else {
    return this->revcomp & this->oligomask;
  }
}
#endif

bool
Block_donep (T this) {
  if (this->last_state == DONE) {
    return true;
  } else {
    return false;
  }
}


extern void
Block_save (T this) {
  /* save->reader = this->reader; -- not necessary */
  /* save->cdnaend = this->cdnaend; -- not necessary */

  this->last_querypos_save = this->last_querypos;
  this->last_state_save = this->last_state;
#ifdef PMAP
  this->aaindex_save = this->aaindex;
  debug(printf("Saving block at last_querypos %d, aaindex %u\n",this->last_querypos,this->aaindex));
#else
  this->forward_save = this->forward;
  this->revcomp_save = this->revcomp;
  debug(printf("Saving block at last_querypos %d, forward %u, revcomp %u\n",
	       this->last_querypos,this->forward,this->revcomp));
#endif

  return;
}

extern void
Block_restore (T this) {
  /* this->reader = save->reader; -- not necessary */
  /* this->cdnaend = save->cdnaend; -- not necessary */

  if (this->cdnaend == FIVE) {
    Reader_reset_start(this->reader,this->last_querypos_save+this->oligosize);
  } else {
    Reader_reset_end(this->reader,this->last_querypos_save-1);
  }

  this->last_querypos = this->last_querypos_save;
  this->last_state = this->last_state_save;
#ifdef PMAP
  this->aaindex = this->aaindex_save;
#else
  this->forward = this->forward_save;
  this->revcomp = this->revcomp_save;
#endif
  return;
}

extern void
Block_reset_ends (T this) {
  Reader_reset_ends(this->reader);

  if (this->cdnaend == FIVE) {
    this->last_querypos = Reader_querystart(this->reader);
  } else {
    this->last_querypos = Reader_queryend(this->reader);
  }
  this->last_state = INIT;

#ifdef PMAP  
  this->aaindex = 0U;
#else
  this->forward = 0U;
  this->revcomp = 0U;
#endif

  return;
}


T
Block_new (cDNAEnd_T cdnaend, Width_T oligosize,
#ifndef PMAP
	   int leftreadshift,
#endif
	   Reader_T reader, int querylength) {
  T new = (T) MALLOC(sizeof(*new));

  new->reader = reader;
  new->cdnaend = cdnaend;
  new->oligosize = oligosize;
#ifndef PMAP
#ifdef HAVE_64_BIT
  new->oligomask = ~(~0ULL << 2*oligosize);
#else
  new->oligomask = ~(~0U << 2*oligosize);
#endif
  new->leftreadshift = leftreadshift;
#endif

  new->querylength = querylength;
  if (cdnaend == FIVE) {
#ifdef PMAP
    new->last_querypos = Reader_startpos(reader);
#else
    new->last_querypos = Reader_startpos(reader) - oligosize;
#endif
  } else if (cdnaend == THREE) {
    new->last_querypos = Reader_endpos(reader) + 1;
  }
  new->last_state = INIT;

#ifdef PMAP  
  new->aaindex = 0U;
#else
  new->forward = 0U;
  new->revcomp = 0U;
#endif

  new->last_querypos_save = new->last_querypos;
  new->last_state_save = new->last_state;
#ifdef PMAP
  new->aaindex_save = new->aaindex;
#else
  new->forward_save = new->forward;
  new->revcomp_save = new->revcomp;
#endif

  return new;
}

void
Block_free (T *old) {
  if (*old) {
    FREE(*old);
  }
  return;
}


bool
Block_next_5 (T this, int genestrand) {
#ifdef DEBUG
  char *nt_fwd, *nt_rev;
#endif

  if (this->last_state == DONE) {
    return false;
  } else {

#ifdef PMAP
    this->last_state = Oligo_next_5(this->last_state,&this->last_querypos,
				    &this->aaindex,this->reader,genestrand);
    debug(printf("Block has aaindex %u at querypos %d\n",
		 this->aaindex,this->last_querypos));
#else
    this->last_state = Oligo_next_5(this->last_state,&this->last_querypos,
				    &this->forward,&this->revcomp,this->reader,genestrand);
    debug(
	  if (this->cdnaend == THREE) {
	    nt_fwd = Oligo_one_nt(this->forward >> this->leftreadshift,this->oligosize);
	    nt_rev = Oligo_one_nt(this->revcomp,this->oligosize);
	  }
	  if (this->cdnaend == FIVE) {
	    nt_fwd = Oligo_one_nt(this->forward,this->oligosize);
	    nt_rev = Oligo_one_nt(this->revcomp >> this->leftreadshift,this->oligosize);
	  }
	  printf("Block has oligo forward %s, revcomp %s at querypos %d\n",
		 nt_fwd,nt_rev,this->last_querypos);
	  FREE(nt_rev);
	  FREE(nt_fwd)
	  );
#endif

    if (this->last_state == DONE) {
      return false;
    } else {
      return true;
    }
  }
}


bool
Block_next_3 (T this, int genestrand) {
#ifdef DEBUG
  char *nt_fwd, *nt_rev;
#endif

  if (this->last_state == DONE) {
    return false;
  } else {

#ifdef PMAP
    this->last_state = Oligo_next_3(this->last_state,&this->last_querypos,
				    &this->aaindex,this->reader,genestrand);
    debug(printf("Block has aaindex %u at querypos %d\n",
		 this->aaindex,this->last_querypos));
#else
    this->last_state = Oligo_next_3(this->last_state,&this->last_querypos,
				    &this->forward,&this->revcomp,this->reader,genestrand);
    debug(
	  if (this->cdnaend == THREE) {
	    nt_fwd = Oligo_one_nt(this->forward >> this->leftreadshift,this->oligosize);
	    nt_rev = Oligo_one_nt(this->revcomp,this->oligosize);
	  }
	  if (this->cdnaend == FIVE) {
	    nt_fwd = Oligo_one_nt(this->forward,this->oligosize);
	    nt_rev = Oligo_one_nt(this->revcomp >> this->leftreadshift,this->oligosize);
	  }
	  printf("Block has oligo forward %s, revcomp %s at querypos %d\n",
		 nt_fwd,nt_rev,this->last_querypos);
	  FREE(nt_rev);
	  FREE(nt_fwd)
	  );
#endif

    if (this->last_state == DONE) {
      return false;
    } else {
      return true;
    }
  }
}


/* Returns whether skipping was successful (as opposed to Block_next, which returns whether scanning can continue) */
bool
Block_skip_5 (T this, int genestrand, int nskip) {
  int init_querypos;

  if (this->last_state == DONE) {
    return false;
  } else {

    init_querypos = this->last_querypos;
#ifdef PMAP
    this->last_state = Oligo_skip_5(this->last_state,&this->last_querypos,
				    &this->aaindex,this->reader,genestrand,nskip);
    debug(printf("Block has aaindex %u at querypos %d\n",
		 this->aaindex,this->last_querypos));
#else
    this->last_state = Oligo_skip_5(this->last_state,&this->last_querypos,
				    &this->forward,&this->revcomp,
				    this->reader,genestrand,nskip);
    debug(printf("Block has oligo %08X at querypos %d\n",
		 this->forward,this->last_querypos));
#endif

    if (this->last_state != VALID) {
      return false;
    } else if (this->last_querypos - init_querypos != nskip) {
      return false;
    } else {
      return true;
    }
  }
}


bool
Block_skip_3 (T this, int genestrand, int nskip) {
  int init_querypos;

  if (this->last_state == DONE) {
    return false;
  } else {

    init_querypos = this->last_querypos;
#ifdef PMAP
    this->last_state = Oligo_skip_3(this->last_state,&this->last_querypos,
				    &this->aaindex,this->reader,genestrand,nskip);
    debug(printf("Block has aaindex %u at querypos %d\n",
		 this->aaindex,this->last_querypos));
#else
    this->last_state = Oligo_skip_3(this->last_state,&this->last_querypos,
				    &this->forward,&this->revcomp,
				    this->reader,genestrand,nskip);
    debug(printf("Block has oligo %08X at querypos %d\n",
		 this->forward,this->last_querypos));
#endif

    if (this->last_state != VALID) {
      return false;
    } else if (this->last_querypos - init_querypos != nskip) {
      return false;
    } else {
      return true;
    }
  }
}


#if 0
/* Returns whether skipping was successful (as opposed to Block_next, which returns whether scanning can continue) */
bool
Block_skipto (T this, int querypos) {
  int nskip;

  if (this->cdnaend == FIVE) {
    nskip = querypos - this->last_querypos;
  } else {
    nskip = this->last_querypos - querypos;
  }
  if (nskip < 0) {
    return false;
  } else {
    return Block_skip(this,nskip);
  }
}
#endif


/* Returns querypos */
#ifdef PMAP

int
Block_process_oligo (Univcoord_T **fwdpositions, int *nfwdhits, 
		     Univcoord_T **revpositions, int *nrevhits,
		     T this, Indexdb_T indexdb_fwd, Indexdb_T indexdb_rev,
		     int indexdb_sizelimit) {

  /* Note that querylength was already multiplied by 3 */
  *fwdpositions = Indexdb_read_with_diagterm(&(*nfwdhits),indexdb_fwd,this->aaindex,
					     /*diagterm*/this->querylength-3*this->last_querypos);
  *revpositions = Indexdb_read_with_diagterm(&(*nrevhits),indexdb_rev,this->aaindex,
					     /*diagterm*/3*this->last_querypos);
  return this->last_querypos;
}

#else

/* In standard mode, indexdb_rev is indexdb_fwd.  They differ for cmet and atoi modes. */
int
Block_process_oligo_5 (Univcoord_T **fwdpositions, int *nfwdhits,
		       Univcoord_T **revpositions, int *nrevhits,
		       T this, Indexdb_T indexdb_fwd, Indexdb_T indexdb_rev) {
#if 0
  printf("block_process: Querypos %d, oligos are %06X and %06X\n",this->last_querypos,this->forward,this->revcomp >> this->leftreadshift);
#endif
  *fwdpositions = Indexdb_read_with_diagterm(&(*nfwdhits),indexdb_fwd,this->forward & this->oligomask,
					     /*diagterm*/this->querylength-this->last_querypos);
  *revpositions = Indexdb_read_with_diagterm(&(*nrevhits),indexdb_rev,(this->revcomp >> this->leftreadshift) & this->oligomask,
					     /*diagterm*/this->last_querypos);

  return this->last_querypos;
}

/* In standard mode, indexdb_rev is indexdb_fwd.  They differ for cmet and atoi modes. */
int
Block_process_oligo_3 (Univcoord_T **fwdpositions, int *nfwdhits,
		       Univcoord_T **revpositions, int *nrevhits,
		       T this, Indexdb_T indexdb_fwd, Indexdb_T indexdb_rev) {
#if 0
  printf("block_process Querypos %d, oligos are %06X and %06X\n",this->last_querypos,this->forward >> this->leftreadshift,this->revcomp);
#endif
  *fwdpositions = Indexdb_read_with_diagterm(&(*nfwdhits),indexdb_fwd,(this->forward >> this->leftreadshift) & this->oligomask,
					     /*diagterm*/this->querylength-this->last_querypos);
  *revpositions = Indexdb_read_with_diagterm(&(*nrevhits),indexdb_rev,this->revcomp & this->oligomask,
					     /*diagterm*/this->last_querypos);

  return this->last_querypos;
}


#endif
