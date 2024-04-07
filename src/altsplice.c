static char rcsid[] = "$Id: f84fc65597fa261a19fa85670d98cedf7e11fc3c $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "altsplice.h"
#include "genomicpos.h"

#include <stdio.h>
#include <string.h>		/* For memcpy */
#include "assert.h"
#include "sedgesort.h"

static int max_insertlength;


/* Altsplice creation */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif


/* Altsplice_trim */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif


/* Altsplice_resolve procedures */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif



#define T Altsplice_T


void
Altsplice_free (T *old, Pathpool_T pathpool) {
#ifdef DEBUG0
  static int call_i = 0;
#endif

  if (*old) {
    debug0(printf("%d: Altsplice_free of %p\n",++call_i,*old));

    /* Allocated by Vectorpool_T, and not reclaiming due to variable lengths */
    /* FREE_OUT((*old)->coords); */
    /* FREE_OUT((*old)->distal_trimpos); */
    /* FREE_OUT((*old)->nmismatches); */
    /* FREE_OUT((*old)->ref_nmismatches); */
    /* FREE_OUT((*old)->distal_probs); */
    Pathpool_free_altsplice(&(*old),pathpool
			    pathpool_trace(__FILE__,__LINE__)); /* Allocated by Pathpool_new_altsplice */
  }
  return;
}


void
Altsplice_print (T this, Univcoord_T medial_coord) {
  int i;
  Chrpos_T splice_distance;

  if (this != NULL) {
    printf("%p",this);

    for (i = 0; i < this->ncoords; i++) {
      if (this->coords[i] > medial_coord) {
	splice_distance = this->coords[i] - medial_coord;
      } else {
	splice_distance = medial_coord - this->coords[i];
      }

#ifdef LARGE_GENOMES
      /* printf(" %u,%llu(%f-%f)->%d",splice_distance,this->coords[i],this->medial_probs[i],this->distal_probs[i],
	 this->distal_trimpos[i]); */
      printf(" %u,%llu(%d,%d)",splice_distance,this->coords[i],this->medial_nmismatches[i],this->distal_nmismatches[i]);
#else
      /* printf(" %u,%u(%f-%f)->%d",splice_distance,this->coords[i],this->medial_probs[i],this->distal_probs[i],
	 this->distal_trimpos[i]); */
      printf(" %u,%u(%d,%d)",splice_distance,this->coords[i],this->medial_nmismatches[i],this->distal_nmismatches[i]);
#endif
    }
  }
  return;
}


T
Altsplice_copy (T old, Pathpool_T pathpool, Vectorpool_T vectorpool) {
  T new;
  int npartners;

#ifdef DEBUG0
  static int call_i = 0;
#endif

  if (old == NULL) {
    return (T) NULL;
  } else {
    new = (T) Pathpool_new_altsplice(pathpool
				     pathpool_trace(__FILE__,__LINE__));

    debug0(printf("%d: Altsplice_copy of %p => %p\n",++call_i,old,new));

    new->boundedp = old->boundedp;
    new->anchor_qpos = old->anchor_qpos;

    new->best_distal_nmatches = old->best_distal_nmatches;
    new->best_distal_length = old->best_distal_length;
    new->best_medial_prob = old->best_medial_prob;
    new->best_distal_prob = old->best_distal_prob;
    new->best_coord = old->best_coord;

    npartners = new->ncoords = old->ncoords;

    new->coords = Vectorpool_new_univcoordvector(vectorpool,npartners);
    new->distal_lengths = Vectorpool_new_intvector(vectorpool,npartners);
    new->distal_trimpos = Vectorpool_new_intvector(vectorpool,npartners);
    new->medial_nmismatches = Vectorpool_new_intvector(vectorpool,npartners);
    new->distal_nmismatches = Vectorpool_new_intvector(vectorpool,npartners);
    new->medial_probs = Vectorpool_new_doublevector(vectorpool,npartners);
    new->distal_probs = Vectorpool_new_doublevector(vectorpool,npartners);

    memcpy(new->coords,old->coords,npartners*sizeof(Univcoord_T));
    memcpy(new->distal_lengths,old->distal_lengths,npartners*sizeof(int));
    memcpy(new->distal_trimpos,old->distal_trimpos,npartners*sizeof(int));
    memcpy(new->medial_nmismatches,old->medial_nmismatches,npartners*sizeof(int));
    memcpy(new->distal_nmismatches,old->distal_nmismatches,npartners*sizeof(int));
    memcpy(new->medial_probs,old->medial_probs,npartners*sizeof(double));
    memcpy(new->distal_probs,old->distal_probs,npartners*sizeof(double));

    return new;
  }
}


static bool
check_ascending (Univcoord_T *coords, int n) {
  Univcoord_T prev_coord;
  int i;

  prev_coord = coords[0];
  for (i = 1; i < n; i++) {
    if (coords[i] < prev_coord) {
      return false;
    }
    prev_coord = coords[i];
  }
 
  return true;
}

static bool
check_descending (Univcoord_T *coords, int n) {
  Univcoord_T prev_coord;
  int i;

  prev_coord = coords[0];
  for (i = 1; i < n; i++) {
    if (coords[i] > prev_coord) {
      return false;
    }
    prev_coord = coords[i];
  }
 
  return true;
}


T
Altsplice_qstart_new (int *best_distal_length, bool boundedp, int anchor_qpos,
		      Univcoord_T *distal_positions,
		      int *distal_lengths, int *distal_trimpos,
		      int *medial_nmismatches, int *distal_nmismatches,
		      double *medial_probs, double *distal_probs, int npartners,
		      Pathpool_T pathpool, Vectorpool_T vectorpool,
		      bool sort_bydistal_p) {
  T new = Pathpool_new_altsplice(pathpool
				 pathpool_trace(__FILE__,__LINE__));
  int *order, i, oldi, k;

#ifdef DEBUG0
  static int call_i = 0;
#endif
  
  debug0(printf("%d: Creating Altsplice_qstart_new %p\n",++call_i,new));

  new->boundedp = boundedp;
  new->anchor_qpos = anchor_qpos;

  /* Find reasonable values for Path_eval_nmatches based on (1) best
     medial prob and (2) best distal prob for that */

  new->best_distal_length = distal_lengths[0];
  new->best_distal_nmatches = distal_lengths[0] - distal_nmismatches[0];
  new->best_medial_prob = medial_probs[0];
  new->best_distal_prob = distal_probs[0];
  new->best_coord = distal_positions[0];

  for (k = 1; k < npartners; k++) {
    if (medial_probs[k] > new->best_medial_prob) {
      new->best_distal_length = distal_lengths[k];
      new->best_distal_nmatches = distal_lengths[k] - distal_nmismatches[k];
      new->best_medial_prob = medial_probs[k];
      new->best_distal_prob = distal_probs[k];
      new->best_coord = distal_positions[k];

    } else if (medial_probs[k] == new->best_medial_prob &&
	       distal_probs[k] > new->best_distal_prob) {
      new->best_distal_length = distal_lengths[k];
      new->best_distal_nmatches = distal_lengths[k] - distal_nmismatches[k];
      /* new->best_medial_prob = medial_probs[k]; */
      new->best_distal_prob = distal_probs[k];
      new->best_coord = distal_positions[k];
    }
  }

  new->ncoords = npartners;

  new->coords = Vectorpool_new_univcoordvector(vectorpool,npartners);
  new->distal_lengths = Vectorpool_new_intvector(vectorpool,npartners);
  new->distal_trimpos = Vectorpool_new_intvector(vectorpool,npartners);
  new->medial_nmismatches = Vectorpool_new_intvector(vectorpool,npartners);
  new->distal_nmismatches = Vectorpool_new_intvector(vectorpool,npartners);
  new->medial_probs = Vectorpool_new_doublevector(vectorpool,npartners);
  new->distal_probs = Vectorpool_new_doublevector(vectorpool,npartners);

  if (sort_bydistal_p == false || npartners == 1 ||
      check_descending(distal_positions,npartners) == true) {
    memcpy(new->coords,distal_positions,npartners*sizeof(Univcoord_T));
    memcpy(new->distal_lengths,distal_lengths,npartners*sizeof(int));
    memcpy(new->distal_trimpos,distal_trimpos,npartners*sizeof(int));
    memcpy(new->medial_nmismatches,medial_nmismatches,npartners*sizeof(int));
    memcpy(new->distal_nmismatches,distal_nmismatches,npartners*sizeof(int));
    memcpy(new->medial_probs,medial_probs,npartners*sizeof(double));
    memcpy(new->distal_probs,distal_probs,npartners*sizeof(double));

  } else {
#ifdef LARGE_GENOMES
    order = Sedgesort_order_uint8(distal_positions,npartners);
#else
    order = Sedgesort_order_uint4(distal_positions,npartners);
#endif
    
    k = 0;
    for (i = npartners - 1; i >= 0; i--) {
      oldi = order[i];

      new->coords[k] = distal_positions[oldi];
      new->distal_lengths[k] = distal_lengths[oldi];
      new->distal_trimpos[k] = distal_trimpos[oldi];
      new->medial_nmismatches[k] = medial_nmismatches[oldi];
      new->distal_nmismatches[k] = distal_nmismatches[oldi];
      new->medial_probs[k] = medial_probs[oldi];
      new->distal_probs[k] = distal_probs[oldi];
      k++;
    }

    FREE(order);
  }

  *best_distal_length = new->best_distal_length;

  return new;
}


T
Altsplice_qend_new (int *best_distal_length, bool boundedp, int anchor_qpos,
		    Univcoord_T *distal_positions,
		    int *distal_lengths, int *distal_trimpos,
		    int *medial_nmismatches, int *distal_nmismatches,
		    double *medial_probs, double *distal_probs, int npartners,
		    Pathpool_T pathpool, Vectorpool_T vectorpool,
		    bool sort_bydistal_p) {
  T new = Pathpool_new_altsplice(pathpool
				 pathpool_trace(__FILE__,__LINE__));
  int *order, i, oldi, k;

#ifdef DEBUG0
  static int call_i = 0;
#endif
  
  debug0(printf("%d: Creating Altsplice_qend_new %p\n",++call_i,new));


  new->boundedp = boundedp;
  new->anchor_qpos = anchor_qpos;

  /* Find reasonable values for Path_eval_nmatches based on (1) best
     medial prob and (2) best distal prob for that */

  new->best_distal_length = distal_lengths[0];
  new->best_distal_nmatches = distal_lengths[0] - distal_nmismatches[0];
  new->best_medial_prob = medial_probs[0];
  new->best_distal_prob = distal_probs[0];
  new->best_coord = distal_positions[0];

  for (k = 1; k < npartners; k++) {
    if (medial_probs[k] > new->best_medial_prob) {
      new->best_distal_length = distal_lengths[k];
      new->best_distal_nmatches = distal_lengths[k] - distal_nmismatches[k];
      new->best_medial_prob = medial_probs[k];
      new->best_distal_prob = distal_probs[k];
      new->best_coord = distal_positions[k];

    } else if (medial_probs[k] == new->best_medial_prob &&
	       distal_probs[k] > new->best_distal_prob) {
      new->best_distal_length = distal_lengths[k];
      new->best_distal_nmatches = distal_lengths[k] - distal_nmismatches[k];
      /* new->best_medial_prob = medial_probs[k]; */
      new->best_distal_prob = distal_probs[k];
      new->best_coord = distal_positions[k];
    }
  }

  new->ncoords = npartners;

  new->coords = Vectorpool_new_univcoordvector(vectorpool,npartners);
  new->distal_lengths = Vectorpool_new_intvector(vectorpool,npartners);
  new->distal_trimpos = Vectorpool_new_intvector(vectorpool,npartners);
  new->medial_nmismatches = Vectorpool_new_intvector(vectorpool,npartners);
  new->distal_nmismatches = Vectorpool_new_intvector(vectorpool,npartners);
  new->medial_probs = Vectorpool_new_doublevector(vectorpool,npartners);
  new->distal_probs = Vectorpool_new_doublevector(vectorpool,npartners);

  if (sort_bydistal_p == false || npartners == 1 ||
      check_ascending(distal_positions,npartners) == true) {
    memcpy(new->coords,distal_positions,npartners*sizeof(Univcoord_T));
    memcpy(new->distal_lengths,distal_lengths,npartners*sizeof(int));
    memcpy(new->distal_trimpos,distal_trimpos,npartners*sizeof(int));
    memcpy(new->medial_nmismatches,medial_nmismatches,npartners*sizeof(int));
    memcpy(new->distal_nmismatches,distal_nmismatches,npartners*sizeof(int));
    memcpy(new->medial_probs,medial_probs,npartners*sizeof(double));
    memcpy(new->distal_probs,distal_probs,npartners*sizeof(double));

  } else {
#ifdef LARGE_GENOMES
    order = Sedgesort_order_uint8(distal_positions,npartners);
#else
    order = Sedgesort_order_uint4(distal_positions,npartners);
#endif
    
    k = 0;
    for (i = 0; i < npartners; i++) {
      oldi = order[i];

      new->coords[k] = distal_positions[oldi];
      new->distal_lengths[k] = distal_lengths[oldi];
      new->distal_trimpos[k] = distal_trimpos[oldi];
      new->medial_nmismatches[k] = medial_nmismatches[oldi];
      new->distal_nmismatches[k] = distal_nmismatches[oldi];
      new->medial_probs[k] = medial_probs[oldi];
      new->distal_probs[k] = distal_probs[oldi];
      k++;
    }

    FREE(order);
  }

  *best_distal_length = new->best_distal_length;

  return new;
}


Univcoord_T
Altsplice_firstcoord (T this) {
  return this->coords[0];
}

Univcoord_T
Altsplice_lastcoord (T this) {
  return this->coords[this->ncoords - 1];
}


/* Returns whether anything is left */
bool
Altsplice_trim_qstart_chrbounds (T this, Univcoord_T chroffset) {
  int endi;
  
  /* Traverse from genome start and go upstream */
  endi = this->ncoords - 1;
  debug3(printf("Entered Altsplice_trim_qstart_chrbounds from %d to 0\n",endi));
  while (endi >= 0 && this->coords[endi] < chroffset + this->distal_lengths[endi]) {
    debug3(printf("Distal coord %u < chroffset %u, so trimming\n",
		  this->coords[endi] - this->distal_lengths[endi],chroffset));
    endi--;
  }

  if (endi < 0) {
    return false;
  } else {
    this->ncoords = endi + 1;
    debug3(printf("Resetting ncoords to be %d\n",this->ncoords));
    return true;
  }
}


/* Returns whether anything is left */
bool
Altsplice_trim_qend_chrbounds (T this, Univcoord_T chrhigh) {
  int endj;
  
  /* Traverse from genome end and go downstream */
  endj = this->ncoords - 1;
  debug3(printf("Entered Altsplice_trim_qstart_chrbounds from %d to 0\n",endj));
  while (endj >= 0 && this->coords[endj] + this->distal_lengths[endj] >= chrhigh) {
    debug3(printf("Distal coord %u >= chrhigh %u, so trimming\n",
		  this->coords[endj] + this->distal_lengths[endj],chrhigh));
    endj--;
  }

  if (endj < 0) {
    return false;
  } else {
    this->ncoords = endj + 1;
    debug3(printf("Resetting ncoords to be %d\n",this->ncoords));
    return true;
  }
}





/* Should look for insertlengths within max_insertlength (continuous
   exon), and choose among those with the best nmismatches and first
   one with a sufficiently high prob.  Could look for sufficiently
   high prob.  Otherwise (intervening exon), insertlength doesn't
   matter, and we should look for the best nmismatches and prob */

/* Resolves the qend of pathL */
/* Diagonals are in order from closest to farthest from pathL, so need to check array from n-1 to 0 */
bool
Altsplice_resolve_qend (Univcoord_T *univdiagonal_L, int *splice_qpos_L, int *distal_trimpos_L,
			int *medial_nmismatches_L, int *distal_nmismatches_L, 
			double *medial_prob_L, double *distal_prob_L,
			T this, int anchor_qpos_L, int querylengthL, int querylengthH,
			Univcoord_T genomicstartH) {

  /* Candidates within max_insertlength are in starti..endi */
  int besti = -1, starti, endi, i;
  int ncandidates;
  int insertlength;
  Univcoord_T genomicendL;
  int best_nmismatches_continuous, best_nmismatches_intervening, nmismatches;
  double best_prob, prob;


  debug9(printf("Altsplice_resolve_qend\n"));

  check_ascending(this->coords,this->ncoords);

  best_nmismatches_continuous = querylengthL;
  best_nmismatches_intervening = querylengthL;

  endi = this->ncoords - 1;
  while (endi >= 0 &&
	 /*genomicendL*/(this->coords[endi] + this->distal_lengths[endi]) >= 
	 genomicstartH /*allow overlap*/+ querylengthH) {

    debug9(genomicendL = this->coords[endi] + this->distal_lengths[endi]);
    debug9(insertlength = genomicstartH - genomicendL + querylengthL + querylengthH);
    debug9(printf(" Overreach: (insertlength %d, nmismatches %d+%d, prob %f+%f)\n",
		  insertlength,this->medial_nmismatches[endi],this->distal_nmismatches[endi],
		  this->medial_probs[endi],this->distal_probs[endi]));
    endi--;
  }


  i = starti = endi;
  while (i >= 0) {
    genomicendL = this->coords[i] + this->distal_lengths[endi];
    insertlength = genomicstartH - genomicendL + querylengthL + querylengthH;
    if (insertlength <= max_insertlength) {
      debug9(printf(" Continuous: (insertlength %d, nmismatches %d+%d, prob %f+%f)\n",
		    insertlength,this->medial_nmismatches[endi],this->distal_nmismatches[endi],
		    this->medial_probs[endi],this->distal_probs[endi]));
      if ((nmismatches = this->medial_nmismatches[i] + this->distal_nmismatches[i]) < best_nmismatches_continuous) {
	best_nmismatches_continuous = nmismatches;
      }
      starti = i;
    } else {
      debug9(printf(" Intervening: (insertlength %d, nmismatches %d+%d, prob %f+%f)\n",
		    insertlength,this->medial_nmismatches[endi],this->distal_nmismatches[endi],
		    this->medial_probs[endi],this->distal_probs[endi]));
      if ((nmismatches = this->medial_nmismatches[i] + this->distal_nmismatches[i]) < best_nmismatches_intervening) {
	best_nmismatches_intervening = nmismatches;
      }
    }
    i--;
  }
  
  debug9(printf("best_nmismatches_continuous %d\n",best_nmismatches_continuous));
  debug9(printf("best_nmismatches_intervening %d\n",best_nmismatches_intervening));

  if (best_nmismatches_continuous < querylengthL) {
    /* Find shortest (first) good continuous solution */
    ncandidates = 0;
    for (i = endi; i >= starti; i--) {
      if (this->medial_nmismatches[i] + this->distal_nmismatches[i] == best_nmismatches_continuous &&
	  (prob = this->medial_probs[i] + this->distal_probs[i]) > 0.9 + 0.9) {
	if (ncandidates++ == 0) {
	  besti = i;
	}
      }
    }
  }

  if (besti < 0 && best_nmismatches_intervening < querylengthL) {
    /* Find best intervening solution */
    best_prob = 0.9;		/* threshold */
    for (i = endi - 1; i >= 0; i--) {
      if (this->medial_nmismatches[i] + this->distal_nmismatches[i] == best_nmismatches_intervening &&
	  (prob = this->medial_probs[i] + this->distal_probs[i]) > best_prob) {
	besti = i;
	best_prob = prob;
      }
    }
  }

  debug9(printf("besti %d\n",besti));

  if (besti < 0) {
    return false;

  } else {
    /* Single intervening exon or same exon */
    *univdiagonal_L = this->coords[besti] + this->distal_lengths[besti];
    *splice_qpos_L = querylengthL - this->distal_lengths[besti];
    *distal_trimpos_L = this->distal_trimpos[besti];
    if (anchor_qpos_L != this->anchor_qpos) {
      /* Anchor qpos was changed, probably by combine_leftright_paths */
      *medial_nmismatches_L = -1;
    } else {
      *medial_nmismatches_L = this->medial_nmismatches[besti];
    }
    *distal_nmismatches_L = this->distal_nmismatches[besti];
    *medial_prob_L = this->medial_probs[besti];
    *distal_prob_L = this->distal_probs[besti];
    return true;
  }
}


/* Diagonals are in order from closest to farthest from pathH, so need to check array from n-1 to 0 */
bool
Altsplice_resolve_qstart (Univcoord_T *univdiagonal_H, int *splice_qpos_H, int *distal_trimpos_H,
			  int *medial_nmismatches_H, int *distal_nmismatches_H,
			  double *medial_prob_H, double *distal_prob_H,
			  T this, int anchor_qpos_H, int querylengthL, int querylengthH,
			  Univcoord_T genomicendL) {

  /* Candidates within max_insertlength are in startj..endj */
  int bestj = -1, startj, endj, j;
  int ncandidates;
  int insertlength;
  Univcoord_T genomicstartH;
  int best_nmismatches_continuous, best_nmismatches_intervening, nmismatches;
  double best_prob, prob;


  debug9(printf("Altsplice_resolve_qstart\n"));

  check_descending(this->coords,this->ncoords);

  best_nmismatches_continuous = querylengthH;
  best_nmismatches_intervening = querylengthH;

  endj = this->ncoords - 1;
  while (endj >= 0 &&
	 /*genomicstartH*/(this->coords[endj] - this->distal_lengths[endj]) /*allow_overlap*/+ querylengthH <= genomicendL) {

    debug9(genomicstartH = this->coords[endj] - this->distal_lengths[endj]);
    debug9(insertlength = genomicstartH - genomicendL + querylengthL + querylengthH);
    debug9(printf(" Overreach: (insertlength %d, nmismatches %d+%d, prob %f+%f)\n",
		  insertlength,this->medial_nmismatches[endj],this->distal_nmismatches[endj],
		  this->medial_probs[endj],this->distal_probs[endj]));
    endj--;
  }

  
  j = startj = endj;
  while (j >= 0) {
    genomicstartH = this->coords[j] - this->distal_lengths[j];
    insertlength = genomicstartH - genomicendL + querylengthL + querylengthH;
    if (insertlength <= max_insertlength) {
      debug9(printf(" Continuous: (insertlength %d, nmismatches %d+%d, prob %f+%f)\n",
		    insertlength,this->medial_nmismatches[endj],this->distal_nmismatches[endj],
		    this->medial_probs[endj],this->distal_probs[endj]));
      if ((nmismatches = this->medial_nmismatches[j] + this->distal_nmismatches[j]) < best_nmismatches_continuous) {
	best_nmismatches_continuous = nmismatches;
      }
      startj = j;
    } else {
      debug9(printf(" Intervening: (insertlength %d, nmismatches %d+%d, prob %f+%f)\n",
		    insertlength,this->medial_nmismatches[endj],this->distal_nmismatches[endj],
		    this->medial_probs[endj],this->distal_probs[endj]));
      if ((nmismatches = this->medial_nmismatches[j] + this->distal_nmismatches[j]) < best_nmismatches_intervening) {
	best_nmismatches_intervening = nmismatches;
      }
    }
    j--;
  }

  debug9(printf("best_nmismatches_continuous %d\n",best_nmismatches_continuous));
  debug9(printf("best_nmismatches_intervening %d\n",best_nmismatches_intervening));

  if (best_nmismatches_continuous < querylengthH) {
    /* Find shortest (first) good continuous solution */
    ncandidates = 0;
    for (j = endj; j >= startj; j--) {
      if (this->medial_nmismatches[j] + this->distal_nmismatches[j] == best_nmismatches_continuous &&
	  (prob = this->medial_probs[j] + this->distal_probs[j]) > 0.9 + 0.9) {
	if (ncandidates++ == 0) {
	  bestj = j;
	}
      }
    }
  }

  if (bestj < 0 && best_nmismatches_intervening < querylengthH) {
    /* Find best intervening solution */
    best_prob = 0.9 + 0.9;		/* threshold */
    for (j = endj - 1; j >= 0; j--) {
      if (this->medial_nmismatches[j] + this->distal_nmismatches[j] == best_nmismatches_intervening &&
	  (prob = this->medial_probs[j] + this->distal_probs[j]) > best_prob) {
	bestj = j;
	best_prob = prob;
      }
    }
  }
	  
  debug9(printf("bestj %d\n",bestj));

  if (bestj < 0) {
    return false;

  } else {
    /* Single intervening exon or same exon */
    *univdiagonal_H = this->coords[bestj] - this->distal_lengths[bestj] + querylengthH;
    *splice_qpos_H = this->distal_lengths[bestj];
    *distal_trimpos_H = this->distal_trimpos[bestj];
    if (anchor_qpos_H != this->anchor_qpos) {
      /* Anchor qpos was changed, probably by combine_leftright_paths */
      *medial_nmismatches_H = -1;
    } else {
      *medial_nmismatches_H = this->medial_nmismatches[bestj];
    }
    *distal_nmismatches_H = this->distal_nmismatches[bestj];
    *medial_prob_H = this->medial_probs[bestj];
    *distal_prob_H = this->distal_probs[bestj];
    return true;
  }
}


bool
Altsplice_resolve_both (Univcoord_T *univdiagonal_L, int *splice_qpos_L, int *distal_trimpos_L,
			int *medial_nmismatches_L, int *distal_nmismatches_L,
			double *medial_prob_L, double *distal_prob_L,

			Univcoord_T *univdiagonal_H, int *splice_qpos_H, int *distal_trimpos_H,
			int *medial_nmismatches_H, int *distal_nmismatches_H,
			double *medial_prob_H, double *distal_prob_H,

			T thisL, int anchor_qpos_L, T thisH, int anchor_qpos_H,
			int querylengthL, int querylengthH) {
  int besti = -1, bestj = -1, besti_intervening = -1, bestj_intervening = -1, i, j;
  int ncandidates;
  int best_insertlength, insertlength;
  Univcoord_T genomicendL, genomicstartH;
  int best_nmismatches_continuous, best_nmismatches_intervening, nmismatches;
  double best_prob, prob;


  debug9(printf("Altsplice_resolve_both\n"));

  best_nmismatches_continuous = best_nmismatches_intervening = querylengthL + querylengthH;

  for (i = 0; i < thisL->ncoords; i++) {
    genomicendL = thisL->coords[i] + thisL->distal_lengths[i];
    for (j = 0; j < thisH->ncoords; j++) {
      genomicstartH = thisH->coords[j] - thisH->distal_lengths[j];
      if (genomicendL < genomicstartH /*allow overlap*/+ querylengthH) {
	/* Look for valid insertlength */
	insertlength = genomicstartH - genomicendL + querylengthL + querylengthH;
	debug9(printf(" (insertlength %u, nmismatches %d+%d + %d+%d, prob %f+%f + %f+%f)\n",
		      insertlength,thisL->medial_nmismatches[i],thisL->distal_nmismatches[i],
		      thisH->medial_nmismatches[j],thisH->distal_nmismatches[j],
		      thisL->medial_probs[i],thisL->distal_probs[i],
		      thisH->medial_probs[j],thisH->distal_probs[j]));
	if (insertlength <= max_insertlength) {
	  if ((nmismatches = thisL->medial_nmismatches[i] + thisL->distal_nmismatches[i] +
	       thisH->medial_nmismatches[j] + thisH->distal_nmismatches[j]) < best_nmismatches_continuous) {
	    best_nmismatches_continuous = nmismatches;
	  }
	} else {
	  if ((nmismatches = thisL->medial_nmismatches[i] + thisL->distal_nmismatches[i] +
	       thisH->medial_nmismatches[j] + thisH->distal_nmismatches[j]) < best_nmismatches_intervening) {
	    best_nmismatches_intervening = nmismatches;
	  }
	}
      }
    }
  }

  debug9(printf("best_nmismatches_continuous %d\n",best_nmismatches_continuous));
  debug9(printf("best_nmismatches_intervening %d\n",best_nmismatches_intervening));

  if (best_nmismatches_continuous < querylengthL + querylengthH ||
      best_nmismatches_intervening < querylengthL + querylengthH) {
    /* Find shortest good continuous solution and best intervening solution */
    ncandidates = 0;
    best_prob = 0.9 + 0.9 + 0.9 + 0.9;	/* threshold */

    for (i = 0; i < thisL->ncoords; i++) {
      genomicendL = thisL->coords[i] + thisL->distal_lengths[i];
      for (j = 0; j < thisH->ncoords; j++) {
	genomicstartH = thisH->coords[j] - thisH->distal_lengths[j];
	if (genomicendL < genomicstartH /*allow overlap*/+ querylengthH) {
	  /* Look for valid insertlength */
	  insertlength = genomicstartH - genomicendL + querylengthL + querylengthH;
	  if (insertlength <= max_insertlength) {
	    if (thisL->medial_nmismatches[i] + thisL->distal_nmismatches[i] +
		thisH->medial_nmismatches[j] + thisH->distal_nmismatches[j] == best_nmismatches_continuous &&
		thisL->distal_probs[i] > 0.9 && thisH->distal_probs[j] > 0.9) {
	      if (ncandidates++ == 0 || insertlength < best_insertlength) {
		besti = i;
		bestj = j;
		best_insertlength = insertlength;
	      }
	    }

	  } else {
	    if (thisL->medial_nmismatches[i] + thisL->distal_nmismatches[i] +
		thisH->medial_nmismatches[j] + thisH->distal_nmismatches[j] == best_nmismatches_intervening &&
		(prob = thisL->medial_probs[i] + thisL->distal_probs[i] +
		 thisH->medial_probs[j] + thisH->distal_probs[j]) > best_prob) {
	      besti_intervening = i;
	      bestj_intervening = j;
	      best_prob = prob;
	    }
	  }
	}
      }
    }
  }


  if (besti < 0) {
    besti = besti_intervening;
    bestj = bestj_intervening;
  }

  debug9(printf("besti and bestj %d,%d\n",besti,bestj));

  if (besti < 0) {
    return false;

  } else {
    *univdiagonal_L = thisL->coords[besti] + thisL->distal_lengths[besti];
    *splice_qpos_L = querylengthL - thisL->distal_lengths[besti];
    *distal_trimpos_L = thisL->distal_trimpos[besti];
    if (anchor_qpos_L != thisL->anchor_qpos) {
      /* Anchor qpos was changed, probably by combine_leftright_paths */
      *medial_nmismatches_L = -1;
    } else {
      *medial_nmismatches_L = thisL->medial_nmismatches[besti];
    }
    *distal_nmismatches_L = thisL->distal_nmismatches[besti];
    *medial_prob_L = thisL->medial_probs[besti];
    *distal_prob_L = thisL->distal_probs[besti];

    *univdiagonal_H = thisH->coords[bestj] - thisH->distal_lengths[bestj] + querylengthH;
    *splice_qpos_H = thisH->distal_lengths[bestj];
    *distal_trimpos_H = thisH->distal_trimpos[bestj];
    if (anchor_qpos_H != thisH->anchor_qpos) {
      /* Anchor qpos was changed, probably by combine_leftright_paths */
      *medial_nmismatches_H = -1;
    } else {
      *medial_nmismatches_H = thisH->medial_nmismatches[bestj];
    }
    *distal_nmismatches_H = thisH->distal_nmismatches[bestj];
    *medial_prob_H = thisH->medial_probs[bestj];
    *distal_prob_H = thisH->distal_probs[bestj];

    return true;
  }
}
  

void
Altsplice_setup (int max_insertlength_in) {

  max_insertlength = max_insertlength_in;

  return;
}
