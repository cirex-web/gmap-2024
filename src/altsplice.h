/* $Id: d8d39fe7fb2786a748e91e693f6e91d8784df8ba $ */
#ifndef ALTSPLICE_INCLUDED
#define ALTSPLICE_INCLUDED

typedef struct Altsplice_T *Altsplice_T;

#include "bool.h"
#include "univcoord.h"
#include "pathpool.h"
#include "vectorpool.h"


#define T Altsplice_T
struct T {

  bool boundedp;

  /* Because of combine_leftright_paths, the medial_nmismatches might not be valid for the combined path */
  int anchor_qpos;		/* The qpos more medial than splice_qpos; beginning of the segment */

  /* Needed by Path_eval_nmatches */
  int best_distal_length;
  int best_distal_nmatches;
  double best_medial_prob;
  double best_distal_prob;

  int ncoords;
  Univcoord_T *coords;
  Univcoord_T best_coord;	/* Needed by Path_genomiclow and Path_genomichigh */

#if 0
  /* For a single medial splice_qpos */
  int distal_length;
  int *distal_trimpos;
  int *nmismatches;
  double medial_prob;
  double *distal_probs;

#else
  /* For handling multiple medial splice_qpos */
  int *distal_lengths;
  int *distal_trimpos;
  int *medial_nmismatches;
  int *distal_nmismatches;
  double *medial_probs;
  double *distal_probs;
#endif

  /* int *ref_nmismatches; */
};


static inline bool
Altsplice_boundedp (T this) {
  return this->boundedp;
}

extern void
Altsplice_free (T *old, Pathpool_T pathpool);

extern void
Altsplice_print (T this, Univcoord_T medial_coord);

extern T
Altsplice_copy (T old, Pathpool_T pathpool, Vectorpool_T vectorpool);

extern T
Altsplice_qstart_new (int *best_distal_length, bool boundedp, int anchor_qpos,
		      Univcoord_T *distal_positions,
		      int *distal_lengths, int *distal_trimpos,
		      int *medial_nmismatches, int *distal_nmismatches,
		      double *medial_probs, double *distal_probs, int npartners,
		      Pathpool_T pathpool, Vectorpool_T vectorpool,
		      bool sort_bydistal_p);

extern T
Altsplice_qend_new (int *best_distal_length, bool boundedp, int anchor_qpos,
		    Univcoord_T *distal_positions,
		    int *distal_lengths, int *distal_trimpos,
		    int *medial_nmismatches, int *distal_nmismatches,
		    double *medial_probs, double *distal_probs, int npartners,
		    Pathpool_T pathpool, Vectorpool_T vectorpool,
		    bool sort_bydistal_p);

extern Univcoord_T
Altsplice_firstcoord (T this);

extern Univcoord_T
Altsplice_lastcoord (T this);

extern bool
Altsplice_trim_qstart_chrbounds (T this, Univcoord_T chroffset);

extern bool
Altsplice_trim_qend_chrbounds (T this, Univcoord_T chrhigh);

extern bool
Altsplice_resolve_qend (Univcoord_T *univdiagonal_L, int *splice_qpos_L, int *distal_trimpos_L,
			int *medial_nmismatches_L, int *distal_nmismatches_L, 
			double *medial_prob_L, double *distal_prob_L,
			T this, int anchor_qpos_L, int querylengthL, int querylengthH,
			Univcoord_T genomicstartH);
extern bool
Altsplice_resolve_qstart (Univcoord_T *univdiagonal_H, int *splice_qpos_H, int *distal_trimpos_H,
			  int *medial_nmismatches_L, int *distal_nmismatches_H,
			  double *medial_prob_H, double *distal_prob_H,
			  T this, int anchor_qpos_H, int querylengthL, int querylengthH,
			  Univcoord_T genomicendL);
extern bool
Altsplice_resolve_both (Univcoord_T *univdiagonal_L, int *splice_qpos_L, int *distal_trimpos_L,
			int *medial_nmismatches_L, int *distal_nmismatches_L,
			double *medial_prob_L, double *distal_prob_L,

			Univcoord_T *univdiagonal_H, int *splice_qpos_H, int *distal_trimpos_H,
			int *medial_nmismatches_H, int *distal_nmismatches_H,
			double *medial_prob_H, double *distal_prob_H,

			T thisL, int anchor_qpos_L, T thisH, int anchor_qpos_H,
			int querylengthL, int querylengthH);

extern void
Altsplice_setup (int max_insertlength_in);

#undef T
#endif


