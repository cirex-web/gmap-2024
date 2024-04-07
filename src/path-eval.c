static char rcsid[] = "$Id: ea0c59dfa9b67bf1967c8df15f8be26907c47aa0 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "path-eval.h"
#include "path-solve.h"
#include "path-trim.h"
#include "transcript-remap.h"

#include <stdio.h>
#include <string.h>		/* For strcpy */
#include <math.h>		/* For rint */
#include <ctype.h>		/* For islower */
#include "fastlog.h"		/* For fasterexp */

#include "assert.h"
#include "list.h"
#include "genomebits_count.h"
#include "junction.h"
#include "mapq.h"

static Genomebits_T genomebits;
static Genomebits_T genomebits_alt;
static Transcriptome_T transcriptome;

static bool *circularp;
static bool *chrsubsetp;
static bool *altlocp;

static int index1part;
static int index1interval;

static Outputtype_T output_type;
static bool md_report_snps_p;
static bool want_random_p;
static bool allow_soft_clips_p;


#ifdef CHECK_ASSERTIONS
#define CHECK_NMISMATCHES 1
#endif

#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Path_eval_nmatches */
#ifdef DEBUG7
#define debug7(x) x
#else
#define debug7(x)
#endif

/* Path_eval_and_sort */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* Path_consolidate */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif


#define T Path_T


/* Previously, some procedures removed a read if path->nmatches <
   nmismatches_allowed.  However, since this procedure sets
   found_score in such cases, found_score is no longer accurate.  One
   solution would be for caller to revise found_score.  Another
   solution is not to compare path->nmatches against
   nmismatches_allowed.  Also, could consider using found_score to
   constrain subsequent procedures. */

/* Sets found_score, score_within_trims, nmatches, ref_nmatches, junction_splice_prob, total_splice_prob */
int
Path_eval_nmatches (int *found_score, T this, Compress_T query_compress_fwd, Compress_T query_compress_rev) {
  int qstart, qend, ninserts;
  Univcoord_T univdiagonal;
  Intlist_T r, x, y;
  Univcoordlist_T q;
  Junction_T junction;
  List_T j;
  /* bool insertionp = false; */
  /* int adj0; deletions - insertions */
  int total_ninserts = 0, nmismatches, ref_nmismatches;


  debug7(printf("\nEntering Path_eval_nmatches on path %p\n",this));
  debug7(Path_print(this));

  Path_expect_fwd(this);

  this->found_score = 0;
  this->score_within_trims = 0;
  this->nmatches = 0;
  this->ref_nmatches = 0;
  this->junction_splice_prob = 0.0;
  this->total_splice_prob = 0.0;


  assert(Univcoordlist_length(this->univdiagonals) == Intlist_length(this->endpoints) - 1);
  assert(Intlist_length(this->nmismatches) == Intlist_length(this->endpoints) - 1);
  assert(Intlist_length(this->ref_nmismatches) == Intlist_length(this->endpoints) - 1);
  assert(List_length(this->junctions) == Intlist_length(this->endpoints) - 2);


  qstart = Intlist_head(this->endpoints);
  nmismatches = Intlist_head(this->nmismatches);
  ref_nmismatches = Intlist_head(this->ref_nmismatches);
  ninserts = 0;
  
  if (this->plusp == true) {
    /* plus */
    j = this->junctions;		/* Put here before we handle querystart_alts */
    if (this->ambig_prob_5 > 0.0) {
      /* Count ambig splice site as matches, but counts against found_score */
      this->total_splice_prob += this->ambig_prob_5;
      /* this->nmatches += qstart; (inferred) */
      /* this->ref_nmatches += qstart; (inferred) */
      this->found_score += qstart;
      /* Does not affect score_within_trims */

    } else if (this->qstart_alts != NULL) {
      this->nmatches += this->qstart_alts->best_distal_nmatches; /* Not nmatches_to_trims, which is 0 for alts_substring */
      debug7(printf("Adding %d matches for qstart_alts => total %d\n",this->qstart_alts->best_distal_nmatches,this->nmatches));
      this->ref_nmatches += this->qstart_alts->best_distal_nmatches; /* Not nmatches_to_trims, which is 0 for alts_substring */
      this->junction_splice_prob += this->qstart_alts->best_medial_prob + this->qstart_alts->best_distal_prob;
      this->total_splice_prob += this->qstart_alts->best_medial_prob + this->qstart_alts->best_distal_prob;
      this->found_score += this->qstart_alts->best_distal_length - this->qstart_alts->best_distal_nmatches;
      this->score_within_trims += this->qstart_alts->best_distal_length - this->qstart_alts->best_distal_nmatches;

    } else if (this->fusion_querystart_junction != NULL) {
      /* Rest of fusions handled below */

    } else {
      this->found_score += qstart;
      /* Does not affect score_within_trims */
    }

    /* Add qpos to get alignstart/alignend */
    for (q = this->univdiagonals, x = this->nmismatches, y = this->ref_nmismatches, r = Intlist_next(this->endpoints); q != NULL;
	 q = Univcoordlist_next(q), x = Intlist_next(x), y = Intlist_next(y), r = Intlist_next(r), j = List_next(j)) {
      qstart += ninserts;
      qend = Intlist_head(r);
#if 0
      if (insertionp == true) {
	nmismatches = ref_nmismatches = -1; /* Recompute nmismatches */
      } else {
	nmismatches = Intlist_head(x);
	ref_nmismatches = Intlist_head(y);
      }
#else
      nmismatches = Intlist_head(x);
      ref_nmismatches = Intlist_head(y);
#endif

      univdiagonal = Univcoordlist_head(q);
      /* left = univdiagonal - (Univcoord_T) this->querylength; */
      debug7(printf("Path_eval_nmatches: ninserts %d, qstart %d..qend %d at univdiagonal %u [%u]\n",
		    ninserts,qstart,qend,univdiagonal,univdiagonal - this->chroffset));

      if (nmismatches >= 0 && ref_nmismatches >= 0) {
	debug7(printf("Checking mismatches at %u from querystart %d to queryend %d\n",univdiagonal - this->chroffset,qstart,qend));
	debug7(printf("%d mismatches expected vs %d measured\n",
		      nmismatches,
		      Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_fwd,
							    univdiagonal,this->querylength,
							    /*pos5*/qstart,/*pos3*/qend,/*plusp*/true,this->genestrand)));
#ifdef CHECK_NMISMATCHES
	assert(nmismatches == Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_fwd,
								    univdiagonal,this->querylength,
								    /*pos5*/qstart,/*pos3*/qend,/*plusp*/true,this->genestrand));
#endif
      } else {
	nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_fwd,
							    univdiagonal,this->querylength,
							    /*pos5*/qstart,/*pos3*/qend,/*plusp*/true,this->genestrand);
	Intlist_head_set(x,nmismatches);		/* Save for Stage3end_new_substrings */
	Intlist_head_set(y,ref_nmismatches);		/* Save for Stage3end_new_substrings */
	debug7(printf("%d (%d ref) mismatches from genome over querypos %d..%d\n",
		      nmismatches,ref_nmismatches,qstart,qend));
      }

      /* Could potentially check here if qstart < qend, but relying upon caller to use endpoints_acceptable_p */
      this->nmatches += (qend - qstart) - nmismatches;
      debug7(printf("(1) Shortcut adds matches of %d = (%d - %d) - nmismatches %d => total %d\n",
		    (qend-qstart)-nmismatches,qend,qstart,nmismatches,this->nmatches));
      this->ref_nmatches += (qend - qstart) - ref_nmismatches;
      this->found_score += nmismatches;
      this->score_within_trims += nmismatches;

      /* Prepare for next iteration */
      qstart = qend;

      if (j == NULL) {
	ninserts = 0;
      } else if ((junction = (Junction_T) List_head(j)) == NULL) {
	/* qstart_junction */
	ninserts = 0;
      } else {
	debug7(printf("Junction: ")); debug7(Junction_print(junction)); debug7(printf("\n"));
	ninserts = Junction_ninserts(junction);
	this->junction_splice_prob += Junction_prob(junction);
	this->total_splice_prob += Junction_prob(junction);
      }
    }

  } else {
    /* minus */

    j = this->junctions;		/* Put here before we handle querystart_alts */
    if (this->ambig_prob_5 > 0.0) {
      /* Count ambig splice site as matches, but counts against found_score */
      this->total_splice_prob += this->ambig_prob_5;
      /* this->nmatches += qstart; (inferred) */
      /* this->ref_nmatches += qstart; (inferred) */
      this->found_score += qstart;
      /* Does not affect score_within_trims */

    } else if (this->qstart_alts != NULL) {
      this->nmatches += this->qstart_alts->best_distal_nmatches; /* Not nmatches_to_trims, which is 0 for alts_substring */
      debug7(printf("Adding %d matches for qstart_alts => total %d\n",this->qstart_alts->best_distal_nmatches,this->nmatches));
      this->ref_nmatches += this->qstart_alts->best_distal_nmatches; /* Not nmatches_to_trims, which is 0 for alts_substring */
      this->junction_splice_prob += this->qstart_alts->best_medial_prob + this->qstart_alts->best_distal_prob;
      this->total_splice_prob += this->qstart_alts->best_medial_prob + this->qstart_alts->best_distal_prob;
      this->found_score += this->qstart_alts->best_distal_length - this->qstart_alts->best_distal_nmatches;
      this->score_within_trims += this->qstart_alts->best_distal_length - this->qstart_alts->best_distal_nmatches;

    } else if (this->fusion_querystart_junction != NULL) {
      /* Fusions handled below */

    } else {
      this->found_score += qstart;
      /* Does not affect score_within_trims */
    }

    /* Subtract qpos to get alignstart/alignend */
    for (q = this->univdiagonals, x = this->nmismatches, y = this->ref_nmismatches, r = Intlist_next(this->endpoints); q != NULL;
	 q = Univcoordlist_next(q), x = Intlist_next(x), y = Intlist_next(y), r = Intlist_next(r), j = List_next(j)) {
      qstart += ninserts;
      qend = Intlist_head(r);
#if 0
      if (insertionp == true) {
	nmismatches = ref_nmismatches = -1; /* Recompute nmismatches */
      } else {
	nmismatches = Intlist_head(x);
	ref_nmismatches = Intlist_head(y);
      }
#else
      nmismatches = Intlist_head(x);
      ref_nmismatches = Intlist_head(y);
#endif
      univdiagonal = Univcoordlist_head(q);
      /* left = univdiagonal - (Univcoord_T) this->querylength; */
      debug7(printf("Path_eval_nmatches: ninserts %d, qstart %d..qend %d at univdiagonal %u [%u]\n",
		    ninserts,qstart,qend,univdiagonal,univdiagonal - this->chroffset));

      if (nmismatches >= 0 && ref_nmismatches >= 0) {
	debug7(printf("Checking mismatches at %u from querystart %d to queryend %d\n",univdiagonal - this->chroffset,qstart,qend));
	debug7(printf("%d mismatches expected vs %d measured\n",
		      nmismatches,Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_rev,
									univdiagonal,this->querylength,
									/*pos5*/qstart,/*pos3*/qend,/*plusp*/false,this->genestrand)));
#ifdef CHECK_NMISMATCHES
	assert(nmismatches == Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_rev,
								    univdiagonal,this->querylength,
								    /*pos5*/qstart,/*pos3*/qend,/*plusp*/false,this->genestrand));
#endif
      } else {
	nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_rev,
							    univdiagonal,this->querylength,
							    /*pos5*/qstart,/*pos3*/qend,/*plusp*/false,this->genestrand);
	Intlist_head_set(x,nmismatches);		/* Save for Stage3end_new_substrings */
	Intlist_head_set(y,ref_nmismatches);		/* Save for Stage3end_new_substrings */
	debug7(printf("%d (%d ref) mismatches from genome over querypos %d..%d\n",
		      nmismatches,ref_nmismatches,this->querylength - qend,this->querylength - qstart));
      }

      /* Could potentially check here if qstart < qend, but relying upon caller to use endpoints_acceptable_p */
      this->nmatches += (qend - qstart) - nmismatches;
      debug7(printf("(2) Shortcut adds matches of %d = (%d - %d) - nmismatches %d => total %d\n",
		    (qend-qstart)-nmismatches,this->querylength - qstart,this->querylength - qend,nmismatches,this->nmatches));
      this->ref_nmatches += (qend - qstart) - ref_nmismatches;
      this->found_score += nmismatches;
      this->score_within_trims += nmismatches;

      /* Prepare for next iteration */
      qstart = qend;

      if (j == NULL) {
	ninserts = 0;
      } else if ((junction = (Junction_T) List_head(j)) == NULL) {
	/* qstart_junction */
	ninserts = 0;
      } else {
	debug7(printf("Junction: ")); debug7(Junction_print(junction)); debug7(printf("\n"));
	ninserts = Junction_ninserts(junction);
	this->junction_splice_prob += Junction_prob(junction);
	this->total_splice_prob += Junction_prob(junction);
      }
    }
  }

  if (this->ambig_prob_3 > 0.0) {
    /* Count ambig splice site as matches, but counts against found_score */
    this->total_splice_prob += this->ambig_prob_3;
    /* this->nmatches += this->querylength - qend; (inferred) */
    /* this->ref_nmatches += this->querylength - qend; (inferred) */
    this->found_score += this->querylength - qend;
    /* Does not affect score_within_trims */

  } else if (this->qend_alts != NULL) {
    this->nmatches += this->qend_alts->best_distal_nmatches; /* Not nmatches_to_trims, which is 0 for alts_substring */
    debug7(printf("Adding %d matches for qend_alts => total %d\n",this->qend_alts->best_distal_nmatches,this->nmatches));
    this->ref_nmatches += this->qend_alts->best_distal_nmatches; /* Not nmatches_to_trims, which is 0 for alts_substring */
    this->junction_splice_prob += this->qend_alts->best_medial_prob + this->qend_alts->best_distal_prob;
    this->total_splice_prob += this->qend_alts->best_medial_prob + this->qend_alts->best_distal_prob;
    this->found_score += this->qend_alts->best_distal_length - this->qend_alts->best_distal_nmatches;
    this->score_within_trims += this->qend_alts->best_distal_length - this->qend_alts->best_distal_nmatches;

  } else if (this->fusion_queryend_junction != NULL) {
    /* Rest of fusions handled below */

  } else {
    this->found_score += this->querylength - qend;
    /* Does not affect score_within_trims */
  }



  /* Fusion */
  if (this->fusion_querystart_junction != NULL || this->fusion_queryend_junction != NULL) {
    this->junction_splice_prob += Junction_prob(this->fusion_querystart_junction);
    this->junction_splice_prob += Junction_prob(this->fusion_queryend_junction);
    this->total_splice_prob += Junction_prob(this->fusion_querystart_junction);
    this->total_splice_prob += Junction_prob(this->fusion_queryend_junction);

    qstart = Intlist_head(this->fusion_endpoints);
    nmismatches = Intlist_head(this->fusion_nmismatches);
    ref_nmismatches = Intlist_head(this->fusion_ref_nmismatches);
    ninserts = 0;

    if (this->fusion_querystart_junction != NULL && this->fusion_plusp == this->plusp) {
      this->found_score += qstart;
    } else if (this->fusion_queryend_junction != NULL && this->fusion_plusp != this->plusp) {
      this->found_score += qstart;
    }

    if (this->fusion_plusp == true) {
      j = this->fusion_junctions;
      /* insertionp = false; */

      /* Add qpos to get alignstart/alignend */
      for (q = this->fusion_univdiagonals, x = this->fusion_nmismatches,
	     y = this->fusion_ref_nmismatches, r = Intlist_next(this->fusion_endpoints); q != NULL;
	   q = Univcoordlist_next(q), x = Intlist_next(x), y = Intlist_next(y), r = Intlist_next(r), j = List_next(j)) {
	qstart += ninserts;
	qend = Intlist_head(r);
#if 0
	if (insertionp == true) {
	  nmismatches = ref_nmismatches = -1; /* Recompute nmismatches */
	} else {
	  nmismatches = Intlist_head(x);
	  ref_nmismatches = Intlist_head(y);
	}
#else
	nmismatches = Intlist_head(x);
	ref_nmismatches = Intlist_head(y);
#endif
	univdiagonal = Univcoordlist_head(q);
	/* left = univdiagonal - (Univcoord_T) this->querylength; */
	debug7(printf("Path_eval_nmatches: ninserts %d, qstart %d..qend %d at univdiagonal %u [%u]\n",
		      ninserts,qstart,qend,univdiagonal,univdiagonal - this->chroffset));

	if (nmismatches >= 0 && ref_nmismatches >= 0) {
	  debug7(printf("Checking fusion mismatches, plus, at %u from querystart %d to queryend %d\n",
			univdiagonal - this->fusion_chroffset,qstart,qend));
	  debug7(printf("%d mismatches expected vs %d measured\n",
			nmismatches,
			Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_fwd,
							      univdiagonal,this->querylength,
							      /*pos5*/qstart,/*pos3*/qend,/*plusp*/true,this->genestrand)));
#ifdef CHECK_NMISMATCHES
	  assert(nmismatches == Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_fwd,
								      univdiagonal,this->querylength,
								      /*pos5*/qstart,/*pos3*/qend,/*plusp*/true,this->genestrand));
#endif
	} else {
	  nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_fwd,
							      univdiagonal,this->querylength,
							      /*pos5*/qstart,/*pos3*/qend,/*plusp*/true,this->genestrand);
	  Intlist_head_set(x,nmismatches);		/* Save for Stage3end_new_substrings */
	  Intlist_head_set(y,ref_nmismatches);		/* Save for Stage3end_new_substrings */
	  debug7(printf("%d (%d ref) mismatches from genome over querypos %d..%d\n",
			nmismatches,ref_nmismatches,qstart,qend));
	}

	/* Could potentially check here if qstart < qend, but relying upon caller to use endpoints_acceptable_p */
	this->nmatches += (qend - qstart) - nmismatches;
	this->ref_nmatches += (qend - qstart) - ref_nmismatches;
	debug7(printf("(3) Shortcut adds matches of %d = (%d - %d) - nmismatches %d => total %d\n",
		      (qend-qstart)-nmismatches,qend,qstart,nmismatches,this->nmatches));
	this->found_score += nmismatches;
	this->score_within_trims += nmismatches;

	/* Prepare for next iteration */
	qstart = qend;

	if (j == NULL) {
	  ninserts = 0;
	} else if ((junction = (Junction_T) List_head(j)) == NULL) {
	  /* qstart_junction */
	  ninserts = 0;
	} else {
	  debug7(printf("Junction: ")); debug7(Junction_print(junction)); debug7(printf("\n"));
	  ninserts = Junction_ninserts(junction);
	  this->junction_splice_prob += Junction_prob(junction);
	  this->total_splice_prob += Junction_prob(junction);
	}
      }

    } else {
      j = this->fusion_junctions; /* Put here before we handle querystart_alts */

      /* Subtract qpos to get alignstart/alignend */
      for (q = this->fusion_univdiagonals, x = this->fusion_nmismatches,
	     y = this->fusion_ref_nmismatches, r = Intlist_next(this->fusion_endpoints); q != NULL;
	   q = Univcoordlist_next(q), x = Intlist_next(x), y = Intlist_next(y), r = Intlist_next(r), j = List_next(j)) {
	qstart += ninserts;
	qend = Intlist_head(r);
#if 0
	if (insertionp == true) {
	  nmismatches = ref_nmismatches = -1; /* Recompute nmismatches */
	} else {
	  nmismatches = Intlist_head(x);
	  ref_nmismatches = Intlist_head(y);
	}
#else
	nmismatches = Intlist_head(x);
	ref_nmismatches = Intlist_head(y);
#endif

	univdiagonal = Univcoordlist_head(q);
	/* left = univdiagonal - (Univcoord_T) this->querylength; */
	debug7(printf("Path_eval_nmatches: ninserts %d, qstart %d..qend %d at univdiagonal %u [%u]\n",
		      ninserts, qstart,qend,univdiagonal,univdiagonal - this->chroffset));

	if (nmismatches >= 0 && ref_nmismatches >= 0) {
	  debug7(printf("Checking fusion mismatches at %u, minus, from querystart %d to queryend %d\n",
			univdiagonal - this->fusion_chroffset,qstart,qend));
	  debug7(printf("%d mismatches expected vs %d measured\n",
			nmismatches,Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_rev,
									  univdiagonal,this->querylength,
									  /*pos5*/qstart,/*pos3*/qend,/*plusp*/false,this->genestrand)));
#ifdef CHECK_NMISMATCHES
	  assert(nmismatches == Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_rev,
								      univdiagonal,this->querylength,
								      /*pos5*/qstart,/*pos3*/qend,/*plusp*/false,this->genestrand));
#endif
	} else {
	  nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_rev,
							      univdiagonal,this->querylength,
							      /*pos5*/qstart,/*pos3*/qend,/*plusp*/false,this->genestrand);
	  Intlist_head_set(x,nmismatches);		/* Save for Stage3end_new_substrings */
	  Intlist_head_set(y,ref_nmismatches);		/* Save for Stage3end_new_substrings */
	  debug7(printf("%d (%d ref) mismatches from genome over querypos %d..%d\n",
			nmismatches,ref_nmismatches,this->querylength - qend,this->querylength - qstart));
	}

	/* Could potentially check here if qstart < qend, but relying upon caller to use endpoints_acceptable_p */
	this->nmatches += (qend - qstart) - nmismatches;
	this->ref_nmatches += (qend - qstart) - ref_nmismatches;
	debug7(printf("(4) Shortcut adds matches of %d = (%d - %d) - nmismatches %d => total %d\n",
		      (qend-qstart)-nmismatches,this->querylength - qstart,this->querylength - qend,nmismatches,this->nmatches));
	this->found_score += nmismatches;
	this->score_within_trims += nmismatches;

	/* Prepare for next iteration */
	qstart = qend;

	if (j == NULL) {
	  ninserts = 0;
	} else if ((junction = (Junction_T) List_head(j)) == NULL) {
	  /* qstart_junction */
	  ninserts = 0;
	} else {
	  debug7(printf("Junction: ")); debug7(Junction_print(junction)); debug7(printf("\n"));
	  ninserts = Junction_ninserts(junction);
	  this->junction_splice_prob += Junction_prob(junction);
	  this->total_splice_prob += Junction_prob(junction);
	}
      }
    }

    if (this->fusion_querystart_junction != NULL && this->fusion_plusp != this->plusp) {
      this->found_score += this->querylength - qend;
    } else if (this->fusion_queryend_junction != NULL && this->fusion_plusp == this->plusp) {
      this->found_score += this->querylength - qend;
    }
  }


  this->nmatches += total_ninserts;

  if (this->found_score < *found_score) {
    *found_score = this->found_score;
  }

  debug7(printf("Path_eval_nmatches returning %d matches in %s for score of %d within trims\n",
		this->nmatches,Intlist_to_string(this->endpoints),this->score_within_trims));

  assert(this->nmatches <= this->querylength);
  return this->nmatches;
}



void
Path_mark_alignment (T path, Compress_T query_compress_fwd, char *queryuc_ptr,
		     Compress_T query_compress_rev, char *queryrc, Pathpool_T pathpool) {

  Compress_T query_compress;
  Univcoord_T univdiagonal;
  Intlist_T q, n;
  Univcoordlist_T p;
  List_T j;
  int nmatches_exonic;
  int pos5, pos3, ninserts;

  if (path->genomic_diff != NULL) {
    /* Already marked */
    return;

  } else {
    path->genomic_diff = Pathpool_new_string(pathpool,path->querylength+1);
    /* genomic_bothdiff = (char *) MALLOC((path->querylength+1) * sizeof(char)); */
    /* Genome_fill_buffer(left,querylength,genomic_diff); */
  }

  if (path->plusp == true) {
    query_compress = query_compress_fwd;
    strcpy(path->genomic_diff,queryuc_ptr); /* Start with query sequence on genomic plus strand */
  } else {
    query_compress = query_compress_rev;
    strcpy(path->genomic_diff,queryrc); /* Start with query sequence on genomic plus strand */
  }

  ninserts = 0;
  for (p = path->univdiagonals, q = path->endpoints, n = path->nmismatches, j = path->junctions; p != NULL;
       p = Univcoordlist_next(p), q = Intlist_next(q), n = Intlist_next(n), j = List_next(j)) {
    assert(Intlist_head(n) >= 0);
    if (Intlist_head(n) > 0 || Compress_non_acgt(query_compress) == true || md_report_snps_p == true) {
      univdiagonal = Univcoordlist_head(p);
      pos5 = Intlist_head(q) + ninserts;
      pos3 = Intlist_head(Intlist_next(q));
      Genomebits_mark_mismatches(&nmatches_exonic,path->genomic_diff,
				 query_compress,univdiagonal,path->querylength,
				 pos5,pos3,/*segment_plusp*/path->plusp,
				 /*query_plusp*/path->plusp,path->genestrand);
    }

    if (j != NULL) {
      ninserts = Junction_ninserts((Junction_T) List_head(j));
    }
  }


  /* Fusion */
  if (path->fusion_querystart_junction != NULL || path->fusion_queryend_junction != NULL) {
    if (path->fusion_plusp == true) {
      query_compress = query_compress_fwd;
    } else {
      query_compress = query_compress_rev;
    }

    ninserts = 0;
    for (p = path->fusion_univdiagonals, q = path->fusion_endpoints,
	   n = path->fusion_nmismatches, j = path->fusion_junctions; p != NULL;
	 p = Univcoordlist_next(p), q = Intlist_next(q), n = Intlist_next(n), j = List_next(j)) {
      if (Intlist_head(n) > 0 || Compress_non_acgt(query_compress) == true || md_report_snps_p == true) {
	univdiagonal = Univcoordlist_head(p);
	pos5 = Intlist_head(q) + ninserts;
	pos3 = Intlist_head(Intlist_next(q));
	Genomebits_mark_mismatches(&nmatches_exonic,path->genomic_diff,
				   query_compress,univdiagonal,path->querylength,
				   pos5,pos3,/*segment_plusp*/path->fusion_plusp,
				   /*query_plusp*/path->plusp,path->genestrand);
      }
    }

    if (j != NULL) {
      ninserts = Junction_ninserts((Junction_T) List_head(j));
    }
  }

  return;
}


bool
Path_eval_perfect_ends_p (T this, Compress_T query_compress_fwd, char *queryuc_ptr,
			  Compress_T query_compress_rev, char *queryrc,
			  int querystart, int queryend, Pathpool_T pathpool) {
  int mod;
  int querypos5, querypos3;

  if (this->found_score == 0) {
    return true;
  } else {
    Path_mark_alignment(this,query_compress_fwd,queryuc_ptr,
			query_compress_rev,queryrc,pathpool);
    for (mod = 0; mod < index1interval; mod++) {
      querypos5 = querystart + mod;
      querypos3 = queryend - mod; /* Typically query_lastpos - mod */
      if (islower(this->genomic_diff[querypos5])) {
	return false;
      } else if (islower(this->genomic_diff[querypos5 + index1part - 1])) {
	return false;
      } else if (islower(this->genomic_diff[querypos3])) {
	return false;
      } else if (islower(this->genomic_diff[querypos3 + index1part - 1])) {
	return false;
      }
    }

    return true;
  }
}


static int
Path_method_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  if (a->method > b->method) {
    return -1;
  } else if (b->method > a->method) {
    return +1;
  } else {
    return 0;
  }
}


/* Used to identify best paths for auxinfo */
int
Path_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  int coverage_a, coverage_b;
  int nbadsplices_a, nbadsplices_b;

  coverage_a = Path_coverage(a);
  coverage_b = Path_coverage(b);

  nbadsplices_a = Path_nbadsplices(a);
  nbadsplices_b = Path_nbadsplices(b);

  if (coverage_a > coverage_b + 20) {
    return -1;
  } else if (coverage_b > coverage_a + 20) {
    return +1;

  } else if (nbadsplices_a < nbadsplices_b) {
    return -1;
  } else if (nbadsplices_b < nbadsplices_a) {
    return +1;

  } else if (a->nmatches > b->nmatches) {
    return -1;
  } else if (b->nmatches > a->nmatches) {
    return +1;

  } else {
    return 0;
  }
}



/* Duplicates with respect to method have already been taken care of */
/* Ignore sensedir, so we keep both sensedirs if they have equivalent matches */
static int
Path_local_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  int coverage_a, coverage_b;
  int nbadsplices_a, nbadsplices_b;
  int nalts_a = 0, nalts_b = 0;

  coverage_a = Path_coverage(a);
  coverage_b = Path_coverage(b);

  nbadsplices_a = Path_nbadsplices(a);
  nbadsplices_b = Path_nbadsplices(b);

  if (a->qstart_alts != NULL) {
    nalts_a++;
  }
  if (a->qend_alts != NULL) {
    nalts_a++;
  }
  if (b->qstart_alts != NULL) {
    nalts_b++;
  }
  if (b->qend_alts != NULL) {
    nalts_b++;
  }


  if (coverage_a > coverage_b + 20) {
    return -1;
  } else if (coverage_b > coverage_a + 20) {
    return +1;

  } else if (nbadsplices_a < nbadsplices_b) {
    return -1;
  } else if (nbadsplices_b < nbadsplices_a) {
    return +1;

  } else if (a->nmatches > b->nmatches) {
    return -1;
  } else if (b->nmatches > a->nmatches) {
    return +1;

  } else if (a->transcripts != NULL && b->transcripts == NULL) {
    return -1;
  } else if (b->transcripts != NULL && a->transcripts == NULL) {
    return +1;

  } else if (nalts_a < nalts_b) {
    return -1;
  } else if (nalts_b < nalts_a) {
    return +1;


#if 0
  /* We want to keep all possible transcript results at a locus */
  } else if (a->junction_splice_prob > b->junction_splice_prob) {
    return -1;
  } else if (b->junction_splice_prob > a->junction_splice_prob) {
    return +1;
  } else if (a->total_splice_prob > b->total_splice_prob) {
    return -1;
  } else if (b->total_splice_prob > a->total_splice_prob) {
    return +1;
#endif

  } else {
    return 0;
  }
}


static int
Path_global_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  int coverage_a, coverage_b;
  int nbadsplices_a, nbadsplices_b;

  coverage_a = Path_coverage(a);
  coverage_b = Path_coverage(b);

  nbadsplices_a = Path_nbadsplices(a);
  nbadsplices_b = Path_nbadsplices(b);

  if (a->transcripts != NULL && b->transcripts == NULL) {
    return -1;
  } else if (b->transcripts != NULL && a->transcripts == NULL) {
    return +1;

  } else if (coverage_a > coverage_b + 20) {
    return -1;
  } else if (coverage_b > coverage_a + 20) {
    return +1;

  } else if (nbadsplices_a < nbadsplices_b) {
    return -1;
  } else if (nbadsplices_b < nbadsplices_a) {
    return +1;

  } else if (a->nmatches > b->nmatches) {
    return -1;
  } else if (b->nmatches > a->nmatches) {
    return +1;
  } else {
    return 0;
  }
}
  

static bool
Path_main_fusion_equal (T a, T b) {

  if (a->fusion_endpoints == NULL) {
    return false;
  } else if (b->fusion_endpoints == NULL) {
    return false;
  } else if (Intlist_equal(a->endpoints,b->fusion_endpoints) == false) {
    return false;
  } else if (Intlist_equal(a->fusion_endpoints,b->endpoints) == false) {
    return false;
  } else if (Univcoordlist_equal(a->univdiagonals,b->fusion_univdiagonals) == false) {
    return false;
  } else if (Univcoordlist_equal(a->fusion_univdiagonals,b->univdiagonals) == false) {
    return false;
  } else {
    return true;
  }
}


T *
Path_eval_and_sort (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
		    T *patharray, int npaths, 
		    Compress_T query_compress_fwd, Compress_T query_compress_rev,
		    char *queryuc_ptr, char *queryrc, char *quality_string,
		    int nmismatches_filter, int mincoverage_filter,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		    Hitlistpool_T hitlistpool, bool filterp) {
  int pathi;
  float maxlik, loglik;
  float total, qual;		/* For Bayesian mapq calculation */

  int randomi, i, j, k, l;
  T temp, path;
  double best_sense_prob, best_antisense_prob;


  debug8(printf("Entered Path_eval_and_sort with patharray %p and npaths %d\n",patharray,npaths));
  if (npaths == 0) {
    /* Skip */
    *npaths_primary = *npaths_altloc = 0;
    *first_absmq = 0;
    *second_absmq = 0;

  } else {
    /* Have already called Path_extend */
    
    /* 0.  Unalias circular alignments and trim chrbounds */
    for (i = 0; i < npaths; i++) {
      path = patharray[i];
      if (circularp[path->chrnum] == true) {
	Path_trim_circular_unalias(path);
      }
      Path_trim_chrbounds(path,query_compress_fwd,query_compress_rev,
			  intlistpool,univcoordlistpool,listpool,pathpool);
    }


    /* 1.  Sort by structure to remove duplicates */
    if (npaths > 1) {
      qsort(patharray,npaths,sizeof(T),Path_structure_cmp);

      k = 0;
      i = 0;
      while (i < npaths) {
	j = i + 1;
	while (j < npaths && Path_structure_cmp(&(patharray[j]),&(patharray[i])) == 0) {
	  j++;
	}
	debug8(printf("Found an identical group by structure (except sensedir) of %d paths => Re-sorting by method_cmp\n",j - i));
	
	qsort(&(patharray[i]),j - i,sizeof(T),Path_method_cmp);
	debug8(printf("(1) Keeping by method_cmp ")); debug8(Path_print(patharray[i]));
	patharray[k++] = patharray[i];
	
	for (l = i + 1; l < j; l++) {
	  debug8(printf("(1) Eliminating by method_cmp ")); debug8(Path_print(patharray[l]));
	  path = patharray[l];
	  Path_free(&path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
	}
	
	i = j;
      }
      npaths = k;
    }
    

    /* 2.  Sort by intervals to find best sensedir at a given locus */
    if (npaths > 1) {
      qsort(patharray,npaths,sizeof(T),Path_interval_cmp);
    
      k = 0;
      i = 0;
      while (i < npaths) {
	j = i + 1;
	while (j < npaths && Path_overlap_p(patharray[j],patharray[i]) == true) {
	  j++;
	}
	debug8(printf("Found an overlapping group of %d => choosing sense vs antisense\n",j - i));

	best_sense_prob = best_antisense_prob = 0.0;
	for (l = i; l < j; l++) {
	  path = patharray[l];
	  if (Path_nbadsplices(path) > 0) {
	    /* Ignore */
	  } else if (path->sensedir == SENSE_FORWARD) {
	    if (path->total_splice_prob > best_sense_prob) {
	      best_sense_prob = path->total_splice_prob;
	    }
	  } else {
	    if (path->total_splice_prob > best_antisense_prob) {
	      best_antisense_prob = path->total_splice_prob;
	    }
	  }
	}

	if (best_sense_prob > best_antisense_prob) {
	  /* Keep only sense */
	  for (l = i; l < j; l++) {
	    path = patharray[l];
	    if (path->sensedir == SENSE_FORWARD) {
	      debug8(printf("(2 sensedir) Keeping sense ")); debug8(Path_print(path));
	      patharray[k++] = path;
	    } else {
	      debug8(printf("(2 sensedir) Eliminating antisense ")); debug8(Path_print(path));
	      Path_free(&path,intlistpool,univcoordlistpool,
			listpool,pathpool,transcriptpool,hitlistpool);
	    }
	  }
	} else if (best_antisense_prob > best_sense_prob) {
	  /* Keep only antisense */
	  for (l = i; l < j; l++) {
	    path = patharray[l];
	    if (path->sensedir == SENSE_ANTI) {
	      debug8(printf("(2 sensedir) Keeping anti ")); debug8(Path_print(path));
	      patharray[k++] = path;
	    } else {
	      debug8(printf("(2 sensedir) Eliminating sense ")); debug8(Path_print(path));
	      Path_free(&path,intlistpool,univcoordlistpool,
			listpool,pathpool,transcriptpool,hitlistpool);
	    }
	  }
	} else {
	  /* Keep both */
	  for (l = i; l < j; l++) {
	    path = patharray[l];
	    patharray[k++] = path;
	  }
	}

	i = j;
      }
      npaths = k;
    }


    /* 3.  Sort by intervals to best paths at a given locus */
    if (npaths > 1) {
      /* Should already be sorted */
      /* qsort(patharray,npaths,sizeof(T),Path_interval_cmp); */
    
      k = 0;
      i = 0;
      while (i < npaths) {
	j = i + 1;
	while (j < npaths && Path_overlap_p(patharray[j],patharray[i]) == true) {
	  j++;
	}
	debug8(printf("Found an overlapping group of %d => re-sorting by Path_local_cmp\n",j - i));

	/* Keep the best ones in the overlapping group */
	qsort(&(patharray[i]),j - i,sizeof(T),Path_local_cmp);
	debug8(printf("(3 best) Keeping by local_cmp ")); debug8(Path_print(patharray[i]));
	patharray[k++] = patharray[i];
      
	for (l = i + 1; l < j; l++) {
	  if (Path_local_cmp(&(patharray[l]),&(patharray[i])) == 0) {
	    debug8(printf("(3 tie) Keeping by local_cmp ")); debug8(Path_print(patharray[l]));
	    patharray[k++] = patharray[l];
	  } else {
	    debug8(printf("(3 worse) Eliminating by local_cmp ")); debug8(Path_print(patharray[l]));
	    path = patharray[l];
	    Path_free(&path,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);
	  }
	}

	i = j;
      }
      npaths = k;
    }


    /* 4.  Find best solution globally */
    /* TODO: Can make this O(n) rather than O(n*log n) */
    if (npaths > 1) {
      qsort(patharray,npaths,sizeof(T),Path_global_cmp);
    }
    debug8(printf("Found the global best solution with %d nmatches\n",patharray[0]->nmatches));

    /* 5.  Check if we should be filtering the result */
    path = patharray[0];
    if (filterp == true && path->score_within_trims > nmismatches_filter) {
      debug8(printf("(4 filter) Best solution has too many %d nmismatches, so eliminating all\n",
		    path->score_within_trims));
      for (k = 0; k < npaths; k++) {
	debug8(printf("(4 filter) Eliminating ")); debug8(Path_print(patharray[k]));
	path = patharray[k];
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
      }
      FREE_OUT(patharray);
      *npaths_primary = *npaths_altloc = 0;
      return (T *) NULL;

    } else if (filterp == true && Path_coverage(path) < mincoverage_filter) {
      debug8(printf("(4 filter) Best solution has too little %d coverage, so eliminating all\n",
		    Path_coverage(path)));
      for (k = 0; k < npaths; k++) {
	debug8(printf("(4 filter) Eliminating ")); debug8(Path_print(patharray[k]));
	path = patharray[k];
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
      }
      FREE_OUT(patharray);
      *npaths_primary = *npaths_altloc = 0;
      return (T *) NULL;

    } else if (filterp == true && allow_soft_clips_p == false &&
	       (Intlist_head(path->endpoints) != 0 || Intlist_last_value(path->endpoints) != path->querylength)) {
      debug8(printf("(4 filter) Best solution has soft clips, so eliminating all\n",
		    Path_coverage(path)));
      for (k = 0; k < npaths; k++) {
	debug8(printf("(4 filter) Eliminating ")); debug8(Path_print(patharray[k]));
	path = patharray[k];
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
      }
      FREE_OUT(patharray);
      *npaths_primary = *npaths_altloc = 0;
      return (T *) NULL;

    } else {
      /* Otherwise, keep all paths equivalent to the best path */
      debug8(printf("(4 best) Keeping best ")); debug8(Path_print(patharray[0]));
      i = 1;			/* Skip the best path */

      for (k = 1; k < npaths; k++) {
	path = patharray[k];
	if (Path_global_cmp(&path,&(patharray[0])) == 0) {
	  debug8(printf("(4 tie) Keeping ")); debug8(Path_print(path));
	  patharray[i++] = path;
	} else {
	  debug8(printf("(4 worse) Eliminating ")); debug8(Path_print(path));
	  Path_free(&path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
	}
      }

      npaths = i;			
    }


    /* 4.  Check for main <=> fusion identity (can occur in single-end
       reads without a paired-end as anchor) */
    if (npaths > 1) {
      for (i = 0; i < npaths; i++) {
	if (patharray[i] != NULL) {
	  for (j = i + 1; j < npaths; j++) {
	    if (patharray[j] != NULL) {
	      if (Path_main_fusion_equal(patharray[j],patharray[i]) == true) {
		debug8(printf("(6) Eliminating ")); debug8(Path_print(patharray[j]));
		path = patharray[j];
		Path_free(&path,intlistpool,univcoordlistpool,
			  listpool,pathpool,transcriptpool,hitlistpool);
		patharray[j] = (Path_T) NULL;
	      }
	    }
	  }
	}
      }

      k = 0;
      for (i = 0; i < npaths; i++) {
	if (patharray[i] != NULL) {
	  patharray[k++] = patharray[i];
	}
      }
      npaths = k;
    }

    if (want_random_p && npaths > 1) {
      /* Randomize among best alignments */
      /* randomi = (int) ((double) i * rand()/((double) RAND_MAX + 1.0)); */
      randomi = (int) (rand() / (((double) RAND_MAX + 1.0) / (double) npaths));
      /* fprintf(stderr,"%d dups => random %d\n",i,randomi); */
      temp = patharray[0];
      patharray[0] = patharray[randomi];
      patharray[randomi] = temp;
    }

    /* Trim alignments at chromosomal bounds */
    for (i = 0; i < npaths; i++) {
      debug8(printf("(8) Trimming ")); debug8(Path_print(patharray[i]));
      path = patharray[i];
      if (circularp[path->chrnum] == false) {
#if 0
	/* Now done at the start of this procedure */
	debug8(printf("Chrnum %d is not circular\n",path->chrnum));
	Path_trim_qstart_trimbounds(path,intlistpool,univcoordlistpool,listpool,pathpool,
				    /*trimbounds*/path->chroffset);
	Path_trim_qend_trimbounds(path,intlistpool,univcoordlistpool,listpool,pathpool,
				  /*trimbounds*/path->chrhigh);
#endif

      } else if (output_type == STD_OUTPUT || output_type == M8_OUTPUT) {
	/* If output type is alignment or m8, then don't want to split up the parts */

      } else if (path->plusp == true) {
	Path_trim_circular(path,query_compress_fwd,intlistpool,univcoordlistpool,listpool);

      } else {
	Path_trim_circular(path,query_compress_rev,intlistpool,univcoordlistpool,listpool);
      }
    }

    /* Compute mapq_loglik */
    if (npaths == 1) {
      path = patharray[0];
      Path_mark_alignment(path,query_compress_fwd,queryuc_ptr,query_compress_rev,queryrc,
			  pathpool);
      path->mapq_loglik = MAPQ_MAXIMUM_SCORE;
      path->mapq_score = MAPQ_max_quality_score(quality_string,path->querylength);
      path->absmq_score = MAPQ_MAXIMUM_SCORE;

      *first_absmq = path->absmq_score;
      *second_absmq = 0;

    } else {
      for (i = 0; i < npaths; i++) {
	path = patharray[i];
	Path_mark_alignment(path,query_compress_fwd,queryuc_ptr,query_compress_rev,queryrc,
			    pathpool);
	path->mapq_loglik =
	  MAPQ_loglik_string(path->genomic_diff,quality_string,path->querylength,path->plusp);
      }
    }
    

    /* Enforce monotonicity */
    for (i = npaths - 1; i > 0; i--) {
      if (patharray[i-1]->mapq_loglik < patharray[i]->mapq_loglik) {
	patharray[i-1]->mapq_loglik = patharray[i]->mapq_loglik;
      }
    }
    maxlik = patharray[0]->mapq_loglik;
    
    /* Subtract maxlik to avoid underflow */
    for (i = 0; i < npaths; i++) {
      patharray[i]->mapq_loglik -= maxlik;
    }

    /* Compute absolute mapq */
    for (i = 0; i < npaths; i++) {
      loglik = patharray[i]->mapq_loglik + MAPQ_MAXIMUM_SCORE;
      if (loglik < 0.0) {
	loglik = 0.0;
      }
      patharray[i]->absmq_score = rint(loglik);
    }
    *first_absmq = patharray[0]->absmq_score;
    if (npaths == 1) {
      *second_absmq = 0;
    } else {
      *second_absmq = patharray[1]->absmq_score;
    }

    /* Compute Bayesian mapq */
    total = 0.0;
    for (i = 0; i < npaths; i++) {
      total += (patharray[i]->mapq_loglik = fasterexp(patharray[i]->mapq_loglik));
    }

    /* Obtain posterior probabilities of being true */
    for (i = 0; i < npaths; i++) {
      patharray[i]->mapq_loglik /= total;
    }

    /* Convert to Phred scores */
    for (i = 0; i < npaths; i++) {
      if ((qual = 1.0 - patharray[i]->mapq_loglik) < 2.5e-10 /* 10^-9.6 */) {
	patharray[i]->mapq_score = 40;
      } else {
	patharray[i]->mapq_score = rint(-10.0 * log10(qual));
      }
    }
  }

  debug8(printf("Exiting Path_eval_and_sort with %d paths\n",npaths));
  
  if (transcriptome != NULL) {
    for (i = 0; i < npaths; i++) {
      path = patharray[i];
      Transcript_velocity_single(path);
    }
  }

  /* Filter for chrsubset */
  /* Want to allow other alignments to be found before filtering */
  *npaths_primary = *npaths_altloc = 0;
  if (chrsubsetp == NULL) {
    for (pathi = 0; pathi < npaths; pathi++) {
      path = patharray[pathi];
      if (altlocp[path->chrnum] == true) {
	(*npaths_altloc) += 1;
      } else {
	(*npaths_primary) += 1;
      }
    }
    return patharray;

  } else {
    k = 0;
    for (pathi = 0; pathi < npaths; pathi++) {
      path = patharray[pathi];
      if (chrsubsetp[path->chrnum] == false) {
	/* Do not save this hit */
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	
      } else {
	/* Save this hit.  Re-use existing array */
	if (altlocp[path->chrnum] == true) {
	  (*npaths_altloc) += 1;
	} else {
	  (*npaths_primary) += 1;
	}
	patharray[k++] = path;
      }
    }

    if ((*npaths_primary) + (*npaths_altloc) == 0) {
      FREE_OUT(patharray);
      return (T *) NULL;
    } else {
      return patharray;
    }
  }
}


#if 0
static int
Path_nmatches_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  if (a->nmatches > b->nmatches) {
    return -1;
  } else if (b->nmatches > a->nmatches) {
    return +1;
  } else if (a->transcriptome_method_p > b->transcriptome_method_p) {
    return -1;
  } else if (b->transcriptome_method_p > a->transcriptome_method_p) {
    return +1;
  } else {
    return 0;
  }
}
#endif


static bool
trims_compatible_p (int common_trim_qstart, int common_trim_qend,
		    int qstart_trim, int qend_trim) {
  if (common_trim_qstart == 0 && common_trim_qend == 0) {
    return true;
  } else if (common_trim_qstart > 0 && common_trim_qend > 0) {
    /* Trimming on both ends by a single path */
    return true;
  } else if (qend_trim > 0 && common_trim_qstart > 0) {
    /* Trimming on a different end from the equivalence class */
    return false;
  } else if (qstart_trim > 0 && common_trim_qend > 0) {
    /* Trimming on a different end from the equivalence class */
    return false;
  } else {
    return true;
  }
}


/* Allowing invalid transcripts but remapping them to the trimmed
   transcript and moving them to valid if appropriate */
/* We have to allow invalid transcripts because Trpath_convert
   procedures call Path_solve_junctions on unsolved paths, resulting
   in invalid transcripts */
#define ALLOW_INVALID_TRANSCRIPTS 1


List_T
Path_consolidate (List_T paths, Shortread_T queryseq,
		  Compress_T query_compress_fwd, Compress_T query_compress_rev,
		  Uintlistpool_T uintlistpool, Intlistpool_T intlistpool,
		  Univcoordlistpool_T univcoordlistpool,
		  Listpool_T listpool, Pathpool_T pathpool,
		  Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  T *patharray;
  List_T class, c;
  List_T equiv_classes, p;
  Intlist_T qstart_trims, qend_trims, q, r;
  T bestpath, path;
  int qstart_trim, qend_trim, common_trim_qstart, common_trim_qend;
  int npaths, i, j, k;
  int best_nmatches;
  int nbest;

  List_T invalid_transcripts, t;
  Transcript_T transcript;
  bool assignedp;


  debug9(printf("Entered Path_consolidate with %d paths\n",List_length(paths)));
  if ((npaths = List_length(paths)) > 1) {
    patharray = (T *) List_to_array(paths,NULL);
    paths = (List_T) NULL;

    qsort(patharray,npaths,sizeof(T),Path_interval_cmp);
    
    i = 0;
    while (i < npaths) {
      /* Consider only paths with the most nmatches.  Use a linear pass instead of sorting */
      best_nmatches = patharray[i]->nmatches;
      bestpath = patharray[i];
      nbest = 1;

      j = i + 1;
      while (j < npaths && Path_overlap_p(patharray[j],patharray[i]) == true) {
	if (patharray[j]->nmatches > best_nmatches) {
	  bestpath = patharray[j];
	  nbest = 1;
	} else if (patharray[j]->nmatches == best_nmatches) {
	  nbest++;
	}
	j++;
      }
      debug9(printf("Found an overlapping group of %d starting at %d, best nmatches %d with %d paths\n",
		    j - i,i,best_nmatches,nbest));
      
      /* Handle the best paths in this overlapping group */
      if (nbest == 1) {
	/* Only a single best path */
	debug9(printf("Keeping by best_nmatches ")); debug9(Path_print(bestpath));
	paths = Hitlist_push(paths,hitlistpool,(void *) bestpath
			     hitlistpool_trace(__FILE__,__LINE__));

	/* Remap path here if transcriptome_method_p is false? */

	for (k = i; k < j; k++) {
	  if (patharray[k]->nmatches < best_nmatches) {
	    path = patharray[k];
	    debug9(printf("(1) Eliminating by best_nmatches ")); debug9(Path_print(path));
	    Path_free(&path,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);
	  }
	}

      } else {
	/* Among the nbest paths, find equivalence classes based on Path_common_structure_p */
	equiv_classes = (List_T) NULL;
	qstart_trims = (Intlist_T) NULL;
	qend_trims = (Intlist_T) NULL;

	for (k = i; k < j; k++) {
	  path = patharray[k];
	  if (path->nmatches < best_nmatches) {
	    debug9(printf("(2) Eliminating by best_nmatches ")); debug9(Path_print(path));
	    Path_free(&path,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);

	  } else {
	    p = equiv_classes;
	    q = qstart_trims;
	    r = qend_trims;
	    while (p != NULL &&
		   (path->transcriptome_method_p != ((T) List_head((List_T) List_head(p)))->transcriptome_method_p ||
		    Path_common_structure_p(&common_trim_qstart,&common_trim_qend,path,
					    (T) List_head((List_T) List_head(p))) == false ||
		    trims_compatible_p(common_trim_qstart,common_trim_qend,
				       Intlist_head(q),Intlist_head(r)) == false)) {
	      p = List_next(p);
	      q = Intlist_next(q);
	      r = Intlist_next(r);
	    }

	    if (p == NULL) {
	      /* Create new equivalence class */
	      debug9(printf("Creating new equivalence class with initial trims 0 and 0\n"));
	      equiv_classes = Listpool_push(equiv_classes,listpool,
					    (void *) Listpool_push(NULL,listpool,(void *) path
								   listpool_trace(__FILE__,__LINE__))
					    listpool_trace(__FILE__,__LINE__));
	      qstart_trims = Intlistpool_push(qstart_trims,intlistpool,0
					      intlistpool_trace(__FILE__,__LINE__));
	      qend_trims = Intlistpool_push(qend_trims,intlistpool,0
					    intlistpool_trace(__FILE__,__LINE__));

	    } else {
	      /* Insert path into equivalence class and revise trims */
	      debug9(printf("Putting path with trims %d and %d into existing equivalence class\n",
			    common_trim_qstart,common_trim_qend));

	      class = (List_T) List_head(p);
	      class = Listpool_push(class,listpool,(void *) path
				    listpool_trace(__FILE__,__LINE__));
	      List_head_set(p,(void *) class);
	      debug9(printf("Revising trims from %d and %d",Intlist_head(q),Intlist_head(r)));
	      if (common_trim_qstart > Intlist_head(q)) {
		Intlist_head_set(q,common_trim_qstart);
	      }
	      if (common_trim_qend > Intlist_head(r)) {
		Intlist_head_set(r,common_trim_qend);
	      }
	      debug9(printf(" to %d and %d\n",Intlist_head(q),Intlist_head(r)));
	    }
	  }
	}
	      
	debug9(printf("Found %d equivalence classes\n",List_length(equiv_classes)));
	for (p = equiv_classes, q = qstart_trims, r = qend_trims;
	     p != NULL; p = List_next(p), q = Intlist_next(q), r = Intlist_next(r)) {
	  /* Process each equivalence class */
	  class = (List_T) List_head(p);

#if 0
	  /* Not needed, since we have separate equivalence classes for transcript and non-transcript methods */
	  /* Find a representative from a transcriptome method */
	  /* non_transcriptomep = false; */
	  bestpath = (T) NULL;
	  for (c = class; c != NULL; c = List_next(c)) {
	    path = (T) List_head(c);
	    if (path->transcriptome_method_p == true) {
	      bestpath = path;
	    } else {
	      /* non_transcriptomep = true; */
	    }
	  }
	  if (bestpath == NULL) {
	    bestpath = (T) List_head(class);
	  }
#else
	  bestpath = (T) List_head(class);
	  debug9(printf("Best path: ")); debug9(Path_print(bestpath));
#endif
	  for (c = class; c != NULL; c = List_next(c)) {
	    path = (T) List_head(c);
	    if (path != bestpath) {
	      debug9(printf("Eliminating ")); debug9(Path_print(path));
	      bestpath->transcripts = List_append(path->transcripts,bestpath->transcripts);
	      path->transcripts = (List_T) NULL;
#ifdef ALLOW_INVALID_TRANSCRIPTS
	      /* If we allow invalid transcripts from transcriptome methods */
	      bestpath->invalid_transcripts = List_append(path->invalid_transcripts,bestpath->invalid_transcripts);
	      path->invalid_transcripts = (List_T) NULL;
#endif
	      Path_free(&path,intlistpool,univcoordlistpool,
			listpool,pathpool,transcriptpool,hitlistpool);
	    }
	  }
	  
	  qstart_trim = Intlist_head(q);
	  qend_trim = Intlist_head(r);
	  debug9(printf("qstart_trim %d, qend_trim %d\n",qstart_trim,qend_trim));
	  if (qstart_trim > 0) {
	    Path_trim_qstart_n(qstart_trim,bestpath,
			       query_compress_fwd,query_compress_rev,
			       intlistpool,univcoordlistpool,listpool,pathpool);
	  }
	  if (qend_trim > 0) {
	    Path_trim_qend_n(qend_trim,bestpath,
			     query_compress_fwd,query_compress_rev,
			     intlistpool,univcoordlistpool,listpool,pathpool);
	  }
	  if (bestpath->genestrand > 0) {
	    Transcript_list_trim(bestpath->transcripts,/*trim5*/qstart_trim,/*trim3*/qend_trim,
				 transcriptome,listpool,transcriptpool);
	  } else {
	    Transcript_list_trim(bestpath->transcripts,/*trim5*/qend_trim,/*trim3*/qstart_trim,
				 transcriptome,listpool,transcriptpool);
	  }

#ifdef ALLOW_INVALID_TRANSCRIPTS
	  /* Hard to trim invalid transcripts, so remap them */
	  invalid_transcripts = bestpath->invalid_transcripts;
	  bestpath->invalid_transcripts = (List_T) NULL;
	  for (t = invalid_transcripts; t != NULL; t = List_next(t)) {
	    transcript = (Transcript_T) List_head(t);
	    if ((assignedp = Transcript_remap_invalid(transcript,bestpath,transcriptome,queryseq,
						      uintlistpool,listpool,transcriptpool)) == true) {
	      Transcript_free(&transcript,listpool,transcriptpool);
	    } else {
	      bestpath->invalid_transcripts = Listpool_push(bestpath->invalid_transcripts,listpool,(void *) transcript
							    listpool_trace(__FILE__,__LINE__));
	    }
	  }
	  Listpool_free_list(&invalid_transcripts,listpool
			     listpool_trace(__FILE__,__LINE__));
#endif

	  bestpath->transcripts = Transcript_list_sort(bestpath->transcripts);
	  bestpath->invalid_transcripts = Transcript_list_sort(bestpath->invalid_transcripts);

	  debug9(printf("Keeping representative of class ")); debug9(Path_print(bestpath));
	  paths = Hitlist_push(paths,hitlistpool,(void *) bestpath
			       hitlistpool_trace(__FILE__,__LINE__));
	  Listpool_free_list(&class,listpool
			     listpool_trace(__FILE__,__LINE__));
	}

	Intlistpool_free_list(&qend_trims,intlistpool
			      intlistpool_trace(__FILE__,__LINE__));
	Intlistpool_free_list(&qstart_trims,intlistpool
			      intlistpool_trace(__FILE__,__LINE__));			      
	Listpool_free_list(&equiv_classes,listpool
			   listpool_trace(__FILE__,__LINE__));
      }

      i = j;
    }
	  
    FREE(patharray);
  }

#ifdef DEBUG9
  printf("Exiting Path_consolidate with %d paths:\n",List_length(paths));
  for (p = paths; p != NULL; p = List_next(p)) {
    Path_print((T) List_head(p));
  }
#endif

  return paths;
}


void
Path_eval_setup (Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in,
		 Transcriptome_T transcriptome_in,
		 bool *circularp_in, bool *chrsubsetp_in, bool *altlocp_in,
		 int index1part_in, int index1interval_in,
		 Outputtype_T output_type_in, bool md_report_snps_p_in,
		 bool want_random_p_in, bool allow_soft_clips_p_in) {

  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;
  transcriptome = transcriptome_in;

  circularp = circularp_in;
  chrsubsetp = chrsubsetp_in;
  altlocp = altlocp_in;

  index1part = index1part_in;
  index1interval = index1interval_in;

  output_type = output_type_in;
  md_report_snps_p = md_report_snps_p_in;
  want_random_p = want_random_p_in;
  allow_soft_clips_p = allow_soft_clips_p_in;

  return;
}

