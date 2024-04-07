static char rcsid[] = "$Id: 7702732fd19dd7818d56e7b5cb360d901fffa7dd $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "path-trim.h"
#include "path-eval.h"

#include <stdio.h>
#include <math.h>		/* For rint */

#include "assert.h"
#include "intlist.h"
#include "univcoord.h"
#include "list.h"
#include "genomebits_count.h"

#include "altsplice.h"
#include "junction.h"


static Genomebits_T genomebits;
static Genomebits_T genomebits_alt;


#define MIN_END 4

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Path_trim_circular */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


/* TODO: Handle alt splices that go outside of chromosome bounds */
/* TODO: Trim fusion paths */

#define T Path_T

static int
count_outside_qstart_byposition (int *nsegments, int *remainder,
				 Intlist_T endpoints, Univcoordlist_T univdiagonals,
				 Univcoord_T trimbounds, int querylength) {
  int noutside = 0;
  Univcoord_T left, univdiagonal;
  int qstart, qend;

  *nsegments = 0;
  *remainder = 0;

  debug(printf("Entered count_outside_qstart_byposition\n"));
  while (univdiagonals != NULL) {
    univdiagonal = Univcoordlist_head(univdiagonals);
    left = univdiagonal - querylength;
    qstart = Intlist_head(endpoints);
    qend = Intlist_second_value(endpoints);
    debug(printf("univdiagonal %u, left %u (%u..%u) vs trimbounds %u\n",
		 univdiagonal,left,left+qstart,left+qend,trimbounds));

    if (left + qstart >= trimbounds) {
      *remainder = 0;
      return noutside;

    } else if (left + qend >= trimbounds) {
      *remainder = (int) (trimbounds - (left + qstart));
      return noutside + (*remainder);

    } else {
      *nsegments += 1;
      noutside += qend - qstart;
      endpoints = Intlist_next(endpoints);
      univdiagonals = Univcoordlist_next(univdiagonals);
    }
  }

  return noutside;
}


static int
count_outside_qstart_bydiagonal (Univcoordlist_T univdiagonals, Univcoord_T trimdiag) {
  int nsegments = 0;
  Univcoord_T univdiagonal;

  debug(printf("Entered count_outside_qstart_bydiagonal\n"));
  while (univdiagonals != NULL) {
    univdiagonal = Univcoordlist_head(univdiagonals);
    /* Allow the same univdiagonal */
    if (univdiagonal >= trimdiag) {
      return nsegments;

    } else {
      nsegments += 1;
      univdiagonals = Univcoordlist_next(univdiagonals);
    }
  }

  return nsegments;
}

static void
trim_qstart (bool *evalp, int nsegments, int remainder, T this,
	     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool) {
  int i;
  int qstart, qend;
  int ignore_int;
  Univcoord_T ignore_univcoord;
  Junction_T junction;
  bool trimmedp = false;
    
  debug(printf("Path needs trimming at qstart (%d segments, %d remainder)\n",
	       nsegments,remainder));
  debug(Path_print(this));

  *evalp = false;
  for (i = 0; i < nsegments; i++) {
    this->endpoints = Intlistpool_pop(this->endpoints,intlistpool,&ignore_int
				      intlistpool_trace(__FILE__,__LINE__));
    this->nmismatches = Intlistpool_pop(this->nmismatches,intlistpool,&ignore_int
					intlistpool_trace(__FILE__,__LINE__));
    this->ref_nmismatches = Intlistpool_pop(this->ref_nmismatches,intlistpool,&ignore_int
					    intlistpool_trace(__FILE__,__LINE__));
    this->univdiagonals = Univcoordlistpool_pop(this->univdiagonals,univcoordlistpool,&ignore_univcoord
						univcoordlistpool_trace(__FILE__,__LINE__));
    this->junctions = Listpool_pop(this->junctions,listpool,(void **) &junction
				   listpool_trace(__FILE__,__LINE__));
    Junction_free(&junction,pathpool);
    trimmedp = true;
  }

  if (remainder == 0) {
    /* Keep endpoint and nmismatches */
  } else {
    Intlist_head_set(this->endpoints,Intlist_head(this->endpoints) + remainder);
    if (Intlist_head(this->nmismatches) == 0) {
      /* Keep nmismatches */
    } else {
      Intlist_head_set(this->nmismatches,-1);
      Intlist_head_set(this->ref_nmismatches,-1);
      *evalp = true;
    }
    trimmedp = true;
  }
    
  /* See if end segment is too short near the chromosome bound */
  qstart = Intlist_head(this->endpoints);
  qend = Intlist_second_value(this->endpoints);
  if (qend - qstart < MIN_END && this->junctions != NULL) {
    /* If so, then trim one more segment */
    debug(printf("Trimming one short segment at qstart\n"));
    this->endpoints = Intlistpool_pop(this->endpoints,intlistpool,&ignore_int
				      intlistpool_trace(__FILE__,__LINE__));
    this->nmismatches = Intlistpool_pop(this->nmismatches,intlistpool,&ignore_int
					intlistpool_trace(__FILE__,__LINE__));
    this->ref_nmismatches = Intlistpool_pop(this->ref_nmismatches,intlistpool,&ignore_int
					    intlistpool_trace(__FILE__,__LINE__));
    this->univdiagonals = Univcoordlistpool_pop(this->univdiagonals,univcoordlistpool,&ignore_univcoord
						univcoordlistpool_trace(__FILE__,__LINE__));
    this->junctions = Listpool_pop(this->junctions,listpool,(void **) &junction
				   listpool_trace(__FILE__,__LINE__));
    Junction_free(&junction,pathpool);
    trimmedp = true;
      
    /* Re-align to the chromosome start */
    Intlist_head_set(this->endpoints,qstart);
    Intlist_head_set(this->nmismatches,-1);
    Intlist_head_set(this->ref_nmismatches,-1);
    *evalp = true;
  }

  if (this->qstart_alts != NULL && trimmedp == true) {
    Altsplice_free(&this->qstart_alts,pathpool);
    this->qstart_alts = (Altsplice_T) NULL;
    this->splice5p = false;
    this->splicetype5 = NO_SPLICE;
    this->ambig_prob_5 = 0.0;
  }
      
  assert(Univcoordlist_length(this->univdiagonals) == Intlist_length(this->endpoints) - 1);
  assert(Intlist_length(this->nmismatches) == Intlist_length(this->endpoints) - 1);
  assert(Intlist_length(this->ref_nmismatches) == Intlist_length(this->endpoints) - 1);
  assert(List_length(this->junctions) == Intlist_length(this->endpoints) - 2);

  return;
}


void
Path_trim_qstart_n (int noutside, T this,
		    Compress_T query_compress_fwd, Compress_T query_compress_rev,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Pathpool_T pathpool) {
  int nsegments = 0, remainder;
  Intlist_T p;
  int overall_qstart;
  int ignore_found_score = this->querylength;
  bool evalp;
  
  debug(printf("Entering Path_trim_qstart_n with noutside %d against %s\n",
	       noutside,Intlist_to_string(this->endpoints)));
  debug(Path_print(this));

  p = this->endpoints;
  overall_qstart = Intlist_head(p);
  while (Intlist_next(p) != NULL && Intlist_second_value(p) - overall_qstart <= noutside) {
    nsegments++;
    p = Intlist_next(p);
  }
  
  debug(printf("Now endpoints is %s\n",Intlist_to_string(p)));
  /* Previously checked qstart = Intlist_head(p)) - overall_qstart >= noutside */
  remainder = noutside - (Intlist_head(p) - overall_qstart);
  debug(printf("Computed %d segments, remainder %d\n",nsegments,remainder));

  /* Call even if nsegments == 0 and remainder == 0 */
  trim_qstart(&evalp,nsegments,remainder,this,
	      intlistpool,univcoordlistpool,listpool,pathpool);

  if (evalp == true) {
    Path_eval_nmatches(&ignore_found_score,this,query_compress_fwd,query_compress_rev);
  }

  /* Used to always call Path_eval_nmatches */
  debug(printf("After trimming: "));
  debug(Path_print(this));
  debug(printf("\n"));

  return;
}


#if 0
/* Will need to call Path_eval_nmatches again after this, so don't
   worry about evalp */
/* Replaced by Path_trim_chrbounds */
void
Path_trim_qstart_trimbounds (T this,
			     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			     Listpool_T listpool, Pathpool_T pathpool, Univcoord_T trimbounds) {
  int noutside, nsegments;
  int remainder;
  bool evalp;

  debug(printf("Path_trim_qstart_trimbounds with trimbounds %u\n",trimbounds));
  debug(Path_print(this));

  if ((noutside = count_outside_qstart_byposition(&nsegments,&remainder,
						  this->endpoints,this->univdiagonals,
						  trimbounds,this->querylength)) == 0) {
    /* Skip */
    debug(printf("No trimming needed\n"));

  } else if (nsegments >= Univcoordlist_length(this->univdiagonals)) {
    fprintf(stderr,"Attempt to trim entire path\n");
    abort();

  } else {
    /* Call even if nsegments == 0 and remainder == 0 */
    trim_qstart(&evalp,nsegments,remainder,this,
		intlistpool,univcoordlistpool,listpool,pathpool);
  }

  return;
}
#endif


/* Returns false is result is not valid, usually because of the mate
   being inside of the path */
bool
Path_trim_qstart_trimdiag (T this,
			   Compress_T query_compress_fwd, Compress_T query_compress_rev,
			   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Pathpool_T pathpool, Univcoord_T trimdiag) {
  int nsegments;
  int ignore_found_score = this->querylength;
  bool evalp;

  debug(printf("Entering Path_trim_qstart_trimdiag with trimdiag %u\n",trimdiag));
  debug(Path_print(this));

  if ((nsegments = count_outside_qstart_bydiagonal(this->univdiagonals,trimdiag)) == 0) {
    /* Skip */
    debug(printf("No trimming needed\n"));
    return true;

  } else if (nsegments >= Univcoordlist_length(this->univdiagonals)) {
    /* fprintf(stderr,"Attempt to trim entire path\n"); */
    return false;

  } else {
    /* Call even if nsegments == 0 and remainder == 0 */
    trim_qstart(&evalp,nsegments,/*remainder*/0,this,
		intlistpool,univcoordlistpool,listpool,pathpool);
    if (evalp == true) {
      Path_eval_nmatches(&ignore_found_score,this,query_compress_fwd,query_compress_rev);
    }
    return true;
  }
}


/* Expecting path to be reversed */
static int
count_outside_qend_byposition (int *nsegments, int *remainder,
			       Intlist_T endpoints, Univcoordlist_T univdiagonals,
			       Univcoord_T trimbounds, int querylength) {
  int noutside = 0;
  Univcoord_T left, univdiagonal;
  int qstart, qend;

  *nsegments = 0;
  *remainder = 0;

  debug(printf("Entered count_outside_qend_byposition\n"));
  while (univdiagonals != NULL) {
    univdiagonal = Univcoordlist_head(univdiagonals);
    left = univdiagonal - querylength;
    qend = Intlist_head(endpoints);
    qstart = Intlist_second_value(endpoints);
    debug(printf("univdiagonal %u, left %u (%u..%u) vs trimbounds %u\n",
		 univdiagonal,left,left+qstart,left+qend,trimbounds));

    if (left + qend < trimbounds) {
      *remainder = 0;
      return noutside;

    } else if (left + qstart < trimbounds) {
      *remainder = (int) ((left + qend) - trimbounds);
      return noutside + (*remainder);

    } else {
      *nsegments += 1;
      noutside += qend - qstart;
      endpoints = Intlist_next(endpoints);
      univdiagonals = Univcoordlist_next(univdiagonals);
    }
  }

  return noutside;
}


static int
count_outside_qend_bydiagonal (Univcoordlist_T univdiagonals, int querylength,
			       Univcoord_T trimdiag) {
  int nsegments = 0;
  Univcoord_T univdiagonal;

  debug(printf("Entered count_outside_qend_bydiagonal\n"));
  while (univdiagonals != NULL) {
    univdiagonal = Univcoordlist_head(univdiagonals);

    /* Allow the same univdiagonal */
    if (univdiagonal - querylength <= trimdiag) {
      return nsegments;

    } else {
      nsegments += 1;
      univdiagonals = Univcoordlist_next(univdiagonals);
    }
  }

  return nsegments;
}


/* Assumes path has been reversed */
static void
trim_qend (bool *evalp, int nsegments, int remainder, T this,
	   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	   Listpool_T listpool, Pathpool_T pathpool) {
  int i;
  int qstart, qend;
  int ignore_int = 0;
  Univcoord_T ignore_univcoord;
  Junction_T junction;
  bool trimmedp = false;

  debug(printf("Path needs trimming at qend (%d segments, %d remainder)\n",
	       nsegments,remainder));
  debug(Path_print(this));

  *evalp = false;
  for (i = 0; i < nsegments; i++) {
    this->endpoints = Intlistpool_pop(this->endpoints,intlistpool,&ignore_int
				      intlistpool_trace(__FILE__,__LINE__));
    this->nmismatches = Intlistpool_pop(this->nmismatches,intlistpool,&ignore_int
					intlistpool_trace(__FILE__,__LINE__));
    this->ref_nmismatches = Intlistpool_pop(this->ref_nmismatches,intlistpool,&ignore_int
					    intlistpool_trace(__FILE__,__LINE__));
    this->univdiagonals = Univcoordlistpool_pop(this->univdiagonals,univcoordlistpool,&ignore_univcoord
						univcoordlistpool_trace(__FILE__,__LINE__));
    this->junctions = Listpool_pop(this->junctions,listpool,(void **) &junction
				   listpool_trace(__FILE__,__LINE__));
    Junction_free(&junction,pathpool);
    trimmedp = true;
  }

  if (remainder == 0) {
    /* Keep endpoint and nmismatches */
  } else {
    Intlist_head_set(this->endpoints,Intlist_head(this->endpoints) - remainder);
    if (Intlist_head(this->nmismatches) == 0) {
      /* Keep nmismatches */
    } else {
      Intlist_head_set(this->nmismatches,-1);
      Intlist_head_set(this->ref_nmismatches,-1);
      *evalp = true;
    }
    trimmedp = true;
  }
    
  /* See if end segment is too short near the chromosome bound */
  qend = Intlist_head(this->endpoints);
  qstart = Intlist_second_value(this->endpoints);
  if (qend - qstart < MIN_END && this->junctions != NULL) {
    /* If so, then trim one more segment */
    debug(printf("Trimming one short segment at qend\n"));
    this->endpoints = Intlistpool_pop(this->endpoints,intlistpool,&ignore_int
				      intlistpool_trace(__FILE__,__LINE__));
    this->nmismatches = Intlistpool_pop(this->nmismatches,intlistpool,&ignore_int
					intlistpool_trace(__FILE__,__LINE__));
    this->ref_nmismatches = Intlistpool_pop(this->ref_nmismatches,intlistpool,&ignore_int
					    intlistpool_trace(__FILE__,__LINE__));
    this->univdiagonals = Univcoordlistpool_pop(this->univdiagonals,univcoordlistpool,&ignore_univcoord
						univcoordlistpool_trace(__FILE__,__LINE__));
    this->junctions = Listpool_pop(this->junctions,listpool,(void **) &junction
				   listpool_trace(__FILE__,__LINE__));
    Junction_free(&junction,pathpool);
    trimmedp = true;
      
    /* Re-align to the chromosome end */
    Intlist_head_set(this->endpoints,qend);
    Intlist_head_set(this->nmismatches,-1);
    Intlist_head_set(this->ref_nmismatches,-1);
    *evalp = true;
  }

  if (this->qend_alts != NULL && trimmedp == true) {
    Altsplice_free(&this->qend_alts,pathpool);
    this->qend_alts = (Altsplice_T) NULL;
    this->splice3p = false;
    this->splicetype3 = NO_SPLICE;
    this->ambig_prob_3 = 0.0;
  }
    
  assert(Univcoordlist_length(this->univdiagonals) == Intlist_length(this->endpoints) - 1);
  assert(Intlist_length(this->nmismatches) == Intlist_length(this->endpoints) - 1);
  assert(Intlist_length(this->ref_nmismatches) == Intlist_length(this->endpoints) - 1);
  assert(List_length(this->junctions) == Intlist_length(this->endpoints) - 2);

  return;
}


void
Path_trim_qend_n (int noutside, T this,
		  Compress_T query_compress_fwd, Compress_T query_compress_rev,
		  Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		  Listpool_T listpool, Pathpool_T pathpool) {
  int nsegments = 0, remainder;
  Intlist_T p;
  int overall_qend;
  int ignore_found_score = this->querylength;
  bool evalp;

  debug(printf("Entering Path_trim_qend_n with noutside %d against %s\n",
	       noutside,Intlist_to_string(this->endpoints)));
  debug(Path_print(this));

  Path_reverse(this,/*expect_fwd_p*/false);

  p = this->endpoints;
  overall_qend = Intlist_head(p);
  while (Intlist_next(p) != NULL && overall_qend - Intlist_second_value(p) <= noutside) {
    nsegments++;
    p = Intlist_next(p);
  }

  debug(printf("Now endpoints is %s\n",Intlist_to_string(p)));
  /* Previously checked overall_qend - (qend = Intlist_head(p)) >= noutside */
  remainder = noutside - (overall_qend - Intlist_head(p));
  debug(printf("Computed %d segments, remainder %d\n",nsegments,remainder));

  /* Call even if nsegments == 0 and remainder == 0 */
  trim_qend(&evalp,nsegments,remainder,this,
	    intlistpool,univcoordlistpool,listpool,pathpool);
  Path_reverse(this,/*expect_fwd_p*/true);
  if (evalp == true) {
    Path_eval_nmatches(&ignore_found_score,this,query_compress_fwd,query_compress_rev);
  }

  /* Used to always call Path_eval_nmatches */

  debug(printf("After trimming: "));
  debug(Path_print(this));
  debug(printf("\n"));

  return;
}



#if 0
/* Replaced by Path_trim_chrbounds */
/* Will need to call Path_eval_nmatches again after this, so don't
   worry about evalp */
void
Path_trim_qend_trimbounds (T this, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Pathpool_T pathpool, Univcoord_T trimbounds) {
  int noutside, nsegments;
  int remainder;
  bool evalp;

  debug(printf("Path_trim_qend_trimbounds with trimbounds %u\n",trimbounds));
  debug(Path_print(this));

  Path_reverse(this,/*expect_fwd_p*/false);

  if ((noutside = count_outside_qend_byposition(&nsegments,&remainder,
						this->endpoints,this->univdiagonals,
						trimbounds,this->querylength)) == 0) {
    /* Skip */
    debug(printf("No trimming needed\n"));

  } else if (nsegments >= Univcoordlist_length(this->univdiagonals)) {
    fprintf(stderr,"Attempt to trim entire path\n");
    abort();

  } else {
    /* Call even if nsegments == 0 and remainder == 0 */
    trim_qend(&evalp,nsegments,remainder,this,
	      intlistpool,univcoordlistpool,listpool,pathpool);
  }

  Path_reverse(this,/*expect_fwd_p*/true);

  return;
}
#endif


/* Returns false is result is not valid, usually because of the mate
   being inside of the path */
bool
Path_trim_qend_trimdiag (T this,
			 Compress_T query_compress_fwd, Compress_T query_compress_rev,
			 Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			 Listpool_T listpool, Pathpool_T pathpool, int querylength,
			 Univcoord_T trimdiag) {
  int nsegments;
  int ignore_found_score = this->querylength;
  bool evalp;

  debug(printf("Entering Path_trim_qend_trimdiag with trimdiag %u\n",trimdiag));
  debug(Path_print(this));

  Path_reverse(this,/*expect_fwd_p*/false);

  if ((nsegments = count_outside_qend_bydiagonal(this->univdiagonals,querylength,trimdiag)) == 0) {
    /* Skip */
    debug(printf("No trimming needed\n"));
    Path_reverse(this,/*expect_fwd_p*/true);
    return true;

  } else if (nsegments >= Univcoordlist_length(this->univdiagonals)) {
    /* fprintf(stderr,"Attempt to trim entire path\n"); */
    Path_reverse(this,/*expect_fwd_p*/true);
    return false;

  } else {
    /* Call even if nsegments == 0 and remainder == 0 */
    trim_qend(&evalp,nsegments,/*remainder*/0,this,
	      intlistpool,univcoordlistpool,listpool,pathpool);
    Path_reverse(this,/*expect_fwd_p*/true);
    if (evalp == true) {
      Path_eval_nmatches(&ignore_found_score,this,
			 query_compress_fwd,query_compress_rev);
    }
    return true;
  }
}


void
Path_trim_chrbounds (T this,
		     Compress_T query_compress_fwd, Compress_T query_compress_rev,
		     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		     Listpool_T listpool, Pathpool_T pathpool) {
  int nsegments, remainder, qstart, qend;
  Univcoord_T univdiagonal;
  Univcoordlist_T p;
  int ignore_found_score = this->querylength;
  bool qstart_evalp, qend_evalp;
  
  debug(printf("Entering Path_trim_chrbounds with chroffset %u and chrhigh %u\n",
	       this->chroffset,this->chrhigh));
  debug(Path_print(this));

  /* qstart */
  if (this->qstart_alts != NULL &&
      Altsplice_trim_qstart_chrbounds(this->qstart_alts,this->chroffset) == false) {
    Altsplice_free(&this->qstart_alts,pathpool);
    this->qstart_alts = (Altsplice_T) NULL;
    this->splice5p = false;
    this->splicetype5 = NO_SPLICE;
    this->ambig_prob_5 = 0.0;
  }

  nsegments = 0;
  p = this->univdiagonals;
  while (p != NULL && Univcoordlist_head(p) < this->chroffset) {
    nsegments++;
    p = Univcoordlist_next(p);
  }
  
  univdiagonal = Univcoordlist_head(p);
  qstart = Intlist_head(this->endpoints);
  debug(printf("Now univdiagonals is %s, qstart %d\n",Univcoordlist_to_string(p),qstart));
  /* remainder = chroffset - (univdiagonal - querylength + endpoint) */
  if (univdiagonal + qstart >= this->chroffset + this->querylength) {
    remainder = 0;
  } else {
    remainder = this->chroffset + this->querylength - univdiagonal - qstart;
  }
  debug(printf("Computed %d segments, remainder %d\n",nsegments,remainder));
  /* Call even if nsegments == 0 and remainder == 0 */
  trim_qstart(&qstart_evalp,nsegments,remainder,this,
	      intlistpool,univcoordlistpool,listpool,pathpool);


  /* qend */
  if (this->qend_alts != NULL &&
      Altsplice_trim_qend_chrbounds(this->qend_alts,this->chrhigh) == false) {
    debug(printf("Removing qend_alts\n"));
    Altsplice_free(&this->qend_alts,pathpool);
    this->qend_alts = (Altsplice_T) NULL;
    this->splice3p = false;
    this->splicetype3 = NO_SPLICE;
    this->ambig_prob_3 = 0.0;
  }

  Path_reverse(this,/*expect_fwd_p*/false);

  nsegments = 0;
  p = this->univdiagonals;
  while (p != NULL && Univcoordlist_head(p) >= this->chrhigh + this->querylength) {
    nsegments++;
    p = Univcoordlist_next(p);
  }

  univdiagonal = Univcoordlist_head(p);
  qend = Intlist_head(this->endpoints);
  debug(printf("Now univdiagonals is %s, qend %d\n",Univcoordlist_to_string(p),qend));
  if (univdiagonal + qend < this->chrhigh + this->querylength) {
    remainder = 0;
  } else {
    remainder = univdiagonal + qend - this->querylength - this->chrhigh;
  }
  debug(printf("Computed %d segments, remainder %d\n",nsegments,remainder));
  /* Call even if nsegments == 0 and remainder == 0 */
  trim_qend(&qend_evalp,nsegments,remainder,this,
	    intlistpool,univcoordlistpool,listpool,pathpool);

  Path_reverse(this,/*expect_fwd_p*/true);

  if (qstart_evalp == true || qend_evalp == true) {
    Path_eval_nmatches(&ignore_found_score,this,query_compress_fwd,query_compress_rev);
  }
  debug(printf("After trimming: "));
  debug(Path_print(this));
  debug(printf("\n"));

  return;
}



static int
count_total (Intlist_T endpoints, Univcoordlist_T univdiagonals) {
  int ntotal = 0;
  int qstart, qend;

  while (univdiagonals != NULL) {
    qstart = Intlist_head(endpoints);
    qend = Intlist_second_value(endpoints);

    ntotal += qend - qstart;
    endpoints = Intlist_next(endpoints);
    univdiagonals = Univcoordlist_next(univdiagonals);
  }

  return ntotal;
}


static int
divide_path (int *nhigh, int *qend_nsegments, int *qend_remainder,
	     Univcoord_T *linear_chrhigh, Univcoord_T *linear_chrlength, T this) {
  int ntotal;

  /* Determine the point of division */
  *linear_chrlength = (this->chrhigh - this->chroffset)/2;
  *linear_chrhigh = this->chroffset + (*linear_chrlength);
  debug2(printf("divide_path\n"));
  debug2(Path_print(this));
  debug2(printf("chrhigh %u, chroffset %u, chrlength %u\n",
		this->chrhigh,this->chroffset,*linear_chrlength));

  ntotal = count_total(this->endpoints,this->univdiagonals);

  Path_reverse(this,/*expect_fwd_p*/false);
  *nhigh = count_outside_qend_byposition(&(*qend_nsegments),&(*qend_remainder),
					 this->endpoints,this->univdiagonals,
					 /*trimbounds*/*linear_chrhigh,this->querylength);
  Path_reverse(this,/*expect_fwd_p*/true);
  debug2(printf("nlow %d, nhigh %d\n",ntotal - (*nhigh),*nhigh));

  return (ntotal - (*nhigh));
}


static void
unalias_path (T this, Univcoord_T linear_chrhigh, Univcoord_T linear_chrlength) {
  Univcoordlist_T u;
  List_T j;
  Altsplice_T altsplice;
  Junction_T junction;
  int i;

  if (this->main_univdiagonal < linear_chrlength) {
    this->main_univdiagonal = 0;
  } else {
    this->main_univdiagonal -= linear_chrlength;
  }

  for (u = this->univdiagonals; u != NULL; u = Univcoordlist_next(u)) {
    Univcoordlist_head_set(u,Univcoordlist_head(u) - linear_chrlength);
  }
  for (j = this->junctions; j != NULL; j = List_next(j)) {
    junction = (Junction_T) List_head(this->junctions);
    if (Junction_type(junction) == DEL_JUNCTION) {
      Junction_set_deletionpos(junction,
			       Junction_deletionpos(junction) - linear_chrlength);
    }
  }
  
  /* It is possible that some qstart alts_coords are above linear_chrhigh and some below */
  /* However, if we don't unalias all of them, we could get a negative splice distance */
  if ((altsplice = this->qstart_alts) != NULL) {
    for (i = 0; i < altsplice->ncoords; i++) {
      /* if (altsplice->coords[i] >= linear_chrhigh) { */
      altsplice->coords[i] -= linear_chrlength;
      /* } */
    }
  }

  if ((altsplice = this->qend_alts) != NULL) {
    for (i = 0; i < altsplice->ncoords; i++) {
      altsplice->coords[i] -= linear_chrlength;
    }
  }

  return;
}


/* Returns true if the entire path is in the upper half, meaning that
   it must have an alias in the lower half */
void
Path_trim_circular_unalias (T this) {
  int nlow, nhigh, qend_nsegments, qend_remainder;
  Univcoord_T linear_chrhigh, linear_chrlength;

  if ((nlow = divide_path(&nhigh,&qend_nsegments,&qend_remainder,
			  &linear_chrhigh,&linear_chrlength,this)) > 0) {
    /* Do not alter the path */
  } else {
    unalias_path(this,linear_chrhigh,linear_chrlength);
  }

  return;
}


void
Path_trim_circular_unalias_pair (T pathL, T pathH) {
  int nlow, nhigh, qend_nsegments, qend_remainder;
  Univcoord_T linear_chrhigh, linear_chrlength;

  if ((nlow = divide_path(&nhigh,&qend_nsegments,&qend_remainder,
			  &linear_chrhigh,&linear_chrlength,pathL)) > 0) {
    /* Do not alter either path */
  } else {
    /* Alter both paths */
    unalias_path(pathL,linear_chrhigh,linear_chrlength);
    unalias_path(pathH,linear_chrhigh,linear_chrlength);
  }

  return;
}


/* For circular alignments: (1) if the entire path is in the lower
   path, we keep it; (2) if the entire path is in the upper half,
   then we unalias it.  Otherwise, if it straddles the circular
   origin, we do not touch the main endpoints and univdiagonals, but
   create circlow and circhigh endpoints */

void
Path_trim_circular (T this, Compress_T query_compress,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool) {
  int nlow, nhigh, qend_nsegments, qend_remainder, i;
  Univcoord_T linear_chrhigh, linear_chrlength, univdiagonal;
  int nmismatches, ref_nmismatches;

  Intlist_T high_endpoints = NULL, high_nmismatches = NULL, high_ref_nmismatches = NULL;
  Univcoordlist_T high_univdiagonals = NULL;
  List_T high_junctions = NULL;
  Junction_T junction;
  Altsplice_T altsplice;

  debug2(printf("Entering Path_trim_circular\n"));

  if ((nlow = divide_path(&nhigh,&qend_nsegments,&qend_remainder,
			  &linear_chrhigh,&linear_chrlength,this)) == 0) {
    /* Entire alignment is in the top half of the double chromosome, so unalias all univdiagonals */
    debug2(printf("Entire alignment is in the top half\n"));
    Path_trim_circular_unalias(this);

  } else if (nhigh == 0) {
    /* Entire alignment is in the bottom half of the double chromosome, so nothing to do */
    debug2(printf("Entire alignment is in the lower half\n"));

  } else {
    /* Divide up the components */
    debug2(printf("Dividing alignment: qend_nsegments %d, qend_remainder %d\n",
		  qend_nsegments,qend_remainder));
    Path_reverse(this,/*expect_fwd_p*/false);

    for (i = 0; i < qend_nsegments; i++) {
      debug2(printf("Pushing %d onto high_endpoints\n",Intlist_head(this->endpoints)));
      high_endpoints = Intlistpool_push(high_endpoints,intlistpool,Intlist_head(this->endpoints)
					intlistpool_trace(__FILE__,__LINE__));
      high_nmismatches = Intlistpool_push(high_nmismatches,intlistpool,Intlist_head(this->nmismatches)
					  intlistpool_trace(__FILE__,__LINE__));
      high_ref_nmismatches = Intlistpool_push(high_ref_nmismatches,intlistpool,Intlist_head(this->ref_nmismatches)
					      intlistpool_trace(__FILE__,__LINE__));
      high_univdiagonals = Univcoordlistpool_push(high_univdiagonals,univcoordlistpool,
						  Univcoordlist_head(this->univdiagonals) - linear_chrlength
						  univcoordlistpool_trace(__FILE__,__LINE__));

      junction = (Junction_T) List_head(this->junctions);
      if (Junction_type(junction) == DEL_JUNCTION) {
	Junction_set_deletionpos(junction,
				 Junction_deletionpos(junction) - linear_chrlength);
      }
      high_junctions = Listpool_push(high_junctions,listpool,(void *) junction
				     listpool_trace(__FILE__,__LINE__));

      this->endpoints = Intlist_next(this->endpoints);
      this->nmismatches = Intlist_next(this->nmismatches);
      this->ref_nmismatches = Intlist_next(this->ref_nmismatches);
      this->univdiagonals = Univcoordlist_next(this->univdiagonals);
      this->junctions = List_next(this->junctions);
    }
			    
    debug2(printf("Pushing %d onto high_endpoints\n",Intlist_head(this->endpoints)));
    high_endpoints = Intlistpool_push(high_endpoints,intlistpool,Intlist_head(this->endpoints)
				      intlistpool_trace(__FILE__,__LINE__));

    high_endpoints = Intlistpool_push(high_endpoints,intlistpool,Intlist_head(this->endpoints) - qend_remainder
				      intlistpool_trace(__FILE__,__LINE__));
    Intlist_head_set(this->endpoints,Intlist_head(this->endpoints) - qend_remainder);

    univdiagonal = Univcoordlist_head(this->univdiagonals);
    if (univdiagonal < this->chroffset + linear_chrlength) {
      /* Can happen if there is splicing across the circular origin */
      high_univdiagonals = Univcoordlistpool_push(high_univdiagonals,univcoordlistpool,univdiagonal
						  univcoordlistpool_trace(__FILE__,__LINE__));
    } else {
      high_univdiagonals = Univcoordlistpool_push(high_univdiagonals,univcoordlistpool,
						  univdiagonal - linear_chrlength
						  univcoordlistpool_trace(__FILE__,__LINE__));
    }

    if (Intlist_head(this->nmismatches) == 0 && Intlist_head(this->ref_nmismatches) == 0) {
      high_nmismatches = Intlistpool_push(high_nmismatches,intlistpool,0
					  intlistpool_trace(__FILE__,__LINE__));
      high_ref_nmismatches = Intlistpool_push(high_ref_nmismatches,intlistpool,0
					      intlistpool_trace(__FILE__,__LINE__));
      
    } else {
      nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress,
							  univdiagonal,this->querylength,
							  /*pos5*/Intlist_second_value(this->endpoints),
							  /*pos3*/Intlist_head(this->endpoints),
							  this->plusp,this->genestrand);
      Intlist_head_set(this->nmismatches,nmismatches);
      Intlist_head_set(this->ref_nmismatches,ref_nmismatches);
      
      nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress,
							  univdiagonal,this->querylength,
							  /*pos5*/Intlist_head(high_endpoints),
							  /*pos3*/Intlist_second_value(high_endpoints),
							  this->plusp,this->genestrand);
      high_nmismatches = Intlistpool_push(high_nmismatches,intlistpool,nmismatches
					  intlistpool_trace(__FILE__,__LINE__));
      high_ref_nmismatches = Intlistpool_push(high_ref_nmismatches,intlistpool,nmismatches
					      intlistpool_trace(__FILE__,__LINE__));
    }
    

    /* Unalias qend_alts, which are above chrlength */
    if ((altsplice = this->qend_alts) != NULL) {
      for (i = 0; i < altsplice->ncoords; i++) {
	altsplice->coords[i] -= linear_chrlength;
      }
    }

    /* Restore low path */
    Path_reverse(this,/*expect_fwd_p*/true);

    if (nlow >= nhigh) {
      /* Make high path the circular part */
      this->circular_high_p = true;
      this->circular_endpoints = high_endpoints;
      this->circular_nmismatches = high_nmismatches;
      this->circular_ref_nmismatches = high_ref_nmismatches;
      this->circular_univdiagonals = high_univdiagonals;
      this->circular_junctions = high_junctions;

    } else {
      /* Make low path the circular part */
      this->circular_high_p = false;
      this->circular_endpoints = this->endpoints;
      this->circular_nmismatches = this->nmismatches;
      this->circular_ref_nmismatches = this->ref_nmismatches;
      this->circular_univdiagonals = this->univdiagonals;
      this->circular_junctions = this->junctions;
      
      this->endpoints = high_endpoints;
      this->nmismatches = high_nmismatches;
      this->ref_nmismatches = high_ref_nmismatches;
      this->univdiagonals = high_univdiagonals;
      this->junctions = high_junctions;
    }

    /* TODO: Fix transcripts (or print them only with the primary entry) */
  }

  debug2(printf("Result of dividing path\n"));
  debug2(Path_print(this));

  return;
}


void
Path_trim_setup (Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in) {
  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;
  return;
}
