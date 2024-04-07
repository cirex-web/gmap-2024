static char rcsid[] = "$Id: diag.c 223519 2020-11-15 02:30:52Z twu $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "diag.h"
#include "diagdef.h"
#include "mem.h"
#include <stdio.h>
#include <stdlib.h>


#define DOMINANCE_END_EQUIV 20

#define EXTRA_BOUNDS 20	/* Was 120. Number of extra diaglines to make active */

#define MAX_DIAGONALS 50
#define MIN_SCORE 10.0


#define T Diag_T

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* dominant diagonals */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* minactive, maxactive */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


#if 0
/* Oligoindex procedures use diagpool_push instead */
T
Diag_new (int querystart, int queryend, Diag_T diagonal) {
  T new = (T) MALLOC(sizeof(*new));

  new->diagonal = diagonal;
  new->querystart = querystart;
  new->queryend = queryend;
  new->nconsecutive = queryend - querystart + 1;
  new->dominatedp = false;
  new->score = 0.0;

  return new;
}
#endif


void
Diag_free (T *old) {
  FREE(*old);
  return;
}

void
Diag_gc (List_T *list) {
  T diagonal;
  List_T p;

  for (p = *list; p != NULL; p = List_next(p)) {
    diagonal = (T) List_head(p);
    FREE(diagonal);
  }
  List_free(&(*list));
  return;
}


Chrpos_T
Diag_diagonal (T this) {
  return this->diagonal;
}

int
Diag_querystart (T this) {
  return this->querystart;
}

int
Diag_queryend (T this) {
  return this->queryend;
}

int
Diag_nconsecutive (T this) {
  return this->nconsecutive;
}

bool
Diag_dominatedp (T this) {
  return this->dominatedp;
}

void
Diag_set_dominatedp (T this) {
  this->dominatedp = true;
}


int
Diag_ascending_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->querystart < y->querystart) {
    return -1;
  } else if (y->querystart < x->querystart) {
    return +1;
  } else if (x->queryend < y->queryend) {
    return -1;
  } else if (y->queryend < x->queryend) {
    return +1;
  } else if (x->diagonal < y->diagonal) {
    return -1;
  } else if (y->diagonal < x->diagonal) {
    return +1;
  } else {
    return 0;
  }
}


int
Diag_descending_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->querystart > y->querystart) {
    return -1;
  } else if (y->querystart > x->querystart) {
    return +1;
  } else if (x->queryend > y->queryend) {
    return -1;
  } else if (y->queryend > x->queryend) {
    return +1;
  } else if (x->diagonal > y->diagonal) {
    return -1;
  } else if (y->diagonal > x->diagonal) {
    return +1;
  } else {
    return 0;
  }
}



static int
nconsecutive_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  if (a->nconsecutive > b->nconsecutive) {
    return -1;
  } else if (b->nconsecutive > a->nconsecutive) {
    return 1;
  } else {
    return 0;
  }
}

static int
diagonal_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  if (a->diagonal < b->diagonal) {
    return -1;
  } else if (b->diagonal < a->diagonal) {
    return +1;
  } else {
    return 0;
  }
}


#if 0
static int
abs_compare (const void *x, const void *y) {
  int absa, absb;
  int a = * (int *) x;
  int b = * (int *) y;

  absa = (a >= 0) ? a : -a;
  absb = (b >= 0) ? b : -b;

  if (absa < absb) {
    return -1;
  } else if (absa > absb) {
    return +1;
  } else if (a >= 0 && b < 0) {
    /* positives before negatives */
    return -1;
  } else if (a < 0 && b >= 0) {
    return +1;
  } else {
    return 0;
  }
}
#endif


/* Note: support (querypos_max - querypos_min) is not a useful metric here */
double
Diag_update_coverage (bool *coveredp, int *ncovered, List_T diagonals, int querylength) {
  List_T p;
  T diag;
  int *scores, querypos, count;

  debug(printf("Entered Diag_update_coverage with %d diagonals\n",List_length(diagonals)));

  *ncovered = 0;
  if (querylength < MAX_QUERYLENGTH_STACK) {
    scores = (int *) CALLOCA(querylength,sizeof(int));
  } else {
    scores = (int *) CALLOC(querylength,sizeof(int));
  }

  for (p = diagonals; p != NULL; p = List_next(p)) {
    diag = (T) List_head(p);
    debug(printf("diagonal at %u from querypos %d to %d\n",diag->diagonal,diag->querystart,diag->queryend));
    scores[diag->querystart] += 1;
    scores[diag->queryend] -= 1;
  }

  count = 0;
  for (querypos = 0; querypos < querylength; querypos++) {
    count += scores[querypos];
    if (count > 0) {
      coveredp[querypos] = true;
    }
    (*ncovered) += (int) coveredp[querypos];
  }

  if (querylength < MAX_QUERYLENGTH_STACK) {
    FREEA(scores);
  } else {
    FREE(scores);
  }

  debug(printf("coverage = %f\n\n",(double) (*ncovered)/(double) querylength));

  return (double) (*ncovered)/(double) querylength;
}


#if 0
static int
diagonal_coverage (int *clear_coverage, T *array, int nunique) {
  int coverage = 0, regionstart, clear_regionstart, i, j = 0, status;
  int *events, nevents;
  T diag;
  
  *clear_coverage = 0;
  nevents = nunique+nunique;

  events = (int *) MALLOCA(nevents * sizeof(int));
  for (i = 0; i < nunique; i++) {
    diag = array[i];
    events[j++] = +diag->querystart;
    events[j++] = -diag->queryend;
  }
  qsort(events,nevents,sizeof(int),abs_compare);

  status = 0;
  for (j = 0; j < nevents; j++) {
    if (events[j] >= 0) {
      status++;
      if (status == 1) {
	/* Start of region.  Now in state 1..MAX_DIAGONALS */
	clear_regionstart = events[j];
	regionstart = events[j];
      } else if (status == MAX_DIAGONALS + 1) {
	/* End of region.  Now in state MAX_DIAGONALS+1 .. Inf */
	*clear_coverage += (events[j]) - clear_regionstart + 1;
      }
    } else {
      status--;
      if (status == MAX_DIAGONALS) {
	/* Start of region.  Now in state 1..MAX_DIAGONALS */
	clear_regionstart = -events[j];
      } else if (status == 0) {
	/* End of region.  Now in state 0 */
	coverage += (-events[j]) - regionstart + 1;
	*clear_coverage += (-events[j]) - clear_regionstart + 1;
      }
    }
  }

  FREEA(events);
  
  return coverage;
}
#endif


int
_querystart (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  if (a->querystart < b->querystart) {
    return -1;
  } else if (b->querystart < a->querystart) {
    return +1;
  } else {
    return 0;
  }
}

static void
print_segment (T this, Chrpos_T chrinit, char *queryseq_ptr, char *genomicseg_ptr) {
  Chrpos_T genomicstart, genomicend, genomicpos;
  int querypos;

#ifdef PMAP
  genomicstart = chrinit + this->diagonal + 3*this->querystart;
  genomicend = chrinit + this->diagonal + 3*this->queryend;
#else
  genomicstart = chrinit + this->diagonal + this->querystart;
  genomicend = chrinit + this->diagonal + this->queryend;
#endif  

  printf(">%d..%d %u..%u diagonal:%u score:%f\n",
	 this->querystart,this->queryend,genomicstart,genomicend,this->diagonal,this->score);
  if (queryseq_ptr != NULL) {
    for (querypos = this->querystart; querypos <= this->queryend; querypos++) {
      printf("%c",queryseq_ptr[querypos]);
    }
    printf("\n");
  }

  if (genomicseg_ptr != NULL) {
    for (genomicpos = genomicstart; genomicpos <= genomicend; genomicpos++) {
      printf("%c",genomicseg_ptr[genomicpos]);
    }
    printf("\n");
  }

  return;
}

/* Called only as extern */
void
Diag_print_segments (List_T diagonals, char *queryseq_ptr, char *genomicseg_ptr) {
  T *array;
  int i, n;

  if ((n = List_length(diagonals)) > 0) {
    array = (T *) MALLOCA(n * sizeof(T));
    List_fill_array((void **) array,diagonals);
    qsort(array,List_length(diagonals),sizeof(T),Diag_ascending_cmp);
    for (i = 0; i < List_length(diagonals); i++) {
      print_segment(array[i],/*chrinit*/0U,queryseq_ptr,genomicseg_ptr);
    }
    FREEA(array);
  }
  return;
}


#if 0
static void
print_segments_for_R_list (List_T list, Chrpos_T chrinit, char *color) {
  List_T p;
  T diag;

  for (p = list; p != NULL; p = List_next(p)) {
    diag = (T) List_head(p);
#ifdef PMAP
    printf("segments(%d,%d+%d,%d,%d+%d,col=\"%s\")  # nconsecutive = %d, score = %.2f\n",
	   diag->querystart,chrinit + diag->diagonal,3*diag->querystart,
	   diag->queryend,chrinit + diag->diagonal,3*diag->queryend,
	   color,diag->nconsecutive,diag->score);
#else
    printf("segments(%d,%d+%d,%d,%d+%d,col=\"%s\")  # nconsecutive = %d, score = %.2f\n",
	   diag->querystart,chrinit + diag->diagonal,diag->querystart,
	   diag->queryend,chrinit + diag->diagonal,diag->queryend,
	   color,diag->nconsecutive,diag->score);
#endif
  }

  return;
}
#endif


#if 0
static void
print_segments_for_R_array (T *array, int n, Chrpos_T chrinit, char *color) {
  int i;
  T diag;

  for (i = 0; i < n; i++) {
    diag = array[i];
#ifdef PMAP
    printf("segments(%d,%d+%d,%d,%d+%d,col=\"%s\")  # nconsecutive = %d, score = %.2f\n",
	   diag->querystart,chrinit + diag->diagonal,3*diag->querystart,
	   diag->queryend,chrinit + diag->diagonal,3*diag->queryend,
	   color,diag->nconsecutive,diag->score);
#else
    printf("segments(%d,%d+%d,%d,%d+%d,col=\"%s\")  # nconsecutive = %d, score = %.2f\n",
	   diag->querystart,chrinit + diag->diagonal,diag->querystart,
	   diag->queryend,chrinit + diag->diagonal,diag->queryend,
	   color,diag->nconsecutive,diag->score);
#endif
  }

  return;
}
#endif


/* Performs sort internally */
static T *
compute_dominance (int *nunique, T *array, int ndiagonals) {
  int superstart, superend, expected_nconsecutive, threshold;
  int i, j, k;
  T super, sub;

  qsort(array,ndiagonals,sizeof(T),nconsecutive_cmp);

  *nunique = ndiagonals;
  i = 0;
  while (i < *nunique) {
    super = array[i];
    superstart = super->querystart;
    superend = super->queryend;

    expected_nconsecutive = superend + 1 - superstart;
    debug(printf("expected_nconsecutive is %d.  actual is %d\n",expected_nconsecutive,super->nconsecutive));
    if (expected_nconsecutive < 100 && super->nconsecutive > expected_nconsecutive - 10) {
      /* Short diagonal.  Rely on arithmetic differences */
      threshold = super->nconsecutive - DOMINANCE_END_EQUIV;
      for (j = i+1; j < *nunique; j++) {
	sub = array[j];
	if (sub->querystart >= superstart && sub->queryend <= superend && sub->nconsecutive < threshold) {
	  sub->dominatedp = true;
	}
      }

      /* Shift array to contain non-dominated diagonals */
      k = i+1;
      for (j = i+1; j < *nunique; j++) {
	sub = array[j];
	if (sub->dominatedp == false) {
	  array[k++] = array[j];
	}
      }
      *nunique = k;

    } else if (expected_nconsecutive >= 100 && super->nconsecutive > expected_nconsecutive * 0.90) {
      /* Long diagonal.  Rely on factor differences */
      /* Needed for very long diagonals, such as in genome-genome alignment */
      threshold = super->nconsecutive * 0.80;
      for (j = i+1; j < *nunique; j++) {
	sub = array[j];
	if (sub->querystart >= superstart && sub->queryend <= superend && sub->nconsecutive < threshold) {
	  sub->dominatedp = true;
	}
      }

      /* Shift array to contain non-dominated diagonals */
      k = i+1;
      for (j = i+1; j < *nunique; j++) {
	sub = array[j];
	if (sub->dominatedp == false) {
	  array[k++] = array[j];
	}
      }
      *nunique = k;
    }

    i++;
  }

  return array;
}
  

static T *
keep_center_diagonal (int *nunique, T *array, Chrpos_T center_diagonal) {
  int nkeep = 0;
  int i, k;


  for (i = 0; i < *nunique; i++) {
    if (array[i]->diagonal + 10000 < center_diagonal || array[i]->diagonal > center_diagonal + 10000) {
      /* Skip */
    } else {
      nkeep++;
    }
  }

  k = 0;
  for (i = 0; i < *nunique; i++) {
    if (array[i]->diagonal + 10000 < center_diagonal || array[i]->diagonal > center_diagonal + 10000) {
      /* Skip */
    } else {
      array[k++] = array[i];
    }
  }

  *nunique = nkeep;
  return array;
}


static void
assign_scores (List_T diagonals, int querylength) {
  int querypos;
  double *cumscores, cum, count;
  List_T p;
  T diag;

  if (querylength < MAX_QUERYLENGTH_STACK) {
    cumscores = (double *) CALLOCA(querylength,sizeof(double));
  } else {
    cumscores = (double *) CALLOC(querylength,sizeof(double));
  }

  /* Record endpoints: +1, 0, -1 */
  for (p = diagonals; p != NULL; p = List_next(p)) {
    diag = (T) List_head(p);
    cumscores[diag->querystart] += 1.0;
    cumscores[diag->queryend] -= 1.0;
  }

  /* Record scores: 1.0/cumsum(endpoints) */
  count = 0.0;
  for (querypos = 0; querypos < querylength; querypos++) {
    count += cumscores[querypos];
    if (count > 0.0) {
      cumscores[querypos] = 1.0/(double) count;
    } else {
      cumscores[querypos] = 0.0;
    }
    /* printf("%d count=%d score=%f\n",querypos,count,cumscores[querypos]); */
  }


  /* Record cumulative scores so we can just subtract to get overall score */
  cum = 0.0;
  for (querypos = 0; querypos < querylength; querypos++) {
    cum += cumscores[querypos];
    /* printf("%d %f %f\n",querypos,cumscores[querypos],cum); */
    cumscores[querypos] = cum;
  }

  for (p = diagonals; p != NULL; p = List_next(p)) {
    diag = (T) List_head(p);
    diag->score = cumscores[diag->queryend] - cumscores[diag->querystart];
  }

  if (querylength < MAX_QUERYLENGTH_STACK) {
    FREEA(cumscores);
  } else {
    FREE(cumscores);
  }

  return;
}

void
Diag_range (int *start, int *end, List_T diagonals, int querylength) {
  List_T p;
  T diag;

  *start = querylength;
  *end = 0;

  for (p = diagonals; p != NULL; p = List_next(p)) {
    diag = (T) List_head(p);
    if (diag->querystart < *start) {
      *start = diag->querystart;
    }
    if (diag->queryend > *end) {
      *end = diag->queryend;
    }
  }
  return;
}


int
Diag_compute_bounds (int *diag_querystart, int *diag_queryend,
		     Chrpos_T *minactive, Chrpos_T *maxactive, List_T diagonals,
		     int querylength, Chrpos_T chrstart, Chrpos_T chrend,
		     Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp) {
  int nunique, ndiagonals, ngooddiagonals, i, j;
  Chrpos_T diagonal, min_diagonal, max_diagonal, center_diagonal;
  int *bins, maxcount;
  int nbins, bini;
  Chrpos_T genomiclength, position, chrinit, chrterm;
  int activestart, activeend;
  int querypos;
  T *array, diag;
  List_T gooddiagonals = NULL, p;
  
  genomiclength = chrend - chrstart;

  if (plusp == true) {
    chrinit = chrstart;
    chrterm = chrend;
  } else {
    chrinit = (chrhigh - chroffset) - chrend;
    chrterm = (chrhigh - chroffset) - chrstart;
  }

  /* Sort diagonals */
  if ((ndiagonals = List_length(diagonals)) == 0) {
    nunique = 0;
    for (querypos = 0; querypos < querylength; querypos++) {
      minactive[querypos] = chrinit;
      maxactive[querypos] = chrterm;
    }
    *diag_querystart = 0;
    *diag_queryend = querylength - 1;

  } else {
    assign_scores(diagonals,querylength);

#if 0
    if (diagnosticp == true) {
      Diag_print_segments(diagonals,queryuc_ptr,genomicuc_ptr);
    }
#endif

    for (p = diagonals; p != NULL; p = List_next(p)) {
      diag = (T) List_head(p);
      if (diag->score >= MIN_SCORE) {
	gooddiagonals = List_push(gooddiagonals,diag);
      }
    }

    if ((ngooddiagonals = List_length(gooddiagonals)) == 0) {
      /* Not possible for List_length(diagonals) == 0) */
      array = (T *) MALLOCA(ndiagonals * sizeof(T));
      List_fill_array((void **) array,diagonals);
      array = compute_dominance(&nunique,array,ndiagonals);

    } else {
      /* Not possible for List_length(gooddiagonals) to be 0 */
      array = (T *) MALLOCA(ngooddiagonals * sizeof(T));
      List_fill_array_and_free((void **) array,&gooddiagonals);
      array = compute_dominance(&nunique,array,ngooddiagonals);
    }

    qsort(array,nunique,sizeof(T),diagonal_cmp);

#ifdef DEBUG0
    printf("Start of diagonals:\n");
    for (i = 0; i < nunique; i++) {
      print_segment(array[i],chrinit,/*queryseq_ptr*/NULL,/*genomicseg_ptr*/NULL);
    }
    printf("End of diagonals\n\n");
#endif

    if (nunique > 100) {
      min_diagonal = array[0]->diagonal;
      max_diagonal = array[nunique-1]->diagonal;
      nbins = (max_diagonal - min_diagonal)/10000 + 1;
      bins = (int *) CALLOC(nbins,sizeof(int));
      for (i = 0; i < nunique; i++) {
	bini = (array[i]->diagonal - min_diagonal)/10000;
	bins[bini] += 1;
      }
      maxcount = 0;
      /* max_bini = -1; */
      diagonal = min_diagonal;
      for (bini = 0; bini < nbins; bini++) {
	if (bins[bini] > maxcount) {
	  maxcount = bins[bini];
	  /* max_bini = bini; */
	  center_diagonal = diagonal;
	}
	diagonal += 10000;
      }
      center_diagonal += 5000;
      debug0(printf("max bin is %d with count %d => center diagonal %u\n",
		    max_bini,maxcount,center_diagonal));

      array = keep_center_diagonal(&nunique,array,center_diagonal);
      FREE(bins);

#ifdef DEBUG0
      printf("Start of center diagonals:\n");
      for (i = 0; i < nunique; i++) {
	print_segment(array[i],chrinit,/*queryseq_ptr*/NULL,/*genomicseg_ptr*/NULL);
      }
      printf("End of center diagonals\n\n");
#endif
    }

    /* Find end regions */
#ifdef ACTIVE_BUFFER
    /* Allow buffer on 5' end to make sure we identify the best initial exon */
    if ((activestart = array[0]->querystart) < ACTIVE_BUFFER) {
      activestart = ACTIVE_BUFFER;
      if (activestart > querylength) {
	activestart = querylength;
      }
    }

    /* Allow buffer on 3' end to make sure we identify the best final exon */
    if ((activeend = array[nunique-1]->queryend) > querylength-ACTIVE_BUFFER) {
      activeend = querylength-ACTIVE_BUFFER;
      if (activeend < 0) {
	activeend = 0;
      }
    }
#else
    /* For setting minactive and maxactive boundaries */
    activestart = array[0]->querystart;
    activeend = array[nunique-1]->queryend;
    debug1(printf("activestart = %d, activeend = %d\n",activestart,activeend));

    /* For query bounds for middle path */
    *diag_querystart = querylength - 1;
    *diag_queryend = 0;
    for (i = 0; i < nunique; i++) {
      if (array[i]->querystart < *diag_querystart) {
	*diag_querystart = array[i]->querystart;
      }
      if (array[i]->queryend > *diag_queryend) {
	*diag_queryend = array[i]->queryend;
      }
    }
    debug1(printf("diag_querystart = %d, diag_queryend = %d\n",*diag_querystart,*diag_queryend));
#endif


    /* Set minactive */
    for (querypos = 0; querypos < activestart; querypos++) {
      minactive[querypos] = 0U;
      debug1(printf("active end: minactive at querypos %d is %u\n",querypos,minactive[querypos]));
    }

    diagonal = array[0]->diagonal;
    for ( ; querypos <= array[0]->queryend; querypos++) {
#ifdef PMAP
      if (diagonal + 3*querypos < EXTRA_BOUNDS) {
	minactive[querypos] = chrinit;
      } else {
	minactive[querypos] = (Chrpos_T) (chrinit + diagonal + 3*querypos - EXTRA_BOUNDS);
      }
#else
      if (diagonal + querypos < EXTRA_BOUNDS) {
	minactive[querypos] = chrinit;
      } else {
	minactive[querypos] = (Chrpos_T) (chrinit + diagonal + querypos - EXTRA_BOUNDS);
      }
#endif
      debug1(printf("to first diagonal: minactive at querypos %d is %u\n",querypos,minactive[querypos]));
    }

    i = 0;
    while (i < nunique) {
      j = i+1;
      while (j < nunique && array[j]->queryend <= array[i]->queryend) {
	j++;
      }
      if (j < nunique) {
	diagonal = array[i]->diagonal; /* Use this diagonal for end of query only if j < nunique */
	for ( ; querypos <= array[j]->queryend; querypos++) {
#ifdef PMAP
	  if (diagonal + 3*querypos < EXTRA_BOUNDS) {
	    minactive[querypos] = chrinit;
	  } else {
	    minactive[querypos] = (Chrpos_T) (chrinit + diagonal + 3*querypos - EXTRA_BOUNDS);
	  }
#else
	  if (diagonal + querypos < EXTRA_BOUNDS) {
	    minactive[querypos] = chrinit;
	  } else {
	    minactive[querypos] = (Chrpos_T) (chrinit + diagonal + querypos - EXTRA_BOUNDS);
	  }
#endif
	  debug1(printf("middle diagonals: minactive at querypos %d is %u\n",querypos,minactive[querypos]));
	}
      }
      i = j;
    }

    debug1(printf("diagonal for region to end of query is %u\n",diagonal));
    for ( ; querypos < querylength; querypos++) {
#ifdef PMAP
      if (diagonal + 3*querypos < EXTRA_BOUNDS) {
	minactive[querypos] = chrinit;
      } else {
	minactive[querypos] = (Chrpos_T) (chrinit + diagonal + 3*querypos - EXTRA_BOUNDS);
      }
#else
      if (diagonal + querypos < EXTRA_BOUNDS) {
	minactive[querypos] = chrinit;
      } else {
	minactive[querypos] = (Chrpos_T) (chrinit + querypos - EXTRA_BOUNDS);
      }
#endif
      debug1(printf("to end of query: minactive at querypos %d is %u\n",querypos,minactive[querypos]));
    }


    /* Set maxactive */
    for (querypos = querylength-1; querypos > activeend; --querypos) {
      maxactive[querypos] = chrterm;
      debug1(printf("active end: maxactive at querypos %d is %u\n",querypos,maxactive[querypos]));
    }

    diagonal = array[nunique-1]->diagonal;
    for ( ; querypos >= array[nunique-1]->querystart; --querypos) {
#ifdef PMAP
      if ((position = diagonal + 3*querypos + EXTRA_BOUNDS) > genomiclength) {
	maxactive[querypos] = chrterm;
      } else {
	maxactive[querypos] = (Chrpos_T) chrinit + position;
      }
#else
      if ((position = diagonal + querypos + EXTRA_BOUNDS) > genomiclength) {
	maxactive[querypos] = chrterm;
      } else {
	maxactive[querypos] = (Chrpos_T) chrinit + position;
      }
#endif
      debug1(printf("to last diagonal: maxactive at querypos %d is %u\n",querypos,maxactive[querypos]));
    }

    i = nunique-1;
    while (i >= 0) {
      j = i-1;
      while (j >= 0 && array[j]->querystart > Diag_querystart(array[i])) {
	--j;
      }
      if (j >= 0) {
	diagonal = array[i]->diagonal; /* Use this diagonal to beginning of query only if j >= 0 */
	for ( ; querypos >= array[j]->querystart; --querypos) {
#ifdef PMAP
	  if ((position = diagonal + 3*querypos + EXTRA_BOUNDS) > genomiclength) {
	    maxactive[querypos] = chrterm;
	  } else {
	    maxactive[querypos] = (Chrpos_T) chrinit + position;
	  }
#else
	  if ((position = diagonal + querypos + EXTRA_BOUNDS) > genomiclength) {
	    maxactive[querypos] = chrterm;
	  } else {
	    maxactive[querypos] = (Chrpos_T) chrinit + position;
	  }
#endif
	  debug1(printf("middle diagonals: maxactive at querypos %d is %u\n",querypos,maxactive[querypos]));
	}
      }
      i = j;
    }

    debug1(printf("diagonal for region to beginning of query is %u\n",diagonal));
    for ( ; querypos >= 0; --querypos) {
#ifdef PMAP
      if ((position = diagonal + 3*querypos + EXTRA_BOUNDS) > genomiclength) {
	maxactive[querypos] = chrterm;
      } else {
	maxactive[querypos] = (Chrpos_T) chrinit + position;
      }
#else
      if ((position = diagonal + querypos + EXTRA_BOUNDS) > genomiclength) {
	maxactive[querypos] = chrterm;
      } else {
	maxactive[querypos] = (Chrpos_T) chrinit + position;
      }
#endif
      debug1(printf("to beginning of query: maxactive at querypos %d is %u\n",querypos,maxactive[querypos]));
    }

    FREEA(array);
  }

  return nunique;
}



void
Diag_max_bounds (Chrpos_T *minactive, Chrpos_T *maxactive,
		 int querylength, Chrpos_T chrstart, Chrpos_T chrend,
		 Univcoord_T chroffset, Univcoord_T chrhigh, bool plusp) {
  Chrpos_T chrinit, chrterm;
  int querypos;
  
  if (plusp == true) {
    chrinit = chrstart;
    chrterm = chrend;
  } else {
    chrinit = (chrhigh - chroffset) - chrend;
    chrterm = (chrhigh - chroffset) - chrstart;
  }

  for (querypos = 0; querypos < querylength; querypos++) {
    minactive[querypos] = chrinit;
    maxactive[querypos] = chrterm;
  }

  return;
}
