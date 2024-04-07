/* $Id$ */
#ifndef AUXINFO_INCLUDED
#define AUXINFO_INCLUDED

typedef struct Auxinfo_T *Auxinfo_T;

#include "bool.h"
#include "path.h"
#include "method.h"
#include "list.h"

#include "ef64.h"
#include "path.h"
#include "chrnum.h"
#include "univdiagpool.h"
#include "intlistpool.h"
#include "univcoord.h"
#include "listpool.h"
#include "pathpool.h"
#include "auxinfopool.h"
#include "transcriptpool.h"
#include "hitlistpool.h"


#define T Auxinfo_T
struct T {
  Chrnum_T chrnum;
  Univcoord_T chroffset;
  Univcoord_T chrhigh;

  bool solvedp;
  bool complete_sense_p;		/* tells whether best_sense_paths came from complete (as opposed to unextended) */
  bool complete_antisense_p;		/* tells whether best_antisense_paths came from complete (as opposed to unextended) */

  Intlist_T best_sense_partners;	/* index into other all_univdiagonals/auxinfo array */
  Intlist_T best_antisense_partners;	/* index into other all_univdiagonals/auxinfo array */

  List_T best_sense_paths;
  List_T best_antisense_paths;

  List_T unextended_sense_paths;
  List_T unextended_antisense_paths;
  List_T complete_sense_paths;
  List_T complete_antisense_paths;

  Method_T method;
  
  int qstart;
  int qend;
  int nmismatches;

  List_T right_univdiags;
  List_T left_univdiags;

  int mate_count;
  int mate_bestk;

  struct Auxinfo_T *rest;
};


extern void
Auxinfo_free_wpaths (T *old, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool,
		     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		     Listpool_T listpool, Pathpool_T pathpool,
		     Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);

extern void
Auxinfo_free (T *old, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool,
	      Intlistpool_T intlistpool, Hitlistpool_T hitlistpool);

extern void
Auxinfo_gc_wpaths (T *array, int n, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool,
		   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		   Listpool_T listpool, Pathpool_T pathpool,
		   Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);

extern void
Auxinfo_gc (T *array, int n, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool,
	    Intlistpool_T intlistpool, Hitlistpool_T hitlistpool);

extern void
Auxinfo_set_best_paths (T this, Hitlistpool_T hitlistpool);
extern void
Auxinfo_set_best_sense_paths (T this, Hitlistpool_T hitlistpool, bool only_complete_p);
extern void
Auxinfo_set_best_antisense_paths (T this, Hitlistpool_T hitlistpool, bool only_complete_p);

extern T
Auxinfo_new (Method_T method, int qstart, int qend, Auxinfopool_T auxinfopool);

extern T
Auxinfo_new_tr (Auxinfopool_T auxinfopool, Path_T path);

extern T
Auxinfo_new_univdiags (Method_T method, int qstart, int qend, int nmismatches,
		       List_T right_univdiags, List_T left_univdiags, Auxinfopool_T auxinfopool);

extern T
Auxinfo_new_sense_path (Path_T path, Hitlistpool_T hitlistpool, Auxinfopool_T auxinfopool);

extern T
Auxinfo_new_antisense_path (Path_T path, Hitlistpool_T hitlistpool, Auxinfopool_T auxinfopool);

extern T
Auxinfo_append (T new, T old);

extern void
Auxinfo_assign_chrinfo (Univcoord_T * univdiagonals, T *auxinfo, int n, int querylength);

extern void
Auxinfo_collect_paths (bool *foundp, List_T *sense_paths, List_T *antisense_paths,
		       T *auxinfo_array, int nunivdiagonals,
		       Hitlistpool_T hitlistpool);

extern void
Auxinfo_setup (EF64_T chromosome_ef64_in);

#undef T
#endif


