/* $Id: 67fbec1ccc58f9b656216306fe14b654b64e1acf $ */
#ifndef EXTENSION_SEARCH_INCLUDED
#define EXTENSION_SEARCH_INCLUDED

#include "method.h"
#include "types.h"
#include "genomicpos.h"
#include "mode.h"
#include "pass.h"

#include "auxinfo.h"

#include "list.h"
#include "iit-read-univ.h"
#include "ef64.h"
#include "genomebits.h"
#include "compress.h"
#include "shortread.h"

#include "stage1hr.h"
#include "indexdb.h"
#include "knownsplicing.h"
#include "knownindels.h"

#include "univdiagpool.h"
#include "auxinfopool.h"
#include "intlistpool.h"
#include "univcoord.h"
#include "listpool.h"


#define T Elt_T
typedef struct T *T;
struct T {
  int min_qstart;
  int max_qend;

  int nmatches;

  Univdiag_T *all_univdiags;
  int n_all_univdiags;

  Univdiag_T *univdiags;	/* Filtered by binary search to generate lowi and highi */
  int nunivdiags;

  int lowi;
  int highi;
};


extern void
Extension_search_setup (Mode_T mode,
			Univcoord_T genomelength_in, int circular_typeint_in, bool *circularp_in, EF64_T chromosome_ef64_in,
			Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in,
			Indexdb_T indexdb_fwd_in, Indexdb_T indexdb_rev_in,
			int index1part_in, int index1interval_in, int maxpaths_search_in,
			int max_insertionlen, int max_deletionlen, Chrpos_T shortsplicedist);

extern void
Elt_gc (List_T *set, Listpool_T listpool, Univdiagpool_T univdiagpool);

extern void
Extension_search (Univcoord_T **_univdiagonals_gplus, Auxinfo_T **auxinfo_gplus, int *nunivdiagonals_gplus,
		  Univcoord_T **_univdiagonals_gminus, Auxinfo_T **auxinfo_gminus, int *nunivdiagonals_gminus,

		  Stage1_T stage1, Compress_T query_compress_fwd, Compress_T query_compress_rev,

		  int querylength, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool,
		  Univcoordlistpool_T univcoordlistpool, Listpool_T listpool);

#undef T
#endif


