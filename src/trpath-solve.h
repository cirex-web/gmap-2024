/* $Id: 87e9c3985d63da93b48245d441d411f3909e32cd $ */
#ifndef TRPATH_SOLVE_INCLUDED
#define TRPATH_SOLVE_INCLUDED

#include "bool.h"
#include "types.h"
#include "trpath.h"
#include "path.h"
#include "method.h"

#include "trdiagdef.h"
#include "trdiag.h"

#include "indel.h"
#include "compress.h"
#include "shortread.h"
#include "genomebits.h"
#include "ef64.h"

#include "intlistpool.h"
#include "uintlistpool.h"
#include "listpool.h"
#include "trpathpool.h"
#include "pathpool.h"
#include "vectorpool.h"



#define T Trpath_T

extern T
Trpath_solve_from_trdiagonal (int *found_score, Trcoord_T trdiagonal, int tstart, int tend,
			      
			      Trnum_T trnum, Trcoord_T troffset, Trcoord_T trhigh,
			      Compress_T query_compress_tr, bool tplusp, int querylength,
			      int *mismatch_positions_alloc, bool want_lowest_coordinate_p,
			      Indelinfo_T indelinfo,

			      Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			      Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
			      Method_T method);

extern T
Trpath_solve_from_diagonals (int *found_score, Trcoord_T middle_trdiagonal,
			     int middle_trdiagonal_qstart, int middle_trdiagonal_qend,
			     int middle_nmismatches,
			     Trdiag_T qstart_trdiag, Trdiag_T qend_trdiag,
			     bool tplusp, int querylength, Compress_T query_compress_tr,
			     int *mismatch_positions_alloc,
			     Trnum_T trnum, Trcoord_T troffset, Trcoord_T trhigh,
			     bool want_lowest_coordinate_p, Indelinfo_T indelinfo,
			     Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			     Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
			     Method_T method);

extern T
Trpath_solve_from_trstart (Trcoord_T trdiagonal,
			   bool tplusp, int querylength, Compress_T query_compress_tr,
			   Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			   Trpathpool_T trpathpool, Method_T method);

extern T
Trpath_solve_from_trend (Trcoord_T trdiagonal,
			 bool tplusp, int querylength, Compress_T query_compress_tr,
			 Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			 Trpathpool_T trpathpool, Method_T method);
extern T
Trpath_solve_from_ends (int *found_score,
			Trcoord_T trdiagonal_i, int pos5_0, int pos3_0,
			Trcoord_T trdiagonal_j, int pos5_1, int pos3_1,
			bool tplusp, int querylength, Compress_T query_compress_tr,
			Trnum_T trnum, Trcoord_T troffset, Trcoord_T trhigh,
			bool want_lowest_coordinate_p, Indelinfo_T indelinfo,
			Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
			Method_T method);


extern void
Trpath_solve_setup (Genomebits_T transcriptomebits_in, EF64_T transcript_ef64_in,
		    int max_insertionlen_in, int max_deletionlen_in);

#undef T
#endif


