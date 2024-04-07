/* $Id: 001d2bb99aac2affc5864188771452c178c9e83e $ */
#ifndef SPLICEENDS_INCLUDED
#define SPLICEENDS_INCLUDED

typedef struct Spliceends_T *Spliceends_T;

#include "bool.h"
#include "types.h"	/* For Splicetype_T */
#include "genomicpos.h"
#include "univcoord.h"
#include "chrnum.h"

#include "compress.h"
#include "genomebits.h"
#include "knownsplicing.h"
#include "indexdb.h"
#include "localdb-read.h"
#include "mergeinfo.h"
#include "stage1hr.h"
#include "spliceendsgen.h"
#include "vectorpool.h"


#define T Spliceends_T
struct T {
  int id;
  bool checkedout_p;

  bool boundedp;

  int nspliceends;
  Splicetype_T splicetype;
  int sensedir;

  int *mismatch_positions_left;	/* allocated memory for computations */
  int *mismatch_positions_right; /* allocated memory for computations */

  int *matchlengths;		

  int *splice_qpos;		/* splice_qpos */
  int *distal_lengths;		/* for qstart, splice_qpos; for qend, querylength - splice_qpos */
  int *distal_trimpos;		/* distal endpoint */

  Univcoord_T *partners;	/* positions, not univdiagonals */

  int *medial_nmismatches;
  int *distal_nmismatches;

  double *medial_probs;
  double *distal_probs;
};


extern void
Spliceends_setup (bool *circularp_in,
		  Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in,
		  Univcoord_T genomelength_in, Indexdb_T indexdb_in,
		  int index1part_in, int index1interval_in,
		  int max_insertionlen_in, int max_deletionlen_in, Chrpos_T shortsplicedist,
		  Localdb_T localdb_in, bool allow_soft_clips_p_in,
		  bool novelsplicingp, bool knownsplicingp);

/* Called only by Spliceendsgen_free_memory */
extern void
Spliceends_free (T *old);

extern T
Spliceends_new (int id, int querylength, Vectorpool_T vectorpool);

extern int
Spliceends_middle_plus (Univcoord_T **diagonals,
			Stage1_T stage1, int qstart, int qend, int querylength,
			Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
			Compress_T query_compress, char *queryptr,
			Univcoord_T *diagonals_alloc, unsigned short *localdb_alloc,
			Localdb_T localdb, int localdb_nmismatches_allowed);

extern int
Spliceends_middle_minus (Univcoord_T **diagonals,
			 Stage1_T stage1, int qstart, int qend, int querylength,
			 Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
			 Compress_T query_compress, char *queryptr,
			 Univcoord_T *diagonals_alloc, unsigned short *localdb_alloc,
			 Localdb_T localdb, int localdb_nmismatches_allowed);

extern int
Spliceends_trim_qstart_nosplice (int *nmismatches_to_trimpos, int *mismatch_positions, int total_nmismatches, 
				 int pos5, int pos3);

extern Univcoord_T
Spliceends_indel_qstart (int nosplice_trimpos, 
			 Univcoord_T univdiagonal, int querylength,
			 Univcoord_T chroffset, Univcoord_T chrhigh,
			 bool plusp, int genestrand,
			 int localdb_nmismatches_allowed, Univcoord_T *novel_diagonals_alloc,
			 unsigned short *localdb_alloc, Stage1_T stage1,
			 Compress_T query_compress, char *queryptr);

extern int
Spliceends_trimmed_qstarts (T *new, int *nosplice_trimpos, int *farsplice_trimpos,
			    int *nosplice_nmismatches, int *farsplice_nmismatches,
			    bool *splice5p, Splicetype_T *splicetype5, double *ambig_prob_5,
			    int try_sensedir, Univcoord_T univdiagonal, int querylength,
			    int qend, int exon_origin, Chrnum_T chrnum, Univcoord_T chroffset,
			    bool plusp, int genestrand, int localdb_nmismatches_allowed, bool innerp, bool salvagep,
			    int *mismatch_positions_alloc, Univcoord_T *novel_diagonals_alloc,
			    unsigned short *localdb_alloc, Stage1_T stage1,
			    Knownsplicing_T knownsplicing, Vectorpool_T vectorpool,
			    Spliceendsgen_T spliceendsgen, Compress_T query_compress, char *queryptr,
			    Genomebits_T genomebits, Genomebits_T genomebits_alt,
			    bool find_splices_p);

extern bool
Spliceends_qstart_trim (int *trimpos, int *nmismatches_to_trimpos,
			int *found_sensedir, Splicetype_T *splicetype, double *ambig_prob_qstart,
			Knownsplicing_T knownsplicing, int try_sensedir,
			Univcoord_T univdiagonal, int querylength,
			int pos3, int exon_origin, Chrnum_T chrnum, Univcoord_T chroffset,
			bool plusp, int genestrand, int *mismatch_positions_alloc,
			Vectorpool_T vectorpool, Spliceendsgen_T spliceendsgen,
			Compress_T query_compress, char *queryptr,
			Genomebits_T genomebits, Genomebits_T genomebits_alt, bool find_splices_p);

extern int
Spliceends_trim_qend_nosplice (int *nmismatches_to_trimpos, int *mismatch_positions, int total_nmismatches,
			       int pos5, int pos3, int querylength);

extern Univcoord_T
Spliceends_indel_qend (int nosplice_trimpos, 
		       Univcoord_T univdiagonal, int querylength,
		       Univcoord_T chroffset, Univcoord_T chrhigh,
		       bool plusp, int genestrand,
		       int localdb_nmismatches_allowed, Univcoord_T *novel_diagonals_alloc,
		       unsigned short *localdb_alloc, Stage1_T stage1,
		       Compress_T query_compress, char *queryptr);

extern int
Spliceends_trimmed_qends (T *new, int *nosplice_trimpos, int *farsplice_trimpos,
			  int *nosplice_nmismatches, int *farsplice_nmismatches,
			  bool *splice3p, Splicetype_T *splicetype3, double *ambig_prob_3,
			  int try_sensedir, Univcoord_T univdiagonal, int querylength,
			  int qstart, int exon_origin, Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			  bool plusp, int genestrand, int localdb_nmismatches_allowed, bool innerp, bool salvagep,
			  int *mismatch_positions_alloc, Univcoord_T *novel_diagonals_alloc,
			  unsigned short *localdb_alloc, Stage1_T stage1,
			  Knownsplicing_T knownsplicing, Vectorpool_T vectorpool,
			  Spliceendsgen_T spliceendsgen, Compress_T query_compress, char *queryptr,
			  Genomebits_T genomebits, Genomebits_T genomebits_alt,
			  bool find_splices_p);

extern bool
Spliceends_qend_trim (int *trimpos, int *nmismatches_to_trimpos,
		      int *found_sensedir, Splicetype_T *splicetype, double *ambig_prob_qend,
		      Knownsplicing_T knownsplicing, int try_sensedir,
		      Univcoord_T univdiagonal, int querylength,
		      int pos5, int exon_origin, Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		      bool plusp, int genestrand, int *mismatch_positions_alloc, 
		      Vectorpool_T vectorpool, Spliceendsgen_T spliceendsgen,
		      Compress_T query_compress, char *queryptr,
		      Genomebits_T genomebits, Genomebits_T genomebits_alt,
		      bool find_splices_p);

extern Univcoord_T *
Spliceends_qstart_resolve (int *ndiagonals, int *local_nmismatches, int pos3, int querylength,
			   Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
			   Compress_T query_compress, char *queryptr, bool plusp, int genestrand,
			   Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			   Stage1_T stage1, int localdb_nmismatches_allowed);

extern Univcoord_T *
Spliceends_qend_resolve (int *ndiagonals, int *local_nmismatches, int pos5, int querylength,
			 Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
			 Compress_T query_compress, char *queryptr, bool plusp, int genestrand,
			 Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			 Stage1_T stage1, int localdb_nmismatches_allowed);

#undef T
#endif

