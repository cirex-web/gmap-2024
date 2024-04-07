/* $Id: 16fa7de374866862d560d01996466711aa0a9919 $ */
#ifndef SPLICE_INCLUDED
#define SPLICE_INCLUDED

typedef struct Spliceinfo_T *Spliceinfo_T;

#include "bool.h"
#include "types.h"
#include "genomicpos.h"
#include "compress.h"
#include "genomebits.h"
#include "knownsplicing.h"
#include "univcoord.h"
#include "indel.h"
#include "stage1hr.h"
#include "localdb-read.h"


#define SPLICE_PROB_HIGH 0.9
#define SPLICE_PROB_LOW 0.2

#define T Spliceinfo_T
struct T {
  int *mismatch_positions_left1;
  int *mismatch_positions_right1;
  int *mismatch_positions_left2;
  int *mismatch_positions_right2;

  int *segmenti_sites_alloc1;
  int *segmenti_knowni_alloc1;
  int *segmenti_sites_alloc2;
  int *segmenti_knowni_alloc2;

  int *segmentk1_sites_alloc1;
  int *segmentk1_knowni_alloc1;
  int *segmentk1_sites_alloc2;
  int *segmentk1_knowni_alloc2;

  int *segmentk2_sites_alloc1;
  int *segmentk2_knowni_alloc1;
  int *segmentk2_sites_alloc2;
  int *segmentk2_knowni_alloc2;

  int *segmentj_sites_alloc1;
  int *segmentj_knowni_alloc1;
  int *segmentj_sites_alloc2;
  int *segmentj_knowni_alloc2;
};


extern void
Splice_setup (Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in,
	      int max_insertionlen_in, int max_deletionlen_in,
	      bool novelsplicingp_in);

extern void
Spliceinfo_free (T *old);

extern T
Spliceinfo_new (int querylength);


extern int
Splice_resolve (int *trimpos5, int *trimpos3,
		Univcoord_T *middle_univdiagonal, int *splice_qpos_i, int *splice_qpos_j,
		int *best_nindels, int *best_indel_pos,
		int *best_nmismatches_i, int *nmismatches_middle, int *best_nmismatches_j,
		int *best_nmismatches_indel,

		int *best_ref_nmismatches_i, int *ref_nmismatches_middle, int *best_ref_nmismatches_j,
		int *best_ref_nmismatches_indel,
		double *best_donor1_prob, double *best_acceptor1_prob,
		double *best_donor2_prob, double *best_acceptor2_prob,

		Univcoord_T univdiagonal_i, Univcoord_T univdiagonal_j,
		Stage1_T stage1, Compress_T query_compress, char *queryptr,
		bool plusp, Univcoord_T chroffset, Univcoord_T chrhigh,

		Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc, Localdb_T localdb,
		int localdb_nmismatches_allowed,

		int pos5, int pos3, int querylength,
		Indelinfo_T indelinfo, Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
		bool sense_forward_p, int genestrand,
		bool check_support_p, bool trim5p, bool trim3p);

extern int
Splice_nomiddle (int *trimpos5, int *trimpos3,
		 int *best_nindels, int *best_indel_pos,
		 int *best_nmismatches_i, int *best_nmismatches_j,
		 int *best_nmismatches_indel,

		 int *best_ref_nmismatches_i, int *best_ref_nmismatches_j,
		 int *best_ref_nmismatches_indel,
		 double *best_donor_prob, double *best_acceptor_prob,

		 Univcoord_T univdiagonal_i, Univcoord_T univdiagonal_j,
		 Compress_T query_compress,
		 bool plusp, Univcoord_T chroffset, Univcoord_T chrhigh,

		 int pos5, int pos3, int querylength,
		 Indelinfo_T indelinfo, Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
		 bool sense_forward_p, int genestrand,
		 bool check_support_p, bool trim5p, bool trim3p);

extern int
Splice_fusion_sense (char *donor1, char *donor2, char *acceptor1, char *acceptor2,
		     double *best_donor_prob, double *best_acceptor_prob,

		     int *best_nmismatches_D, int *best_nmismatches_A,
		     int *best_ref_nmismatches_D, int *best_ref_nmismatches_A,
		       
		     Univcoord_T univdiagonalD, Univcoord_T univdiagonalA,
		     Compress_T query_compress_D, bool plusDp, Univcoord_T chroffset_D,
		     Compress_T query_compress_A, bool plusAp, Univcoord_T chroffset_A,
		     
		     int queryposD, int queryposA, int querylength,
		     Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
		     int genestrand);


extern int
Splice_fusion_antisense (char *donor1, char *donor2, char *acceptor1, char *acceptor2,
			 double *best_donor_prob, double *best_acceptor_prob,

			 int *best_nmismatches_A, int *best_nmismatches_D,
			 int *best_ref_nmismatches_A, int *best_ref_nmismatches_D,
		       
			 Univcoord_T univdiagonalA, Univcoord_T univdiagonalD,
			 Compress_T query_compress_A, bool plusAp, Univcoord_T chroffset_A,
			 Compress_T query_compress_D, bool plusDp, Univcoord_T chroffset_D,
			 
			 int queryposA, int queryposD, int querylength,
			 Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
			 int genestrand);

#undef T
#endif

