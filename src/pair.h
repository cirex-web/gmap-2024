/* $Id: pair.h 224673 2021-10-25 14:58:40Z twu $ */
#ifndef PAIR_INCLUDED
#define PAIR_INCLUDED

typedef struct Pair_T *Pair_T;

#include "bool.h"
#include "genomicpos.h"
#include "chrnum.h"
#include "list.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "sequence.h"
#include "reader.h"		/* For cDNAEnd_T */
#include "uintlist.h"
#include "genome.h"
#include "chimera.h"
#include "filestring.h"
#include "pairpool.h"


#define MATCHESPERGAP 3

typedef enum {CDS_CDNA, CDS_GENOMIC} CDStype_T;
typedef enum {CIGAR_ACTION_IGNORE, CIGAR_ACTION_WARNING, CIGAR_ACTION_NOPRINT, CIGAR_ACTION_ABORT} Cigar_action_T;


#define T Pair_T

extern void
Pair_setup (bool novelsplicingp_in, IIT_T splicesites_iit_in, int trim_indel_score_in,
	    bool print_margin_p_in, bool gff3_separators_p_in,
	    bool sam_insert_0M_p_in, bool force_xs_direction_p_in,
	    bool md_lowercase_variant_p_in, bool snps_p_in,
	    bool gff3_phase_swap_p_in, CDStype_T cdstype_in,
	    bool cigar_extended_p_in, Cigar_action_T cigar_action_in);
extern int
Pair_querypos (T this);
extern Chrpos_T
Pair_genomepos (T this);
extern char
Pair_cdna (T this);
extern char
Pair_comp (T this);
extern char
Pair_genome (T this);
extern char
Pair_genomealt (T this);
extern bool
Pair_gapp (T this);
extern bool
Pair_shortexonp (T this);
extern void
Pair_print_ends (List_T pairs);

extern void
Pair_set_genomepos (struct Pair_T *pairarray, int npairs, Univcoord_T chroffset,
		    Univcoord_T chrhigh, bool watsonp);
extern void
Pair_subtract_genomepos (struct T *pairs, int npairs, Chrpos_T adjustment);

#if 0
extern void
Pair_set_genomepos_list (List_T pairs, Univcoord_T chroffset, Univcoord_T chrhigh,
			 bool watsonp);
#endif
extern List_T
Pair_clip_bounded_list_5 (List_T source, int minpos, int maxpos);
extern List_T
Pair_clip_bounded_list_3 (List_T source, int minpos, int maxpos);
extern int
Pair_clip_bounded_array (struct T *source, int npairs, int minpos, int maxpos);

extern List_T
Pair_protect_end5 (List_T pairs);
extern List_T
Pair_protect_end3 (List_T pairs);
extern void
Pair_protect_list (List_T pairs);

extern T
Pair_new_out (int querypos, Chrpos_T genomepos, char cdna, char comp, char genome);
extern void
Pair_free_out (T *old);

extern int
Pair_translation_length (struct T *pairs, int npairs);
extern void
Pair_print_continuous (Filestring_T fp, struct T *pairs, int npairs, bool watsonp,
		       bool genomefirstp, int invertmode, bool nointronlenp);

extern void
Pair_print_continuous_byexon (Filestring_T fp, struct T *pairs, int npairs, bool watsonp, int invertmode);
extern void
Pair_print_alignment (Filestring_T fp, struct T *pairs, int npairs, Chrnum_T chrnum,
		      Univcoord_T chroffset, Univ_IIT_T chromosome_iit, bool watsonp,
		      int invertmode, bool nointronlenp, int wraplength);

extern void
Pair_print_pathsummary (Filestring_T fp, int pathnum, T start, T end, Chrnum_T chrnum,
			Univcoord_T chroffset, Univ_IIT_T chromosome_iit, bool referencealignp, 
			IIT_T altstrain_iit, char *strain, Univ_IIT_T contig_iit, char *dbversion,
			int querylength_given, int skiplength, int trim_start, int trim_end,
			int nexons, int matches, int unknowns, int mismatches, 
			int qopens, int qindels, int topens, int tindels,
			bool watsonp, int cdna_direction,
			int translation_start, int translation_end, int translation_length,
			int relaastart, int relaaend);

extern void
Pair_print_coordinates (Filestring_T fp, struct T *pairs, int npairs, Chrnum_T chrnum,
			Univcoord_T chroffset, Univ_IIT_T chromosome_iit,
			bool watsonp, int invertmode);

extern int
Pair_cmp (const void *a, const void *b);

extern void
Pair_dump_one (T this, bool zerobasedp);
extern void
Pair_dump_list (List_T pairs, bool zerobasedp);
extern void
Pair_dump_array (struct T *pairs, int npairs, bool zerobasedp);
extern void
Pair_dump_array_stderr (struct T *pairs, int npairs, bool zerobasedp);
extern void
Pair_dump_genome_array (struct T *pairs, int npairs);
extern void
Pair_dump_comp_array (struct T *pairs, int npairs);
extern Chrpos_T
Pair_genomicpos (struct T *pairs, int npairs, int querypos, bool headp);
extern int
Pair_codon_changepos (struct T *pairs, int npairs, int aapos, int cdna_direction);

extern bool
Pair_identical_p (List_T pairs1, List_T pairs2);
extern void
Pair_check_list_pairs (List_T pairs);
extern void
Pair_check_list_path (List_T path);
extern bool
Pair_check_array_pairs (struct T *pairs, int npairs);
extern bool
Pair_check_array_path (struct T *path, int npairs);

extern void
Pair_print_exonsummary (Filestring_T fp, struct T *pairs, int npairs, Chrnum_T chrnum,
			Univcoord_T chroffset, Genome_T genome, Genome_T genomealt,
			Univ_IIT_T chromosome_iit, bool watsonp, int cdna_direction,
			bool genomefirstp, int invertmode);

extern int
Pair_cigar_length (List_T tokens);
extern void
Pair_print_tokens (Filestring_T fp, List_T tokens);
extern void
Pair_tokens_free (List_T *tokens);
extern List_T
Pair_tokens_copy (List_T old);

extern void
Pair_print_gff3 (Filestring_T fp, struct T *pairs, int npairs, int pathnum, char *accession, char *restofheader,
		 T start, T end, Genome_T genome, Chrnum_T chrnum, Univ_IIT_T chromosome_iit,
		 int translation_end, int querylength_given, int skiplength, int matches, int mismatches, 
		 int qindels, int tindels, int unknowns, bool watsonp, int cdna_direction,
		 bool gff_gene_format_p, bool gff_estmatch_format_p, char *sourcename);

#ifdef GSNAP
extern void
Pair_print_m8 (Filestring_T fp, struct T *pairs_querydir, int npairs, bool invertedp,
	       Chrnum_T chrnum, Shortread_T queryseq, Shortread_T headerseq,
	       char *acc_suffix, Univ_IIT_T chromosome_iit);
#endif

#ifndef PMAP
extern void
Pair_print_bedpe (Filestring_T fp, struct T *pairs_querydir, int npairs,
		  Chrnum_T chrnum, bool watsonp, Univ_IIT_T chromosome_iit);
#endif

extern void
Pair_fix_cdna_direction_array (struct T *pairs_querydir, int npairs, int cdna_direction);
extern int
Pair_guess_cdna_direction_array (int *sensedir, struct T *pairs_querydir, int npairs, bool invertedp,
				 Genome_T genome, Genome_T genomealt, Univcoord_T chroffset, bool watsonp);
extern int
Pair_gsnap_nsegments (int *total_nmismatches, int *total_nindels, int *nintrons,
		      int *nindelbreaks, struct T *pairs, int npairs, int querylength);
extern int
Pair_tokens_cigarlength (List_T tokens);


extern int
Pair_circularpos (int *alias, struct T *pairs, int npairs, Chrpos_T chrlength, bool plusp, int querylength);
extern void
Pair_alias_circular (struct T *pairs, int npairs, Chrpos_T chrlength);
extern void
Pair_unalias_circular (struct T *pairs, int npairs, Chrpos_T chrlength);

extern void
Pair_print_sam_nomapping (Filestring_T fp, char *abbrev, char *acc1, char *acc2, char *queryseq_ptr,
			  char *quality_string, int querylength, int quality_shift,
			  bool first_read_p, bool sam_paired_p, char *sam_read_group_id);
extern void
Pair_print_sam_nomapping_flipped (Filestring_T fp, char *abbrev, char *chrstring,
				  char *genomicseg, int genomiclength,
				  bool first_read_p, bool sam_paired_p, char *sam_read_group_id);


extern struct T *
Pair_hardclip (int *clipped_npairs, int hardclip_start, int hardclip_end,
	       struct T *pairs, int npairs, int querylength);

extern List_T
Pair_clean_cigar (List_T tokens, bool watsonp);
extern List_T
Pair_compute_cigar (bool *intronp, int *hardclip_start, int *hardclip_end, struct T *pairs, int npairs, int querylength_given,
		    bool watsonp, bool flippedp, int chimera_part);

extern void
Pair_print_sam (Filestring_T fp, char *abbrev, struct Pair_T *pairarray, int npairs,
		char *acc1, char *acc2, Genome_T genome, Chrnum_T chrnum, Univ_IIT_T chromosome_iit,
		char *queryseq_ptr, char *quality_string,
		int hardclip_low, int hardclip_high, int querylength_given,
		bool watsonp, int sensedir, int chimera_part, Chimera_T chimera,
		int quality_shift, bool first_read_p, int pathnum, int npaths_primary, int npaths_altloc,
		int absmq_score, int second_absmq, Chrpos_T chrpos, Chrpos_T chrlength,
		int mapq_score, bool sam_paired_p, char *sam_read_group_id);

extern void
Pair_print_sam_flipped (Filestring_T fp, char *abbrev, struct T *pairarray, int npairs,
			char *acc1, char *acc2, Genome_T genome, Chrnum_T chrnum, Univ_IIT_T chromosome_iit,
			char *genomicseg, int genomiclength, int hardclip_low, int hardclip_high, int querylength_given,
			bool watsonp, int sensedir, int chimera_part, Chimera_T chimera,
			int quality_shift, bool first_read_p, int pathnum, int npaths_primary, int npaths_altloc,
			int absmq_score, int second_absmq, int querypos,
			int mapq_score, bool sam_paired_p, char *sam_read_group_id);

extern List_T
Pair_compute_md_string (int *nmismatches_refdiff, int *nmismatches_bothdiff, int *nindels,
			struct T *pairs, int npairs, bool watsonp, List_T cigar_tokens);

extern Uintlist_T
Pair_exonbounds (struct T *pairs, int npairs);

extern void
Pair_print_pslformat_nt (Filestring_T fp, struct T *pairs, int npairs, T start, T end,
			 Sequence_T queryseq, Genome_T genome, Chrnum_T chrnum,
			 Univ_IIT_T chromosome_iit, int matches, int unknowns, int mismatches, 
			 bool watsonp);


extern void
Pair_print_pslformat_pro (Filestring_T fp, struct T *pairs, int npairs, T start, T end,
			  Sequence_T queryseq, Genome_T genome, Chrnum_T chrnum,
			  Univ_IIT_T chromosome_iit, bool watsonp, int cdna_direction);

extern void
Pair_print_exons (Filestring_T fp, struct T *pairs, int npairs, int wraplength, int ngap, bool cdnap);

extern void
Pair_print_protein_genomic (Filestring_T fp, struct T *ptr, int npairs, int wraplength, bool forwardp);
#ifdef PMAP
extern void
Pair_print_nucleotide_cdna (Filestring_T fp, struct T *ptr, int npairs, int wraplength);
#else
extern void
Pair_print_protein_cdna (Filestring_T fp, struct T *ptr, int npairs, int wraplength, bool forwardp);
#endif

extern void
Pair_print_iit_map (Filestring_T fp, Sequence_T queryseq, char *accession,
		    T start, T end, Chrnum_T chrnum, Univ_IIT_T chromosome_iit);
extern void
Pair_print_iit_exon_map (Filestring_T fp, struct T *pairs, int npairs, Sequence_T queryseq, char *accession,
			 T start, T end, Chrnum_T chrnum, Univ_IIT_T chromosome_iit);
extern void
Pair_print_splicesites (Filestring_T fp, struct T *pairs, int npairs, char *accession,
			int nexons, Chrnum_T chrnum, Univ_IIT_T chromosome_iit, bool watsonp);
extern void
Pair_print_introns (Filestring_T fp, struct T *pairs, int npairs, char *accession,
		    int nexons, Chrnum_T chrnum, Univ_IIT_T chromosome_iit);
extern void
Pair_print_mask_introns (Filestring_T fp, struct T *pairs, int npairs,
			 Chrpos_T chrlength, int wraplength, bool include_utr_p);


extern int
Pair_nmatches_posttrim (int *max_match_length, List_T pairs, int pos5, int pos3);
extern int
Pair_array_nmatches_posttrim (struct T *pairs, int npairs, int pos5, int pos3);
extern int
Pair_nmismatches_region (int *nindelbreaks, int *nbadintrons, struct T *pairs, int npairs,
			 int trim_left, int trim_right, int start_amb_nmatches, int end_amb_nmatches,
			 int querylength);

extern int
Pair_goodness_simple (List_T pairs);
extern void
Pair_fracidentity_simple (int *matches, int *unknowns, int *mismatches, List_T pairs);
extern void
Pair_fracidentity (int *matches, int *unknowns, int *mismatches, 
		   int *qopens, int *qindels, int *topens, int *tindels,
		   int *ncanonical, int *nsemicanonical, int *nnoncanonical,
		   double *min_splice_prob, List_T pairs, int cdna_direction);
extern int
Pair_fracidentity_array (int *matches, int *unknowns, int *mismatches, int *qopens, int *qindels, 
			 int *topens, int *tindels, int *ncanonical, int *nsemicanonical, int *nnoncanonical,
			 double *min_splice_prob, struct T *ptr, int npairs, int cdna_direction);
extern int
Pair_fracidentity_score (List_T pairs);

extern double
Pair_frac_error (List_T pairs, int cdna_direction);

extern void
Pair_fracidentity_bounded (int *matches, int *unknowns, int *mismatches, 
			   int *qopens, int *qindels, int *topens, int *tindels,
			   int *ncanonical, int *nsemicanonical, int *nnoncanonical,
			   struct T *pairs, int npairs, 
			   int cdna_direction, int minpos, int maxpos);
extern void
Pair_matchscores (int *matchscores, struct T *ptr, int npairs);
extern int
Pair_maxnegscore (List_T pairs);


extern void
Pair_pathscores (bool *gapp, int *pathscores, struct T *ptr, int npairs, 
		 int cdna_direction, int querylength, cDNAEnd_T cdnaend, int pre_extension_slop);

extern int
Pair_cdna_direction (List_T pairs);
extern int
Pair_nexons_approx (List_T pairs);
extern int
Pair_nexons (struct T *pairs, int npairs);
extern bool
Pair_consistentp (int *ncanonical, struct T *pairs, int npairs, int cdna_direction);

#ifndef PMAP
extern void
Pairarray_chrpos_bounds (Chrpos_T *chrpos_start, Chrpos_T *chrpos_end,
			 struct T *pairarray, int npairs);
#endif

extern Chrpos_T
Pairarray_genomicbound_from_start (struct T *pairarray, int npairs, int overlap);
extern Chrpos_T
Pairarray_genomicbound_from_end (struct T *pairarray, int npairs, int overlap);
extern char *
Pairarray_genomic_sequence (int *seqlength, struct T *pairarray, int npairs);


extern T
Pair_start_bound (int *cdna_direction, List_T pairs, int breakpoint);
extern T
Pair_end_bound (int *cdna_direction, List_T pairs, int breakpoint);


extern void
Pair_trim_distances (int *trim5, int *trim3, List_T pairs);

extern List_T
Pair_trim_ends (bool *trim5p, bool *trim3p, List_T pairs, int ambig_end_length_5, int ambig_end_length_3);
extern void
Pair_split_circular (List_T *pairs_below, List_T *pairs_above, List_T pairs,
		     Chrpos_T chrlength, Pairpool_T pairpool, bool plusp);

extern void
Pair_flip (T this);

#undef T
#endif
