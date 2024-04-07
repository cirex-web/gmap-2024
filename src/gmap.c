static char rcsid[] = "$Id: 50ceb981edd1a50960517cebc97adb53ac1e129d $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* Needed to define pthread_t on Solaris */
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For strcpy */
#include <strings.h>		/* For rindex */
#include <ctype.h>
#include <math.h>		/* For rint */

#include "simd.h"

#if 0
#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#ifdef HAVE_SSE4_1
#include <smmintrin.h>
#endif
#if !defined(HAVE_SSE4_2)
/* Skip popcnt */
#elif defined(HAVE_POPCNT)
#include <immintrin.h>
#endif

#if !defined(HAVE_SSE4_2)
/* Skip mm_popcnt */
#elif defined(HAVE_MM_POPCNT)
#include <nmmintrin.h>
#endif

#if !defined(HAVE_SSE4_2)
/* Skip lzcnt/tzcnt */
#elif defined(HAVE_LZCNT) || defined(HAVE_TZCNT)
#include <immintrin.h>
#endif
#endif


#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#include <signal.h>

#include "except.h"
#include "mem.h"
#include "bool.h"
#include "fopen.h"
#include "access.h"
#include "filesuffix.h"

#include "sequence.h"
#include "match.h"
#include "matchpool.h"
#include "pairpool.h"
#include "diagpool.h"
#include "cellpool.h"
#include "stopwatch.h"
#include "translation.h"	/* For Translation_setup */
#include "genome.h"
#include "genome-write.h"
#include "compress-write.h"
#include "stage1.h"
#include "gregion.h"
#ifdef PMAP
#include "oligoindex_pmap.h"
#else
#include "oligoindex_hr.h"	/* For Oligoindex_hr_setup */
#endif
#include "stage2.h"
#include "splicestringpool.h"
#include "splicetrie.h"
#include "splicetrie_build.h"
#include "dynprog.h"
#include "dynprog_single.h"
#include "dynprog_genome.h"
#include "dynprog_end.h"
#include "pair.h"
#include "stage3.h"
#include "comp.h"
#include "chimera.h"
#ifdef PMAP
#include "oligop.h"		/* For Oligop_setup */
#include "backtranslation.h"
#else
#include "oligo.h"		/* For Oligo_setup */
#endif
#include "indexdb.h"
#include "result.h"
#include "request.h"
#include "intlist.h"
#include "list.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "datadir.h"

#include "filestring.h"
#include "outputtype.h"
#include "output.h"
#include "inbuffer.h"
#include "samheader.h"
#include "outbuffer.h"

#include "getopt.h"


#define MAX_QUERYLENGTH_FOR_ALLOC    100000
#define MAX_GENOMICLENGTH_FOR_ALLOC 1000000


#define STAGE1_FIRSTPAIR_SIZELIMIT 10000
#define STAGE1_STUTTER_SIZELIMIT 100
#define STAGE1_FILLIN_SIZELIMIT 10


#define POSSIBLE_OLIGOS 65536	/* 4^8 */
#define MAX_OLIGODEPTH 3.0
#define MAX_BADOLIGOS 0.30	/* Setting to 1.0 effectively turns this check off */
#define MAX_REPOLIGOS 0.40	/* Setting to 1.0 effectively turns this check off */

/* Value of 1 can miss end exons, but values larger than 1 can lead to
   very long (or infinite?) run times when combined with
   --intronlength */
#define MAX_CHIMERA_ITER 3

#define CHIMERA_PENALTY 30	/* A small value for chimera_margin will reduce this  */
#define CHIMERA_IDENTITY 0.98
#define CHIMERA_PVALUE 0.01
#define CHIMERA_FVALUE 6.634897	/* qnorm(CHIMERA_PVALUE/2)^2 */
#define CHIMERA_SLOP 90	/* in nucleotides */
#define CHIMERA_EXTEND 20	/* Was previously 8, but this missed exon-exon boundaries */

#define MIN_MATCHES 20


#define MAX_NALIGNMENTS 10


/* #define EXTRACT_GENOMICSEG 1 */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Chimera detection */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Chimera detection, details */
#ifdef DEBUG2A
#define debug2a(x) x
#else
#define debug2a(x)
#endif

/* stage3list_remove_duplicates */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif



/************************************************************************
 *   Global variables
 ************************************************************************/

static int translation_code = 1;
static bool alt_initiation_codons_p = false;

static Univ_IIT_T chromosome_iit = NULL;
static Univ_IIT_T altscaffold_iit = NULL;
static Univcoord_T genomelength;
static int circular_typeint = -1;
static int nchromosomes;
static bool *circularp = NULL;
static bool any_circular_p;

static bool *altlocp = NULL;
static Univcoord_T *alias_starts = NULL;
static Univcoord_T *alias_ends = NULL;
static Univ_IIT_T contig_iit = NULL;

/* For the usual case when aligning against a genome index or the -g flag */
/* Not used for --cmdline, --selfalign, --pairalign, or user_genomicsegs*/
static Genome_T global_genome = NULL; /* Set for everything except selfalign, pairalign, or user_genomicsegs */
static Genome_T global_genomealt = NULL;
#if 0
static Genomecomp_T *genomecomp_blocks = NULL;
#endif

#ifdef PMAP
static Alphabet_T required_alphabet = AA0;
static Alphabet_T alphabet = AA20; /* Initialize in case we have a usersegment */
static int alphabet_size = 20;	   /* Initialize in case we have a usersegment */
static Width_T index1part_aa = 7;
#else
static Width_T index1part;
#endif

static Indexdb_T indexdb_fwd = NULL;
static Indexdb_T indexdb_rev = NULL;

/* static Localdb_T localdb = NULL; */

static Width_T required_index1part = 0;
static Width_T index1interval;
static Width_T required_index1interval = 0;

/* static Width_T local1part = 8; */
/* static Width_T required_local1part = 0; */
/* static Width_T local1interval; */
/* static Width_T required_local1interval = 0; */

static IIT_T altstrain_iit = NULL;

/* Cmet and AtoI */
static char *user_modedir = NULL; /* user_cmetdir, user_atoidir */
static Mode_T mode = STANDARD;


static char *user_snpsdir = NULL;
static char *snps_root = (char *) NULL;
static IIT_T map_iit = NULL;
static int *map_divint_crosstable = NULL;

#ifdef PMAP
#if 0
static Width_T minindexsize = 3;	/* In stage 2; in aa */
static Width_T maxindexsize = 6;	/* In stage 2; in aa */
#endif
/* Now controlled by defect_rate */
static int maxpeelback = 20;	/* Needs to be at least indexsize
				   because stage 2 jumps by indexsize.
				   Also should exceed length of
				   repeated nucleotides (e.g., a
				   string of consecutive T's) */
#else
/* Making minindexsize too small can lead to spurious exons in stage 2 */
/* FOOBAR */
#if 0
static Width_T minindexsize = 8;	/* In stage 2; in nt.  Used if sampling required in stage 1. */
static Width_T maxindexsize = 8;	/* In stage 2; in nt */
#endif
/* was 20, but that is not sufficient in some cases */
static int maxpeelback = 60;	/* Needs to be at least indexsize
				   because stage 2 jumps by indexsize.
				   Also should exceed length of
				   repeated nucleotides (e.g., a
				   string of consecutive T's) */
#endif
static int maxpeelback_distalmedial = 100; /* Needs to be longer to fix bad end exons */

/* static int stuttercycles = 2; */
static int stutterhits = 3;
static int sufflookback = 60;
static int nsufflookback = 5;

#if 0
static int maxoligohits = 400; /* Must be smaller than ALLOC in oligoindex.c */
#endif
static int nullgap = 600;
static int extramaterial_end = 10;
static int extramaterial_paired = 8; /* Should be at least indexsize in nt */
static int extraband_single = 6; /* This is in addition to length2 -
				    length1.  If onesidegap is true in
				    dynprog.c, then this is equivalent
				    to extraband_single of 0.  Needs
				    to be > 0 to handle default
				    close_indels_mode. */
static int extraband_end = 6; /* Was 6.  Shouldn't differ from 0, since onesidegapp is true?
				 This is only on both sides of main diagonal */
static int extraband_paired = 14; /* This is in addition to length2 - length1 */

static int user_open = 0;
static int user_extend = 0;
static bool user_dynprog_p = false;

static Stopwatch_T stopwatch = NULL;


/************************************************************************
 *   Program options
 ************************************************************************/

/* Input options */
static char *user_genomedir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;

static char *user_genomicsegs = NULL; /* For the -g flag */

static char *user_cmdline = NULL;
static bool user_selfalign_p = false;
static bool user_pairalign_p = false;

static unsigned int part_modulus = 0;
static unsigned int part_interval = 1;

static char *read_files_command = NULL;


/* Compute options */
bool multiple_sequences_p = false;

static bool sharedp = false;
static bool preload_shared_memory_p = false;
static bool unload_shared_memory_p = false;
static bool expand_offsets_p = false;

#ifdef HAVE_MMAP
static Access_mode_T offsetsstrm_access = USE_ALLOCATE;
static Access_mode_T positions_access = USE_MMAP_PRELOAD;
static Access_mode_T locoffsetsstrm_access = USE_ALLOCATE;
static Access_mode_T locpositions_access = USE_ALLOCATE;

static Access_mode_T genome_access = USE_MMAP_PRELOAD;
#else
static Access_mode_T offsetsstrm_access = USE_ALLOCATE;
static Access_mode_T positions_access = USE_ALLOCATE;
static Access_mode_T locoffsetsstrm_access = USE_ALLOCATE;
static Access_mode_T locpositions_access = USE_ALLOCATE;

static Access_mode_T genome_access = USE_ALLOCATE;
#endif

static int min_intronlength = 9;
static int max_deletionlength = 30;
static int maxtotallen_bound = 2400000;

static bool split_large_introns_p = false;

/* Need to set higher than 200,000 for many human genes, such as ALK */
static int maxintronlen = 500000; /* Was used previously in stage 1.  Now used only in stage 2 and Stage3_mergeable. */
static int maxintronlen_ends = 10000; /* Used in stage 3 */

static int end_trimming_score = 0;
static int minendexon = 12;
static int maxextension = 1000000; /* Used in stage 1.  Not adjustable by user */
static int chimera_margin = 30;	/* Useful for finding readthroughs */
static int index1interval = 3; /* Stage 1 interval if user provides a genomic segment */
/* static char *referencefile = NULL; */

#if 0
#ifndef PMAP
static bool literalrefp = false;
#endif
#endif


/* static bool altstrainp = false; */
#ifdef HAVE_PTHREAD
static pthread_t output_thread_id, *worker_thread_ids;
static pthread_key_t global_request_key;
static int nworkers = 1;	/* (int) sysconf(_SC_NPROCESSORS_ONLN) */
#else
static int nworkers = 0;	/* (int) sysconf(_SC_NPROCESSORS_ONLN) */
#endif
#ifndef PMAP
static bool prune_poor_p = false;
static bool prune_repetitive_p = false;
#endif
static int canonical_mode = 1;
static bool cross_species_p = false;
static int homopolymerp = false;

static char *user_chrsubsetname = NULL;
static Univcoord_T chrsubset_start = 0;
static Univcoord_T chrsubset_end = -1;

static int close_indels_mode = +1;
static double microexon_spliceprob = 0.95;
static int suboptimal_score_start = -1; /* Determined by simulations to have minimal effect */
static int suboptimal_score_end = 3; /* Determined by simulations to have diminishing returns above 3 */

static int trim_indel_score = -2; /* was -4 */


/* Output options */
static unsigned int output_buffer_size = 1000;
static Printtype_T printtype = SIMPLE;
static bool exception_raise_p = true;
static bool debug_graphic_p = false;
static bool stage1debug = false;
static bool diag_debug = false;
static Stage3debug_T stage3debug = NO_STAGE3DEBUG;
static bool timingp = false;
static bool checkp = false;
static int maxpaths_report = 5;	/* 0 means 1 if nonchimeric, 2 if chimeric */
static bool quiet_if_excessive_p = false;
static double suboptimal_score_float = 0.50;
static bool require_splicedir_p = false;


/* GFF3 */
static bool gff3_separators_p = true;
static bool gff3_phase_swap_p = false;
static GFF3_fasta_annotation_T gff3_fasta_annotation_type = NO_ANNOTATION;
static CDStype_T cdstype = CDS_CDNA;

/* SAM */
/* Applicable to PMAP? */
static bool sam_paired_p = false;
static bool user_quality_shift = false;
static int quality_shift = 0;
static bool sam_headers_p = true;
static char *sam_read_group_id = NULL;
static char *sam_read_group_name = NULL;
static char *sam_read_group_library = NULL;
static char *sam_read_group_platform = NULL;
static bool sam_insert_0M_p = false;
static bool sam_cigar_extended_p = false;
static bool sam_flippedp = false;
static Cigar_action_T cigar_action = CIGAR_ACTION_WARNING;


static bool orderedp = false;
static bool failsonlyp = false;
static bool nofailsp = false;
static bool checksump = false;
static int chimera_overlap = 0;
static bool force_xs_direction_p = false;
static bool md_lowercase_variant_p = false;

/* Map file options */
static char *user_mapdir = NULL;
static char *map_iitfile = NULL;
static bool map_exons_p = false;
static bool map_bothstrands_p = false;
static bool print_comment_p = false;
static int nflanking = 0;

/* Alignment options */
static bool fulllengthp = false;
static int cds_startpos = -1;
static bool truncatep = false;
static int sense_try = 0;		/* both */
static int sense_filter = 0;		/* both */
static int cdna_direction_try = 0;	/* both */
static double min_trimmed_coverage = 0.0;
static double min_identity = 0.0;
static bool strictp = true;
/* static int proteinmode = 1; */
static bool nointronlenp = false;
static int print_margin_p = true;
static int invertmode = 0;
static int ngap = 3;
static int wraplength = 50;


/* Splicing IIT */
static bool novelsplicingp = true; /* Can be disabled with --nosplicing flag */
static bool knownsplicingp = false;
static bool distances_observed_p = false;
static Chrpos_T shortsplicedist = 2000000;
static char *user_splicingdir = (char *) NULL;
static char *splicing_file = (char *) NULL;
static IIT_T splicing_iit = NULL;
static bool amb_closest_p = false;

static int donor_typeint = -1;		/* for splicing_iit */
static int acceptor_typeint = -1;	/* for splicing_iit */

static int *splicing_divint_crosstable = NULL;
static Univcoord_T *splicesites = NULL;
static Splicetype_T *splicetypes = NULL;
static Chrpos_T *splicedists = NULL; /* maximum observed splice distance for given splice site */
static List_T *splicestrings = NULL;
static Genomecomp_T *splicefrags_ref = NULL;
static Genomecomp_T *splicefrags_alt = NULL;
static int nsplicesites = 0;

/* Splicing via splicesites */
static int *nsplicepartners_skip = NULL;
static int *nsplicepartners_obs = NULL;
static int *nsplicepartners_max = NULL;

static bool splicetrie_precompute_p = true;
static Trieoffset_T *trieoffsets_obs = NULL;
static Triecontent_T *triecontents_obs = NULL;
static Trieoffset_T *trieoffsets_max = NULL;
static Triecontent_T *triecontents_max = NULL;


/* Input/output */
static char *split_output_root = NULL;
static char *failedinput_root = NULL;
static bool appendp = false;
static Inbuffer_T inbuffer = NULL;
static Outbuffer_T outbuffer = NULL;
static unsigned int input_buffer_size = 1000; /* previously inbuffer_nspaces */


#ifdef PMAP
/* Used alphabetically: 01235789ABbCcDdEefGgHIiKkLlMmNnOoPQRSstuVvwXxYZ */
#else
/* Used alphabetically: 01235789AaBbCcDdEeFfGgHIijKkLlMmNnOoPpQRSsTtuVvwXxYZ */
#endif

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */
#ifdef PMAP
  {"alphabet", required_argument, 0, 'a'}, /* required_alphabet */
#endif
  {"kmer", required_argument, 0, 'k'}, /* required_index1part, index1part */
  {"sampling", required_argument, 0, 0}, /* required_nterval, index1interval */
  {"gseg", required_argument, 0, 'g'}, /* user_genomicsegs */
  {"selfalign", no_argument, 0, '1'}, /* user_selfalign_p */
  {"pairalign", no_argument, 0, '2'}, /* user_pairalign_p */
  {"cmdline", required_argument, 0, 0}, /* user_cmdline */
  {"part", required_argument, 0, 'q'}, /* part_modulus, part_interval */
  {"input-buffer-size", required_argument, 0, 0}, /* input_buffer_size */

  {"read-files-command", required_argument, 0, 0}, /* read_files_command */
  

  /* Compute options */
  {"use-shared-memory", required_argument, 0, 0}, /* sharedp */
  {"preload-shared-memory", no_argument, 0, 0},	  /* preload_shared_memory_p */
  {"unload-shared-memory", no_argument, 0, 0},	  /* unload_shared_memory_p */
#ifdef HAVE_MMAP
  {"batch", required_argument, 0, 'B'}, /* offsetsstrm_access, positions_access, genome_access */
#endif
  {"expand-offsets", required_argument, 0, 0}, /* expand_offsets_p */
  {"max-deletionlength",required_argument, 0, 0}, /* max_deletionlength */
  {"min-intronlength", required_argument, 0, 0}, /* min_intronlength */

  {"intronlength", required_argument, 0, 'K'},		/* maxintronlen, maxintronlen_ends */
  {"max-intronlength-middle", required_argument, 0, 0}, /* maxintronlen */
  {"max-intronlength-ends", required_argument, 0, 0}, /* maxintronlen_ends */
  {"split-large-introns", no_argument, 0, 0},	      /* split_large_introns_p */

  {"end-trimming-score", required_argument, 0, 0},      /* end_trimming_score */
  {"trim-end-exons", required_argument, 0, 0}, /* minendexon */
  {"totallength", required_argument, 0, 'L'}, /* maxtotallen_bound */
  {"chimera-margin", required_argument, 0, 'x'}, /* chimera_margin */
  {"no-chimeras", no_argument, 0, 0},		 /* chimera_margin */
#if 0
  {"reference", required_argument, 0, 'w'}, /* referencefile */
#else
  {"localsplicedist", required_argument, 0, 'w'}, /* shortsplicedist */
#endif
  {"translation-code", required_argument, 0, 0}, /* translation_code */
  {"alt-start-codons", no_argument, 0, 0}, /* alt_initiation_codons_p */

  {"nthreads", required_argument, 0, 't'}, /* nworkers */
  {"splicingdir", required_argument, 0, 0}, /* user_splicingdir */
  {"nosplicing", no_argument, 0, 0},	    /* novelsplicingp */
  {"use-splicing", required_argument, 0, 's'}, /* splicing_iit, knownsplicingp (was previously altstrainp) */
  {"chrsubset", required_argument, 0, 'c'}, /* user_chrsubsetname */
  {"canonical-mode", required_argument, 0, 0}, /* canonical_mode */
  {"cross-species", no_argument, 0, 0}, /* cross_species_p */
  {"indel-open", required_argument, 0, 0}, /* user_open, user_dynprog_p */
  {"indel-extend", required_argument, 0, 0}, /* user_extend, user_dynprog_p */
  {"homopolymer", no_argument, 0, 0},	/* homopolymerp */
#ifndef PMAP
  {"prunelevel", required_argument, 0, 'p'}, /* prune_poor_p, prune_repetitive_p */
#endif
  {"allow-close-indels", required_argument, 0, 0}, /* close_indels_mode, extraband_single */
  {"microexon-spliceprob", required_argument, 0, 0}, /* microexon_spliceprob */
  {"stage2-start", required_argument, 0, 0},	     /* suboptimal_score_start */
  {"stage2-end", required_argument, 0, 0},	     /* suboptimal_score_end */

  {"cmetdir", required_argument, 0, 0}, /* user_modedir */
  {"atoidir", required_argument, 0, 0}, /* user_modedir */
  {"mode", required_argument, 0, 0}, /* mode */

  /* Output options */
  {"output-buffer-size", required_argument, 0, 0}, /* output_buffer_size */
  {"summary", no_argument, 0, 'S'}, /* printtype */
  {"align", no_argument, 0, 'A'}, /* printtype */
  {"continuous", no_argument, 0, '3'}, /* printtype */
  {"continuous-by-exon", no_argument, 0, '4'}, /* printtype */
  {"noexceptions", no_argument, 0, '0'}, /* exception_raise_p */
  {"graphic", no_argument, 0, '6'}, /* debug_graphic_p */
  {"stage3debug", required_argument, 0, '8'}, /* stage3debug */
  {"diagnostic", no_argument, 0, '9'}, /* checkp */
  {"npaths", required_argument, 0, 'n'}, /* maxpaths_report */
#if 0
  {"quiet-if-excessive", no_argument, 0, 0}, /* quiet_if_excessive_p */
#endif
  {"format", required_argument, 0, 'f'}, /* printtype */
  {"failsonly", no_argument, 0, 0}, /* failsonlyp */
  {"nofails", no_argument, 0, 0}, /* nofailsp */
  {"split-output", required_argument, 0, 0}, /* split_output_root */
  {"failed-input", required_argument, 0, 0}, /* failedinput_root */
  {"append-output", no_argument, 0, 0},	     /* appendp */
  {"suboptimal-score", required_argument, 0, 0}, /* suboptimal_score_float */
  {"require-splicedir", no_argument, 0, 0}, /* require_splicedir_p */

  {"gff3-add-separators", required_argument, 0, 0}, /* gff3_separators_p */
  {"gff3-swap-phase", required_argument, 0, 0}, /* gff3_phase_swap_p */
  {"gff3-fasta-annotation", required_argument, 0, 0}, /* gff3_fasta_annotation_type */
  {"gff3-cds", required_argument, 0, 0}, /* cdstype */

#ifndef PMAP
  {"quality-protocol", required_argument, 0, 0}, /* quality_shift */
  {"quality-print-shift", required_argument, 0, 'j'}, /* quality_shift */
  {"no-sam-headers", no_argument, 0, 0},	/* sam_headers_p */
  {"sam-use-0M", no_argument, 0, 0},		/* sam_insert_0M_p */
  {"sam-extended-cigar", no_argument, 0, 0},	/* sam_cigar_extended_p */
  {"sam-flipped", no_argument, 0, 0},	/* sam_flippedp */
  {"read-group-id", required_argument, 0, 0},	/* sam_read_group_id */
  {"read-group-name", required_argument, 0, 0},	/* sam_read_group_name */
  {"read-group-library", required_argument, 0, 0}, /* sam_read_group_library */
  {"read-group-platform", required_argument, 0, 0}, /* sam_read_group_platform */
  {"force-xs-dir", no_argument, 0, 0},		    /* force_xs_direction_p */
  {"md-lowercase-snp", no_argument, 0, 0},	    /* md_lowercase_variant_p */
  {"action-if-cigar-error", required_argument, 0, 0},	/* cigar_action */
#endif

  {"ordered", no_argument, 0, 'O'}, /* orderedp */
  {"md5", no_argument, 0, '5'}, /* checksump */
  {"chimera-overlap", required_argument, 0, 'o'}, /* chimera_overlap */
  {"snpsdir", required_argument, 0, 'V'},   /* user_snpsdir */
  {"use-snps", required_argument, 0, 'v'}, /* snps_root */

  /* Map file options */
  {"mapdir", required_argument, 0, 'M'}, /* user_mapdir */
  {"map", required_argument, 0, 'm'},	/* map_iitfile */
  {"mapexons", no_argument, 0, 'e'}, /* map_exons_p */
  {"mapboth", no_argument, 0, 'b'}, /* map_bothstrands_p */
  {"nflanking", required_argument, 0, 'u'}, /* nflanking */
  {"print-comment", no_argument, 0, 0},	    /* print_comment_p */

  /* Alignment options */
  {"exons", required_argument, 0, 'E'}, /* printtype */
#ifdef PMAP
  {"protein_gen", no_argument, 0, 'P'}, /* printtype */
  {"nucleotide", no_argument, 0, 'Q'}, /* printtype */
#else
  {"protein_dna", no_argument, 0, 'P'}, /* printtype */
  {"protein_gen", no_argument, 0, 'Q'}, /* printtype */
  {"fulllength", no_argument, 0, 'F'}, /* fulllengthp */
  {"cdsstart", required_argument, 0, 'a'}, /* cds_startpos */
  {"truncate", no_argument, 0, 'T'}, /* truncatep */
  {"strand", required_argument, 0, 0}, /* cdna_direction_try */
  {"direction", required_argument, 0, 'z'}, /* sense_try, sense_filter */
#endif
  {"tolerant", no_argument, 0, 'Y'}, /* strictp */
  {"nolengths", no_argument, 0, 0},	/* nointronlenp */
  {"nomargin", no_argument, 0, 0},	     /* print_margin_p */
  {"invertmode", required_argument, 0, 'I'}, /* invertmode */
  {"introngap", required_argument, 0, 'i'}, /* ngap */
  {"wraplength", required_argument, 0, 'l'}, /* wraplength */
  
  /* Filtering options */
  {"min-trimmed-coverage", required_argument, 0, 0}, /* min_trimmed_coverage */
  {"min-identity", required_argument, 0, 0},	/* min_identity */

  /* Diagnostic options */
  {"time", no_argument, 0, 0},	/* timingp */

  /* Help options */
  {"check", no_argument, 0, 0}, /* check_compiler_assumptions */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  char *genomedir;

  fprintf(stdout,"\n");
#ifdef PMAP
  fprintf(stdout,"PMAP: Protein Mapping and Alignment Program\n");
#else
  fprintf(stdout,"GMAP: Genomic Mapping and Alignment Program\n");
#endif
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Build target: %s\n",TARGET);
  fprintf(stdout,"Features: ");
#ifdef HAVE_PTHREAD
  fprintf(stdout,"pthreads enabled, ");
#else
  fprintf(stdout,"no pthreads, ");
#endif
#ifdef HAVE_ALLOCA
  fprintf(stdout,"alloca available, ");
#else
  fprintf(stdout,"no alloca, ");
#endif
#ifdef HAVE_ZLIB
  fprintf(stdout,"zlib available, ");
#else
  fprintf(stdout,"no zlib, ");
#endif
#ifdef HAVE_MMAP
  fprintf(stdout,"mmap available, ");
#else
  fprintf(stdout,"no mmap, ");
#endif
#ifdef WORDS_BIGENDIAN
  fprintf(stdout,"bigendian, ");
#else
  fprintf(stdout,"littleendian, ");
#endif
#ifdef HAVE_SIGACTION
  fprintf(stdout,"sigaction available, ");
#else
  fprintf(stdout,"no sigaction, ");
#endif
#ifdef HAVE_64_BIT
  fprintf(stdout,"64 bits available");
#else
  fprintf(stdout,"64 bits not available");
#endif
  fprintf(stdout,"\n");

  fprintf(stdout,"Popcnt:");
#ifdef HAVE_POPCNT
  fprintf(stdout," popcnt/lzcnt/tzcnt");
#endif
#ifdef HAVE_MM_POPCNT
  fprintf(stdout," mm_popcnt");
#endif
#ifdef HAVE_BUILTIN_POPCOUNT
  fprintf(stdout," builtin_popcount");
#endif
  fprintf(stdout,"\n");

  fprintf(stdout,"Builtin functions:");
#ifdef HAVE_BUILTIN_CLZ
  fprintf(stdout," builtin_clz");
#endif
#ifdef HAVE_BUILTIN_CTZ
  fprintf(stdout," builtin_ctz");
#endif
#ifdef HAVE_BUILTIN_POPCOUNT
  fprintf(stdout," builtin_popcount");
#endif
  fprintf(stdout,"\n");


  fprintf(stdout,"SIMD functions compiled:");
#ifdef HAVE_ALTIVEC
  fprintf(stdout," Altivec");
#endif
#ifdef HAVE_MMX
  fprintf(stdout," MMX");
#endif
#ifdef HAVE_SSE
  fprintf(stdout," SSE");
#endif
#ifdef HAVE_SSE2
  fprintf(stdout," SSE2");
#endif
#ifdef HAVE_SSE3
  fprintf(stdout," SSE3");
#endif
#ifdef HAVE_SSSE3
  fprintf(stdout," SSSE3");
#endif
#ifdef HAVE_SSE4_1
  fprintf(stdout," SSE4.1");
#endif
#ifdef HAVE_SSE4_2
  fprintf(stdout," SSE4.2");
#endif
#ifdef HAVE_AVX2
  fprintf(stdout," AVX2");
#endif
#ifdef HAVE_AVX512
  fprintf(stdout," AVX512");
#endif
#ifdef HAVE_AVX512BW
  fprintf(stdout," AVX512BW");
#endif
  fprintf(stdout,"\n");


#ifdef PMAP
  fprintf(stdout,"Stage 1 index size: %d aa\n",index1part_aa);
#endif
  fprintf(stdout,"Sizes: off_t (%d), size_t (%d), unsigned int (%d), long int (%d), long long int (%d)\n",
	  (int) sizeof(off_t),(int) sizeof(size_t),(int) sizeof(unsigned int),(int) sizeof(long int),(int) sizeof(long long int));
  fprintf(stdout,"Default gmap directory (compiled): %s\n",GMAPDB);
  genomedir = Datadir_find_genomedir(/*user_genomedir*/NULL);
  fprintf(stdout,"Default gmap directory (environment): %s\n",genomedir);
  FREE(genomedir);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

/* This flag is not well-supported, and therefore hidden, but
   kept for backwards compatibility */
/*  -R, --rel=STRING               Release\n\ */

static void
print_program_usage ();


static void
check_compiler_assumptions () {
  unsigned int x = rand(), y = rand();
#ifdef HAVE_SSE2
  int z;
  __m128i a;
#ifdef HAVE_SSE4_1
  char negx, negy;
#endif
#endif


#ifdef HAVE_SSE2
  fprintf(stderr,"Checking compiler assumptions for SSE2: ");
  fprintf(stderr,"%08X %08X",x,y);
  a = _mm_xor_si128(_mm_set1_epi32(x),_mm_set1_epi32(y));
  z = _mm_cvtsi128_si32(a);
  fprintf(stderr," xor=%08X\n",z);
#endif

#ifdef HAVE_SSE4_1
  if ((negx = (char) x) > 0) {
    negx = -negx;
  }
  if ((negy = (char) y) > 0) {
    negy = -negy;
  }

  fprintf(stderr,"Checking compiler assumptions for SSE4.1: ");
  fprintf(stderr,"%d %d",negx,negy);
  a = _mm_max_epi8(_mm_set1_epi8(negx),_mm_set1_epi8(negy));
  z = _mm_extract_epi8(a,0);
  fprintf(stderr," max=%d => ",z);
  if (negx > negy) {
    if (z == (int) negx) {
      fprintf(stderr,"compiler sign extends\n"); /* technically incorrect, but SIMD procedures behave properly */
    } else {
      fprintf(stderr,"compiler zero extends\n");
    }
  } else {
    if (z == (int) negy) {
      fprintf(stderr,"compiler sign extends\n"); /* technically incorrect, but SIMD procedures behave properly */
    } else {
      fprintf(stderr,"compiler zero extends\n");
    }
  }
#endif

#ifdef HAVE_SSE4_2
  fprintf(stderr,"Checking compiler options for SSE4.2: ");
  fprintf(stderr,"%08X ",x);
#ifdef HAVE_LZCNT
  fprintf(stderr,"_lzcnt_u32=%d ",_lzcnt_u32(x));
#endif
#ifdef HAVE_BUILTIN_CLZ
  fprintf(stderr,"__builtin_clz=%d ",__builtin_clz(x));
#endif
#ifdef HAVE_TZCNT
  fprintf(stderr,"_tzcnt_u32=%d ",_tzcnt_u32(x));
#endif
#ifdef HAVE_BUILTIN_CTZ
  fprintf(stderr,"__builtin_ctz=%d ",__builtin_ctz(x));
#endif

#ifdef HAVE_POPCNT
  fprintf(stderr,"_popcnt32=%d ",_popcnt32(x));
#endif
#if defined(HAVE_MM_POPCNT)
  fprintf(stderr,"_mm_popcnt_u32=%d ",_mm_popcnt_u32(x));
#endif
#if defined(HAVE_BUILTIN_POPCOUNT)
  fprintf(stderr,"__builtin_popcount=%d ",__builtin_popcount(x));
#endif
  fprintf(stderr,"\n");

#endif

  fprintf(stderr,"Finished checking compiler assumptions\n");

  return;
}


/************************************************************************/


/* Call before Stage1_compute */
static Diagnostic_T
evaluate_query (bool *poorp, bool *repetitivep, char *queryuc_ptr, int querylength,
		Oligoindex_T oligoindex) {
  Diagnostic_T diagnostic;

  diagnostic = Diagnostic_new();

#ifdef PMAP
  Oligoindex_set_inquery(&diagnostic->query_badoligos,&diagnostic->query_repoligos,
			 &diagnostic->query_trimoligos,&diagnostic->query_trim_start,
			 &diagnostic->query_trim_end,oligoindex,queryuc_ptr,
			 /*querystart*/0,/*queryend*/querylength);
  *poorp = false;
  *repetitivep = false;
#else
  diagnostic->query_oligodepth = 
    Oligoindex_set_inquery(&diagnostic->query_badoligos,&diagnostic->query_repoligos,
			   &diagnostic->query_trimoligos,&diagnostic->query_trim_start,
			   &diagnostic->query_trim_end,oligoindex,queryuc_ptr,
			   /*querystart*/0,/*queryend*/querylength,/*trimp*/true);

  debug2(printf("query_trimoligos %d, fraction badoligos %f = %d/%d, oligodepth %f, fraction repoligos %f = %d/%d\n",
		diagnostic->query_trimoligos,
		(double) diagnostic->query_badoligos/(double) diagnostic->query_trimoligos,
		diagnostic->query_badoligos,diagnostic->query_trimoligos,
		diagnostic->query_oligodepth,
		(double) diagnostic->query_repoligos/(double) diagnostic->query_trimoligos,
		diagnostic->query_repoligos,diagnostic->query_trimoligos));

  if (diagnostic->query_trimoligos == 0) {
    *poorp = true;
  } else if (((double) diagnostic->query_badoligos/(double) diagnostic->query_trimoligos > MAX_BADOLIGOS) ||
	     (diagnostic->query_trim_end - diagnostic->query_trim_start < 80 && diagnostic->query_badoligos > 0)) {
    *poorp = true;
  } else {
    *poorp = false;
  }

  if (diagnostic->query_trimoligos == 0) {
    *repetitivep = false;
  } else if (diagnostic->query_oligodepth > MAX_OLIGODEPTH || 
	     (double) diagnostic->query_repoligos/(double) diagnostic->query_trimoligos > MAX_REPOLIGOS) {
    *repetitivep = true;
  } else {
    *repetitivep = false;
  }
#endif

  return diagnostic;
}




static Stage3_T *
stage3array_from_list (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
		       List_T stage3list, bool chimerap, bool remove_overlaps_p) {
  Stage3_T *array1, *array0, x, y;
  bool *eliminate;
  int norig_primary, norig_altloc, i_primary, i_altloc, i, j;
  int threshold_score;

  Univcoord_T alias_start, alias_end;

  debug(printf("Entering stage3array_from_list with %d entries\n",List_length(stage3list)));

  /* Stage3_recompute_goodness(stage3list); -- No longer necessary */
  Stage3_compute_mapq(stage3list);

  if (stage3list == NULL) {
    *first_absmq = 0;
    *second_absmq = 0;
    *npaths_primary = *npaths_altloc = 0;
    return (Stage3_T *) NULL;

#if 0
  } else if (mergedp == true) {
    debug(printf("mergedp is true\n"));
    Stage3_count_paths(&norig_primary,&norig_altloc,stage3list);
    array0 = (Stage3_T *) List_to_array_out(stage3list,NULL);
    List_free(&stage3list);
    *first_absmq = 0;
    *second_absmq = 0;
    *npaths_primary = norig_primary;
    *npaths_altloc = norig_altloc;
    return array0;
#endif

  } else if (chimerap == true) {
    debug(printf("chimerap is true\n"));
    Stage3_count_paths(&norig_primary,&norig_altloc,stage3list);
    array0 = (Stage3_T *) List_to_array_out(stage3list,NULL);
    List_free(&stage3list);
    *first_absmq = Stage3_absmq_score(array0[0]);
    if (norig_primary + norig_altloc <= 2) {
      *second_absmq = 0;
    } else {
      qsort(&(array0[2]),norig_primary + norig_altloc - 2,sizeof(Stage3_T),Stage3_cmp);
      *second_absmq = Stage3_absmq_score(array0[2]);
    }
    *npaths_primary = norig_primary;
    *npaths_altloc = norig_altloc;
    return array0;

  } else if (remove_overlaps_p == false) {
    debug(printf("remove_overlaps_p is false\n"));
    Stage3_count_paths(&norig_primary,&norig_altloc,stage3list);
    array0 = (Stage3_T *) List_to_array_out(stage3list,NULL);
    List_free(&stage3list);
    qsort(array0,norig_primary + norig_altloc,sizeof(Stage3_T),Stage3_cmp);

    if (suboptimal_score_float < 1.0) {
      threshold_score = Stage3_goodness(array0[0]) * suboptimal_score_float;
      debug(printf("threshold score %d = goodness %d * suboptimal score_float %f\n",
		   threshold_score,Stage3_goodness(array0[0]),suboptimal_score_float));
    } else {
      threshold_score = Stage3_goodness(array0[0]) - (int) suboptimal_score_float;
      debug(printf("threshold score %d = goodness %d - suboptimal score %d\n",
		   threshold_score,Stage3_goodness(array0[0]),(int) suboptimal_score_float));
    }

    if (Stage3_altloc_chr(&alias_start,&alias_end,array0[0]) == false) {
      i_primary = 1;
      i_altloc = 0;
    } else {
      i_primary = 0;
      i_altloc = 1;
    }
    i = 1;
    while (i < norig_primary + norig_altloc && Stage3_goodness(array0[i]) >= threshold_score) {
      if (Stage3_altloc_chr(&alias_start,&alias_end,array0[i]) == false) {
	i_primary++;
      } else {
	i_altloc++;
      }
      i++;
    }
    while (i < norig_primary + norig_altloc) {
      Stage3_free(&(array0[i]));
      i++;
    }

    *npaths_primary = i_primary;
    *npaths_altloc = i_altloc;
    *first_absmq = Stage3_absmq_score(array0[0]);
    if ((*npaths_primary) + (*npaths_altloc) < 2) {
      *second_absmq = 0;
    } else {
      *second_absmq = Stage3_absmq_score(array0[1]);
    }

    return array0;

  } else {
    debug(printf("remove_overlaps_p is true\n"));
    Stage3_count_paths(&norig_primary,&norig_altloc,stage3list);
    eliminate = (bool *) CALLOCA(norig_primary + norig_altloc,sizeof(bool));

    /* Initial sort to remove subsumed alignments */
    array0 = (Stage3_T *) MALLOCA((norig_primary + norig_altloc) * sizeof(Stage3_T));
    List_fill_array_and_free((void **) array0,&stage3list);
    qsort(array0,norig_primary + norig_altloc,sizeof(Stage3_T),Stage3_cmp);

    for (i = 0; i < norig_primary + norig_altloc; i++) {
      x = array0[i];
      debug(printf("%d: chr %d:%u..%u, goodness %d, matches %d, npairs %d\n",
		   i,Stage3_chrnum(x),Stage3_chrstart(x),Stage3_chrend(x),Stage3_goodness(x),Stage3_matches(x),Stage3_npairs(x)));
      for (j = i+1; j < norig_primary + norig_altloc; j++) {
	y = array0[j];
	if (Stage3_overlap(x,y)) {
	  eliminate[j] = true;
	}
      }
    }


    *npaths_primary = *npaths_altloc = 0;
    for (i = 0; i < norig_primary + norig_altloc; i++) {
      if (eliminate[i] == false) {
	if (Stage3_altloc_chr(&alias_start,&alias_end,array0[i]) == false) {
	  (*npaths_primary)++;
	} else {
	  (*npaths_altloc)++;
	}
      }
    }

    array1 = (Stage3_T *) MALLOC_OUT(((*npaths_primary) + (*npaths_altloc)) * sizeof(Stage3_T)); /* Return value */
    j = 0;
    for (i = 0; i < norig_primary + norig_altloc; i++) {
      x = array0[i];
      if (eliminate[i] == true) {
	Stage3_free(&x);
      } else {
	array1[j++] = x;
      }
    }
    FREEA(array0);
    FREEA(eliminate);

    if (suboptimal_score_float < 1.0) {
      threshold_score = Stage3_goodness(array1[0]) * suboptimal_score_float;
      debug(printf("threshold score %d = goodness %d * suboptimal score %f\n",
		   threshold_score,Stage3_goodness(array1[0]),suboptimal_score_float));
    } else {
      threshold_score = Stage3_goodness(array1[0]) - (int) suboptimal_score_float;
      debug(printf("threshold score %d = goodness %d - suboptimal score %d\n",
		   threshold_score,Stage3_goodness(array1[0]),(int) suboptimal_score_float));
    }

    if (Stage3_altloc_chr(&alias_start,&alias_end,array1[0]) == false) {
      i_primary = 1;
      i_altloc = 0;
    } else {
      i_primary = 0;
      i_altloc = 1;
    }
    i = 1;
    while (i < (*npaths_primary) + (*npaths_altloc) && Stage3_goodness(array1[i]) >= threshold_score) {
      if (Stage3_altloc_chr(&alias_start,&alias_end,array1[i]) == false) {
	i_primary++;
      } else {
	i_altloc++;
      }
      i++;
    }
    while (i < (*npaths_primary) + (*npaths_altloc)) {
      Stage3_free(&(array1[i]));
      i++;
    }

    *npaths_primary = i_primary;
    *npaths_altloc = i_altloc;
    *first_absmq = Stage3_absmq_score(array1[0]);
    if ((*npaths_primary) + (*npaths_altloc) < 2) {
      *second_absmq = 0;
    } else {
      *second_absmq = Stage3_absmq_score(array1[1]);
    }
    return array1;
  }
}


static List_T
update_stage3middle_list (List_T stage3middle_list, Sequence_T queryseq,
#ifdef PMAP
			  Sequence_T queryntseq,
#endif
			  Sequence_T queryuc, Stage2_alloc_T stage2_alloc,
			  Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
			  Genome_T genome, Genome_T genomealt,
			  Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
			  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
			  Chrpos_T chrstart, Chrpos_T chrend, bool watsonp, int genestrand,
			  Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
			  Stopwatch_T worker_stopwatch) {
  /* int stage2_source, stage2_indexsize; */
  /* double stage3_runtime; */

#ifdef PMAP
  Sequence_T genomicuc = NULL;
  char *genomicseg_ptr = NULL, *genomicuc_ptr = NULL;
#elif defined(EXTRACT_GENOMICSEG)
  Sequence_T genomicuc = NULL;
#endif

  List_T all_stage2results, all_stage3middle_results = NULL, p;
  Stage2_T stage2;
  Stage3middle_T stage3middle;
#ifdef PMAP
  int subseq_offset;
#endif


#ifdef PMAP_OLD
  /* Previously used for PMAP */
  if (user_genomicsegs == NULL && uncompressedp == false && straintype == 0) {
    genomicuc = Sequence_alias(genomicseg);
  } else {
    genomicuc = Sequence_uppercase(genomicseg);
  }
  genomicseg_ptr = Sequence_fullpointer(genomicseg);
  genomicuc_ptr = Sequence_fullpointer(genomicuc);
#elif defined(EXTRACT_GENOMICSEG)
  if (user_genomicsegs == NULL && uncompressedp == false && straintype == 0) {
    genomicuc = Sequence_alias(genomicseg);
  } else {
    genomicuc = Sequence_uppercase(genomicseg);
  }
  genomicseg_ptr = Sequence_fullpointer(genomicseg);
  genomicuc_ptr = Sequence_fullpointer(genomicuc);
#endif

#if 0
  if (canonical_mode == 0) {
    do_final_p = false;
  } else if (canonical_mode == 1) {
    do_final_p = true;
  } else if (lowidentityp == false) {
    do_final_p = false;
  } else {
    do_final_p = true;
  }
#endif

  debug(printf("Entering update_stage3middle_list with %d results\n",List_length(stage3middle_list)));
  debug2(printf("Beginning Stage2_compute with chrstart %u and chrend %u and query_subseq_offset %d\n",
		chrstart,chrend,Sequence_subseq_offset(queryseq)));
  all_stage2results = Stage2_compute(Sequence_trimpointer(queryseq),Sequence_trimpointer(queryuc),
				     Sequence_trimlength(queryseq),/*query_offset*/0,
				     chrstart,chrend,chroffset,chrhigh,/*plusp*/watsonp,genestrand,
				     stage2_alloc,/*proceed_pctcoverage*/0.3,oligoindices_major,
				     genome,genomealt,pairpool,diagpool,cellpool,
				     /*localp*/true,/*skip_repetitive_p*/true,
				     /*favor_right_p*/false,/*max_nalignments*/MAX_NALIGNMENTS,
				     worker_stopwatch,diag_debug);
  debug(printf("End of Stage2_compute\n"));


  for (p = all_stage2results; p != NULL; p = List_next(p)) {
    stage2 = (Stage2_T) List_head(p);
    stage3middle = Stage3_compute_middle(Stage2_middle(stage2),Stage2_all_starts(stage2),Stage2_all_ends(stage2),
#ifdef PMAP
					 /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
					 /*queryseq_ptr*/Sequence_subseq_pointer(queryntseq,subseq_offset),
					 /*queryuc_ptr*/Sequence_subseq_pointer(queryntseq,subseq_offset),
					 /*querylength*/Sequence_subseq_length(queryntseq,subseq_offset),
#else
					 /*queryseq_ptr*/Sequence_fullpointer(queryseq),
					 /*queryuc_ptr*/Sequence_fullpointer(queryuc),
					 /*querylength*/Sequence_fulllength(queryseq),
#endif
					 chrnum,chroffset,chrhigh,chrlength,
					 watsonp,genestrand,/*jump_late_p*/watsonp ? false : true,maxpeelback,
					 oligoindices_minor,diagpool,cellpool,
					 genome,genomealt,pairpool,dynprogL,dynprogM,dynprogR,sense_try);
    Stage2_free(&stage2);
    all_stage3middle_results = List_push(all_stage3middle_results,(void *) stage3middle);
  }
  List_free(&all_stage2results);

  return List_append(all_stage3middle_results,stage3middle_list);
}



/* Combination of update_stage3middle_list and Stage3_compute_ends,
   Needed for solving middle segments of chimeras */
static List_T
update_stage3list (List_T stage3list, Sequence_T queryseq,
#ifdef PMAP
		   Sequence_T queryntseq,
#endif
		   Sequence_T queryuc, Stage2_alloc_T stage2_alloc,
		   Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
		   Genome_T genome, Genome_T genomealt, Pairpool_T pairpool,
		   Diagpool_T diagpool, Cellpool_T cellpool, int straintype, char *strain,
		   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, Chrpos_T chrlength,
		   Chrpos_T chrstart, Chrpos_T chrend, bool watsonp, int genestrand,
		   Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		   int min_matches, Stopwatch_T worker_stopwatch) {
  /* int stage2_source, stage2_indexsize; */
  /* double stage3_runtime; */

#ifdef PMAP
  Sequence_T genomicuc = NULL;
  char *genomicseg_ptr = NULL, *genomicuc_ptr = NULL;
#elif defined(EXTRACT_GENOMICSEG)
  Sequence_T genomicuc = NULL;
#endif
  List_T all_stage2results, p;
  Stage2_T stage2;
  Stage3_T stage3;

  struct Pair_T *pairarray;
  List_T pairs;
  int goodness;
  int npairs, cdna_direction, matches, unknowns, mismatches, qopens, qindels, topens, tindels,
    ncanonical, nsemicanonical, nnoncanonical;
  int sensedir;
  int nmatches_posttrim, max_match_length, ambig_end_length_5, ambig_end_length_3;
  Splicetype_T ambig_splicetype_5, ambig_splicetype_3;
  double ambig_prob_5, ambig_prob_3;
  double min_splice_prob;
#ifdef PMAP
  int subseq_offset;
#endif


#ifdef PMAP_OLD
  /* Previously used for PMAP */
  if (user_genomicsegs == NULL && uncompressedp == false && straintype == 0) {
    genomicuc = Sequence_alias(genomicseg);
  } else {
    genomicuc = Sequence_uppercase(genomicseg);
  }
  genomicseg_ptr = Sequence_fullpointer(genomicseg);
  genomicuc_ptr = Sequence_fullpointer(genomicuc);
#elif defined(EXTRACT_GENOMICSEG)
  if (user_genomicsegs == NULL && uncompressedp == false && straintype == 0) {
    genomicuc = Sequence_alias(genomicseg);
  } else {
    genomicuc = Sequence_uppercase(genomicseg);
  }
  genomicseg_ptr = Sequence_fullpointer(genomicseg);
  genomicuc_ptr = Sequence_fullpointer(genomicuc);
#endif

#if 0
  if (canonical_mode == 0) {
    do_final_p = false;
  } else if (canonical_mode == 1) {
    do_final_p = true;
  } else if (lowidentityp == false) {
    do_final_p = false;
  } else {
    do_final_p = true;
  }
#endif

  debug(printf("Entering update_stage3list with %d results\n",List_length(stage3list)));
  debug2(printf("Beginning Stage2_compute with chrstart %u and chrend %u and query_subseq_offset %d\n",
		chrstart,chrend,Sequence_subseq_offset(queryseq)));
  all_stage2results = Stage2_compute(Sequence_trimpointer(queryseq),Sequence_trimpointer(queryuc),
				     Sequence_trimlength(queryseq),/*query_offset*/0,
				     chrstart,chrend,chroffset,chrhigh,/*plusp*/watsonp,genestrand,
				     stage2_alloc,/*proceed_pctcoverage*/0.3,oligoindices_major,
				     genome,genomealt,pairpool,diagpool,cellpool,
				     /*localp*/true,/*skip_repetitive_p*/true,
				     /*favor_right_p*/false,/*max_nalignments*/MAX_NALIGNMENTS,
				     worker_stopwatch,diag_debug);

  debug(printf("End of Stage2_compute\n"));

  for (p = all_stage2results; p != NULL; p = List_next(p)) {
    stage2 = (Stage2_T) List_head(p);

    /* Stopwatch_start(worker_stopwatch); */
#ifdef PMAP
    subseq_offset = Sequence_subseq_offset(queryseq); /* in nucleotides */
#endif
    pairarray = Stage3_compute_one(&cdna_direction,&sensedir,&pairs,&npairs,&goodness,
				   &matches,&nmatches_posttrim,&max_match_length,
				   &ambig_end_length_5,&ambig_end_length_3,
				   &ambig_splicetype_5,&ambig_splicetype_3,
				   &ambig_prob_5,&ambig_prob_3,
				   &unknowns,&mismatches,&qopens,&qindels,&topens,&tindels,
				   &ncanonical,&nsemicanonical,&nnoncanonical,&min_splice_prob,
				   Stage2_middle(stage2),Stage2_all_starts(stage2),Stage2_all_ends(stage2),
#ifdef PMAP
				   /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
				   /*queryseq_ptr*/Sequence_subseq_pointer(queryntseq,subseq_offset),
				   /*queryuc_ptr*/Sequence_subseq_pointer(queryntseq,subseq_offset),
				   /*querylength*/Sequence_subseq_length(queryntseq,subseq_offset),
				   /*skiplength*/Sequence_skiplength(queryntseq),
				   /*query_subseq_offset*/subseq_offset,
#else
				   /*queryseq_ptr*/Sequence_fullpointer(queryseq),
				   /*queryuc_ptr*/Sequence_fullpointer(queryuc),
				   /*querylength*/Sequence_fulllength(queryseq),
				   /*skiplength*/Sequence_skiplength(queryseq),
				   /*query_subseq_offset*/Sequence_subseq_offset(queryseq),
#endif
				   chrnum,chroffset,chrhigh,
				   /*knownsplice_limit_low*/0U,/*knownsplice_limit_high*/-1U,
				   watsonp,genestrand,/*jump_late_p*/watsonp ? false : true,maxpeelback,
				   oligoindices_minor,diagpool,cellpool,
				   genome,genomealt,pairpool,dynprogL,dynprogM,dynprogR,sense_try,sense_filter);
    /* stage3_runtime = Stopwatch_stop(worker_stopwatch); */
    if (pairarray == NULL) {
      /* Skip */
    } else if (matches < min_matches) {
      FREE_OUT(pairarray);
    } else if ((stage3 = Stage3_new(pairarray,pairs,npairs,goodness,cdna_direction,sensedir,
				    matches,unknowns,mismatches,
				    qopens,qindels,topens,tindels,ncanonical,nsemicanonical,nnoncanonical,
				    genome,genomealt,chrnum,chroffset,chrhigh,chrlength,watsonp,genestrand,
				    /*querylength*/Sequence_fulllength(queryseq),
				    /*skiplength*/Sequence_skiplength(queryseq),
				    /*trimlength*/Sequence_trimlength(queryseq),
				    straintype,strain,altstrain_iit)) != NULL) {
      debug(printf("Pushing %p onto stage3list\n",stage3));
      stage3list = List_push(stage3list,(void *) stage3);
    }

    Stage2_free(&stage2);
  }

  List_free(&all_stage2results);

#ifdef PMAP_OLD
  Sequence_free(&genomicuc);
#elif defined(EXTRACT_GENOMICSEG)
  Sequence_free(&genomicuc);
#endif

  return stage3list;
}



#if 0
/* This code is duplicated in get-genome.c */
static int
index_compare (const void *a, const void *b) {
  int index1 = * (int *) a;
  int index2 = * (int *) b;
  int type1, type2;
  Chrpos_T pos1, pos2;

  type1 = Interval_type(IIT_interval(altstrain_iit,index1));
  type2 = Interval_type(IIT_interval(altstrain_iit,index2));
  
  if (type1 < type2) {
    return -1;
  } else if (type1 > type2) {
    return +1;
  } else {
    /* Store in descending genomic position, so right shifting works
       in Genome_patch_strain */
    pos1 = Interval_low(IIT_interval(altstrain_iit,index1));
    pos2 = Interval_low(IIT_interval(altstrain_iit,index2));

    if (pos1 > pos2) {
      return -1;
    } else if (pos1 < pos2) {
      return +1;
    } else {
      return 0;
    }
  }
}
#endif


static Stage3_T *
stage3_self_align (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
		   Sequence_T queryseq, Pairpool_T pairpool) {
  List_T stage3list = NULL;
  List_T pairs;
  Stage3_T stage3;
  struct Pair_T *pairarray;
  int npairs;
  
  pairarray = Stage3_compute_self_align(&pairs,&npairs,/*queryseq_ptr*/Sequence_fullpointer(queryseq),
					pairpool);

  if ((stage3 = Stage3_new(pairarray,pairs,npairs,/*goodness*/0,/*cdna_direction*/0,/*sensedir*/true,
			   /*matches*/npairs,/*unknowns*/0,/*mismatches*/0,
			   /*qopens*/0,/*qindels*/0,/*topens*/0,/*tindels*/0,
			   /*ncanonical*/0,/*nsemicanonical*/0,/*nnoncanonical*/0,
			   /*genome*/NULL,/*genomealt*/NULL,
			   /*chrnum*/0,/*chroffset*/0,/*chrhigh*/0,/*chrlength*/0,
			   /*watsonp*/0,/*genestrand*/0,/*querylength*/npairs,
			   /*skiplength*/0,/*trimlength*/npairs,
			   /*straintype*/0,/*strain*/NULL,altstrain_iit)) == NULL) {
    *npaths_primary = *npaths_altloc = 0;
    return (Stage3_T *) NULL;
  } else {
    stage3list = List_push(stage3list,(void *) stage3);
    return stage3array_from_list(&(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),
				 stage3list,/*chimerap*/false,/*remove_overlaps_p*/true);
  }
}


/* Not sure how to treat genestrand for this case */
static Stage3_T *
stage3_skip_stage1 (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
		    Sequence_T queryseq, Sequence_T queryuc,
#ifdef PMAP
		    Sequence_T queryntseq,
#endif
		    Stage2_alloc_T stage2_alloc,
		    Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
		    Genome_T genome, Genome_T genomealt, Pairpool_T pairpool,
		    Diagpool_T diagpool, Cellpool_T cellpool,
		    Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		    Stopwatch_T worker_stopwatch) {
  List_T stage3list, stage3middle_list, p;
  Stage3middle_T stage3middle;
  Stage3_T stage3;
  bool watsonp;

  struct Pair_T *pairarray;
  List_T pairs;
  int goodness;
  int npairs, cdna_direction, matches, unknowns, mismatches, qopens, qindels, topens, tindels,
    ncanonical, nsemicanonical, nnoncanonical;
  int sensedir;
  int nmatches_posttrim, max_match_length, ambig_end_length_5, ambig_end_length_3;
  Splicetype_T ambig_splicetype_5, ambig_splicetype_3;
  double ambig_prob_5, ambig_prob_3;
  double min_splice_prob;
  int min_matches;
#ifdef PMAP
  int subseq_offset;
#endif

  Univcoord_T chroffset, chrhigh;
  Chrpos_T chrlength, chrpos;
  Chrnum_T chrnum = 0;

#ifdef PMAP
  Sequence_T revcomp;
#endif
		    
  chroffset = chrpos = 0U;
  chrhigh = chrlength = Genome_genomelength(genome);
  if ((min_matches = Genome_genomelength(genome)/2) > MIN_MATCHES) {
    min_matches = MIN_MATCHES;
  }

  if (cdna_direction_try >= 0) {
    stage3middle_list = update_stage3middle_list(/*stage3middle_list*/NULL,queryseq,
#ifdef PMAP
						 queryntseq,
#endif
						 queryuc,stage2_alloc,oligoindices_major,oligoindices_minor,
						 genome,genomealt,pairpool,diagpool,cellpool,
						 chrnum,chroffset,chrhigh,chrlength,
						 /*chrstart*/0,/*chrend*/chrhigh,/*watsonp*/true,
						 /*genestrand for usersegment*/0,
						 dynprogL,dynprogM,dynprogR,worker_stopwatch);
  }

#ifdef PMAP
  revcomp = Sequence_revcomp(usersegment);
#endif

  if (cdna_direction_try <= 0) {
    stage3middle_list = update_stage3middle_list(stage3middle_list,queryseq,
#ifdef PMAP
						 queryntseq,
#endif
						 queryuc,stage2_alloc,oligoindices_major,oligoindices_minor,
						 genome,genomealt,pairpool,diagpool,cellpool,
						 chrnum,chroffset,chrhigh,chrlength,
						 /*chrstart*/0,/*chrend*/chrhigh,/*watsonp*/false,
						 /*genestrand for usersegment*/0,
						 dynprogL,dynprogM,dynprogR,worker_stopwatch);
  }

#ifdef PMAP
  Sequence_free(&revcomp);
#endif

  if (stage3middle_list == NULL) {
    *npaths_primary = *npaths_altloc = 0;
    return (Stage3_T *) NULL;

  } else {
    stage3list = (List_T) NULL;
    for (p = stage3middle_list; p != NULL; p = List_next(p)) {
      stage3middle = (Stage3middle_T) List_head(p);
      watsonp = Stage3middle_watsonp(stage3middle);

#ifdef PMAP
      subseq_offset = Sequence_subseq_offset(queryseq); /* in nucleotides */
#endif
      pairarray = Stage3_compute_ends(&cdna_direction,&sensedir,&pairs,&npairs,&goodness,
				      &matches,&nmatches_posttrim,&max_match_length,
				      &ambig_end_length_5,&ambig_end_length_3,
				      &ambig_splicetype_5,&ambig_splicetype_3,
				      &ambig_prob_5,&ambig_prob_3,
				      &unknowns,&mismatches,&qopens,&qindels,&topens,&tindels,
				      &ncanonical,&nsemicanonical,&nnoncanonical,&min_splice_prob,
				      stage3middle,
#ifdef PMAP
				      /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
				      /*queryseq_ptr*/Sequence_subseq_pointer(queryntseq,subseq_offset),
				      /*queryuc_ptr*/Sequence_subseq_pointer(queryntseq,subseq_offset),
				      /*querylength*/Sequence_subseq_length(queryntseq,subseq_offset),
				      /*skiplength*/Sequence_skiplength(queryntseq),
				      /*query_subseq_offset*/subseq_offset,
#else
				      /*queryseq_ptr*/Sequence_fullpointer(queryseq),
				      /*queryuc_ptr*/Sequence_fullpointer(queryuc),
				      /*querylength*/Sequence_fulllength(queryseq),
				      /*skiplength*/Sequence_skiplength(queryseq),
				      /*query_subseq_offset*/Sequence_subseq_offset(queryseq),
#endif
				      /*knownsplice_limit_low*/0U,/*knownsplice_limit_high*/-1U,
				      maxpeelback,genome,genomealt,pairpool,dynprogL,dynprogM,dynprogR,
				      sense_filter,oligoindices_minor,diagpool,cellpool);
      /* stage3_runtime = Stopwatch_stop(worker_stopwatch); */
      if (pairarray == NULL) {
	/* Skip */
      } else if (matches < min_matches) {
	FREE_OUT(pairarray);
      } else if ((stage3 = Stage3_new(pairarray,pairs,npairs,goodness,cdna_direction,sensedir,
				      matches,unknowns,mismatches,
				      qopens,qindels,topens,tindels,ncanonical,nsemicanonical,nnoncanonical,
				      genome,genomealt,chrnum,chroffset,chrhigh,chrlength,
				      watsonp,/*genestrand for usersegment*/0,
				      /*querylength*/Sequence_fulllength(queryseq),
				      /*skiplength*/Sequence_skiplength(queryseq),
				      /*trimlength*/Sequence_trimlength(queryseq),
				      /*straintype*/0,/*strain*/NULL,altstrain_iit)) != NULL) {
	debug(printf("Pushing %p onto stage3list\n",stage3));
	stage3list = List_push(stage3list,(void *) stage3);
      }
      Stage3middle_free(&stage3middle);
    }
    List_free(&stage3middle_list);

    return stage3array_from_list(&(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),
				 stage3list,/*chimerap*/false,/*remove_overlaps_p*/true);
  }
}


#if 0
static List_T
stage3list_remove_duplicates (List_T stage3list) {
  List_T unique = NULL;
  Stage3_T *array;
  int best_score;
  Chrpos_T shortest_genomiclength;
  int n, besti, i, j, k;
  
  if ((n = List_length(stage3list)) == 0) {
    return (List_T) NULL;
  } else if (n == 1) {
    return stage3list;
  } else {
    array = (Stage3_T *) List_to_array(stage3list,NULL);
    List_free(&stage3list);
    qsort(array,n,sizeof(Stage3_T),Stage3_position_cmp);

    i = 0;
    while (i < n) {
      best_score = Stage3_goodness(array[i]);
      shortest_genomiclength = Stage3_genomiclength(array[i]);
      besti = i;
      debug3(printf("i = %d, score %d, genomiclength %u\n",
		    i,best_score,shortest_genomiclength));

      j = i + 1;
      while (j < n && Stage3_position_cmp(&(array[i]),&(array[j])) == 0) {
	debug3(printf("  j = %d, score %d, genomiclength %u\n",
		      j,Stage3_goodness(array[j]),Stage3_genomiclength(array[j])));

	if (Stage3_goodness(array[j]) < best_score) {
	  best_score = Stage3_goodness(array[j]);
	  shortest_genomiclength = Stage3_genomiclength(array[j]);
	  besti = j;

	} else if (Stage3_goodness(array[j]) == best_score &&
		   Stage3_genomiclength(array[j]) < shortest_genomiclength) {
	  best_score = Stage3_goodness(array[j]);
	  shortest_genomiclength = Stage3_genomiclength(array[j]);
	  besti = j;
	}

	j++;
      }
      debug3(printf("  => besti = %d, score %d, genomiclength %u\n",
		    besti,best_score,shortest_genomiclength));

      for (k = i; k < j; k++) {
	if (k == besti) {
	  unique = List_push(unique,(void *) array[besti]);
	} else {
	  Stage3_free(&(array[k]));
	}
      }

      i = j;
    }
    
    FREE(array);

    return unique;
  }
}
#endif


#if 0
static List_T
stage3list_remove_empties (List_T stage3list) {
  List_T nonempty = NULL, p;
  Stage3_T stage3;
  
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    if (Stage3_pairs == NULL) {
      debug2(printf("Removing empty stage3 %p\n",stage3));
      Stage3_free(&stage3);
    } else {
      nonempty = List_push(nonempty,(void *) stage3);
    }
  }

  return nonempty;
}
#endif


static List_T
stage3list_sort (List_T stage3list) {
  List_T sorted = NULL;
  Stage3_T *array;
  int n, i;

  if ((n = List_length(stage3list)) == 0) {
    return (List_T) NULL;
  } else if (n == 1) {
    return stage3list;
  } else {
    array = (Stage3_T *) List_to_array(stage3list,NULL);
    List_free(&stage3list);
    qsort(array,n,sizeof(Stage3_T),Stage3_cmp);
    for (i = n-1; i >= 0; i--) {
      sorted = List_push(sorted,(void *) array[i]);
    }
    FREE(array);

    return sorted;
  }
}


static List_T
stage3list_filter_and_sort (Chimera_T *chimera, List_T stage3list) {
  List_T sorted = NULL;
  Stage3_T *array, stage3;
  int n, i;

  if ((n = List_length(stage3list)) == 0) {
    return (List_T) NULL;

  } else if (n == 1) {
    stage3 = (Stage3_T) List_head(stage3list);
    if (Stage3_passes_filter(stage3,min_trimmed_coverage,min_identity) == false) {
      Stage3_free(&stage3);
      List_free(&stage3list);
      return (List_T) NULL;
    } else {
      return stage3list;
    }

  } else if (*chimera == NULL) {
    array = (Stage3_T *) List_to_array(stage3list,NULL);
    List_free(&stage3list);
    qsort(array,n,sizeof(Stage3_T),Stage3_cmp);
    for (i = n-1; i >= 0; i--) {
      if (Stage3_passes_filter(array[i],min_trimmed_coverage,min_identity) == false) {
	Stage3_free(&(array[i]));
      } else {
	sorted = List_push(sorted,(void *) array[i]);
      }
    }
    FREE(array);
    return sorted;

  } else if (Stage3_passes_filter_chimera(*chimera,min_trimmed_coverage,min_identity) == true) {
    array = (Stage3_T *) List_to_array(stage3list,NULL);
    List_free(&stage3list);
    qsort(array,n,sizeof(Stage3_T),Stage3_cmp);
    for (i = n-1; i >= 0; i--) {
      if (Stage3_chimera_left_p(array[i]) == true) {
	sorted = List_push(sorted,(void *) array[i]);
      } else if (Stage3_chimera_right_p(array[i]) == true) {
	sorted = List_push(sorted,(void *) array[i]);
      } else if (Stage3_passes_filter(array[i],min_trimmed_coverage,min_identity) == false) {
	Stage3_free(&(array[i]));
      } else {
	sorted = List_push(sorted,(void *) array[i]);
      }
    }
    FREE(array);
    return sorted;

  } else {
    array = (Stage3_T *) List_to_array(stage3list,NULL);
    List_free(&stage3list);
    qsort(array,n,sizeof(Stage3_T),Stage3_cmp);
    for (i = n-1; i >= 0; i--) {
      if (Stage3_chimera_left_p(array[i]) == true) {
	Stage3_free(&(array[i]));
      } else if (Stage3_chimera_right_p(array[i]) == true) {
	Stage3_free(&(array[i]));
      } else if (Stage3_passes_filter(array[i],min_trimmed_coverage,min_identity) == false) {
	Stage3_free(&(array[i]));
      } else {
	sorted = List_push(sorted,(void *) array[i]);
      }
    }
    FREE(array);

    Chimera_free(&(*chimera));
    *chimera = (Chimera_T) NULL;

    return sorted;
  }
}


/* Each gregion has its own genestrand */
static List_T
stage3_from_gregions (List_T stage3list, List_T gregions,
		      Sequence_T queryseq, Sequence_T queryuc,
#ifdef PMAP
		      Sequence_T queryntseq,
#endif
		      Stage2_alloc_T stage2_alloc,
		      Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
		      Genome_T genome, Genome_T genomealt, Pairpool_T pairpool,
		      Diagpool_T diagpool, Cellpool_T cellpool,
		      Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		      int min_matches, Stopwatch_T worker_stopwatch) {
  Gregion_T gregion, *gregion_array;
  int ngregions, ncovered, max_ncovered, stage2_source;
  int n, i;

  List_T stage3middle_list = NULL;
  Stage3middle_T stage3middle, *stage3middle_array;
  Stage3_T stage3;
  bool watsonp;
  int genestrand;

  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh;
  Chrpos_T chrlength;

  struct Pair_T *pairarray;
  List_T pairs;
  int goodness, best_score;
  int npairs, cdna_direction, matches, unknowns, mismatches, qopens, qindels, topens, tindels,
    ncanonical, nsemicanonical, nnoncanonical;
  int sensedir;
  int nmatches_posttrim, max_match_length, ambig_end_length_5, ambig_end_length_3;
  Splicetype_T ambig_splicetype_5, ambig_splicetype_3;
  double ambig_prob_5, ambig_prob_3;
  double min_splice_prob;
#ifdef PMAP
  int subseq_offset;
#endif

#if 0
  int *indexarray, nindices, straintype, j;
#endif
  void *item;

#ifdef EXTRACT_GENOMICSEG
  genomicuc_ptr = Sequence_fullpointer(genomicuc);
  Sequence_T genomicseg = NULL, genomicuc = NULL;
#endif
		    
  if ((ngregions = List_length(gregions)) > 0) {
    gregion_array = (Gregion_T *) List_to_array(gregions,NULL);
    List_free(&gregions);

    for (i = 0; i < ngregions; i++) {
      gregion = gregion_array[i];

#if defined(EXTRACT_GENOMICSEG)
      genomicseg = Genome_get_segment(genome,Gregion_genomicstart(gregion),Gregion_genomiclength(gregion),
				      /*chromosome_iit*/NULL,Gregion_revcompp(gregion));
      genomicuc = Sequence_uppercase(genomicseg);
      genomicuc_ptr = Sequence_fullpointer(genomicuc);
#endif
      ncovered = Stage2_scan(&stage2_source,Sequence_trimpointer(queryuc),Sequence_trimlength(queryseq),
			     genome,Gregion_chrstart(gregion),Gregion_chrend(gregion),
			     Gregion_chroffset(gregion),Gregion_chrhigh(gregion),
			     /*plusp*/Gregion_revcompp(gregion) ? false : true,Gregion_genestrand(gregion),
			     stage2_alloc,oligoindices_major,diagpool);
      Gregion_set_ncovered(gregion,ncovered,stage2_source);
#if defined(EXTRACT_GENOMICSEG)
      Sequence_free(&genomicuc);
      Sequence_free(&genomicseg);
#endif
    }
    qsort(gregion_array,ngregions,sizeof(Gregion_T),Gregion_cmp);
    max_ncovered = Gregion_ncovered(gregion_array[0]);
    debug(printf("max_ncovered of gregion_array[0] = %d\n",max_ncovered));
    if (max_ncovered < 0.10*Sequence_fulllength(queryseq)) {
      debug(printf("coverage is too short, so skipping\n"));
      for (i = 0; i < ngregions; i++) {
	Gregion_free(&(gregion_array[i]));
      }
      FREE(gregion_array);

    } else {
      gregions = (List_T) NULL;
      i = 0;
      while (i < ngregions && Gregion_ncovered(gregion_array[i]) > 0.25*max_ncovered) {
	debug(printf("Keeping %d ncovered relative to %d\n",Gregion_ncovered(gregion_array[i]),max_ncovered));
	gregions = List_push(gregions,(void *) gregion_array[i]);
	i++;
      }
      while (i < ngregions) {
	debug(printf("Discarding array %d with ncovered = %d\n",i,Gregion_ncovered(gregion_array[i])));
	Gregion_free(&(gregion_array[i]));
	i++;
      }
      FREE(gregion_array);
    }

    while (gregions != NULL) {
      gregions = List_pop(gregions,&item);
      gregion = (Gregion_T) item;

      /* if (Match_usep(match) == true) { */
      stage3middle_list = update_stage3middle_list(stage3middle_list,queryseq,
#ifdef PMAP
						   queryntseq,
#endif
						   queryuc,stage2_alloc,oligoindices_major,oligoindices_minor,
						   genome,genomealt,pairpool,diagpool,cellpool,Gregion_chrnum(gregion),
						   Gregion_chroffset(gregion),Gregion_chrhigh(gregion),Gregion_chrlength(gregion),
						   Gregion_chrstart(gregion),Gregion_chrend(gregion),
						   Gregion_plusp(gregion),Gregion_genestrand(gregion),
						   dynprogL,dynprogM,dynprogR,worker_stopwatch);
      Gregion_free(&gregion);
    }

    if (stage3middle_list != NULL) {
      stage3middle_array = (Stage3middle_T *) List_to_array_n(&n,stage3middle_list);
      qsort(stage3middle_array,n,sizeof(Stage3middle_T),Stage3middle_cmp);
      List_free(&stage3middle_list);

      best_score = Stage3middle_goodness(stage3middle_array[0]);
      i = 0;
    
      while (i < n && Stage3middle_goodness(stage3middle_array[i]) > best_score - 20) {
	stage3middle = stage3middle_array[i];
	debug(printf("Processing stage3middle %d with goodness %d\n",i,Stage3middle_goodness(stage3middle)));

	chrnum = Stage3middle_chrnum(stage3middle);
	chroffset = Stage3middle_chroffset(stage3middle);
	chrhigh = Stage3middle_chrhigh(stage3middle);
	chrlength = Stage3middle_chrlength(stage3middle);
	watsonp = Stage3middle_watsonp(stage3middle);
	genestrand = Stage3middle_genestrand(stage3middle);

#ifdef PMAP
	subseq_offset = Sequence_subseq_offset(queryseq); /* in nucleotides */
#endif
	pairarray = Stage3_compute_ends(&cdna_direction,&sensedir,&pairs,&npairs,&goodness,
					&matches,&nmatches_posttrim,&max_match_length,
					&ambig_end_length_5,&ambig_end_length_3,
					&ambig_splicetype_5,&ambig_splicetype_3,
					&ambig_prob_5,&ambig_prob_3,
					&unknowns,&mismatches,&qopens,&qindels,&topens,&tindels,
					&ncanonical,&nsemicanonical,&nnoncanonical,&min_splice_prob,
					stage3middle,
#ifdef PMAP
					/*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
					/*queryseq_ptr*/Sequence_subseq_pointer(queryntseq,subseq_offset),
					/*queryuc_ptr*/Sequence_subseq_pointer(queryntseq,subseq_offset),
					/*querylength*/Sequence_subseq_length(queryntseq,subseq_offset),
					/*skiplength*/Sequence_skiplength(queryntseq),
					/*query_subseq_offset*/subseq_offset,
#else
					/*queryseq_ptr*/Sequence_fullpointer(queryseq),
					/*queryuc_ptr*/Sequence_fullpointer(queryuc),
					/*querylength*/Sequence_fulllength(queryseq),
					/*skiplength*/Sequence_skiplength(queryseq),
					/*query_subseq_offset*/Sequence_subseq_offset(queryseq),
#endif
					/*knownsplice_limit_low*/0U,/*knownsplice_limit_high*/-1U,
					maxpeelback,genome,genomealt,pairpool,dynprogL,dynprogM,dynprogR,
					sense_filter,oligoindices_minor,diagpool,cellpool);
	/* stage3_runtime = Stopwatch_stop(worker_stopwatch); */
	if (pairarray == NULL) {
	  /* Skip */
	} else if (matches < min_matches) {
	  FREE_OUT(pairarray);
	} else if ((stage3 = Stage3_new(pairarray,pairs,npairs,goodness,cdna_direction,sensedir,
					matches,unknowns,mismatches,
					qopens,qindels,topens,tindels,ncanonical,nsemicanonical,nnoncanonical,
					genome,genomealt,chrnum,chroffset,chrhigh,chrlength,watsonp,genestrand,
					/*querylength*/Sequence_fulllength(queryseq),
					/*skiplength*/Sequence_skiplength(queryseq),
					/*trimlength*/Sequence_trimlength(queryseq),
					/*straintype*/0,/*strain*/NULL,altstrain_iit)) != NULL) {
	  debug(printf("Pushing %p onto stage3list\n",stage3));
	  stage3list = List_push(stage3list,(void *) stage3);
	}
	Stage3middle_free(&stage3middle);
	i++;
      }

      while (i < n) {
	stage3middle = stage3middle_array[i];
	debug(printf("Ignoring stage3middle %d with goodness %d\n",i,Stage3middle_goodness(stage3middle)));
	Stage3middle_free(&stage3middle);
	i++;
      }

      FREE(stage3middle_array);
    }
  }

#ifdef PMAP_OLD
  Sequence_free(&genomicuc);
#elif defined(EXTRACT_GENOMICSEG)
  Sequence_free(&genomicuc);
#endif

  return stage3list;
}


static bool
middle_piece_local_p (int *querystart, int *queryend,
		      Chrpos_T *chrstart, Chrpos_T *chrend,
		      Chrnum_T *chrnum, Univcoord_T *chroffset, Univcoord_T *chrhigh,
		      Chrpos_T *chrlength, bool *plusp, int *genestrand,
		      Stage3_T from, Stage3_T to) {

  debug2(printf("? middle_piece_local_p from [%p] %d..%d (%u..%u) -> to [%p] %d..%d (%u..%u) => ",
		from,Stage3_querystart(from),Stage3_queryend(from),
		Stage3_chrstart(from),Stage3_chrend(from),
		to,Stage3_querystart(to),Stage3_queryend(to),
		Stage3_chrstart(to),Stage3_chrend(to)));

  if (Stage3_chimera_right_p(from) == true) {
    debug2(printf("false, because from is already part of a chimera on its right\n"));
    return false;
    
  } else if (Stage3_chimera_left_p(to) == true) {
    debug2(printf("false, because to is already part of a chimera on its left\n"));
    return false;

  } else if ((*chrnum = Stage3_chrnum(from)) != Stage3_chrnum(to)) {
    /* Different chromosomes */
    debug2(printf("different chromosomes\n"));
    return false;

  } else if (Stage3_watsonp(from) != Stage3_watsonp(to)) {
    /* Different strands */
    debug2(printf("different strands\n"));
    return false;

  } else if (Stage3_genestrand(from) != Stage3_genestrand(to)) {
    /* Different genestrands */
    debug2(printf("different genestrands\n"));
    return false;

  } else if (Stage3_querystart(to) <= Stage3_queryend(from) + CHIMERA_SLOP) {
    /* Already joinable */
    debug2(printf("wrong query order or already joinable\n"));
    return false;

  } else if ((*plusp = Stage3_watsonp(from)) == true) {
    if (Stage3_chrend(from) < Stage3_chrstart(to) &&
	Stage3_chrend(from) + 1000000 > Stage3_chrstart(to)) {
      debug2(printf("true, because %u < %u and %u + %u > %u\n",
		    Stage3_chrend(from),Stage3_chrstart(to),
		    Stage3_chrend(from),1000000,Stage3_chrstart(to)));
      Univ_IIT_interval_bounds(&(*chroffset),&(*chrhigh),&(*chrlength),chromosome_iit,
			       *chrnum,circular_typeint);
      *querystart = Stage3_queryend(from);
      *queryend = Stage3_querystart(to);
      *chrstart = Stage3_chrend(from);
      *chrend = Stage3_chrstart(to);
      *genestrand = Stage3_genestrand(from);
      return true;
    } else {
      debug2(printf("false, watsonp true, from_end %u, to start %u\n",
		    Stage3_chrend(from),Stage3_chrstart(to)));
      return false;
    }

  } else {
    if (Stage3_chrstart(to) < Stage3_chrend(from) &&
	Stage3_chrstart(to) + 1000000 > Stage3_chrend(from)) {
      debug2(printf("true, because %u < %u and %u + %u > %u\n",
		    Stage3_chrstart(to),Stage3_chrend(from),
		    Stage3_chrstart(to),1000000,Stage3_chrend(from)));
      Univ_IIT_interval_bounds(&(*chroffset),&(*chrhigh),&(*chrlength),chromosome_iit,
			       *chrnum,circular_typeint);
      *querystart = Stage3_queryend(from);
      *queryend = Stage3_querystart(to);
      *chrstart = Stage3_chrstart(to);
      *chrend = Stage3_chrend(from);
      *genestrand = Stage3_genestrand(from);
      return true;
    } else {
      debug2(printf("false, watsonp false, from_end %u, to start %u\n",
		    Stage3_chrend(from),Stage3_chrstart(to)));
      return false;
    }
  }
}


static bool
middle_piece_chimera_p (int *querystart, int *queryend, Stage3_T from, Stage3_T to) {

  debug2(printf("? middle_piece_chimera_p from [%p] %d..%d (%u..%u) -> to [%p] %d..%d (%u..%u) => ",
		from,Stage3_querystart(from),Stage3_queryend(from),
		Stage3_chrstart(from),Stage3_chrend(from),
		to,Stage3_querystart(to),Stage3_queryend(to),
		Stage3_chrstart(to),Stage3_chrend(to)));

  if (Stage3_chimera_right_p(from) == true) {
    debug2(printf("false, because from is already part of a chimera on its right\n"));
    return false;
    
  } else if (Stage3_chimera_left_p(to) == true) {
    debug2(printf("false, because to is already part of a chimera on its left\n"));
    return false;

  } else if (Stage3_querystart(to) <= Stage3_queryend(from) + CHIMERA_SLOP) {
    /* Already joinable */
    debug2(printf("wrong query order or already joinable\n"));
    return false;

  } else {
    *querystart = Stage3_queryend(from);
    *queryend = Stage3_querystart(to);
    return true;
  }
}


/* Does not alter stage3list.  Puts Stage3_T objects into stage3array_sub1, stage3array_sub2, or both */
static void
local_separate_paths (Stage3_T **stage3array_sub1, int *npaths_sub1, 
		      Stage3_T **stage3array_sub2, int *npaths_sub2,
		      List_T stage3list) {
  List_T p;
  Stage3_T from, to, stage3;
  Stage3_T *by_queryend, *by_querystart;
  Chrnum_T chrnum;
  int npaths, i, j, k, kstart, kend;
  int queryend;

  debug2(printf("local_separate_paths called with list length %d\n",List_length(stage3list)));
#ifdef DEBUG2
  for (p = stage3list; p != NULL; p = List_next(p)) {
    printf("%p\n",List_head(p));
  }
#endif

  if (stage3list == NULL) {
    *stage3array_sub1 = (Stage3_T *) NULL;
    *npaths_sub1 = 0;
    *stage3array_sub2 = (Stage3_T *) NULL;
    *npaths_sub2 = 0;
    return;

  } else {
    for (p = stage3list; p != NULL; p = List_next(p)) {
      stage3 = (Stage3_T) List_head(p);
      Stage3_clear_joinable(stage3);
    }
  }

  by_queryend = (Stage3_T *) List_to_array_out_n(&npaths,stage3list);
  qsort(by_queryend,npaths,sizeof(Stage3_T),Stage3_chrnum_queryend_cmp);

  by_querystart = (Stage3_T *) List_to_array_out_n(&npaths,stage3list);
  qsort(by_querystart,npaths,sizeof(Stage3_T),Stage3_chrnum_querystart_cmp);

#ifdef DEBUG2
  for (i = 0; i < npaths; i++) {
    stage3 = (Stage3_T) by_queryend[i];
    printf("from: %p query %d..%d, chrnum %d, genomic %u..%u\t",
	   stage3,Stage3_querystart(stage3),Stage3_queryend(stage3),
	   Stage3_chrnum(stage3),Stage3_genomicstart(stage3),Stage3_genomicend(stage3));

    stage3 = (Stage3_T) by_querystart[i];
    printf("to: %p query %d..%d, chrnum %d, genomic %u..%u\n",
	   stage3,Stage3_querystart(stage3),Stage3_queryend(stage3),
	   Stage3_chrnum(stage3),Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
  }
#endif

  kend = 0;
  for (i = 0; i < npaths; i++) {
    debug2(printf("queryend %d:",i));
    from = by_queryend[i];

    /* Find matching chromosomal bounds for querystart */
    chrnum = Stage3_chrnum(from);
    while (kend < npaths && Stage3_chrnum(by_querystart[kend]) == chrnum) {
      kend++;
    }
    kstart = kend - 1;
    while (kstart >= 0 && Stage3_chrnum(by_querystart[kstart]) == chrnum) {
      kstart--;
    }
    kstart++;
    debug2(printf(" querystart bounded by %d..%d:",kstart,kend));


    /* Find matching querystart */
    queryend = Stage3_queryend(from);
    j = kstart;
    while (j < kend && Stage3_querystart(by_querystart[j]) < queryend + CHIMERA_SLOP) {
      j++;
    }
    j--;

    while (j >= kstart && Stage3_querystart(by_querystart[j]) > queryend - CHIMERA_SLOP) {
      j--;
    }
    j++;

    while (j < kend && Stage3_querystart(by_querystart[j]) < queryend + CHIMERA_SLOP) {
      to = by_querystart[j];

      debug2(printf(" %d",j));
      if (Chimera_local_join_p(from,to,CHIMERA_SLOP) == true) {
	debug2(printf("(to %d)",i));
	Stage3_set_joinable_left(from);
	Stage3_set_joinable_right(to);
      }

      j++;
    }
    debug2(printf("\n"));
  }

  FREE(by_querystart);
  FREE(by_queryend);


  *npaths_sub1 = *npaths_sub2 = 0;
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    debug2(printf("Stage3 %p.  joinable_left_p %d, joinable_right_p %d\n",
		  stage3,Stage3_joinable_left_p(stage3),Stage3_joinable_right_p(stage3)));
    if (Stage3_joinable_left_p(stage3) == true) {
      debug2(printf("Putting stage3 %p into local sub1\n",stage3));
      (*npaths_sub1)++;
    }
    if (Stage3_joinable_right_p(stage3) == true) {
      debug2(printf("Putting stage3 %p into local sub2\n",stage3));
      (*npaths_sub2)++;
    }
  }

  if (*npaths_sub1 == 0 || *npaths_sub2 == 0) {
    *stage3array_sub1 = (Stage3_T *) NULL;
    *npaths_sub1 = 0;
    *stage3array_sub2 = (Stage3_T *) NULL;
    *npaths_sub2 = 0;
    
  } else {
    *stage3array_sub1 = (Stage3_T *) MALLOC((*npaths_sub1) * sizeof(Stage3_T)); /* Return value */
    *stage3array_sub2 = (Stage3_T *) MALLOC((*npaths_sub2) * sizeof(Stage3_T)); /* Return value */
    j = k = 0;
    for (p = stage3list; p != NULL; p = List_next(p)) {
      stage3 = (Stage3_T) List_head(p);
      /* Note: it is possible that the same stage3 object gets put into both lists */
      if (Stage3_joinable_left_p(stage3) == true) {
	debug2(printf("Putting %p into sub1\n",stage3));
	(*stage3array_sub1)[j++] = stage3;
      }
      if (Stage3_joinable_right_p(stage3) == true) {
	debug2(printf("Putting %p into sub2\n",stage3));
	(*stage3array_sub2)[k++] = stage3;
      }
    }
  }

  debug2(printf("local_separate_paths returning %d paths\n",List_length(stage3list)));
#ifdef DEBUG2
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    printf("%p %p\n",stage3,Stage3_pairs(stage3));
  }
#endif

  return;
}


/* Does not alter stage3list.  Puts Stage3_T objects into stage3array_sub1, stage3array_sub2, or both */
static void
distant_separate_paths (Stage3_T **stage3array_sub1, int *npaths_sub1, 
			Stage3_T **stage3array_sub2, int *npaths_sub2,
			List_T stage3list) {
  List_T p;
  Stage3_T from, to, stage3;
  Stage3_T *by_queryend, *by_querystart;
  int npaths, i, j, k;
  int queryend;

  debug2(printf("distant_separate_paths called with list length %d\n",List_length(stage3list)));
#ifdef DEBUG2
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    printf("%p %p\n",stage3,Stage3_pairs(stage3));
  }
#endif


  if (stage3list == NULL) {
    *stage3array_sub1 = (Stage3_T *) NULL;
    *npaths_sub1 = 0;
    *stage3array_sub2 = (Stage3_T *) NULL;
    *npaths_sub2 = 0;
    return;

  } else {
    for (p = stage3list; p != NULL; p = List_next(p)) {
      stage3 = (Stage3_T) List_head(p);
      Stage3_clear_joinable(stage3);
    }
  }

  by_queryend = (Stage3_T *) List_to_array_out_n(&npaths,stage3list);
  qsort(by_queryend,npaths,sizeof(Stage3_T),Stage3_queryend_cmp);

  by_querystart = (Stage3_T *) List_to_array_out_n(&npaths,stage3list);
  qsort(by_querystart,npaths,sizeof(Stage3_T),Stage3_querystart_cmp);

  j = 0;
  for (i = 0; i < npaths; i++) {
    from = by_queryend[i];
    queryend = Stage3_queryend(from);

    while (j < npaths && Stage3_querystart(by_querystart[j]) < queryend + CHIMERA_SLOP) {
      j++;
    }
    j--;

    while (j >= 0 && Stage3_querystart(by_querystart[j]) > queryend - CHIMERA_SLOP) {
      j--;
    }
    j++;

    while (j < npaths && Stage3_querystart(by_querystart[j]) < queryend + CHIMERA_SLOP) {
      to = by_querystart[j];

      if (Chimera_distant_join_p(from,to,CHIMERA_SLOP) == true) {
	debug2(printf("Found distant join from %d to %d\n",i,j));
	Stage3_set_joinable_left(from);
	Stage3_set_joinable_right(to);
      }

      j++;
    }
  }

  FREE(by_querystart);
  FREE(by_queryend);


  *npaths_sub1 = *npaths_sub2 = 0;
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    if (Stage3_joinable_left_p(stage3) == true) {
      (*npaths_sub1)++;
    }
    if (Stage3_joinable_right_p(stage3) == true) {
      (*npaths_sub2)++;
    }
  }

  if (*npaths_sub1 == 0 || *npaths_sub2 == 0) {
    *stage3array_sub1 = (Stage3_T *) NULL;
    *npaths_sub1 = 0;
    *stage3array_sub2 = (Stage3_T *) NULL;
    *npaths_sub2 = 0;
  } else {
    *stage3array_sub1 = (Stage3_T *) MALLOC((*npaths_sub1) * sizeof(Stage3_T)); /* Return value */
    *stage3array_sub2 = (Stage3_T *) MALLOC((*npaths_sub2) * sizeof(Stage3_T)); /* Return value */
    j = k = 0;
    for (p = stage3list; p != NULL; p = List_next(p)) {
      stage3 = (Stage3_T) List_head(p);
      /* Note: it is possible that the same stage3 object gets put into both lists */
      if (Stage3_joinable_left_p(stage3) == true) {
	debug2(printf("Putting stage3 %p into distant sub1\n",stage3));
	(*stage3array_sub1)[j++] = stage3;
      }
      if (Stage3_joinable_right_p(stage3) == true) {
	debug2(printf("Putting stage3 %p into distant sub2\n",stage3));
	(*stage3array_sub2)[k++] = stage3;
      }
    }
  }

  debug2(printf("distant_separate_paths returning %d paths\n",List_length(stage3list)));
#ifdef DEBUG2
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    printf("%p %p\n",stage3,Stage3_pairs(stage3));
  }
#endif

  return;
}


static Stage3_T
merge_left_and_right_readthrough (Stage3_T *stage3array_sub1, int bestfrom,
				  Stage3_T *stage3array_sub2, int bestto,
				  int breakpoint, int queryntlength,
#ifdef PMAP
				  char *queryaaseq_ptr,
#endif
				  Sequence_T queryseq, char *queryseq_ptr, char *queryuc_ptr,
				  Genome_T genome, Genome_T genomealt, Pairpool_T pairpool,
				  Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
				  Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool) {
  Stage3_T stage3, best0, best1;

  best0 = stage3array_sub1[bestfrom];
  best1 = stage3array_sub2[bestto];

  debug2(printf("\nEntering merge_left_and_right_readthrough with bestfrom %d: %p, bestto %d: %p\n",
		bestfrom,best0,bestto,best1));
  debug2(printf("Running Stage3_merge_local\n"));

  if ((stage3 = Stage3_merge_local(best0,best1,/*minpos1*/0,/*maxpos1*/breakpoint,
				   /*minpos2*/breakpoint+1,/*maxpos2*/queryntlength,queryseq,
#ifdef PMAP
				   queryaaseq_ptr,
#endif
				   queryseq_ptr,queryuc_ptr,oligoindices_minor,diagpool,cellpool,
				   genome,genomealt,pairpool,dynprogL,dynprogM,dynprogR,maxpeelback)) == NULL) {
    return (Stage3_T) NULL;

  } else {
    debug2(printf("Changing genomicend of merged stage3 from %u to %u\n",Stage3_genomicend(stage3),Stage3_genomicend(best1)));
    Stage3_set_genomicend(stage3,Stage3_genomicend(best1));
#ifndef PMAP
    Stage3_guess_cdna_direction(best0,genome,genomealt);
#endif
    return stage3;
  }
}


#if 0
/* Returns a list with only two Stage3_T objects */
static List_T
merge_left_and_right_transloc (Stage3_T *stage3array_sub1, int npaths_sub1, int bestfrom,
			       Stage3_T *stage3array_sub2, int npaths_sub2, int bestto,
			       List_T stage3list) {
  List_T newstage3list, p;
  Stage3_T best0, best1, stage3, *array;
  int i, k;

  best0 = stage3array_sub1[bestfrom];
  best1 = stage3array_sub2[bestto];

  debug2(printf("\nEntering merge_left_and_right_transloc with bestfrom %d: %p, bestto %d: %p, and stage3list %d\n",
		bestfrom,best0,bestto,best1,List_length(stage3list)));

  debug2(printf("Before Stage3_merge_chimera, best0 is %p, query %d..%d\n",
		best0,Stage3_querystart(best0),Stage3_queryend(best0)));
  debug2(Stage3_print_ends(best0));
  debug2(printf("Before Stage3_merge_chimera, best1 is %p, query %d..%d\n",
		best1,Stage3_querystart(best1),Stage3_queryend(best1)));
  debug2(Stage3_print_ends(best1));

  debug2(printf("Rearranging paths\n"));
  newstage3list = (List_T) NULL;

  debug(printf("Pushing %p onto newstage3list\n",best0));
  debug(printf("Pushing %p onto newstage3list\n",best1));
  if (Stage3_npairs(best0) == 0) {
    Stage3_free(&best0);
    best0 = (Stage3_T) NULL;
  } else {
    newstage3list = List_push(newstage3list,(void *) best0);
    debug2(Stage3_print_ends(best0));
  }
  if (Stage3_npairs(best1) == 0) {
    Stage3_free(&best1);
    best1 = (Stage3_T) NULL;
  } else {
    newstage3list = List_push(newstage3list,(void *) best1);
    debug2(Stage3_print_ends(best1));
  }
  
  if (List_length(stage3list) > 2) {
    /* Push rest of results, taking care not to have duplicates */
    array = (Stage3_T *) MALLOCA((List_length(stage3list) - 2) * sizeof(Stage3_T));
    k = 0;
    for (p = stage3list; p != NULL; p = List_next(p)) {
      stage3 = (Stage3_T) List_head(p);
      if (Stage3_npairs(stage3) == 0) {
	Stage3_free(&stage3);
      } else if (stage3 == best0 || stage3 == best1) {
	/* Skip */
      } else {
	array[k++] = stage3;
	debug(printf("Pushing %p onto newstage3list\n",stage3));
	newstage3list = List_push(newstage3list,(void *) stage3);
      }
    }
    qsort(array,k,sizeof(Stage3_T),Stage3_identity_cmp);
    FREEA(array);
  }
    
  List_free(&stage3list);
  return List_reverse(newstage3list);
}
#endif


static int
find_breakpoint (int *cdna_direction, int *chimerapos, int *chimeraequivpos, int *exonexonpos,
		 char *donor1, char *donor2, char *acceptor2, char *acceptor1,
		 bool *donor_watsonp, bool *acceptor_watsonp, double *donor_prob, double *acceptor_prob,
		 Stage3_T from, Stage3_T to,
#ifdef PMAP
		 Sequence_T queryntseq,
#endif
		 Sequence_T queryseq, Sequence_T queryuc, int queryntlength,
		 Genome_T genome, Genome_T genomealt, Pairpool_T pairpool,
		 Univ_IIT_T chromosome_iit) {
  int breakpoint, rangelow, rangehigh, leftpos, rightpos, midpos;
  int maxpeelback_from, maxpeelback_to;
  int found_cdna_direction, try_cdna_direction;
  char comp;			/* Not really used anywhere */

  int queryjump;
  int genomejump;
  bool max_extend_p;
  Chrpos_T left_chrlength, right_chrlength;
  Univcoord_T chroffset, chrhigh;

  if (Stage3_queryend(from) < Stage3_querystart(to)) {
    /* Gap exists between the two parts */
    if ((leftpos = Stage3_queryend(from) - CHIMERA_EXTEND) < 0) {
      leftpos = 0;
    }
    if ((rightpos = Stage3_querystart(to) + CHIMERA_EXTEND) >= queryntlength) {
      rightpos = queryntlength - 1;
    }
    maxpeelback_from = CHIMERA_EXTEND;
    maxpeelback_to = CHIMERA_EXTEND;
    debug2(printf("overlap: leftpos %d, rightpos %d, queryntlength %d, maxpeelback_from %d, maxpeelback_to %d\n",
		  leftpos,rightpos,queryntlength,maxpeelback_from,maxpeelback_to));

    if (Stage3_watsonp(from) == true && Stage3_watsonp(to) == true) {
      queryjump = Stage3_querystart(to) - Stage3_queryend(from) - 1;
      genomejump = Stage3_genomicstart(to) - Stage3_genomicend(from) - 1U;
      max_extend_p = ((int) genomejump == queryjump) ? false : true;
      debug2(printf("gap exists: genomejump = %u, queryjump = %d, max_extend_p = %d\n",genomejump,queryjump,max_extend_p));
    } else if (Stage3_watsonp(from) == false && Stage3_watsonp(to) == false) {
      queryjump = Stage3_querystart(to) - Stage3_queryend(from) - 1;
      genomejump = Stage3_genomicend(from) - Stage3_genomicstart(to) - 1U;
      max_extend_p = ((int) genomejump == queryjump) ? false : true;
      debug2(printf("gap exists: genomejump = %u, queryjump = %d, max_extend_p = %d\n",genomejump,queryjump,max_extend_p));
    } else {
      max_extend_p = false;
    }
    
  } else {
    /* Two parts overlap */
    if ((leftpos = Stage3_querystart(to) - CHIMERA_EXTEND) < 0) {
      leftpos = 0;
    }
    if ((rightpos = Stage3_queryend(from) + CHIMERA_EXTEND) >= queryntlength) {
      rightpos = queryntlength - 1;
    }
    midpos = (leftpos+rightpos)/2;
    /* maxpeelback_from = rightpos - Stage3_querystart(to); */
    /* maxpeelback_to = Stage3_queryend(from) - leftpos; */
    maxpeelback_from = rightpos - midpos;
    maxpeelback_to = midpos - leftpos;
    debug2(printf("overlap: leftpos %d, rightpos %d, midpos %d, queryntlength %d, maxpeelback_from %d, maxpeelback_to %d\n",
		  leftpos,rightpos,midpos,queryntlength,maxpeelback_from,maxpeelback_to));
#if 0
    if (Stage3_watsonp(from) == true && Stage3_watsonp(to) == true) {
      queryjump = Stage3_queryend(from) - Stage3_querystart(to) - 1;
      genomejump = Stage3_genomicend(from) - Stage3_genomicstart(to) - 1U;
      max_extend_p = (genomejump == queryjump) ? false : true;
    } else if (Stage3_watsonp(from) == false && Stage3_watsonp(to) == false) {
      queryjump = Stage3_queryend(from) - Stage3_querystart(to) - 1;
      genomejump = Stage3_genomicstart(to) - Stage3_genomicend(from) - 1U;
      max_extend_p = (genomejump == queryjump) ? false : true;
    } else {
      max_extend_p = false;
    }
#else
    debug2(printf("parts overlap: max_extend_p is false\n"));
    max_extend_p = false;
#endif
  }

  debug2(printf("Before Stage3_extend_right, bestfrom is %p, query %d..%d, rightpos %d, pairs %p\n",
		from,Stage3_querystart(from),Stage3_queryend(from),rightpos,Stage3_pairs(from)));
  debug2(Stage3_print_ends(from));
  debug2(printf("Before Stage3_extend_left, bestto is %p, query %d..%d, leftpos %d, pairs %p\n",
		to,Stage3_querystart(to),Stage3_queryend(to),leftpos,Stage3_pairs(to)));
  debug2(Stage3_print_ends(to));
  
  Stage3_extend_right(from,/*goal*/rightpos,
#ifdef PMAP
		      /*querylength*/Sequence_fulllength(queryntseq),
		      /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
		      /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
#else
		      /*querylength*/Sequence_fulllength(queryseq),
		      /*queryseq_ptr*/Sequence_fullpointer(queryseq),
		      /*queryuc_ptr*/Sequence_fullpointer(queryuc),
#endif
		      max_extend_p,genome,genomealt,pairpool,Stage3_genestrand(from),maxpeelback_from);

  Stage3_extend_left(to,/*goal*/leftpos,
#ifdef PMAP
		     /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
		     /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
#else
		     /*queryseq_ptr*/Sequence_fullpointer(queryseq),
		     /*queryuc_ptr*/Sequence_fullpointer(queryuc),
#endif
		     max_extend_p,genome,genomealt,pairpool,Stage3_genestrand(to),maxpeelback_to);

  debug2(printf("Before Chimera_find_breakpoint, bestfrom is %p, query %d..%d, pairs %p\n",
                 from,Stage3_querystart(from),Stage3_queryend(from),Stage3_pairs(from)));
  debug2(Stage3_print_ends(from));
  debug2(printf("Before Chimera_find_breakpoint, bestto is %p, query %d..%d, pairs %p\n",
                 to,Stage3_querystart(to),Stage3_queryend(to),Stage3_pairs(to)));
  debug2(Stage3_print_ends(to));

  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&left_chrlength,chromosome_iit,Stage3_chrnum(from),circular_typeint);
  Univ_IIT_interval_bounds(&chroffset,&chrhigh,&right_chrlength,chromosome_iit,Stage3_chrnum(to),circular_typeint);

  if ((*chimerapos = Chimera_find_breakpoint(&(*chimeraequivpos),&rangelow,&rangehigh,
					     &(*donor1),&(*donor2),&(*acceptor2),&(*acceptor1),
					     from,to,queryntlength,genome,genomealt,
					     left_chrlength,right_chrlength)) < 0) {
    /* TODO: Allow finding a breakpoint for DNA-Seq, which needs no donor or acceptor nucleotides */
    debug2(printf("Chimera_find_breakpoint returns no value\n"));
    *donor_prob = *acceptor_prob = 0.0;
    *donor_watsonp = *acceptor_watsonp = true;
    *cdna_direction = 0;
    return -1;

  } else {
    debug2(printf("Chimera_find_breakpoint has chimerapos %d..%d\n",*chimerapos,*chimeraequivpos));

    Stage3_trim_right(from,/*goal*/rangehigh,
		      /*queryseq_ptr*/Sequence_fullpointer(queryseq),
		      genome,genomealt,pairpool);

    Stage3_trim_left(to,/*goal*/rangelow,
		     /*queryseq_ptr*/Sequence_fullpointer(queryseq),
		     genome,genomealt,pairpool);
    
    debug2(printf("Before Chimera_find_exonexon, bestfrom is %p, query %d..%d, pairs %p\n",
                  from,Stage3_querystart(from),Stage3_queryend(from),Stage3_pairs(from)));
    debug2(printf("Before Chimera_find_exonexon, bestto is %p, query %d..%d, pairs %p\n",
                  to,Stage3_querystart(to),Stage3_queryend(to),Stage3_pairs(to)));

    if ((*exonexonpos = Chimera_find_exonexon(&found_cdna_direction,&try_cdna_direction,
					      &(*donor1),&(*donor2),&(*acceptor2),&(*acceptor1),
					      &comp,&(*donor_watsonp),&(*acceptor_watsonp),&(*donor_prob),&(*acceptor_prob),
					      /*left_part*/from,/*right_part*/to,genome,genomealt ? genomealt : genome,
					      chromosome_iit,/*breakpoint_start*/Stage3_querystart(to),
					      /*breakpoint_end*/Stage3_queryend(from))) <= 0) {
      /* Couldn't find a good exon-exon junction, so rely on sequence */
      *donor_prob = *acceptor_prob = 0.0;
      *donor_watsonp = *acceptor_watsonp = true;
      
      debug2(printf("Chimera_find_breakpoint returns boundary at %d..%d (switch can occur at %d..%d)\n",
		    *chimerapos,*chimeraequivpos,(*chimerapos)-1,*chimeraequivpos));
      
      breakpoint = ((*chimerapos) + (*chimeraequivpos))/2;
      *cdna_direction = try_cdna_direction;
      debug2(printf("Exon-exon boundary not found, but setting breakpoint to be %d\n",breakpoint));
      return breakpoint;
      
    } else {
      /* Use the exon-exon solution */
      breakpoint = *chimerapos = *chimeraequivpos = *exonexonpos;
      *cdna_direction = found_cdna_direction;
      debug2(printf("Exon-exon boundary found at %d, which is breakpoint.  Comp = %c\n",
		    *exonexonpos,comp));
      return breakpoint;
    }
  }
}


/* Can potentially include a larger stage3list */
static List_T
check_for_local (bool *mergedp, List_T stage3list, int effective_start, int effective_end,
		 Sequence_T queryseq, Sequence_T queryuc,
#ifdef PMAP
		 Sequence_T queryntseq,
#endif
		 int queryntlength, Stage2_alloc_T stage2_alloc,
		 Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
		 Matchpool_T matchpool, Genome_T genome, Genome_T genomealt, Pairpool_T pairpool,
		 Diagpool_T diagpool, Cellpool_T cellpool,
		 Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, int min_matches) {
  List_T gregions = NULL, p;
  Stage3_T *stage3array_sub1 = NULL, *stage3array_sub2 = NULL, from, to, stage3, stage3_merged;
  Sequence_T querysubseq = NULL, querysubuc = NULL;
  Diagnostic_T diagnostic;
  int bestfrom, bestto;
  int five_margin, three_margin, five_score = 0, three_score = 0;
  int extension;
  int npaths_sub1 = 0, npaths_sub2 = 0;
  bool lowidentityp, poorp, repetitivep;

  int max_single_goodness;
  int breakpoint, chimerapos, chimeraequivpos, exonexonpos;
  int chimera_cdna_direction;
  char donor1, donor2, acceptor2, acceptor1;
  bool donor_watsonp, acceptor_watsonp;
  double donor_prob, acceptor_prob;
  
  int kstart1, kstart2, kend1, kend2;
  Chrnum_T chrnum;
#ifdef DEBUG2
  int k;
#endif


#ifdef PMAP
  five_margin = effective_start - 3*Sequence_trim_start(queryseq);
  three_margin = 3*Sequence_trim_end(queryseq) - effective_end;
  debug2(printf("Margins are %d = %d - %d on the 5' end and %d = %d - %d on the 3' end\n",
		five_margin,effective_start,3*Sequence_trim_start(queryseq),
		three_margin,3*Sequence_trim_end(queryseq),effective_end));
#else
  five_margin = effective_start - Sequence_trim_start(queryseq);
  three_margin = Sequence_trim_end(queryseq) - effective_end;
  debug2(printf("Margins are %d = %d - %d on the 5' end and %d = %d - %d on the 3' end\n",
		five_margin,effective_start,Sequence_trim_start(queryseq),
		three_margin,Sequence_trim_end(queryseq),effective_end));
#endif

#ifdef DEBUG2A
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    Pair_dump_array(Stage3_pairarray(stage3),Stage3_npairs(stage3),/*zerobasedp*/true);
    printf("\n");
  }
#endif

  /* Stage3_recompute_goodness(stage3list); */
  max_single_goodness = 0;
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    if (Stage3_goodness(stage3) > max_single_goodness) {
      max_single_goodness = Stage3_goodness(stage3);
    }
  }
  debug2(printf("max single goodness = %d\n",max_single_goodness));


  debug2(printf("Running local_separate_paths\n"));
  local_separate_paths(&stage3array_sub1,&npaths_sub1,&stage3array_sub2,&npaths_sub2,
		       stage3list);
  debug2(printf("local: npaths_sub1 %d, npaths_sub2 %d, stage3list %d\n",
		npaths_sub1,npaths_sub2,List_length(stage3list)));

  if (npaths_sub1 == 0 && npaths_sub2 == 0) {
    /* Need to compute on margin explicitly */
    if (five_margin < chimera_margin && three_margin < chimera_margin) {
      debug2(printf("Insufficient margins\n"));
    } else if (five_margin > three_margin) {
#if 0
      /* extension makes it harder to find the other alignment.  The merging process will help fill in any gap. */
      extension = CHIMERA_SLOP;
      debug2(printf("Comparing extension %d with %d = (effective_start %d)/2\n",
		    extension,effective_start/2,effective_start));
      if (extension > effective_start/2) {
	/* Extension occupies more than 1/3 of sequence */
	debug2(printf("Proposed extension of %d is too long relative to effective_start %d\n",extension,effective_start));
	extension = effective_start/3;
      }
#else
      extension = 0;
#endif
      if ((querysubseq = Sequence_subsequence(queryseq,0,effective_start+extension)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,0,effective_start+extension)) != NULL) {
	  debug2(printf("5 margin > 3 margin.  "));
	  debug2(printf("Beginning Stage1_compute on 5' margin from effective_start %d (%d..%d)\n",
			effective_start,0,effective_start+extension));
	  debug2a(Sequence_stdout(querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(querysubuc),Sequence_fulllength(querysubuc),
				      Oligoindex_array_elt(oligoindices_major,0));
	  if (poorp == true || repetitivep == true) {
	    debug2(printf("Subsequence is poor or repetitive\n"));
	  } else {
	    if (mode == CMET_NONSTRANDED || mode == ATOI_NONSTRANDED || mode == TTOC_NONSTRANDED) {
	      gregions = Stage1_compute_nonstranded(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,
						    chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
						    stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    } else {
	      gregions = Stage1_compute(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,/*genestrand*/0,
					chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
					stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    }
	    debug2(printf("A.  Performing Stage 3 starting with list length %d\n",List_length(stage3list)));
	    stage3list = stage3_from_gregions(stage3list,gregions,querysubseq,querysubuc,
#ifdef PMAP
					      queryntseq,
#endif
					      stage2_alloc,oligoindices_major,oligoindices_minor,
					      genome,genomealt,pairpool,diagpool,cellpool,
					      dynprogL,dynprogM,dynprogR,min_matches,/*worker_stopwatch*/NULL);
#ifdef DEBUG2
	    for (p = stage3list; p != NULL; p = List_next(p)) {
	      stage3 = (Stage3_T) List_head(p);
	      printf("%d..%d, %u..%u\n",
		     Stage3_querystart(stage3),Stage3_queryend(stage3),
		     Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
	    }
#endif
	  }
	  Diagnostic_free(&diagnostic);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);
	}
	Sequence_free(&querysubseq);
      }

      /* And recompute on original part, just in case stage 1 was led astray by the ends */
      if ((querysubseq = Sequence_subsequence(queryseq,effective_start,queryntlength)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,effective_start,queryntlength)) != NULL) {
	  debug2(printf("Recomputing on original part.  "));
	  debug2(printf("Beginning Stage1_compute on 5' margin from effective_start %d (%d..%d)\n",
			effective_start,effective_start,queryntlength));
	  debug2a(Sequence_stdout(querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(querysubuc),Sequence_fulllength(querysubuc),
				      Oligoindex_array_elt(oligoindices_major,0));
	  if (poorp == true || repetitivep == true) {
	    debug2(printf("Subsequence is poor or repetitive\n"));
	  } else {
	    if (mode == CMET_NONSTRANDED || mode == ATOI_NONSTRANDED || mode == TTOC_NONSTRANDED) {
	      gregions = Stage1_compute_nonstranded(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,
						    chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
						    stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    } else {
	      gregions = Stage1_compute(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,/*genestrand*/0,
					chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
					stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    }
	    debug2(printf("B.  Performing Stage 3 starting with list length %d\n",List_length(stage3list)));
	    stage3list = stage3_from_gregions(stage3list,gregions,querysubseq,querysubuc,
#ifdef PMAP
					      queryntseq,
#endif
					      stage2_alloc,oligoindices_major,oligoindices_minor,
					      genome,genomealt,pairpool,diagpool,cellpool,
					      dynprogL,dynprogM,dynprogR,min_matches,/*worker_stopwatch*/NULL);
#ifdef DEBUG2
	    for (p = stage3list; p != NULL; p = List_next(p)) {
	      stage3 = (Stage3_T) List_head(p);
	      printf("%d..%d, %u..%u\n",
		     Stage3_querystart(stage3),Stage3_queryend(stage3),
		     Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
	    }
#endif
	  }
	  Diagnostic_free(&diagnostic);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);
	}
	Sequence_free(&querysubseq);
      }

      debug2(printf("Running local_separate_paths\n"));
      local_separate_paths(&stage3array_sub1,&npaths_sub1,&stage3array_sub2,&npaths_sub2,
			   stage3list);
      debug2(printf("local: npaths_sub1 %d, npaths_sub2 %d, stage3list %d\n",
		    npaths_sub1,npaths_sub2,List_length(stage3list)));
      
    } else {
#if 0
      /* extension makes it harder to find the other alignment.  The merging process will help fill in any gap. */
      extension = CHIMERA_SLOP;
      debug2(printf("Comparing extension %d with %d = (queryntlength %d - effective_end %d)/2\n",
		    extension,(queryntlength-effective_end)/2,queryntlength,effective_end));
      if (extension > (queryntlength - effective_end)/2) {
	/* Extension occupies more than 1/3 of sequence */
	debug2(printf("Proposed extension of %d is too long relative to queryntlength %d and effective_end %d\n",
		      extension,queryntlength,effective_end));
	extension = (queryntlength - effective_end)/3;
      }
#else
      extension = 0;
#endif
      if ((querysubseq = Sequence_subsequence(queryseq,effective_end-extension,queryntlength)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,effective_end-extension,queryntlength)) != NULL) {
	  debug2(printf("5 margin <= 3 margin.  "));
	  debug2(printf("Beginning Stage1_compute on 3' margin from effective_end %d (%d..%d) (extension %d)\n",
			effective_end,effective_end-extension,queryntlength,extension));
	  debug2(Sequence_stdout(querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(querysubuc),Sequence_fulllength(querysubuc),
				      Oligoindex_array_elt(oligoindices_major,0));
	  if (poorp == true || repetitivep == true) {
	    debug2(printf("Subsequence is poor or repetitive\n"));
	  } else {
	    if (mode == CMET_NONSTRANDED || mode == ATOI_NONSTRANDED || mode == TTOC_NONSTRANDED) {
	      gregions = Stage1_compute_nonstranded(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,
						    chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
						    stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    } else {
	      gregions = Stage1_compute(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,/*genestrand*/0,
					chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
					stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    }
	    debug2(printf("C.  Performing Stage 3 with list length %d\n",List_length(stage3list)));
	    stage3list = stage3_from_gregions(stage3list,gregions,querysubseq,querysubuc,
#ifdef PMAP
					      queryntseq,
#endif
					      stage2_alloc,oligoindices_major,oligoindices_minor,
					      genome,genomealt,pairpool,diagpool,cellpool,
					      dynprogL,dynprogM,dynprogR,min_matches,/*worker_stopwatch*/NULL);
#ifdef DEBUG2
	    for (p = stage3list; p != NULL; p = List_next(p)) {
	      stage3 = (Stage3_T) List_head(p);
	      printf("%d..%d, %u..%u\n",
		     Stage3_querystart(stage3),Stage3_queryend(stage3),
		     Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
	    }
#endif
	  }
	  Diagnostic_free(&diagnostic);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);
	}
	Sequence_free(&querysubseq);
      }

      /* And recompute on original part, just in case stage 1 was led astray by the ends */
      if ((querysubseq = Sequence_subsequence(queryseq,0,effective_end)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,0,effective_end)) != NULL) {
	  debug2(printf("Recomputing on original part.  "));
	  debug2(printf("Beginning Stage1_compute on 3' margin from effective_end %d (%d..%d), extension %d\n",
			effective_end,0,effective_end,extension));
	  debug2(Sequence_stdout(querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(querysubuc),Sequence_fulllength(querysubuc),
				      Oligoindex_array_elt(oligoindices_major,0));
	  if (poorp == true || repetitivep == true) {
	    debug2(printf("Subsequence is poor or repetitive\n"));
	  } else {
	    if (mode == CMET_NONSTRANDED || mode == ATOI_NONSTRANDED || mode == TTOC_NONSTRANDED) {
	      gregions = Stage1_compute_nonstranded(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,
						    chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
						    stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    } else {
	      gregions = Stage1_compute(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,/*genestrand*/0,
					chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
					stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    }
	    debug2(printf("D.  Performing Stage 3 with list length %d\n",List_length(stage3list)));
	    stage3list = stage3_from_gregions(stage3list,gregions,querysubseq,querysubuc,
#ifdef PMAP
					      queryntseq,
#endif
					      stage2_alloc,oligoindices_major,oligoindices_minor,
					      genome,genomealt,pairpool,diagpool,cellpool,
					      dynprogL,dynprogM,dynprogR,min_matches,/*worker_stopwatch*/NULL);
#ifdef DEBUG2
	    for (p = stage3list; p != NULL; p = List_next(p)) {
	      stage3 = (Stage3_T) List_head(p);
	      printf("%d..%d, %u..%u\n",
		     Stage3_querystart(stage3),Stage3_queryend(stage3),
		     Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
	    }
#endif
	  }
	  Diagnostic_free(&diagnostic);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);

	}
	Sequence_free(&querysubseq);
      }

      debug2(printf("Running local_separate_paths\n"));
      local_separate_paths(&stage3array_sub1,&npaths_sub1,&stage3array_sub2,&npaths_sub2,
			   stage3list);
      debug2(printf("local: npaths_sub1 %d, npaths_sub2 %d, stage3list %d\n",
		    npaths_sub1,npaths_sub2,List_length(stage3list)));
    }
  }

  *mergedp = false;
  if (npaths_sub1 == 0 && npaths_sub2 == 0) {
    /* Skip */

  } else if (npaths_sub1 == 0) {
    /* Skip */
    FREE(stage3array_sub2);

  } else if (npaths_sub2 == 0) {
    /* Skip */
    FREE(stage3array_sub1);

  } else {
    /* Iterate for each chromosome */
    qsort(stage3array_sub1,npaths_sub1,sizeof(Stage3_T),Stage3_chrnum_cmp);
    qsort(stage3array_sub2,npaths_sub2,sizeof(Stage3_T),Stage3_chrnum_cmp);


    kend1 = kend2 = 0;
    *mergedp = false;
    /* List_free(&stage3list); */

    while (kend1 < npaths_sub1 && kend2 < npaths_sub2) {
      kstart1 = kend1;
      kstart2 = kend2;
      chrnum = Stage3_chrnum(stage3array_sub1[kstart1]);
      while (kend1 < npaths_sub1 && Stage3_chrnum(stage3array_sub1[kend1]) == chrnum) {
	kend1++;
      }
      while (kend2 < npaths_sub2 && Stage3_chrnum(stage3array_sub2[kend2]) == chrnum) {
	kend2++;
      }

#ifdef DEBUG2
      printf("Chimera_bestpath left\n");
      for (k = kstart1; k < kend1; k++) {
	stage3 = stage3array_sub1[k];
	printf("%d..%d, %d:%u..%u\n",
	       Stage3_querystart(stage3),Stage3_queryend(stage3),
	       Stage3_chrnum(stage3),Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
      }
      printf("Chimera_bestpath right\n");
      for (k = kstart2; k < kend2; k++) {
	stage3 = stage3array_sub2[k];
	printf("%d..%d, %d:%u..%u\n",
	       Stage3_querystart(stage3),Stage3_queryend(stage3),
	       Stage3_chrnum(stage3),Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
      }
#endif

      if (Chimera_bestpath(&five_score,&three_score,&chimerapos,&chimeraequivpos,&bestfrom,&bestto,
			   &(stage3array_sub1[kstart1]),/*npaths1*/kend1-kstart1,
			   &(stage3array_sub2[kstart2]),/*npaths2*/kend2-kstart2,
			   queryntlength,CHIMERA_SLOP,/*circularp*/NULL,/*localp*/true) == false) {
	/* Skip */
	debug2(printf("Chimera_bestpath returns false\n"));

      } else {
	from = stage3array_sub1[kstart1 + bestfrom];
	to = stage3array_sub2[kstart2 + bestto];
	debug2(printf("Chimera_bestpath returns bestfrom %d (%d..%d, %u..%u) to bestto %d (%d..%d, %u..%u)\n",
		      bestfrom,Stage3_querystart(from),Stage3_queryend(from),Stage3_genomicstart(from),Stage3_genomicend(from),
		      bestto,Stage3_querystart(to),Stage3_queryend(to),Stage3_genomicstart(to),Stage3_genomicend(to)));

	if ((breakpoint = find_breakpoint(&chimera_cdna_direction,&chimerapos,&chimeraequivpos,&exonexonpos,
					  &donor1,&donor2,&acceptor2,&acceptor1,
					  &donor_watsonp,&acceptor_watsonp,&donor_prob,&acceptor_prob,from,to,
#ifdef PMAP
					  queryntseq,
#endif
					  queryseq,queryuc,queryntlength,
					  genome,genomealt,pairpool,chromosome_iit)) <= 0) {
	  debug2(printf("Cannot find breakpoint\n"));

	} else {
	  debug2(printf("find_breakpoint returns %d\n",breakpoint));

	  /* Check to see if we can merge chimeric parts */
	  debug2(printf("Before Stage3_mergeable, bestfrom is %p, query %d..%d, pairs %p\n",
			from,Stage3_querystart(from),Stage3_queryend(from),Stage3_pairs(from)));
	  debug2(printf("Before Stage3_mergeable, bestto is %p, query %d..%d, pairs %p\n",
			to,Stage3_querystart(to),Stage3_queryend(to),Stage3_pairs(to)));
	
	  if (Stage3_mergeable(from,to,breakpoint,queryntlength) == true) {
	    debug2(printf("Mergeable! -- Merging left and right as a readthrough\n"));
	    if ((stage3_merged =
		 merge_left_and_right_readthrough(&(stage3array_sub1[kstart1]),/*npaths1:kend1-kstart1,*/bestfrom,
						  &(stage3array_sub2[kstart2]),/*npaths2:kend2-kstart2,*/bestto,
						  breakpoint,queryntlength,queryseq,
#ifdef PMAP
						  /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
						  /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
						  /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
#else
						  /*queryseq_ptr*/Sequence_fullpointer(queryseq),
						  /*queryuc_ptr*/Sequence_fullpointer(queryuc),
#endif
						  genome,genomealt,pairpool,dynprogL,dynprogM,dynprogR,
						  oligoindices_minor,diagpool,cellpool)) != NULL) {
	      stage3list = List_push(stage3list,(void *) stage3_merged);
	    }

	    debug2(printf("After merge_left_and_right_readthrough, bestfrom is %p, query %d..%d, pairs %p\n",
			  from,Stage3_querystart(from),Stage3_queryend(from),Stage3_pairs(from)));
	    debug2(printf("After merge_left_and_right_readthrough, bestto is %p, query %d..%d, pairs %p\n",
			  to,Stage3_querystart(to),Stage3_queryend(to),Stage3_pairs(to)));
	  }
	}
      }
    }

    FREE(stage3array_sub2);
    FREE(stage3array_sub1);

    /* stage3list = List_reverse(stage3list); */
  }

  debug2(printf("check_for_local returning list of length %d\n",List_length(stage3list)));
#ifdef DEBUG2
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    printf("%p %p\n",stage3,Stage3_pairs(stage3));
  }
#endif

  /* stage3list = stage3list_remove_empties(stage3list); */

#if 0
  /* Should be handled by apply_stage3 loop */
  /* Needed after calls to stage3_from_gregions */
  Stage3_recompute_goodness(stage3list);
  stage3list = stage3list_remove_duplicates(stage3list);
#endif

  return stage3list;
}


static List_T
check_for_chimera (bool *mergedp, Chimera_T *chimera, List_T stage3list, int effective_start, int effective_end,
		   Sequence_T queryseq, Sequence_T queryuc,
#ifdef PMAP
		   Sequence_T queryntseq,
#endif
		   int queryntlength, Stage2_alloc_T stage2_alloc,
		   Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
		   Matchpool_T matchpool, Genome_T genome, Genome_T genomealt, Pairpool_T pairpool,
		   Diagpool_T diagpool, Cellpool_T cellpool,
		   Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, int min_matches) {
  List_T gregions = NULL, p;
  Stage3_T new_left, new_right;
  Stage3_T *stage3array_sub1 = NULL, *stage3array_sub2 = NULL, from, to, stage3;
  Sequence_T querysubseq = NULL, querysubuc = NULL;
  Diagnostic_T diagnostic;
  int bestfrom, bestto;
  int five_margin, three_margin, five_score = 0, three_score = 0;
  int extension;
  int npaths_sub1 = 0, npaths_sub2 = 0;
  bool lowidentityp, poorp, repetitivep;

  int max_single_goodness, chimeric_goodness, penalty, matches0, matches1;
  int breakpoint, chimerapos, chimeraequivpos, exonexonpos;
  int chimera_cdna_direction;
  char donor1, donor2, acceptor2, acceptor1;
  bool donor_watsonp, acceptor_watsonp;
  double donor_prob, acceptor_prob;
  

  debug2(printf("check_for_chimera called with %d paths\n",List_length(stage3list)));
#ifdef DEBUG2
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    printf("%p %p\n",stage3,Stage3_pairs(stage3));
  }
#endif



#ifdef PMAP
  five_margin = effective_start - 3*Sequence_trim_start(queryseq);
  three_margin = 3*Sequence_trim_end(queryseq) - effective_end;
  debug2(printf("Margins are %d = %d - %d on the 5' end and %d = %d - %d on the 3' end\n",
		five_margin,effective_start,3*Sequence_trim_start(queryseq),
		three_margin,3*Sequence_trim_end(queryseq),effective_end));
#else
  five_margin = effective_start - Sequence_trim_start(queryseq);
  three_margin = Sequence_trim_end(queryseq) - effective_end;
  debug2(printf("Margins are %d = %d - %d on the 5' end and %d = %d - %d on the 3' end\n",
		five_margin,effective_start,Sequence_trim_start(queryseq),
		three_margin,Sequence_trim_end(queryseq),effective_end));
#endif

#ifdef DEBUG2A
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    Pair_dump_array(Stage3_pairarray(stage3),Stage3_npairs(stage3),/*zerobasedp*/true);
    printf("\n");
  }
#endif

  /* Stage3_recompute_goodness(stage3list); */
  max_single_goodness = 0;
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    if (Stage3_goodness(stage3) > max_single_goodness) {
      max_single_goodness = Stage3_goodness(stage3);
    }
  }
  debug2(printf("max single goodness = %d\n",max_single_goodness));


  debug2(printf("Running distant_separate_paths\n"));
  distant_separate_paths(&stage3array_sub1,&npaths_sub1,&stage3array_sub2,&npaths_sub2,
			 stage3list);
  debug2(printf("chimera: npaths_sub1 %d, npaths_sub2 %d, stage3list %d\n",
		npaths_sub1,npaths_sub2,List_length(stage3list)));

  if (npaths_sub1 == 0 && npaths_sub2 == 0) {
    /* Need to compute on margin explicitly */
    if (five_margin < chimera_margin && three_margin < chimera_margin) {
      debug2(printf("Insufficient margins\n"));
    } else if (five_margin > three_margin) {
      extension = CHIMERA_SLOP;
      debug2(printf("Comparing extension %d with %d = (effective_start %d)/2\n",
		    extension,effective_start/2,effective_start));
      if (extension > effective_start/2) {
	/* Extension occupies more than 1/3 of sequence */
	debug2(printf("Proposed extension of %d is too long relative to effective_start %d\n",extension,effective_start));
	extension = effective_start/3;
      }
      if ((querysubseq = Sequence_subsequence(queryseq,0,effective_start+extension)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,0,effective_start+extension)) != NULL) {
	  debug2(printf("5 margin > 3 margin.  "));
	  debug2(printf("Beginning Stage1_compute on 5' margin from effective_start %d (%d..%d)\n",
			effective_start,0,effective_start+extension));
	  debug2a(Sequence_stdout(querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(querysubuc),Sequence_fulllength(querysubuc),
				      Oligoindex_array_elt(oligoindices_major,0));
	  if (poorp == true || repetitivep == true) {
	    debug2(printf("Subsequence is poor or repetitive\n"));
	  } else {
	    if (mode == CMET_NONSTRANDED || mode == ATOI_NONSTRANDED || mode == TTOC_NONSTRANDED) {
	      gregions = Stage1_compute_nonstranded(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,
						    chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
						    stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    } else {
	      gregions = Stage1_compute(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,/*genestrand*/0,
					chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
					stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    }
	    debug2(printf("A.  Performing Stage 3 starting with list length %d\n",List_length(stage3list)));
	    stage3list = stage3_from_gregions(stage3list,gregions,querysubseq,querysubuc,
#ifdef PMAP
					      queryntseq,
#endif
					      stage2_alloc,oligoindices_major,oligoindices_minor,
					      genome,genomealt,pairpool,diagpool,cellpool,
					      dynprogL,dynprogM,dynprogR,min_matches,/*worker_stopwatch*/NULL);
#ifdef DEBUG2
	    for (p = stage3list; p != NULL; p = List_next(p)) {
	      stage3 = (Stage3_T) List_head(p);
	      printf("%d..%d, %u..%u\n",
		     Stage3_querystart(stage3),Stage3_queryend(stage3),
		     Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
	    }
#endif
	  }
	  Diagnostic_free(&diagnostic);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);
	}
	Sequence_free(&querysubseq);
      }

      /* And recompute on original part, just in case stage 1 was led astray by the ends */
      if ((querysubseq = Sequence_subsequence(queryseq,effective_start,queryntlength)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,effective_start,queryntlength)) != NULL) {
	  debug2(printf("Recomputing on original part.  "));
	  debug2(printf("Beginning Stage1_compute on 5' margin from effective_start %d (%d..%d)\n",
			effective_start,effective_start,queryntlength));
	  debug2a(Sequence_stdout(querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(querysubuc),Sequence_fulllength(querysubuc),
				      Oligoindex_array_elt(oligoindices_major,0));
	  if (poorp == true || repetitivep == true) {
	    debug2(printf("Subsequence is poor or repetitive\n"));
	  } else {
	    if (mode == CMET_NONSTRANDED || mode == ATOI_NONSTRANDED || mode == TTOC_NONSTRANDED) {
	      gregions = Stage1_compute_nonstranded(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,
						    chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
						    stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    } else {
	      gregions = Stage1_compute(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,/*genestrand*/0,
					chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
					stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    }
	    debug2(printf("B.  Performing Stage 3 starting with list length %d\n",List_length(stage3list)));
	    stage3list = stage3_from_gregions(stage3list,gregions,querysubseq,querysubuc,
#ifdef PMAP
					      queryntseq,
#endif
					      stage2_alloc,oligoindices_major,oligoindices_minor,
					      genome,genomealt,pairpool,diagpool,cellpool,
					      dynprogL,dynprogM,dynprogR,min_matches,/*worker_stopwatch*/NULL);
#ifdef DEBUG2
	    for (p = stage3list; p != NULL; p = List_next(p)) {
	      stage3 = (Stage3_T) List_head(p);
	      printf("%d..%d, %u..%u\n",
		     Stage3_querystart(stage3),Stage3_queryend(stage3),
		     Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
	    }
#endif
	  }
	  Diagnostic_free(&diagnostic);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);
	}
	Sequence_free(&querysubseq);
      }

      debug2(printf("Running distant_separate_paths\n"));
      distant_separate_paths(&stage3array_sub1,&npaths_sub1,&stage3array_sub2,&npaths_sub2,
			     stage3list);
      debug2(printf("chimera: npaths_sub1 %d, npaths_sub2 %d, stage3list %d\n",
		    npaths_sub1,npaths_sub2,List_length(stage3list)));

    } else {
      extension = CHIMERA_SLOP;
      debug2(printf("Comparing extension %d with %d = (queryntlength %d - effective_end %d)/2\n",
		    extension,(queryntlength-effective_end)/2,queryntlength,effective_end));
      if (extension > (queryntlength - effective_end)/2) {
	/* Extension occupies more than 1/3 of sequence */
	debug2(printf("Proposed extension of %d is too long relative to queryntlength %d and effective_end %d\n",
		      extension,queryntlength,effective_end));
	extension = (queryntlength - effective_end)/3;
      }
      if ((querysubseq = Sequence_subsequence(queryseq,effective_end-extension,queryntlength)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,effective_end-extension,queryntlength)) != NULL) {
	  debug2(printf("5 margin <= 3 margin.  "));
	  debug2(printf("Beginning Stage1_compute on 3' margin from effective_end %d (%d..%d)\n",
			effective_end,effective_end-extension,queryntlength));
	  debug2(Sequence_stdout(querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(querysubuc),Sequence_fulllength(querysubuc),
				      Oligoindex_array_elt(oligoindices_major,0));
	  if (poorp == true || repetitivep == true) {
	    debug2(printf("Subsequence is poor or repetitive\n"));
	  } else {
	    if (mode == CMET_NONSTRANDED || mode == ATOI_NONSTRANDED || mode == TTOC_NONSTRANDED) {
	      gregions = Stage1_compute_nonstranded(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,
						    chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
						    stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    } else {
	      gregions = Stage1_compute(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,/*genestrand*/0,
					chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
					stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    }
	    debug2(printf("C.  Performing Stage 3 with list length %d\n",List_length(stage3list)));
	    stage3list = stage3_from_gregions(stage3list,gregions,querysubseq,querysubuc,
#ifdef PMAP
					      queryntseq,
#endif
					      stage2_alloc,oligoindices_major,oligoindices_minor,
					      genome,genomealt,pairpool,diagpool,cellpool,
					      dynprogL,dynprogM,dynprogR,min_matches,/*worker_stopwatch*/NULL);
#ifdef DEBUG2
	    for (p = stage3list; p != NULL; p = List_next(p)) {
	      stage3 = (Stage3_T) List_head(p);
	      printf("%d..%d, %u..%u\n",
		     Stage3_querystart(stage3),Stage3_queryend(stage3),
		     Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
	    }
#endif
	  }
	  Diagnostic_free(&diagnostic);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);
	}
	Sequence_free(&querysubseq);
      }

      /* And recompute on original part, just in case stage 1 was led astray by the ends */
      if ((querysubseq = Sequence_subsequence(queryseq,0,effective_end)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,0,effective_end)) != NULL) {
	  debug2(printf("Recomputing on original part.  "));
	  debug2(printf("Beginning Stage1_compute on 3' margin from effective_end %d (%d..%d)\n",
			effective_end,0,effective_end));
	  debug2(Sequence_stdout(querysubseq,/*uppercasep*/true,wraplength,/*trimmedp*/true));

	  diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(querysubuc),Sequence_fulllength(querysubuc),
				      Oligoindex_array_elt(oligoindices_major,0));
	  if (poorp == true || repetitivep == true) {
	    debug2(printf("Subsequence is poor or repetitive\n"));
	  } else {
	    if (mode == CMET_NONSTRANDED || mode == ATOI_NONSTRANDED || mode == TTOC_NONSTRANDED) {
	      gregions = Stage1_compute_nonstranded(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,
						    chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
						    stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    } else {
	      gregions = Stage1_compute(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,/*genestrand*/0,
					chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
					stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    }
	    debug2(printf("D.  Performing Stage 3 with list length %d\n",List_length(stage3list)));
	    stage3list = stage3_from_gregions(stage3list,gregions,querysubseq,querysubuc,
#ifdef PMAP
					      queryntseq,
#endif
					      stage2_alloc,oligoindices_major,oligoindices_minor,
					      genome,genomealt,pairpool,diagpool,cellpool,
					      dynprogL,dynprogM,dynprogR,min_matches,/*worker_stopwatch*/NULL);
#ifdef DEBUG2
	    for (p = stage3list; p != NULL; p = List_next(p)) {
	      stage3 = (Stage3_T) List_head(p);
	      printf("%d..%d, %u..%u\n",
		     Stage3_querystart(stage3),Stage3_queryend(stage3),
		     Stage3_genomicstart(stage3),Stage3_genomicend(stage3));
	    }
#endif
	  }
	  Diagnostic_free(&diagnostic);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);

	}
	Sequence_free(&querysubseq);
      }

      debug2(printf("Running distant_separate_paths\n"));
      distant_separate_paths(&stage3array_sub1,&npaths_sub1,&stage3array_sub2,&npaths_sub2,
			     stage3list);
      debug2(printf("chimera: npaths_sub1 %d, npaths_sub2 %d, stage3list %d\n",
		    npaths_sub1,npaths_sub2,List_length(stage3list)));
    }
  }

  *mergedp = false;
  *chimera = (Chimera_T) NULL;
  if (npaths_sub1 == 0 || npaths_sub2 == 0) {
    /* Skip */

  } else if (Chimera_bestpath(&five_score,&three_score,&chimerapos,&chimeraequivpos,&bestfrom,&bestto,
			      stage3array_sub1,npaths_sub1,stage3array_sub2,npaths_sub2,queryntlength,
			      CHIMERA_SLOP,circularp,/*localp*/false) == false) {
    /* Skip */
    debug2(printf("Chimera_bestpath returns false, so skipping\n"));
    FREE(stage3array_sub2);
    FREE(stage3array_sub1);

  } else {
    from = stage3array_sub1[bestfrom];
    to = stage3array_sub2[bestto];
    debug2(printf("Chimera_bestpath returns bestfrom %d (%d..%d, %u..%u) to bestto %d (%d..%d, %u..%u)\n",
		  bestfrom,Stage3_querystart(from),Stage3_queryend(from),Stage3_genomicstart(from),Stage3_genomicend(from),
		  bestto,Stage3_querystart(to),Stage3_queryend(to),Stage3_genomicstart(to),Stage3_genomicend(to)));

    chimeric_goodness = Stage3_chimeric_goodness(&matches0,&matches1,from,to,chimerapos);
    debug2(printf("chimeric goodness = %d\n",chimeric_goodness));
    
    penalty = CHIMERA_PENALTY;
    if (chimera_margin < penalty) {
      /* User is looking for higher sensitivity */
      penalty = chimera_margin;
    }

    if (chimeric_goodness < max_single_goodness + penalty) {
      debug2(printf("chimeric goodness not good enough relative to max_single_goodness %d and penalty %d\n",
		    max_single_goodness,penalty));

    } else if ((breakpoint = find_breakpoint(&chimera_cdna_direction,&chimerapos,&chimeraequivpos,&exonexonpos,
					     &donor1,&donor2,&acceptor2,&acceptor1,
					     &donor_watsonp,&acceptor_watsonp,&donor_prob,&acceptor_prob,from,to,
#ifdef PMAP
					     queryntseq,
#endif
					     queryseq,queryuc,queryntlength,
					     genome,genomealt,pairpool,chromosome_iit)) <= 0) {
      debug2(printf("find_breakpoint returns no value\n"));

    } else {
      debug2(printf("find_breakpoint returns %d\n",breakpoint));

      /* Check to see if we can merge chimeric parts */
      debug2(printf("Before Stage3_mergeable, bestfrom is %p, query %d..%d, pairs %p\n",
		    from,Stage3_querystart(from),Stage3_queryend(from),Stage3_pairs(from)));
      debug2(printf("Before Stage3_mergeable, bestto is %p, query %d..%d, pairs %p\n",
		    to,Stage3_querystart(to),Stage3_queryend(to),Stage3_pairs(to)));

      if (maxpaths_report != 1 && /* if maxpaths_report == 1, then don't want distant chimeras */
	  Stage3_mergeable(from,to,breakpoint,queryntlength) == false &&
	  Stage3_test_bounds(from,0,chimeraequivpos+chimera_overlap) == true &&
	  Stage3_test_bounds(to,chimerapos+1-chimera_overlap,queryntlength) == true &&
	  Stage3_merge_chimera(&new_left,&new_right,/*best0*/from,/*best1*/to,
			       /*minpos1*/0,/*maxpos1*/breakpoint,
			       /*minpos2*/breakpoint+1,/*maxpos2*/queryntlength,queryseq,
#ifdef PMAP
			       Sequence_fullpointer(queryntseq),Sequence_fullpointer(queryntseq),
#else
			       Sequence_fullpointer(queryseq),Sequence_fullpointer(queryuc),
#endif
			       genome,genomealt,pairpool,dynprogL,dynprogR,maxpeelback) == true) {

	debug2(printf("Not mergeable -- Merging left and right as a transloc\n"));
	*chimera = Chimera_new(new_left,new_right,chimerapos,chimeraequivpos,exonexonpos,chimera_cdna_direction,
			       donor1,donor2,acceptor2,acceptor1,donor_watsonp,acceptor_watsonp,
			       donor_prob,acceptor_prob);

	debug2(printf("Before merge_left_and_right_transloc, bestfrom is %p, query %d..%d\n",
		      from,Stage3_querystart(from),Stage3_queryend(from)));
	debug2(printf("Before merge_left_and_right_transloc, bestto is %p, query %d..%d\n",
		      to,Stage3_querystart(to),Stage3_queryend(to)));
	
	/* Used to call merge_left_and_right_transloc */
	for (p = stage3list; p != NULL; p = List_next(p)) {
	  stage3 = (Stage3_T) List_head(p);
	  Stage3_free(&stage3);
	}
	List_free(&stage3list);

	stage3list = List_push(NULL,(void *) new_right);
	stage3list = List_push(stage3list,(void *) new_left);
      }

      debug2(printf("After Stage3_mergeable, bestfrom is %p, query %d..%d, pairs %p\n",
		    from,Stage3_querystart(from),Stage3_queryend(from),Stage3_pairs(from)));
      debug2(printf("After Stage3_mergeable, bestto is %p, query %d..%d, pairs %p\n",
		    to,Stage3_querystart(to),Stage3_queryend(to),Stage3_pairs(to)));
    }

    FREE(stage3array_sub2);
    FREE(stage3array_sub1);
  }

  debug2(printf("check_for_chimera returning list of length %d\n",List_length(stage3list)));
#ifdef DEBUG2
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    printf("%p %p\n",stage3,Stage3_pairs(stage3));
  }
#endif

#if 0
  /* Should be handled by apply_stage3 loop */
  /* Needed after calls to stage3_from_gregions */
  Stage3_recompute_goodness(stage3list);
  stage3list = stage3list_remove_duplicates(stage3list);
#endif

  return stage3list;
}


/* Needs to guarantee that all elements of stage3list and middlepieces end up in result */
/* The Stage3_T objects from and to come from stage3list */
/* The Stage3_T object middle does not come from stage3list (but from middlepieces in caller) */
static List_T
merge_middlepieces (List_T stage3list, Stage3_T from, Stage3_T to, Stage3_T middle,
		    bool mergeableAp, bool mergeableBp,
		    int breakpointA, int breakpointB, Sequence_T queryseq,
#ifdef PMAP
		    Sequence_T queryntseq,
#endif
		    Sequence_T queryuc, int queryntlength,
		    Genome_T genome, Genome_T genomealt, Pairpool_T pairpool,
		    Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		    Oligoindex_array_T oligoindices_minor, Diagpool_T diagpool, Cellpool_T cellpool) {
  Stage3_T stage3_merged;

  debug2(printf("Entered merge_middlepieces\n"));

  if (mergeableAp == true && mergeableBp == true) {
    if ((stage3_merged =
	 merge_left_and_right_readthrough(/*stage3array_sub1*/&from,/*npaths_sub1:1,*//*bestfrom*/0,
					  /*stage3array_sub2*/&middle,/*npaths_sub2:1,*//*bestto*/0,
					  breakpointA,queryntlength,queryseq,
#ifdef PMAP
					  /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
					  /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
					  /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
#else
					  /*queryseq_ptr*/Sequence_fullpointer(queryseq),
					  /*queryuc_ptr*/Sequence_fullpointer(queryuc),
#endif
					  genome,genomealt,pairpool,dynprogL,dynprogM,dynprogR,
					  oligoindices_minor,diagpool,cellpool)) == NULL) {
      debug2(printf("Merge of from/middle failed, so trying merge of middle/to\n"));
      if ((stage3_merged =
	   merge_left_and_right_readthrough(/*stage3array_sub1*/&middle,/*npaths_sub1:1,*//*bestfrom*/0,
					    /*stage3array_sub2*/&to,/*npaths_sub2:1,*//*bestto*/0,
					    breakpointB,queryntlength,queryseq,
#ifdef PMAP
					    /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
					    /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
					    /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
#else
					    /*queryseq_ptr*/Sequence_fullpointer(queryseq),
					    /*queryuc_ptr*/Sequence_fullpointer(queryuc),
#endif
					    genome,genomealt,pairpool,dynprogL,dynprogM,dynprogR,
					    oligoindices_minor,diagpool,cellpool)) == NULL) {
	debug2(printf("And merge of middle/to also failed\n"));
      } else {
	debug2(printf("But merge of middle/to succeeded\n"));
	stage3list = List_push(stage3list,(void *) stage3_merged);
      }

    } else {
      debug2(printf("Merge of from/middle succeeded.  Now merging result/to\n"));
      stage3list = List_push(stage3list,(void *) stage3_merged);
      if ((stage3_merged =
	   merge_left_and_right_readthrough(/*stage3array_sub1*/&stage3_merged,/*npaths_sub1:1,*//*bestfrom*/0,
					    /*stage3array_sub2*/&to,/*npaths_sub2:1,*//*bestto*/0,
					    breakpointB,queryntlength,queryseq,
#ifdef PMAP
					    /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
					    /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
					    /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
#else
					    /*queryseq_ptr*/Sequence_fullpointer(queryseq),
					    /*queryuc_ptr*/Sequence_fullpointer(queryuc),
#endif
					    genome,genomealt,pairpool,dynprogL,dynprogM,dynprogR,
					    oligoindices_minor,diagpool,cellpool)) == NULL) {
	debug2(printf("But merge of result/to failed\n"));
      } else {
	debug2(printf("And merge of result/to succeeded\n"));
	stage3list = List_push(stage3list,(void *) stage3_merged);
      }
    }

  } else if (mergeableBp == true) {
    if ((stage3_merged =
	 merge_left_and_right_readthrough(/*stage3array_sub1*/&middle,/*npaths_sub1:1,*//*bestfrom*/0,
					  /*stage3array_sub2*/&to,/*npaths_sub2:1,*//*bestto*/0,
					  breakpointB,queryntlength,queryseq,
#ifdef PMAP
					  /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
					  /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
					  /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
#else
					  /*queryseq_ptr*/Sequence_fullpointer(queryseq),
					  /*queryuc_ptr*/Sequence_fullpointer(queryuc),
#endif
					  genome,genomealt,pairpool,dynprogL,dynprogM,dynprogR,
					  oligoindices_minor,diagpool,cellpool)) != NULL) {
      stage3list = List_push(stage3list,(void *) stage3_merged);
    }

  } else if (mergeableAp == true) {
    if ((stage3_merged =
	 merge_left_and_right_readthrough(/*stage3array_sub1*/&from,/*npaths_sub1:1,*//*bestfrom*/0,
					  /*stage3array_sub2*/&middle,/*npaths_sub2:1,*//*bestto*/0,
					  breakpointA,queryntlength,queryseq,
#ifdef PMAP
					  /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
					  /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
					  /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
#else
					  /*queryseq_ptr*/Sequence_fullpointer(queryseq),
					  /*queryuc_ptr*/Sequence_fullpointer(queryuc),
#endif
					  genome,genomealt,pairpool,dynprogL,dynprogM,dynprogR,
					  oligoindices_minor,diagpool,cellpool)) != NULL) {
      stage3list = List_push(stage3list,(void *) stage3_merged);
    }
  }

  return stage3list;
}



/* Returns stage3list with additional merged alignments and middle pieces */
static List_T
check_middle_piece_local (bool *foundp, List_T stage3list, Sequence_T queryseq, Sequence_T queryuc,
#ifdef PMAP
			  Sequence_T queryntseq,
#endif
			  int queryntlength, Stage2_alloc_T stage2_alloc,
			  Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
			  Genome_T genome, Genome_T genomealt, Pairpool_T pairpool,
			  Diagpool_T diagpool, Cellpool_T cellpool,
			  Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, int min_matches) {
  Sequence_T querysubseq = NULL, querysubuc = NULL;
  int npaths, i, j;
  Stage3_T from = NULL, to = NULL, middle = NULL;
  Stage3_T *by_queryend, *by_querystart;
  List_T r;
  bool plusp;
  int genestrand;

  int querystart, queryend;
  Chrpos_T chrstart, chrend, chrlength;
  Univcoord_T chroffset, chrhigh;
  Chrnum_T chrnum;

  int breakpointA = 0, chimeraposA, chimeraequivposA, exonexonposA;
  char donorA1, donorA2, acceptorA2, acceptorA1;
  bool donor_watsonp_A, acceptor_watsonp_A;
  double donor_prob_A, acceptor_prob_A;

  int breakpointB = 0, chimeraposB, chimeraequivposB, exonexonposB;
  char donorB1, donorB2, acceptorB2, acceptorB1;
  bool donor_watsonp_B, acceptor_watsonp_B;
  double donor_prob_B, acceptor_prob_B;

  int chimera_cdna_direction_A, chimera_cdna_direction_B;
  bool mergeableAp, mergeableBp;

  List_T all_middlepieces = NULL, middlepieces;
#ifdef DEBUG2A
  List_T p;
  Stage3_T stage3;
#endif


#ifdef DEBUG2A
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    Pair_dump_array(Stage3_pairarray(stage3),Stage3_npairs(stage3),/*zerobasedp*/true);
    printf("\n");
  }
#endif

  *foundp = false;

  by_queryend = (Stage3_T *) List_to_array_out_n(&npaths,stage3list);
  qsort(by_queryend,npaths,sizeof(Stage3_T),Stage3_queryend_cmp);

  by_querystart = (Stage3_T *) List_to_array_out_n(&npaths,stage3list);
  qsort(by_querystart,npaths,sizeof(Stage3_T),Stage3_querystart_cmp);

  j = 0;
  for (i = 0; i < npaths && *foundp == false; i++) {
    from = by_queryend[i];
    queryend = Stage3_queryend(from);

    while (j < npaths && Stage3_querystart(by_querystart[j]) < queryend) {
      j++;
    }
    j--;

    while (j >= 0 && Stage3_querystart(by_querystart[j]) > queryend) {
      j--;
    }
    j++;

    for ( ; j < npaths && *foundp == false; j++) {
      to = by_querystart[j];

      if (middle_piece_local_p(&querystart,&queryend,&chrstart,&chrend,
			       &chrnum,&chroffset,&chrhigh,&chrlength,&plusp,&genestrand,
			       from,to) == true) {
	debug2(printf("Found middle piece missing from %d to %d\n",i,j));

	if ((querysubseq = Sequence_subsequence(queryseq,querystart,queryend)) != NULL) {
	  if ((querysubuc = Sequence_subsequence(queryuc,querystart,queryend)) != NULL) {
	    debug2(printf("Performing Stage 3 on %d..%d against %u..%u\n",
			  querystart,queryend,chrstart,chrend));
	    if ((middlepieces = update_stage3list(/*stage3list*/NULL,querysubseq,
#ifdef PMAP
						  queryntseq,
#endif
						  querysubuc,stage2_alloc,oligoindices_major,oligoindices_minor,
						  genome,genomealt,pairpool,diagpool,cellpool,
						  /*straintype*/0,/*strain*/NULL,chrnum,
						  chroffset,chrhigh,chrlength,chrstart,chrend,plusp,genestrand,
						  dynprogL,dynprogM,dynprogR,min_matches,/*worker_stopwatch*/NULL)) != NULL) {
	      middlepieces = stage3list_sort(middlepieces);

	      /* 1.  Look first for middle piece that joins locally on both ends */
	      r = middlepieces;
	      mergeableAp = mergeableBp = false;
	      while (r != NULL && (mergeableAp == false || mergeableBp == false)) {
		middle = (Stage3_T) List_head(r);
		if (Chimera_local_join_p(from,middle,CHIMERA_SLOP) == true && Chimera_local_join_p(middle,to,CHIMERA_SLOP) == true) {
		  if ((breakpointA = find_breakpoint(&chimera_cdna_direction_A,&chimeraposA,&chimeraequivposA,&exonexonposA,
						     &donorA1,&donorA2,&acceptorA2,&acceptorA1,
						     &donor_watsonp_A,&acceptor_watsonp_A,&donor_prob_A,&acceptor_prob_A,
						     from,/*to*/middle,
#ifdef PMAP
						     queryntseq,
#endif
						     queryseq,queryuc,queryntlength,
						     genome,genomealt,pairpool,chromosome_iit)) <= 0) {
		    mergeableAp = false;
		  } else {
		    mergeableAp = Stage3_mergeable(from,/*to*/middle,breakpointA,queryntlength);
		  }

		  if ((breakpointB = find_breakpoint(&chimera_cdna_direction_B,&chimeraposB,&chimeraequivposB,&exonexonposB,
						     &donorB1,&donorB2,&acceptorB2,&acceptorB1,
						     &donor_watsonp_B,&acceptor_watsonp_B,&donor_prob_B,&acceptor_prob_B,
						     /*from*/middle,to,
#ifdef PMAP
						     queryntseq,
#endif
						     queryseq,queryuc,queryntlength,
						     genome,genomealt,pairpool,chromosome_iit)) <= 0) {
		    mergeableBp = false;
		  } else {
		    mergeableBp = Stage3_mergeable(/*from*/middle,to,breakpointB,queryntlength);
		  }
		}
		r = List_next(r);
	      }	/* End of while loop looking for dual merge */

	      if (mergeableAp == true && mergeableBp == true) {
		debug2(printf("Middle segment %p found and mergeable locally with both! -- Merging three as a readthrough.\n",middle));
		*foundp = true;
	      } else {
		/* 2.  Look for middle piece that joins locally on one end */
		r = middlepieces;
		mergeableAp = mergeableBp = false;
		while (r != NULL && mergeableAp == false && mergeableBp == false) {
		  middle = (Stage3_T) List_head(r);
		  if (Chimera_local_join_p(from,middle,CHIMERA_SLOP) == true && Chimera_local_join_p(middle,to,CHIMERA_SLOP) == true) {
		    if ((breakpointA = find_breakpoint(&chimera_cdna_direction_A,&chimeraposA,&chimeraequivposA,&exonexonposA,
						       &donorA1,&donorA2,&acceptorA2,&acceptorA1,
						       &donor_watsonp_A,&acceptor_watsonp_A,&donor_prob_A,&acceptor_prob_A,
						       from,/*to*/middle,
#ifdef PMAP
						       queryntseq,
#endif
						       queryseq,queryuc,queryntlength,
						       genome,genomealt,pairpool,chromosome_iit)) <= 0) {
		      mergeableAp = false;
		    } else {
		      mergeableAp = Stage3_mergeable(from,/*to*/middle,breakpointA,queryntlength);
		    }

		    if ((breakpointB = find_breakpoint(&chimera_cdna_direction_B,&chimeraposB,&chimeraequivposB,&exonexonposB,
						       &donorB1,&donorB2,&acceptorB2,&acceptorB1,
						       &donor_watsonp_B,&acceptor_watsonp_B,&donor_prob_B,&acceptor_prob_B,
						       /*from*/middle,to,
#ifdef PMAP
						       queryntseq,
#endif
						       queryseq,queryuc,queryntlength,
						       genome,genomealt,pairpool,chromosome_iit)) <= 0) {
		      mergeableBp = false;
		    } else {
		      mergeableBp = Stage3_mergeable(/*from*/middle,to,breakpointB,queryntlength);
		    }
		  }
		  r = List_next(r);
		} /* End of while loop looking for single merge */

		if (mergeableAp == true || mergeableBp == true) {
		  *foundp = true;
		}
	      }

	      stage3list = merge_middlepieces(stage3list,from,to,middle,mergeableAp,mergeableBp,
					      breakpointA,breakpointB,queryseq,
#ifdef PMAP
					      queryntseq,
#endif
					      queryuc,queryntlength,genome,genomealt,pairpool,
					      dynprogL,dynprogM,dynprogR,
					      oligoindices_minor,diagpool,cellpool);
	      all_middlepieces = List_append(all_middlepieces,middlepieces);
	    }

	    Sequence_free(&querysubuc);
	  }
	  Sequence_free(&querysubseq);
	}
      }
    }
  }

  FREE(by_querystart);
  FREE(by_queryend);

  stage3list = List_append(stage3list,all_middlepieces);

  return stage3list;
}


/* Returns stage3list with additional merged alignments and middle pieces */
static List_T
check_middle_piece_chimera (bool *foundp, List_T stage3list, Sequence_T queryseq, Sequence_T queryuc,
#ifdef PMAP
			    Sequence_T queryntseq,
#endif
			    int queryntlength, Stage2_alloc_T stage2_alloc,
			    Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
			    Matchpool_T matchpool, Genome_T genome, Genome_T genomealt, Pairpool_T pairpool,
			    Diagpool_T diagpool, Cellpool_T cellpool,
			    Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, int min_matches) {
  Sequence_T querysubseq = NULL, querysubuc = NULL;
  int npaths, i, j;
  Stage3_T bestfrom, bestto, from, to, middle, stage3_merged;
  Stage3_T *by_queryend, *by_querystart;
  List_T r;
  int querystart, queryend, maxdist, dist;

  int breakpointA, chimeraposA, chimeraequivposA, exonexonposA;
  char donorA1, donorA2, acceptorA2, acceptorA1;
  bool donor_watsonp_A, acceptor_watsonp_A;
  double donor_prob_A, acceptor_prob_A;

  int breakpointB, chimeraposB, chimeraequivposB, exonexonposB;
  char donorB1, donorB2, acceptorB2, acceptorB1;
  bool donor_watsonp_B, acceptor_watsonp_B;
  double donor_prob_B, acceptor_prob_B;

  int chimera_cdna_direction_A, chimera_cdna_direction_B;
  bool mergeableAp, mergeableBp;

  List_T middlepieces = NULL;
  Diagnostic_T diagnostic;
  List_T gregions;
  bool lowidentityp, poorp, repetitivep;

#ifdef DEBUG2A
  List_T p;
  Stage3_T stage3;
#endif


#ifdef DEBUG2A
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    Pair_dump_array(Stage3_pairarray(stage3),Stage3_npairs(stage3),/*zerobasedp*/true);
    printf("\n");
  }
#endif

  by_queryend = (Stage3_T *) List_to_array_out_n(&npaths,stage3list);
  qsort(by_queryend,npaths,sizeof(Stage3_T),Stage3_queryend_cmp);

  by_querystart = (Stage3_T *) List_to_array_out_n(&npaths,stage3list);
  qsort(by_querystart,npaths,sizeof(Stage3_T),Stage3_querystart_cmp);

  maxdist = 0;
  j = 0;
  for (i = 0; i < npaths; i++) {
    from = by_queryend[i];
    queryend = Stage3_queryend(from);

    while (j < npaths && Stage3_querystart(by_querystart[j]) < queryend) {
      j++;
    }
    j--;

    while (j >= 0 && Stage3_querystart(by_querystart[j]) > queryend) {
      j--;
    }
    j++;

    if (j < npaths) {
      /* Should have the first querystart just after queryend */
      to = by_querystart[j];

      if ((dist = Stage3_queryend(to) - Stage3_querystart(from)) > maxdist) {
	bestfrom = from;
	bestto = to;
	maxdist = dist;
      }
    }
  }

  FREE(by_querystart);
  FREE(by_queryend);


  *foundp = false;
  if (maxdist < CHIMERA_SLOP) {
    debug2(printf("maxdist %d < CHIMERA_SLOP %d\n",maxdist,CHIMERA_SLOP));
  } else {
    if (middle_piece_chimera_p(&querystart,&queryend,bestfrom,bestto) == true) {
      if ((querysubseq = Sequence_subsequence(queryseq,querystart,queryend)) != NULL) {
	if ((querysubuc = Sequence_subsequence(queryuc,querystart,queryend)) != NULL) {
	  debug2(printf("Performing Stage 3 on %d..%d\n",querystart,queryend));

	  diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(querysubuc),
				      Sequence_fulllength(querysubuc),Oligoindex_array_elt(oligoindices_major,0));
	  if (poorp == true || repetitivep == true) {
	    debug2(printf("Subsequence is poor or repetitive\n"));
	  } else {
	    if (mode == CMET_NONSTRANDED || mode == ATOI_NONSTRANDED || mode == TTOC_NONSTRANDED) {
	      gregions = Stage1_compute_nonstranded(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,
						    chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
						    stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    } else {
	      gregions = Stage1_compute(&lowidentityp,querysubuc,indexdb_fwd,indexdb_rev,/*genestrand*/0,
					chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
					stutterhits,diagnostic,/*worker_stopwatch*/NULL,/*nbest*/10);
	    }
	    debug2(printf("Performing Stage 3 starting with list length %d\n",List_length(stage3list)));
	    middlepieces = stage3_from_gregions(/*stage3list*/NULL,gregions,querysubseq,querysubuc,
#ifdef PMAP
						queryntseq,
#endif
						stage2_alloc,oligoindices_major,oligoindices_minor,
						genome,genomealt,pairpool,diagpool,cellpool,
						dynprogL,dynprogM,dynprogR,min_matches,/*worker_stopwatch*/NULL);
	  }
	  Diagnostic_free(&diagnostic);

	  /* Above function frees gregions */
	  Sequence_free(&querysubuc);
	}
	Sequence_free(&querysubseq);
      }
    }

    if (middlepieces != NULL) {
      middlepieces = stage3list_sort(middlepieces);

      r = middlepieces;
      mergeableAp = mergeableBp = false;
      while (r != NULL && mergeableAp == false && mergeableBp == false) {
	middle = (Stage3_T) List_head(r);
	if (middle != bestfrom && middle != bestto) {
	  if (Chimera_local_join_p(bestfrom,middle,CHIMERA_SLOP) == true) {
	    if ((breakpointA = find_breakpoint(&chimera_cdna_direction_A,&chimeraposA,&chimeraequivposA,&exonexonposA,
					       &donorA1,&donorA2,&acceptorA2,&acceptorA1,
					       &donor_watsonp_A,&acceptor_watsonp_A,&donor_prob_A,&acceptor_prob_A,
					       bestfrom,/*to*/middle,
#ifdef PMAP
					       queryntseq,
#endif
					       queryseq,queryuc,queryntlength,
					       genome,genomealt,pairpool,chromosome_iit)) <= 0) {
	      mergeableAp = false;
	    } else {
	      mergeableAp = Stage3_mergeable(bestfrom,/*to*/middle,breakpointA,queryntlength);
	    }
	  }
	  if (Chimera_local_join_p(middle,bestto,CHIMERA_SLOP) == true) {
	    if ((breakpointB = find_breakpoint(&chimera_cdna_direction_B,&chimeraposB,&chimeraequivposB,&exonexonposB,
					       &donorB1,&donorB2,&acceptorB2,&acceptorB1,
					       &donor_watsonp_B,&acceptor_watsonp_B,&donor_prob_B,&acceptor_prob_B,
					       /*from*/middle,to,
#ifdef PMAP
					       queryntseq,
#endif
					       queryseq,queryuc,queryntlength,
					       genome,genomealt,pairpool,chromosome_iit)) <= 0) {
	      mergeableBp = false;
	    } else {
	      mergeableBp = Stage3_mergeable(/*from*/middle,bestto,breakpointB,queryntlength);
	    }
	  }
	}
	r = List_next(r);
      }

      if (mergeableAp == true) {
	debug2(printf("Middle segment %p found and mergeable locally with from! -- Merging as a readthrough.  cdna_direction = %d\n",
 	              middle,chimera_cdna_direction_A));
	if ((stage3_merged =
	     merge_left_and_right_readthrough(/*stage3array_sub1*/&bestfrom,/*npaths_sub1:1,*//*bestfrom*/0,
					      /*stage3array_sub2*/&middle,/*npaths_sub2:1,*//*bestto*/0,
					      breakpointA,queryntlength,queryseq,
#ifdef PMAP
					      /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
					      /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
					      /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
#else
					      /*queryseq_ptr*/Sequence_fullpointer(queryseq),
					      /*queryuc_ptr*/Sequence_fullpointer(queryuc),
#endif
					      genome,genomealt,pairpool,dynprogL,dynprogM,dynprogR,
					      oligoindices_minor,diagpool,cellpool)) != NULL) {
	  stage3list = List_push(stage3list,(void *) stage3_merged);
	  *foundp = true;
	}

      } else if (mergeableBp == true) {
	debug2(printf("Middle segment %p found and mergeable locally with to! -- Merging as a readthrough.  cdna_direction = %d\n",
		      middle,chimera_cdna_direction_B));
	if ((stage3_merged = 
	     merge_left_and_right_readthrough(/*stage3array_sub1*/&middle,/*npaths_sub1:1,*//*bestfrom*/0,
					      /*stage3array_sub2*/&bestto,/*npaths_sub2:1,*//*bestto*/0,
					      breakpointB,queryntlength,queryseq,
#ifdef PMAP
					      /*queryaaseq_ptr*/Sequence_fullpointer(queryseq),
					      /*queryseq_ptr*/Sequence_fullpointer(queryntseq),
					      /*queryuc_ptr*/Sequence_fullpointer(queryntseq),
#else
					      /*queryseq_ptr*/Sequence_fullpointer(queryseq),
					      /*queryuc_ptr*/Sequence_fullpointer(queryuc),
#endif
					      genome,genomealt,pairpool,dynprogL,dynprogM,dynprogR,
					      oligoindices_minor,diagpool,cellpool)) != NULL) {
	  stage3list = List_push(stage3list,(void *) stage3_merged);
	  *foundp = true;
	}

      } else {
	debug2(printf("Middle segment found but notmergeable\n"));
      }

    }
  }

  stage3list = List_append(stage3list,middlepieces);

  return stage3list;
}



static List_T
stage3_with_stage1 (bool *mergedp, Chimera_T *chimera, List_T gregions, Sequence_T queryseq, Sequence_T queryuc,
#ifdef PMAP
		    Sequence_T queryntseq,
#endif
		    Stage2_alloc_T stage2_alloc,
		    Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
		    Matchpool_T matchpool, Genome_T genome, Genome_T genomealt, Pairpool_T pairpool,
		    Diagpool_T diagpool, Cellpool_T cellpool,
		    Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR, int min_matches,
		    Stopwatch_T worker_stopwatch) {
  List_T stage3list, newstage3list, split_objects, p, q;
  Stage3_T nonchimericbest, chimera1, chimera2, stage3, newstage3;
  bool testlocalp, testchimerap, foundp;
  int effective_start, effective_end;
  int queryntlength;
  int iter;

  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh;
  Chrpos_T chrlength;
  List_T pairs_below, pairs_above;
  bool watsonp, chimerap;
  int cdna_direction, genestrand, sensedir;

  
  *mergedp = false;
  *chimera = NULL;

  debug(printf("Calling stage3_from_gregions\n"));
  stage3list = stage3_from_gregions(/*stage3list*/(List_T) NULL,gregions,queryseq,queryuc,
#ifdef PMAP
				    queryntseq,
#endif
				    stage2_alloc,oligoindices_major,oligoindices_minor,
				    genome,genomealt,pairpool,diagpool,cellpool,
				    dynprogL,dynprogM,dynprogR,min_matches,worker_stopwatch);

  debug2(printf("Initial search gives stage3list of length %d\n",List_length(stage3list)));
#ifdef DEBUG2
  for (p = stage3list; p != NULL; p = List_next(p)) {
    Stage3_print_ends(List_head(p));
  }
#endif

  if (diag_debug == true) {
    return stage3list;		/* really diagonals */
  }

  queryntlength = Sequence_ntlength(queryseq);

  if (stage3list != NULL) {
    iter = 0;
    testlocalp = true;
    while (testlocalp == true && iter++ < MAX_CHIMERA_ITER) {
      debug2(printf("\n\n*** Testing for local on %d Stage3_T objects, iter %d ***\n",
		    List_length(stage3list),iter));

      /* Stage3_recompute_goodness(stage3list); */
      /* stage3list = stage3list_remove_duplicates(stage3list); */
      stage3list = stage3list_sort(stage3list);

#ifdef DEBUG2
      for (p = stage3list; p != NULL; p = List_next(p)) {
	Stage3_print_ends(List_head(p));
      }
      printf("\n");
#endif
      nonchimericbest = (Stage3_T) List_head(stage3list);
      debug2(printf("nonchimericbest is %p\n",nonchimericbest));

#if 0
      if (List_length(stage3list) <= 1) {
	debug2(printf("Only 0 or 1 alignments, so won't look for local\n"));
	testlocalp = false;
      }
      else 
#endif

      if (Stage3_domain(nonchimericbest) < chimera_margin) {
	debug2(printf("Existing alignment is too short, so won't look for local\n"));
	testlocalp = false;

#if 0
      } else if (Stage3_fracidentity(nonchimericbest) < CHIMERA_IDENTITY &&
		 Chimera_alignment_break(&effective_start,&effective_end,nonchimericbest,Sequence_ntlength(queryseq),CHIMERA_FVALUE) >= chimera_margin
		 ) {
	debug2(printf("Break in alignment quality at %d..%d detected, so will look for local\n",
		      effective_start,effective_end));
	testlocalp = true;
#endif

      } else if (Stage3_largemargin(&effective_start,&effective_end,nonchimericbest,Sequence_ntlength(queryseq)) >= chimera_margin) {
	debug2(printf("Large margin at %d..%d detected (%d >= %d), so will look for local\n",
		      effective_start,effective_end,Stage3_largemargin(&effective_start,&effective_end,nonchimericbest,Sequence_ntlength(queryseq)),chimera_margin));
	testlocalp = true;
	
      } else {
	debug2(printf("Good alignment already with identity %f, so won't look for local\n",
		      Stage3_fracidentity(nonchimericbest)));
	testlocalp = false;
      }

      if (testlocalp == true) {
	testlocalp = false;
	debug2(printf("Checking for local, starting with list length %d, effective_start %d, effective_end %d\n",
		      List_length(stage3list),effective_start,effective_end));
	stage3list = check_for_local(&(*mergedp),stage3list,effective_start,effective_end,
				     queryseq,queryuc,
#ifdef PMAP
				     queryntseq,
#endif
				     queryntlength,stage2_alloc,oligoindices_major,oligoindices_minor,
				     matchpool,genome,genomealt,pairpool,
				     diagpool,cellpool,dynprogL,dynprogM,dynprogR,min_matches);
	debug2(printf("After check for local, we still have %d paths\n",List_length(stage3list)));
	
#if 0
	/* For some reason, we need to filter out cases where npairs is 0 */
	old = stage3list;
	stage3list = (List_T) NULL;
	for (p = old; p != NULL; p = List_next(p)) {
          stage3 = (Stage3_T) List_head(p);
          if (Stage3_npairs(stage3) == 0) {
            Stage3_free(&stage3);
          } else {
            stage3list = List_push(stage3list,(void *) stage3);
          }
        }
	List_free(&old);
#endif

	if (*mergedp == true) {
	  testlocalp = true;	/* Local merge */
	} else if (iter == 1) {
	  /* Check for middle pieces only on first iteration */
	  debug2(printf("Checking for middle piece local, starting with list length %d\n",List_length(stage3list)));
	  stage3list = check_middle_piece_local(&foundp,stage3list,queryseq,queryuc,
#ifdef PMAP
						queryntseq,
#endif
						queryntlength,stage2_alloc,oligoindices_major,oligoindices_minor,
						genome,genomealt,pairpool,diagpool,cellpool,
						dynprogL,dynprogM,dynprogR,min_matches);
	  if (foundp == true) {
	    /* Iterate */
	    testlocalp = true;
	  }
	} else {
	  testlocalp = false;
	}
      }
    }
  }

  if (stage3list != NULL) {
    iter = 0;
    testchimerap = true;
    while (testchimerap == true && iter++ < MAX_CHIMERA_ITER) {
      debug2(printf("\n\n*** Testing for chimera on %d Stage3_T objects, iter %d ***\n",
		    List_length(stage3list),iter));

      /* Stage3_recompute_goodness(stage3list); */
      /* stage3list = stage3list_remove_duplicates(stage3list); */
      stage3list = stage3list_sort(stage3list);

#ifdef DEBUG2
      for (p = stage3list; p != NULL; p = List_next(p)) {
	Stage3_print_ends(List_head(p));
      }
      printf("\n");
#endif
      nonchimericbest = (Stage3_T) List_head(stage3list);
      debug2(printf("nonchimericbest is %p\n",nonchimericbest));

      if (novelsplicingp == false) {
	testchimerap = false;

      } else if (chimera_margin <= 0) {
	debug2(printf("turned off\n"));
	testchimerap = false;

      } else if (maxpaths_report == 1) {
	debug2(printf("maxpaths set to 1\n"));
	testchimerap = false;

      } else if (Stage3_domain(nonchimericbest) < chimera_margin) {
	debug2(printf("Existing alignment is too short, so won't look for chimera\n"));
	testchimerap = false;

#if 0
      } else if (Stage3_fracidentity(nonchimericbest) < CHIMERA_IDENTITY &&
		 Chimera_alignment_break(&effective_start,&effective_end,nonchimericbest,Sequence_ntlength(queryseq),CHIMERA_FVALUE) >= chimera_margin
		 ) {
	debug2(printf("Break in alignment quality at %d..%d detected, so will look for chimera\n",
		      effective_start,effective_end));
	testchimerap = true;
#endif

      } else if (Stage3_largemargin(&effective_start,&effective_end,nonchimericbest,Sequence_ntlength(queryseq)) >= chimera_margin) {
	debug2(printf("Large margin at %d..%d detected (%d >= %d), so will look for chimera\n",
		      effective_start,effective_end,Stage3_largemargin(&effective_start,&effective_end,nonchimericbest,Sequence_ntlength(queryseq)),chimera_margin));
	testchimerap = true;
	
      } else {
	debug2(printf("Good alignment already with identity %f, so won't look for chimera\n",
		      Stage3_fracidentity(nonchimericbest)));
	testchimerap = false;
      }

      if (testchimerap == true) {
	testchimerap = false;
	debug2(printf("Checking for chimera, starting with list length %d, effective_start %d, effective_end %d\n",
		      List_length(stage3list),effective_start,effective_end));
	stage3list = check_for_chimera(&(*mergedp),&(*chimera),stage3list,effective_start,effective_end,
				       queryseq,queryuc,
#ifdef PMAP
				       queryntseq,
#endif
				       queryntlength,stage2_alloc,oligoindices_major,oligoindices_minor,
				       matchpool,genome,genomealt,pairpool,diagpool,cellpool,
				       dynprogL,dynprogM,dynprogR,min_matches);
	debug2(printf("chimera is %p\n",*chimera));
	if (*chimera != NULL) {
	  testchimerap = false;
	} else {
	  if (*mergedp == true) {
	    testchimerap = true;	/* Local merge */
	  } else if (iter == 1) {
	    /* Check for middle pieces only on first iteration */
	    debug2(printf("Checking for middle piece chimera, starting with list length %d\n",List_length(stage3list)));
	    stage3list = check_middle_piece_chimera(&foundp,stage3list,queryseq,queryuc,
#ifdef PMAP
						    queryntseq,
#endif
						    queryntlength,stage2_alloc,oligoindices_major,oligoindices_minor,
						    matchpool,genome,genomealt,pairpool,diagpool,cellpool,
						    dynprogL,dynprogM,dynprogR,min_matches);
	    if (foundp == true) {
	      /* Iterate */
	      testchimerap = true;
	    } else {
	      testchimerap = false;
	    }
	  } else {
	    testchimerap = false;
	  }
	}
	debug2(printf("testchimerap is %d\n",testchimerap));
      }
    }
  }

  debug2(printf("apply_stage3 returning list of length %d\n",List_length(stage3list)));
#ifdef DEBUG2
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    printf("%p %p\n",stage3,Stage3_pairs(stage3));
  }
#endif

  /* Split on large introns */
  if (split_large_introns_p == true) {
    /* Need to check whether stage3 belongs to a chimera */
    chimerap = false;
    for (p = stage3list; p != NULL; p = List_next(p)) {
      stage3 = (Stage3_T) List_head(p);
      if (Stage3_chimera_left_p(stage3) == true || Stage3_chimera_right_p(stage3) == true) {
	chimerap = true;
      }
    }

    if (chimerap == false) {
      newstage3list = (List_T) NULL;
      for (p = stage3list; p != NULL; p = List_next(p)) {
	stage3 = (Stage3_T) List_head(p);
	if (Stage3_npairs(stage3) == 0) {
	  Stage3_free(&stage3);
	} else if ((split_objects = Stage3_split(stage3,queryseq,genome,genomealt,pairpool)) == NULL) {
	  debug(printf("Pushing %p onto newstage3list\n",stage3));
	  newstage3list = List_push(newstage3list,(void *) stage3);
	} else {
	  for (q = split_objects; q != NULL; q = List_next(q)) {
	    newstage3 = (Stage3_T) List_head(q);
	    debug(printf("Pushing %p onto newstage3list\n",newstage3));
	    newstage3list = List_push(newstage3list,(void *) newstage3);
	  }
	  List_free(&split_objects);
	  Stage3_free(&stage3);
	}
      }
      List_free(&stage3list);
      stage3list = newstage3list;
    }
  }


  /* Split circular alignments (need to guarantee that chimeras do not
     contain alignments to circular chromosomes) */
  newstage3list = (List_T) NULL;
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    chrnum = Stage3_chrnum(stage3);
    if (circularp[chrnum] == false) {
      newstage3list = List_push(newstage3list,(void *) stage3);
    } else {
      chroffset = Stage3_chroffset(stage3);
      chrhigh = Stage3_chrhigh(stage3);
      chrlength = Stage3_chrlength(stage3);
      watsonp = Stage3_watsonp(stage3);

      Pair_split_circular(&pairs_below,&pairs_above,Stage3_pairs(stage3),
			  chrlength,pairpool,watsonp);
#if 0
      printf("PAIRS BELOW\n");
      Pair_dump_list(pairs_below,true);
      printf("PAIRS ABOVE\n");
      Pair_dump_list(pairs_above,true);
#endif
	  
      cdna_direction = Stage3_cdna_direction(stage3);
      genestrand = Stage3_genestrand(stage3);
      sensedir = Stage3_sensedir(stage3);

      if ((newstage3 = Stage3_new_from_pairs(pairs_below,cdna_direction,watsonp,genestrand,sensedir,
					     genome,genomealt,pairpool,queryseq,/*query_subseq_offset*/0,
					     chrnum,chroffset,chrhigh,chrlength)) != NULL) {
	debug(printf("Pushing %p onto stage3list\n",newstage3));
	newstage3list = List_push(newstage3list,(void *) newstage3);
      }
      if ((newstage3 = Stage3_new_from_pairs(pairs_above,cdna_direction,watsonp,genestrand,sensedir,
					     genome,genomealt,pairpool,queryseq,/*query_subseq_offset*/0,
					     chrnum,chroffset,chrhigh,chrlength)) != NULL) {
	debug(printf("Pushing %p onto stage3list\n",newstage3));
	newstage3list = List_push(newstage3list,(void *) newstage3);
      }
      Stage3_free(&stage3);
    }
  }
  List_free(&stage3list);
  stage3list = newstage3list;


  /* Needed after call to stage3_from_gregions */
  /* Stage3_recompute_goodness(stage3list); */

  /* Final call, so do both filtering and sorting */
  Stage3_recompute_coverage(stage3list,queryseq);
  stage3list = stage3list_filter_and_sort(&(*chimera),stage3list);
  debug2(printf("After filter and sort, have %d paths\n",List_length(stage3list)));

  if (*chimera != NULL && List_length(stage3list) > 2) {
    /* Compare chimera against non-chimeric alignments */
    chimera1 = (Stage3_T) List_head(stage3list);
    chimera2 = (Stage3_T) List_head(List_next(stage3list));
    nonchimericbest = (Stage3_T) List_head(List_next(List_next(stage3list)));
    debug2(printf("chimera1 %d, chimera2 %d\n",Stage3_goodness(chimera1),Stage3_goodness(chimera2)));
    debug2(printf("%p non-chimeric %d %d..%d\n",
		  nonchimericbest,Stage3_goodness(nonchimericbest),Stage3_querystart(nonchimericbest),Stage3_queryend(nonchimericbest)));

    if (Stage3_queryend(nonchimericbest) > (Stage3_querystart(chimera2) + Stage3_queryend(chimera2))/2 &&
	Stage3_querystart(nonchimericbest) < (Stage3_querystart(chimera1) + Stage3_queryend(chimera1))/2) {
      stage3list = List_pop(stage3list,(void **) &chimera1);
      stage3list = List_pop(stage3list,(void **) &chimera2);
      Stage3_free(&chimera1);
      Stage3_free(&chimera2);
      Chimera_free(&(*chimera));
      *chimera = (Chimera_T) NULL;
    }
  }

  debug2(printf("apply_stage3 returning %d paths\n",List_length(stage3list)));
#ifdef DEBUG2
  for (p = stage3list; p != NULL; p = List_next(p)) {
    stage3 = (Stage3_T) List_head(p);
    printf("%p %p\n",stage3,Stage3_pairs(stage3));
  }
#endif

  return stage3list;
}


static Filestring_T
process_request (Filestring_T *fp_failedinput, double *worker_runtime, Request_T request,
		 Matchpool_T matchpool, Genome_T genome, Genome_T genomealt,
		 Pairpool_T pairpool, Diagpool_T diagpool, Cellpool_T cellpool,
		 Stage2_alloc_T stage2_alloc, Oligoindex_array_T oligoindices_major, Oligoindex_array_T oligoindices_minor,
		 Dynprog_T dynprogL, Dynprog_T dynprogM, Dynprog_T dynprogR,
		 Stopwatch_T worker_stopwatch) {
  Filestring_T fp;
  Result_T result;
  int jobid;
  Diagnostic_T diagnostic;
  Sequence_T queryseq, queryuc;
  Chimera_T chimera = NULL;
  bool mergedp, lowidentityp;
  bool repetitivep = false, poorp = false;

  List_T gregions = NULL, stage3list;
  Stage3_T *stage3array;
  int npaths_primary, npaths_altloc, first_absmq, second_absmq;
#ifdef PMAP
  Sequence_T queryntseq;
#endif

  jobid = Request_id(request);
  queryseq = Request_queryseq(request);
  Matchpool_reset(matchpool);
  Pairpool_reset(pairpool);
  Diagpool_reset(diagpool);
  Cellpool_reset(cellpool);


  if (worker_stopwatch != NULL) {
    Stopwatch_start(worker_stopwatch);
  }

  if (Sequence_fulllength_given(queryseq) <= 0) {
    result = Result_new(jobid,/*mergedp*/false,(Chimera_T) NULL,(Stage3_T *) NULL,
			/*npaths_primary*/0,/*npaths_altloc*/0,/*first_absmq*/0,/*second_absmq*/0,
			/*diagnostic*/NULL,EMPTY_SEQUENCE);
      
  } else if (Sequence_fulllength_given(queryseq) < 
#ifdef PMAP
	     index1part_aa
#else
	     index1part
#endif
	     ) {
    result = Result_new(jobid,/*mergedp*/false,(Chimera_T) NULL,(Stage3_T *) NULL,
			/*npaths_primary*/0,/*npaths_altloc*/0,/*first_absmq*/0,/*second_absmq*/0,
			/*diagnostic*/NULL,SHORT_SEQUENCE);

  } else {			/* Sequence_fulllength_given(queryseq) > 0 */
    queryuc = Sequence_uppercase(queryseq);
#ifdef PMAP
    queryntseq = Sequence_convert_to_nucleotides(queryseq);
#endif

    diagnostic = evaluate_query(&poorp,&repetitivep,Sequence_fullpointer(queryuc),
				Sequence_fulllength(queryuc),Oligoindex_array_elt(oligoindices_major,0));

#ifndef PMAP
    if (poorp == true && prune_poor_p == true) {
      result = Result_new(jobid,/*mergedp*/false,(Chimera_T) NULL,(Stage3_T *) NULL,
			  /*npaths_primary*/0,/*npaths_altloc*/0,/*first_absmq*/0,/*second_absmq*/0,
			  diagnostic,POOR_SEQUENCE);
    } else if (repetitivep == true && prune_repetitive_p == true) {
      result = Result_new(jobid,/*mergedp*/false,(Chimera_T) NULL,(Stage3_T *) NULL,
			  /*npaths_primary*/0,/*npaths_altloc*/0,/*first_absmq*/0,/*second_absmq*/0,
			  diagnostic,REPETITIVE);
    }
#endif

    if (user_selfalign_p == true) {
      stage3array = stage3_self_align(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,
				      queryseq,pairpool);
      result = Result_new(jobid,/*mergedp*/false,(Chimera_T) NULL,stage3array,npaths_primary,npaths_altloc,
			  first_absmq,second_absmq,diagnostic,NO_FAILURE);
      
    } else if (indexdb_fwd == NULL) {
      /* Skip stage 1 */
      stage3array = stage3_skip_stage1(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,queryseq,queryuc,
#ifdef PMAP
				       queryntseq,
#endif
				       stage2_alloc,oligoindices_major,oligoindices_minor,
				       genome,genomealt,pairpool,diagpool,cellpool,
				       dynprogL,dynprogM,dynprogR,worker_stopwatch);

      result = Result_new(jobid,/*mergedp*/false,(Chimera_T) NULL,stage3array,npaths_primary,npaths_altloc,
			  first_absmq,second_absmq,diagnostic,NO_FAILURE);

    } else {		/* Not user segment and not maponly */
#ifndef PMAP
#if 0
      /* Don't do Sequence_trim, because it affects sequences like NM_018406 */
      Sequence_trim(queryseq,diagnostic->query_trim_start,diagnostic->query_trim_end);
      Sequence_trim(queryuc,diagnostic->query_trim_start,diagnostic->query_trim_end);
#endif
#endif

      debug(printf("Calling stage 1\n"));
      if (mode == CMET_NONSTRANDED || mode == ATOI_NONSTRANDED || mode == TTOC_NONSTRANDED) {
	gregions = Stage1_compute_nonstranded(&lowidentityp,queryuc,indexdb_fwd,indexdb_fwd,
					      chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
					      stutterhits,diagnostic,worker_stopwatch,/*nbest*/10);

      } else {
	gregions = Stage1_compute(&lowidentityp,queryuc,indexdb_fwd,indexdb_rev,/*genestrand*/0,
				  chromosome_iit,chrsubset_start,chrsubset_end,matchpool,
				  stutterhits,diagnostic,worker_stopwatch,/*nbest*/10);
      }
      debug(printf("Got %d gregions\n",List_length(gregions)));

      if (stage1debug == true) {
	/* result = Result_new_stage1debug(jobid,gregions,diagnostic,NO_FAILURE); */
	abort();
      } else {
	debug(printf("Applying stage 3\n"));
	stage3list = stage3_with_stage1(&mergedp,&chimera,gregions,queryseq,queryuc,
#ifdef PMAP
					queryntseq,
#endif
					stage2_alloc,oligoindices_major,oligoindices_minor,
					matchpool,genome,genomealt,pairpool,diagpool,cellpool,
					dynprogL,dynprogM,dynprogR,/*min_matches*/MIN_MATCHES,worker_stopwatch);
	if (diag_debug == true) {
#if 0
	  result = Result_new_diag_debug(jobid,/*diagonals*/stage3list,diagnostic,NO_FAILURE);
#endif
	  abort();
	} else if (stage3list == NULL) {
	  result = Result_new(jobid,mergedp,chimera,/*stage3array*/NULL,/*npaths_primary*/0,/*npaths_altloc*/0,
			      /*first_absmq*/0,/*second_absmq*/0,diagnostic,NO_FAILURE);
	} else if (chimera == NULL) {
	  stage3array = stage3array_from_list(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,
					      stage3list,/*chimerap*/false,/*remove_overlaps_p*/true);
	  debug2(printf("chimera is NULL.  npaths_primary %d, npaths_altloc %d\n",npaths_primary,npaths_altloc));
	  result = Result_new(jobid,mergedp,/*chimera*/NULL,stage3array,npaths_primary,npaths_altloc,
			      first_absmq,second_absmq,diagnostic,NO_FAILURE);
	} else {
	  stage3array = stage3array_from_list(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,
					      stage3list,/*chimerap*/true,/*remove_overlaps_p*/false);
	  debug2(printf("chimera is not NULL.  npaths_primary %d, npaths_altloc %d\n",npaths_primary,npaths_altloc));
	  result = Result_new(jobid,mergedp,chimera,stage3array,npaths_primary,npaths_altloc,
			      first_absmq,second_absmq,diagnostic,NO_FAILURE);
	}
      }

      Oligoindex_clear_inquery(Oligoindex_array_elt(oligoindices_major,0),/*queryuc_ptr*/Sequence_fullpointer(queryuc),
			       /*querystart*/0,/*queryend*/Sequence_fulllength(queryuc));

    } /* Matches not user segment and not maponly */

#ifdef PMAP
    Sequence_free(&queryntseq);
#endif
    Sequence_free(&queryuc);
  } /* Matches sequence length > 0 */

  fp = Output_filestring_fromresult(&(*fp_failedinput),result,request,/*headerseq*/queryseq);
  *worker_runtime = worker_stopwatch == NULL ? 0.00 : Stopwatch_stop(worker_stopwatch);
  Result_free(&result);
  return fp;
}


#ifdef HAVE_SIGACTION
static const Except_T sigfpe_error = {"SIGFPE--arithmetic exception"};
static const Except_T sigsegv_error = {"SIGSEGV--segmentation violation"};
static const Except_T sigtrap_error = {"SIGTRAP--hardware fault"};
static const Except_T misc_signal_error = {"Miscellaneous signal"};

static void
signal_handler (int sig) {
  Request_T request;
  Sequence_T queryseq;

  if (sig == SIGPIPE) {
    /* Allow pipe */
    return;
  }

  switch (sig) {
  case SIGABRT: fprintf(stderr,"Signal received: SIGABRT\n"); break;
  case SIGFPE: fprintf(stderr,"Signal received: SIGFPE\n"); break;
  case SIGHUP: fprintf(stderr,"Signal received: SIGHUP\n"); break;
  case SIGILL:
    fprintf(stderr,"Signal received: SIGILL\n");
    fprintf(stderr,"An illegal instruction means that this program is being run on a computer\n");
    fprintf(stderr,"  with different features than the computer used to compile the program\n");
    fprintf(stderr,"You may need to re-compile the program on the same computer type as the target machine\n");
    fprintf(stderr,"  or re-compile with fewer features by doing something like\n");
    fprintf(stderr,"  ./configure --disable-simd\n");
    break;
  case SIGINT: fprintf(stderr,"Signal received: SIGINT\n"); break;
  case SIGPIPE: fprintf(stderr,"Signal received: SIGPIPE\n"); break;
  case SIGQUIT: fprintf(stderr,"Signal received: SIGQUIT\n"); break;
  case SIGSEGV: fprintf(stderr,"Signal received: SIGSEGV\n"); break;
  case SIGSYS: fprintf(stderr,"Signal received: SIGSYS\n"); break;
  case SIGTERM: fprintf(stderr,"Signal received: SIGTERM\n"); break;
  case SIGTRAP: fprintf(stderr,"Signal received: SIGTRAP\n"); break;
  case SIGXCPU: fprintf(stderr,"Signal received: SIGXCPU\n"); break;
  case SIGXFSZ: fprintf(stderr,"Signal received: SIGXFSZ\n"); break;
  }

  Access_emergency_cleanup();

#ifdef HAVE_PTHREAD
  request = (Request_T) pthread_getspecific(global_request_key);
  if (request == NULL) {
    /* fprintf(stderr,"Unable to retrieve request for thread\n"); */
  } else {
    queryseq = Request_queryseq(request);
    if (queryseq == NULL) {
      fprintf(stderr,"Unable to retrieve queryseq for request\n");
    } else {
      fprintf(stderr,"Problem sequence: ");
      fprintf(stderr,"%s (%d bp)\n",Sequence_accession(queryseq),Sequence_fulllength(queryseq));
    }
  }
#endif

  exit(9);

  return;
}
#endif


#define POOL_FREE_INTERVAL 200

static void
single_thread () {
  Stage2_alloc_T stage2_alloc;
  Oligoindex_array_T oligoindices_major, oligoindices_minor;
  Dynprog_T dynprogL, dynprogM, dynprogR;
  Matchpool_T matchpool;
  Pairpool_T pairpool;
  Diagpool_T diagpool;
  Cellpool_T cellpool;
  Stopwatch_T worker_stopwatch;
  Request_T request;
  Genome_T genome, genomealt;
  Filestring_T fp, fp_failedinput;
  Sequence_T queryseq;
  int jobid = 0;
  double worker_runtime;

#ifdef MEMUSAGE
  long int memusage, memusage_constant = 0;
  char acc[100+1], comma0[20], comma1[20], comma2[20], comma3[20], comma4[20], comma5[20];
#endif

  stage2_alloc = Stage2_alloc_new(MAX_QUERYLENGTH_FOR_ALLOC);
  oligoindices_major = Oligoindex_array_new_major(MAX_QUERYLENGTH_FOR_ALLOC,MAX_GENOMICLENGTH_FOR_ALLOC);
  oligoindices_minor = Oligoindex_array_new_minor(MAX_QUERYLENGTH_FOR_ALLOC,MAX_GENOMICLENGTH_FOR_ALLOC);
  dynprogL = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired,
			 /*doublep*/true);
  dynprogM = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired,
			 /*doublep*/false);
  dynprogR = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired,
			 /*doublep*/true);
  matchpool = Matchpool_new();
  pairpool = Pairpool_new();
  diagpool = Diagpool_new();
  cellpool = Cellpool_new();
  worker_stopwatch = (timingp == true) ? Stopwatch_new() : (Stopwatch_T) NULL;

  /* Except_stack_create(); -- requires pthreads */

#ifdef MEMUSAGE
  memusage_constant += Mem_usage_report_std_heap();
  Genomicpos_commafmt_fill(comma0,memusage_constant);
  Mem_usage_reset_heap_baseline(0);
#endif

  while ((request = Inbuffer_get_request(inbuffer)) != NULL) {
    genome = Request_genome(request);
    genomealt = Request_genomealt(request);

#ifdef MEMUSAGE
    queryseq = Request_queryseq(request);
    fprintf(stderr,"Single thread starting %s\n",Sequence_accession(queryseq));
    Mem_usage_reset_stack_max();
    Mem_usage_reset_heap_max();
#endif

    TRY
      fp = process_request(&fp_failedinput,&worker_runtime,request,
			   matchpool,genome,genomealt,pairpool,diagpool,cellpool,
			   stage2_alloc,oligoindices_major,oligoindices_minor,
			   dynprogL,dynprogM,dynprogR,worker_stopwatch);
      if (timingp == true) {
        queryseq = Request_queryseq(request);
        fprintf(stderr,"%s\t%.6f\n",Sequence_accession(queryseq),worker_runtime);
      }

    ELSE
      queryseq = Request_queryseq(request);
      if (Sequence_accession(queryseq) == NULL) {
        fprintf(stderr,"Problem with unnamed sequence (%d bp)\n",Sequence_fulllength_given(queryseq));
      } else {
        fprintf(stderr,"Problem with sequence %s (%d bp)\n",
  	      Sequence_accession(queryseq),Sequence_fulllength_given(queryseq));
      }
      fprintf(stderr,"To obtain a core dump, re-run program on problem sequence with the -0 [zero] flag\n");
      fprintf(stderr,"Exiting...\n");
      exit(9);
    RERAISE;
    END_TRY;

    /* Only a single thread, but this way Outbuffer_print_filestrings can assume Filestring_stringify has been called */
    Filestring_stringify(fp);
    if (fp_failedinput != NULL) {
      Filestring_stringify(fp_failedinput);
    }
    Outbuffer_print_filestrings(fp,fp_failedinput);

    if (jobid % POOL_FREE_INTERVAL == 0) {
      Pairpool_free_memory(pairpool);
      Diagpool_free_memory(diagpool);
      Cellpool_free_memory(cellpool);
      Matchpool_free_memory(matchpool);
    }

#ifdef MEMUSAGE
    /* Copy acc before we free the request */
    queryseq = Request_queryseq(request);
    strncpy(acc,Sequence_accession(queryseq),100);
    acc[100] = '\0';
#endif

    Request_free(&request);

#ifdef MEMUSAGE
    Genomicpos_commafmt_fill(comma1,Mem_usage_report_std_heap_max());
    Genomicpos_commafmt_fill(comma2,Mem_usage_report_std_heap());
    Genomicpos_commafmt_fill(comma3,Mem_usage_report_keep());
    Genomicpos_commafmt_fill(comma4,Mem_usage_report_in());
    Genomicpos_commafmt_fill(comma5,Mem_usage_report_out());

    fprintf(stderr,"Acc %s: constant %s  max %s  std %s  keep %s  in %s  out %s\n",
	    acc,comma0,comma1,comma2,comma3,comma4,comma5);

    if ((memusage = Mem_usage_report_std_heap()) != 0) {
      fprintf(stderr,"Memory leak in single thread of %ld bytes\n",memusage);
      fflush(stdout);
      exit(9);
    }
#endif
  }

#ifdef MEMUSAGE
  Mem_usage_std_heap_add(memusage_constant);
#endif

  /* Except_stack_destroy(); -- requires pthreads */

  if (worker_stopwatch != NULL) {
    Stopwatch_free(&worker_stopwatch);
  }
  Cellpool_free(&cellpool);
  Diagpool_free(&diagpool);
  Pairpool_free(&pairpool);
  Matchpool_free(&matchpool);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_array_free(&oligoindices_minor);
  Oligoindex_array_free(&oligoindices_major);
  Stage2_alloc_free(&stage2_alloc);

#ifdef MEMUSAGE
  Mem_usage_set_threadname("main");
#endif

  return;
}


#ifdef HAVE_PTHREAD
static void *
worker_thread (void *data) {
  Stage2_alloc_T stage2_alloc;
  Oligoindex_array_T oligoindices_major, oligoindices_minor;
  Dynprog_T dynprogL, dynprogM, dynprogR;
  Matchpool_T matchpool;
  Pairpool_T pairpool;
  Diagpool_T diagpool;
  Cellpool_T cellpool;
  Stopwatch_T worker_stopwatch;
  Request_T request;
  Genome_T genome, genomealt;
  Filestring_T fp, fp_failedinput;
  Sequence_T queryseq;
  int worker_jobid = 0;
  double worker_runtime;
#if defined(DEBUG) || defined(MEMUSAGE)
  long int worker_id = (long int) data;
#endif

#ifdef MEMUSAGE
  long int memusage_constant = 0, memusage, max_memusage;
  char threadname[12];
  char acc[100+1], comma0[20], comma1[20], comma2[20], comma3[20], comma4[20], comma5[20];
  sprintf(threadname,"thread-%ld",worker_id);
  Mem_usage_set_threadname(threadname);
#endif

  /* Thread-specific data and storage */
  stage2_alloc = Stage2_alloc_new(MAX_QUERYLENGTH_FOR_ALLOC);
  oligoindices_major = Oligoindex_array_new_major(MAX_QUERYLENGTH_FOR_ALLOC,MAX_GENOMICLENGTH_FOR_ALLOC);
  oligoindices_minor = Oligoindex_array_new_minor(MAX_QUERYLENGTH_FOR_ALLOC,MAX_GENOMICLENGTH_FOR_ALLOC);
  dynprogL = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired,
			 /*doublep*/true);
  dynprogM = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired,
			 /*doublep*/false);
  dynprogR = Dynprog_new(nullgap,EXTRAQUERYGAP,maxpeelback,extramaterial_end,extramaterial_paired,
			 /*doublep*/true);
  matchpool = Matchpool_new();
  pairpool = Pairpool_new();
  diagpool = Diagpool_new();
  cellpool = Cellpool_new();
  worker_stopwatch = (timingp == true) ? Stopwatch_new() : (Stopwatch_T) NULL;

  Except_stack_create();

#ifdef MEMUSAGE
  memusage_constant += Mem_usage_report_std_heap();
  Genomicpos_commafmt_fill(comma0,memusage_constant);
  Mem_usage_reset_heap_baseline(0);
#endif

  while ((request = Inbuffer_get_request(inbuffer)) != NULL) {
    debug(printf("worker_thread %ld got request %d\n",worker_id,Request_id(request)));
    pthread_setspecific(global_request_key,(void *) request);

    genome = Request_genome(request);
    genomealt = Request_genomealt(request);

#ifdef MEMUSAGE
    queryseq = Request_queryseq(request);
    fprintf(stderr,"Thread %d starting %s\n",worker_id,Sequence_accession(queryseq));
    Mem_usage_reset_stack_max();
    Mem_usage_reset_heap_max();
#endif

    TRY
      fp = process_request(&fp_failedinput,&worker_runtime,request,
			   matchpool,genome,genomealt,pairpool,diagpool,cellpool,
			   stage2_alloc,oligoindices_major,oligoindices_minor,
			   dynprogL,dynprogM,dynprogR,worker_stopwatch);
      if (timingp == true) {
        queryseq = Request_queryseq(request);
        fprintf(stderr,"%s\t%.6f\n",Sequence_accession(queryseq),worker_runtime);
      }

    ELSE
      queryseq = Request_queryseq(request);
      if (queryseq == NULL) {
	fprintf(stderr,"NULL");
      } else if (Sequence_accession(queryseq) == NULL) {
	fprintf(stderr,"unnamed (%d bp)",Sequence_fulllength_given(queryseq));
      } else {
	fprintf(stderr,"%s (%d bp)",Sequence_accession(queryseq),Sequence_fulllength_given(queryseq));
      }
      fprintf(stderr,"\n");
      fprintf(stderr,"To obtain a core dump, re-run program on problem sequence with the -0 [zero] flag\n");

      fprintf(stderr,"Exiting...\n");
      exit(9);
    RERAISE;
    END_TRY;

    debug(printf("worker_thread %ld putting filestring for request %d\n",worker_id,Request_id(request)));

    /* Parallelize the stringify operation by performing by worker thread and not the output thread */
    Filestring_stringify(fp);
    if (fp_failedinput != NULL) {
      Filestring_stringify(fp_failedinput);
    }

    Outbuffer_put_filestrings(outbuffer,Request_id(request),fp,fp_failedinput);

    if (worker_jobid % POOL_FREE_INTERVAL == 0) {
      Pairpool_free_memory(pairpool);
      Diagpool_free_memory(diagpool);
      Cellpool_free_memory(cellpool);
      Matchpool_free_memory(matchpool);
    }

#ifdef MEMUSAGE
    /* Copy acc before we free the request */
    queryseq = Request_queryseq(request);
    strncpy(acc,Sequence_accession(queryseq),100);
    acc[100] = '\0';
#endif

    Request_free(&request);

#ifdef MEMUSAGE
    Genomicpos_commafmt_fill(comma1,Mem_usage_report_std_heap_max());
    Genomicpos_commafmt_fill(comma2,Mem_usage_report_std_heap());
    Genomicpos_commafmt_fill(comma3,Mem_usage_report_keep());
    Genomicpos_commafmt_fill(comma4,Mem_usage_report_in());
    Genomicpos_commafmt_fill(comma5,Mem_usage_report_out());

    fprintf(stderr,"Acc %s, thread %d: constant %s  max %s  std %s  keep %s  in %s  out %s\n",
	    acc,worker_id,comma0,comma1,comma2,comma3,comma4,comma5);

    if ((memusage = Mem_usage_report_std_heap()) != 0) {
      fprintf(stderr,"Memory leak in worker thread %ld of %ld bytes\n",worker_id,memusage);
      fflush(stdout);
      exit(9);
    }
#endif
  }

#ifdef MEMUSAGE
  Mem_usage_std_heap_add(memusage_constant);
#endif

  Except_stack_destroy();

  if (worker_stopwatch != NULL) {
    Stopwatch_free(&worker_stopwatch);
  }
  Cellpool_free(&cellpool);
  Diagpool_free(&diagpool);
  Pairpool_free(&pairpool);
  Matchpool_free(&matchpool);
  Dynprog_free(&dynprogR);
  Dynprog_free(&dynprogM);
  Dynprog_free(&dynprogL);
  Oligoindex_array_free(&oligoindices_minor);
  Oligoindex_array_free(&oligoindices_major);
  Stage2_alloc_free(&stage2_alloc);

#ifdef MEMUSAGE
  Mem_usage_set_threadname("main");
#endif

  return (void *) NULL;
}
#endif


void
check_map_iit (IIT_T map_iit, Univ_IIT_T chromosome_iit) {
  char *typestring, *lookup, *p;
  int type, destranded_len;
  bool errorp = false;

  for (type = 1; type < IIT_ntypes(map_iit); type++) {
    lookup = typestring = IIT_typestring(map_iit,type);
    if ((p = rindex(typestring,'+')) != NULL) {
      destranded_len = (p - typestring)/sizeof(char);
      lookup = (char *) MALLOC((destranded_len+1)*sizeof(char));
      strncpy(lookup,typestring,destranded_len);
      lookup[destranded_len] = '\0';

    } else if ((p = rindex(typestring,'-')) != NULL) {
      destranded_len = (p - typestring)/sizeof(char);
      lookup = (char *) MALLOC((destranded_len+1)*sizeof(char));
      strncpy(lookup,typestring,destranded_len);
      lookup[destranded_len] = '\0';
    }

    if (Univ_IIT_find_one(chromosome_iit,lookup) < 0) {
      if (p != NULL) {
	fprintf(stderr,"Warning: In %s, type %s (without the %s) does not correspond to a known chromosome in %s.\n",
		map_iitfile,typestring,p,dbversion);
      } else {
	fprintf(stderr,"Warning: In %s, type %s does not correspond to a known chromosome in %s.\n",
		map_iitfile,typestring,dbversion);
      }
      errorp = true;
    }

    if (p != NULL) {
      FREE(lookup);
    }
  }
  if (errorp == true) {
    fprintf(stderr,"Known chromosomes: ");
    Univ_IIT_dump_labels(stderr,chromosome_iit);
  }
  return;
}


void
parse_part (unsigned int *part_modulus, unsigned int *part_interval, char *string) {
  char *p = string;

  if (sscanf(p,"%u",&(*part_modulus)) < 1) {
    fprintf(stderr,"Cannot parse first integer from %s\n",string);
    exit(9);
  }

  while (*p != '\0' && isdigit(*p)) {
    p++;
  }
  while (*p != '\0' && !isdigit(*p)) {
    p++;
  }
  if (sscanf(p,"%d",&(*part_interval)) < 1) {
    fprintf(stderr,"Cannot parse first integer from %s\n",string);
    exit(9);
  }
  if ((*part_modulus) >= (*part_interval)) {
    fprintf(stderr,"In %s, batch number %u must be less than the number of batches %u\n",
	    string,*part_modulus,*part_interval);
    exit(9);
  }
  if (*part_interval == 0) {
    fprintf(stderr,"Bad batch specification %s.  Batch interval cannot be 0.\n",string);
    exit(9);
  }

  return;
}


static char *
check_valid_int (char *string) {
  char *p = string;

  if (*p == '+' || *p == '-') {
    p++;
  }

  if (!isdigit(*p)) {
    fprintf(stderr,"value %s is not a valid int\n",string);
    exit(9);
    return NULL;
  }
  while (*p != '\0' && isdigit(*p)) {
    p++;
  }

  if (*p == 'e') {
    p++;
    if (*p == '+') {
      p++;
    }
    if (!isdigit(*p)) {
      return false;
    }
    while (*p != '\0' && isdigit(*p)) {
      p++;
    }
  }

  if (*p == '\0') {
    return string;
  } else {
    fprintf(stderr,"value %s is not a valid int\n",string);
    exit(9);
    return NULL;
  }
}


static double
check_valid_float (char *string, const char *option) {
  double value;
  char *p = string;

  if (*p == '+' || *p == '-') {
    p++;
  }

  while (*p != '\0' && isdigit(*p)) {
    p++;
  }
  if (*p == '\0') {
    if ((value = atof(string)) > 1.0 || value < 0.0) {
      fprintf(stderr,"Value for option %s should be between 0.0 and 1.0\n",option);
      exit(9);
    } else {
      return value;
    }
  }

  if (*p == '.') {
    p++;
  }

  if (!isdigit(*p)) {
    fprintf(stderr,"Value %s for option %s is not a valid float\n",string,option);
    exit(9);
    return 0.0;
  }
  while (*p != '\0' && isdigit(*p)) {
    p++;
  }

  if (*p == 'e') {
    p++;
    if (*p == '+' || *p == '-') {
      p++;
    }
    if (!isdigit(*p)) {
      fprintf(stderr,"Value %s for option %s is not a valid float\n",string,option);
      exit(9);
      return 0.0;
    }
    while (*p != '\0' && isdigit(*p)) {
      p++;
    }
  }

  if (*p == '\0') {
    if ((value = atof(string)) > 1.0 || value < 0.0) {
      fprintf(stderr,"Value for option %s should be between 0.0 and 1.0\n",option);
      exit(9);
    } else {
      return value;
    }
  } else {
    fprintf(stderr,"Value %s for option %s is not a valid float\n",string,option);
    exit(9);
    return 0.0;
  }
}

static char *
check_valid_float_or_int (char *string) {
  char *p = string;

  if (*p == '+' || *p == '-') {
    p++;
  }

  while (*p != '\0' && isdigit(*p)) {
    p++;
  }
  if (*p == '\0') {
    return string;
  }

  if (*p == '.') {
    p++;
  }

  if (!isdigit(*p)) {
    fprintf(stderr,"value %s is not a valid float\n",string);
    exit(9);
    return NULL;
  }
  while (*p != '\0' && isdigit(*p)) {
    p++;
  }

  if (*p == 'e') {
    p++;
    if (*p == '+' || *p == '-') {
      p++;
    }
    if (!isdigit(*p)) {
      fprintf(stderr,"value %s is not a valid float\n",string);
      exit(9);
      return NULL;
    }
    while (*p != '\0' && isdigit(*p)) {
      p++;
    }
  }

  if (*p == '\0') {
    return string;
  } else {
    fprintf(stderr,"value %s is not a valid float\n",string);
    exit(9);
    return NULL;
  }
}


static int
parse_command_line (int argc, char *argv[], int optind) {
  int opt, c;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;
  char **argstart;

  int len;
  int user_ngap = -1;


  fprintf(stderr,"GMAP version %s called with args:",PACKAGE_VERSION);
  argstart = &(argv[-optind]);
  for (c = 1; c < argc + optind; c++) {
    fprintf(stderr," %s",argstart[c]);
  }
  fprintf(stderr,"\n");

  while ((opt = getopt_long(argc,argv,
#ifdef PMAP
			    "q:D:a:d:k:g:2B:K:w:L:x:1t:s:c:SA03468:9n:f:ZO5o:V:v:M:m:ebu:E:PQYI:i:l:",
#else
			    "q:D:d:k:g:2B:K:w:L:x:1t:s:c:p:SA03468:9n:f:ZO5o:V:v:M:m:ebu:E:PQFa:Tz:j:YI:i:l:",
#endif
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 0:
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"version")) {
	print_program_version();
	return 1;
      } else if (!strcmp(long_name,"check")) {
	check_compiler_assumptions();
	return 1;
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	return 1;

      } else if (!strcmp(long_name,"time")) {
	timingp = true;

      } else if (!strcmp(long_name,"use-shared-memory")) {
	if (!strcmp(optarg,"1")) {
	  sharedp = true;
	} else if (!strcmp(optarg,"0")) {
	  sharedp = false;
	} else {
	  fprintf(stderr,"--use-shared-memory flag must be 0 or 1\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"preload-shared-memory")) {
	preload_shared_memory_p = true;

      } else if (!strcmp(long_name,"unload-shared-memory")) {
	unload_shared_memory_p = true;

      } else if (!strcmp(long_name,"expand-offsets")) {
	fprintf(stderr,"Note: --expand-offsets flag is no longer supported.  With the latest algorithms, it doesn't improve speed much.  Ignoring this flag");

      } else if (!strcmp(long_name,"sampling")) {
	required_index1interval = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"cmdline")) {
	user_cmdline = optarg;

      } else if (!strcmp(long_name,"suboptimal-score")) {
	suboptimal_score_float = atof(check_valid_float_or_int(optarg));
	if (suboptimal_score_float > 1.0 && suboptimal_score_float != rint(suboptimal_score_float)) {
	  fprintf(stderr,"Cannot specify fractional value %f for --suboptimal-score except between 0.0 and 1.0\n",
		  suboptimal_score_float);
	  return 9;
	}

      } else if (!strcmp(long_name,"require-splicedir")) {
	require_splicedir_p = true;

      } else if (!strcmp(long_name,"splicingdir")) {
	user_splicingdir = optarg;

      } else if (!strcmp(long_name,"nosplicing")) {
	novelsplicingp = false;

      } else if (!strcmp(long_name,"no-chimeras")) {
	chimera_margin = 0;

      } else if (!strcmp(long_name,"translation-code")) {
	translation_code = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"alt-start-codons")) {
	alt_initiation_codons_p = true;

      } else if (!strcmp(long_name,"max-deletionlength")) {
	max_deletionlength = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"min-intronlength")) {
	min_intronlength = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"max-intronlength-middle")) {
	maxintronlen = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"max-intronlength-ends")) {
	maxintronlen_ends = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"split-large-introns")) {
	split_large_introns_p = true;

      } else if (!strcmp(long_name,"end-trimming-score")) {
	end_trimming_score = atoi(check_valid_int(optarg));
	if (end_trimming_score > 0) {
	  fprintf(stderr,"end-trimming-score should be 0 or negative\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"trim-end-exons")) {
	minendexon = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"allow-close-indels")) {
	if (!strcmp(optarg,"0")) {
	  /* Disallow */
	  close_indels_mode = -1;
	  extraband_single = 0;
	} else if (!strcmp(optarg,"1")) {
	  /* Always allow */
	  close_indels_mode = +1;
	  extraband_single = 3;
	} else if (!strcmp(optarg,"2")) {
	  /* Allow for high-quality alignments */
	  close_indels_mode = 0;
	  extraband_single = 3;
	} else {
	  fprintf(stderr,"allow-close-indels argument %s not recognized.  Only allow 0, 1, or 2.  Run 'gsnap --help' for more information.\n",optarg);
	  return 9;
	}
      } else if (!strcmp(long_name,"microexon-spliceprob")) {
	microexon_spliceprob = check_valid_float(optarg,long_name);
      } else if (!strcmp(long_name,"stage2-start")) {
	suboptimal_score_start = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"stage2-end")) {
	suboptimal_score_end = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"canonical-mode")) {
	if (!strcmp(optarg,"0")) {
	  canonical_mode = 0;
	} else if (!strcmp(optarg,"1")) {
	  canonical_mode = 1;
	} else if (!strcmp(optarg,"2")) {
	  canonical_mode = 2;
	} else {
	  fprintf(stderr,"Canonical level %s not recognized.\n",optarg);
	  fprintf(stderr,"0=low reward for canonical introns, 1=high reward for canonical introns (default)\n");
	  fprintf(stderr,"2=low reward for high-identity seqs, high reward otherwise\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"cross-species")) {
	cross_species_p = true;

      } else if (!strcmp(long_name,"indel-open")) {
	user_open = atoi(check_valid_int(optarg));
	if (user_open > 0) {
	  fprintf(stderr,"Expecting a negative value for --indel-open\n");
	  exit(9);
	} else if (user_open < -127) {
	  /* Limited to -127 because of SIMD dynamic programming on 8-bit values */
	  user_open = -127;
	}
	  
	user_dynprog_p = true;

      } else if (!strcmp(long_name,"indel-extend")) {
	user_extend = atoi(check_valid_int(optarg));
	if (user_extend > 0) {
	  fprintf(stderr,"Expecting a negative value for --indel-extend\n");
	  exit(9);
	} else if (user_extend < -127) {
	  /* Limited to -127 because of SIMD dynamic programming on 8-bit values */
	  user_extend = -127;
	}
	user_dynprog_p = true;

      } else if (!strcmp(long_name,"homopolymer")) {
	homopolymerp = true;

      } else if (!strcmp(long_name,"cmetdir")) {
	user_modedir = optarg;

      } else if (!strcmp(long_name,"atoidir")) {
	user_modedir = optarg;

      } else if (!strcmp(long_name,"mode")) {
	if (!strcmp(optarg,"standard")) {
	  mode = STANDARD;
	} else if (!strcmp(optarg,"cmet-stranded")) {
	  mode = CMET_STRANDED;
	} else if (!strcmp(optarg,"cmet-nonstranded")) {
	  mode = CMET_NONSTRANDED;
	  fprintf(stderr,"Non-stranded mode not yet working properly\n");
	  exit(9);
	} else if (!strcmp(optarg,"atoi-stranded")) {
	  mode = ATOI_STRANDED;
	} else if (!strcmp(optarg,"atoi-nonstranded")) {
	  mode = ATOI_NONSTRANDED;
	  fprintf(stderr,"Non-stranded mode not yet working properly\n");
	  exit(9);
	} else if (!strcmp(optarg,"ttoc-stranded")) {
	  mode = TTOC_STRANDED;
	} else if (!strcmp(optarg,"ttoc-nonstranded")) {
	  mode = TTOC_NONSTRANDED;
	  fprintf(stderr,"Non-stranded mode not yet working properly\n");
	  exit(9);
	} else {
	  fprintf(stderr,"--mode must be standard, cmet-stranded, cmet-nonstranded, atoi-stranded, atoi-nonstranded, ttoc-stranded, or ttoc-nonstranded\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"min-trimmed-coverage")) {
	min_trimmed_coverage = check_valid_float(optarg,long_name);
      } else if (!strcmp(long_name,"min-identity")) {
	min_identity = check_valid_float(optarg,long_name);

      } else if (!strcmp(long_name,"read-files-command")) {
	read_files_command = optarg;

      } else if (!strcmp(long_name,"input-buffer-size")) {
	input_buffer_size = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"output-buffer-size")) {
	output_buffer_size = atoi(check_valid_int(optarg));
      } else if (!strcmp(long_name,"print-comment")) {
	print_comment_p = true;

      } else if (!strcmp(long_name,"strand")) {
	if (!strcmp(optarg,"plus")) {
	  cdna_direction_try = +1;
	} else if (!strcmp(optarg,"minus")) {
	  cdna_direction_try = -1;
	} else if (!strcmp(optarg,"both")) {
	  cdna_direction_try = 0;
	} else {
	  fprintf(stderr,"strand %s not recognized.  Must be plus, minus, or both (default)\n",optarg);
	  return 9;
	}

      } else if (!strcmp(long_name,"nolengths")) {
	nointronlenp = true;
      } else if (!strcmp(long_name,"nomargin")) {
	print_margin_p = false;
      } else if (!strcmp(long_name,"failsonly")) {
	if (nofailsp == true) {
	  fprintf(stderr,"Cannot specify both --nofails and --failsonly\n");
	  return 9;
	} else {
	  failsonlyp = true;
	}
      } else if (!strcmp(long_name,"failed-input")) {
	failedinput_root = optarg;
#if 0
      } else if (!strcmp(long_name,"quiet-if-excessive")) {
	quiet_if_excessive_p = true;
#endif
      } else if (!strcmp(long_name,"nofails")) {
	if (failsonlyp == true) {
	  fprintf(stderr,"Cannot specify both --nofails and --failsonly\n");
	  return 9;
	} else {
	  nofailsp = true;
	}
      } else if (!strcmp(long_name,"split-output")) {
	split_output_root = optarg;
      } else if (!strcmp(long_name,"append-output")) {
	appendp = true;

      } else if (!strcmp(long_name,"gff3-add-separators")) {
	if (!strcmp(optarg,"1")) {
	  gff3_separators_p = true;
	} else if (!strcmp(optarg,"0")) {
	  gff3_separators_p = false;
	} else {
	  fprintf(stderr,"--gff3-add-separators flag must be 0 or 1\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"gff3-swap-phase")) {
	if (!strcmp(optarg,"1")) {
	  gff3_phase_swap_p = true;
	} else if (!strcmp(optarg,"0")) {
	  gff3_phase_swap_p = false;
	} else {
	  fprintf(stderr,"--gff3-swap-phase flag must be 0 or 1\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"gff3-fasta-annotation")) {
	if (!strcmp(optarg,"0")) {
	  gff3_fasta_annotation_type = NO_ANNOTATION;
	} else if (!strcmp(optarg,"1")) {
	  gff3_fasta_annotation_type = INSERT_ANNOTATION;
	} else if (!strcmp(optarg,"2")) {
	  gff3_fasta_annotation_type = KEYVALUE_ANNOTATION;
	} else {
	  fprintf(stderr,"--gff3-fasta-annotation flag must be 0, 1, or 2\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"gff3-cds")) {
	if (!strcmp(optarg,"cdna")) {
	  cdstype = CDS_CDNA;
	} else if (!strcmp(optarg,"genomic")) {
	  cdstype = CDS_GENOMIC;
	} else {
	  fprintf(stderr,"--gff3-cds flag must be cdna or genomic\n");
	  return 9;
	}

#ifndef PMAP
      } else if (!strcmp(long_name,"no-sam-headers")) {
	sam_headers_p = false;
      } else if (!strcmp(long_name,"sam-use-0M")) {
	sam_insert_0M_p = true;
      } else if (!strcmp(long_name,"sam-extended-cigar")) {
	sam_cigar_extended_p = true;
      } else if (!strcmp(long_name,"sam-flipped")) {
	sam_flippedp = true;
      } else if (!strcmp(long_name,"quality-protocol")) {
	if (user_quality_shift == true) {
	  fprintf(stderr,"Cannot specify both -j (--quality-print-shift) and --quality-protocol\n");
	  return 9;
	} else if (!strcmp(optarg,"illumina")) {
	  quality_shift = -31;
	  user_quality_shift = true;
	} else if (!strcmp(optarg,"sanger")) {
	  quality_shift = 0;
	  user_quality_shift = true;
	} else {
	  fprintf(stderr,"The only values allowed for --quality-protocol are illumina or sanger\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"force-xs-dir")) {
	force_xs_direction_p = true;

      } else if (!strcmp(long_name,"md-lowercase-snp")) {
	md_lowercase_variant_p = true;

      } else if (!strcmp(long_name,"action-if-cigar-error")) {
	if (!strcmp(optarg,"ignore")) {
	  cigar_action = CIGAR_ACTION_IGNORE;
	} else if (!strcmp(optarg,"warning")) {
	  cigar_action = CIGAR_ACTION_WARNING;
	} else if (!strcmp(optarg,"noprint")) {
	  cigar_action = CIGAR_ACTION_NOPRINT;
	} else if (!strcmp(optarg,"abort")) {
	  cigar_action = CIGAR_ACTION_ABORT;
	} else {
	  fprintf(stderr,"The only values allowed for --action-if-cigar-error are ignore, warning, noprint, abort\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"read-group-id")) {
	sam_read_group_id = optarg;
      } else if (!strcmp(long_name,"read-group-name")) {
	sam_read_group_name = optarg;
      } else if (!strcmp(long_name,"read-group-library")) {
	sam_read_group_library = optarg;
      } else if (!strcmp(long_name,"read-group-platform")) {
	sam_read_group_platform = optarg;
#endif
      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'gmap --help'",long_name);
	return 9;
      }
      break;

    case 'q': parse_part(&part_modulus,&part_interval,optarg); break;
    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;
#ifdef PMAP
    case 'a':
      if ((required_alphabet = Alphabet_find(optarg)) == AA0) {
	return 9;
      }
      break;
    case 'k': required_index1part = atoi(check_valid_int(optarg)); break;
#else
    case 'k':
      required_index1part = atoi(check_valid_int(optarg));
      if (required_index1part > MAXIMUM_KMER) {
	fprintf(stderr,"The value for k-mer size must be %d or less\n",MAXIMUM_KMER);
	return 9;
      }
      break;
#endif
#if 0
    case 'G': uncompressedp = true; break;
#endif
    case 'g': user_genomicsegs = optarg; break;
    case '1': user_selfalign_p = true; break;
    case '2': user_pairalign_p = true; break;

    case 'B': 
      if (!strcmp(optarg,"5")) {
	fprintf(stderr,"Note: Batch mode 5 is now the same as batch mode 4.\n");
	offsetsstrm_access = USE_ALLOCATE; /* Doesn't matter */
	positions_access = USE_ALLOCATE;
	locoffsetsstrm_access = USE_ALLOCATE; /* Doesn't matter */
	locpositions_access = USE_ALLOCATE;

	genome_access = USE_ALLOCATE;

      } else if (!strcmp(optarg,"4")) {
	offsetsstrm_access = USE_ALLOCATE;
	positions_access = USE_ALLOCATE;
	locoffsetsstrm_access = USE_ALLOCATE;
	locpositions_access = USE_ALLOCATE;

	genome_access = USE_ALLOCATE;
#ifdef HAVE_MMAP

      } else if (!strcmp(optarg,"3")) {
	offsetsstrm_access = USE_ALLOCATE;
	positions_access = USE_ALLOCATE;
	locoffsetsstrm_access = USE_ALLOCATE;
	locpositions_access = USE_ALLOCATE;

	genome_access = USE_MMAP_PRELOAD; /* was batch_genome_p = true */

      } else if (!strcmp(optarg,"2")) {
	offsetsstrm_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	positions_access = USE_MMAP_PRELOAD; /* was batch_positions_p = true */
	locoffsetsstrm_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	locpositions_access = USE_MMAP_PRELOAD; /* was batch_positions_p = true */

	genome_access = USE_MMAP_PRELOAD; /* was batch_genome_p = true */

      } else if (!strcmp(optarg,"1")) {
	offsetsstrm_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	positions_access = USE_MMAP_PRELOAD; /* was batch_positions_p = true */
	locoffsetsstrm_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	locpositions_access = USE_MMAP_PRELOAD; /* was batch_positions_p = true */

	genome_access = USE_MMAP_ONLY; /* was batch_genome_p = false */

      } else if (!strcmp(optarg,"0")) {
	offsetsstrm_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	positions_access = USE_MMAP_ONLY; /* was batch_positions_p = false */
	locoffsetsstrm_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	locpositions_access = USE_MMAP_ONLY; /* was batch_positions_p = false */

	genome_access = USE_MMAP_ONLY; /* was batch_genome_p = false */
#endif

      } else {
#ifdef HAVE_MMAP
	fprintf(stderr,"Batch mode %s not recognized.  Only allow 0-5.  Run 'gmap --help' for more information.\n",optarg);
#else
	fprintf(stderr,"Batch mode %s not recognized.  Only allow 4-5, since mmap is disabled.  Run 'gmap --help' for more information.\n",optarg);
#endif
	return 9;
      }
      break;

    case 'K': maxintronlen = maxintronlen_ends = atoi(check_valid_int(optarg)); break;

    case 'w': shortsplicedist = strtoul(check_valid_int(optarg),NULL,10); break;

    case 'L': maxtotallen_bound = atoi(check_valid_int(optarg)); break;
    case 'x':
#ifdef PMAP
      chimera_margin = atoi(check_valid_int(optarg))/3; 
#else
      chimera_margin = atoi(check_valid_int(optarg)); 
#endif
      if (chimera_margin <= 0) {
	/* Disable finding of chimeras */
#if 0
      } else if (chimera_margin < CHIMERA_SLOP) {
	/* Not sure why chimera_margin should be tied to CHIMERA_SLOP */
	chimera_margin = CHIMERA_SLOP;
#endif
      }
      break;
      /* case 'w': referencefile = optarg; break; */

#ifdef HAVE_PTHREAD
    case 't': nworkers = atoi(check_valid_int(optarg)); break;
#else
    case 't': fprintf(stderr,"This version of GMAP has pthreads disabled, so ignoring the value of %s for -t\n",optarg); break;
#endif

    case 's': splicing_file = optarg; knownsplicingp = true; break;
    case 'c': user_chrsubsetname = optarg; break;

#ifndef PMAP
    case 'p': switch (atoi(check_valid_int(optarg))) {
      case 0: prune_poor_p = false, prune_repetitive_p = false; break;
      case 1: prune_poor_p = true; prune_repetitive_p = false; break;
      case 2: prune_poor_p = false; prune_repetitive_p = true; break;
      case 3: prune_poor_p = true; prune_repetitive_p = true; break;
      default: fprintf(stderr,"Prune level %s not recognized.\n",optarg);
	fprintf(stderr,"0=no pruning, 1=poor seqs, 2=repetitive seqs, 3=both poor and repetitive seqs (default)\n");
	return 9;
      }
      break;
#endif

    case 'S': printtype = SUMMARY; break;
    case 'A': printtype = ALIGNMENT; break;
    case '0': exception_raise_p = false; break; /* Allows signals to pass through */
    case '3': printtype = CONTINUOUS; break;
    case '4': printtype = CONTINUOUS_BY_EXON; break;
    case '6': debug_graphic_p = true; break;
    case '8':
      if (!strcmp(optarg,"stage1")) {
	stage1debug = true;
      } else if (!strcmp(optarg,"diag")) {
	diag_debug = true;
      } else if (!strcmp(optarg,"stage2")) {
	stage3debug = POST_STAGE2;
      } else if (!strcmp(optarg,"singles")) {
	stage3debug = POST_SINGLES;
      } else if (!strcmp(optarg,"introns")) {
	stage3debug = POST_INTRONS;
      } else if (!strcmp(optarg,"hmm")) {
	stage3debug = POST_HMM;
      } else if (!strcmp(optarg,"smoothing")) {
	stage3debug = POST_SMOOTHING;
      } else if (!strcmp(optarg,"dualintrons")) {
	stage3debug = POST_DUAL_INTRONS;
      } else if (!strcmp(optarg,"cycles")) {
	stage3debug = POST_CYCLES;
      } else if (!strcmp(optarg,"dualbreaks")) {
	stage3debug = POST_DUAL_BREAKS;
      } else if (!strcmp(optarg,"middle")) {
	stage3debug = POST_MIDDLE;
      } else if (!strcmp(optarg,"ends")) {
	stage3debug = POST_ENDS;
      } else if (!strcmp(optarg,"canonical")) {
	stage3debug = POST_CANONICAL;
      } else if (!strcmp(optarg,"trim")) {
	stage3debug = POST_CANONICAL;
      } else if (!strcmp(optarg,"changepoint")) {
	stage3debug = POST_CHANGEPOINT;
      } else if (!strcmp(optarg,"distalmedial")) {
	stage3debug = POST_DISTAL_MEDIAL;
      } else {
	fprintf(stderr,"Allowed arguments for -8 flag are stage2, smoothing, singles, introns, hmm, dualbreaks, cycles, canonical, changepoint, distalmedial\n");
	return 9;
      }
      break;
    case '9': checkp = true; break;
    case 'n':
      maxpaths_report = atoi(check_valid_int(optarg));
      if (maxpaths_report == 1) {
	fprintf(stderr,"Note: -n 1 will not report chimeric alignments.  If you want a single alignment plus chimeras, use -n 0 instead.\n");
      }
      break;
    case 'f':
      if (!strcmp(optarg,"1") || !strcmp(optarg,"psl_nt")) {
	printtype = PSL_NT;
#ifdef PMAP
      } else if (!strcmp(optarg,"0") || !strcmp(optarg,"psl_pro")) {
	printtype = PSL_PRO;
#else
      } else if (!strcmp(optarg,"psl")) {
	printtype = PSL_NT;
      } else if (!strcmp(optarg,"6") || !strcmp(optarg,"splicesites")) {
	printtype = SPLICESITES;
      } else if (!strcmp(optarg,"introns")) {
	printtype = INTRONS;
      } else if (!strcmp(optarg,"mask_introns")) {
	printtype = MASK_INTRONS;
      } else if (!strcmp(optarg,"mask_utr_introns")) {
	printtype = MASK_UTR_INTRONS;
      } else if (!strcmp(optarg,"samse")) {
	printtype = SAM;
	sam_paired_p = false;
      } else if (!strcmp(optarg,"sampe")) {
	printtype = SAM;
	sam_paired_p = true;
      } else if (!strcmp(optarg,"bedpe")) {
	printtype = BEDPE;
#endif
      } else if (!strcmp(optarg,"2") || !strcmp(optarg,"gff3_gene")) {
	printtype = GFF3_GENE;
      } else if (!strcmp(optarg,"3") || !strcmp(optarg,"gff3_match_cdna")) {
	printtype = GFF3_MATCH_CDNA;
      } else if (!strcmp(optarg,"4") || !strcmp(optarg,"gff3_match_est")) {
	printtype = GFF3_MATCH_EST;
      } else if (!strcmp(optarg,"7") || !strcmp(optarg,"map_exons")) {
	printtype = MAP_EXONS;
      } else if (!strcmp(optarg,"8") || !strcmp(optarg,"map_ranges")) {
	printtype = MAP_RANGES;
      } else if (!strcmp(optarg,"9") || !strcmp(optarg,"coords")) {
	printtype = COORDS;
      } else {
	fprintf(stderr,"Output format \"%s\" not recognized.  Allowed formats are:\n",optarg);
	fprintf(stderr,"  psl_nt (1)\n");
#ifdef PMAP
	fprintf(stderr,"  psl_pro (0)\n");
#else
	fprintf(stderr,"  psl\n");
	fprintf(stderr,"  splicesites (6)\n");
	fprintf(stderr,"  introns\n");
	fprintf(stderr,"  mask_introns\n");
	fprintf(stderr,"  mask_utr_introns\n");
	fprintf(stderr,"  samse\n");
	fprintf(stderr,"  sampe\n");
	fprintf(stderr,"  bedpe\n");
#endif
	fprintf(stderr,"  gff3_gene (2)\n");
	fprintf(stderr,"  gff3_match_cdna (3)\n");
	fprintf(stderr,"  gff3_match_est (4)\n");
	fprintf(stderr,"  map_exons (7)\n");
	fprintf(stderr,"  map_ranges (8)\n");
	fprintf(stderr,"  coords (9)\n");
	return 9;
      }
      break;
    case 'O': orderedp = true; break;
    case '5': checksump = true; break;
    case 'o': chimera_overlap = atoi(check_valid_int(optarg)); break;

    case 'V': user_snpsdir = optarg; break;
    case 'v': snps_root = optarg; break;

    case 'M': user_mapdir = optarg; break;
    case 'm': 
      map_iitfile = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(map_iitfile,optarg);
      if ((len = strlen(map_iitfile)) > 4 && strcmp(&(map_iitfile[len-4]),".iit") == 0) {
	map_iitfile[len-4] = '\0';
      }
      break;

    case 'e': map_exons_p = true; break;
    case 'b': map_bothstrands_p = true; break;
    case 'u': nflanking = atoi(check_valid_int(optarg)); break;

    case 'E': 
      if (!strcmp(optarg,"cdna")) {
	printtype = EXONS_CDNA;
      } else if (!strcmp(optarg,"genomic")) {
	printtype = EXONS_GENOMIC;
      } else if (!strcmp(optarg,"cdna+introns")) {
	printtype = EXONS_CDNA_WINTRONS;
      } else if (!strcmp(optarg,"genomic+introns")) {
	printtype = EXONS_GENOMIC_WINTRONS;
      } else {
	fprintf(stderr,"Argument to -E flag must be either \"cdna\" or \"genomic\"\n");
	return 9;
      }
      break;

#ifdef PMAP
    case 'P': printtype = PROTEIN_GENOMIC; break;
    case 'Q': printtype = CDNA; break; 
#else
    case 'P': printtype = CDNA; break;
    case 'Q': printtype = PROTEIN_GENOMIC; break;
    case 'F': fulllengthp = true; break;
    case 'a': cds_startpos = atoi(check_valid_int(optarg)); break;
    case 'T': truncatep = true; fulllengthp = true; break;
    case 'z':
      if (!strcmp(optarg,"sense_force")) {
	sense_try = +1;
	sense_filter = 0;
      } else if (!strcmp(optarg,"antisense_force")) {
	sense_try = -1;
	sense_filter = 0;
      } else if (!strcmp(optarg,"sense_filter")) {
	sense_try = 0;
	sense_filter = +1;
      } else if (!strcmp(optarg,"antisense_filter")) {
	sense_try = 0;
	sense_filter = -1;
      } else if (!strcmp(optarg,"auto")) {
	sense_try = 0;
	sense_filter = 0;
      } else {
	fprintf(stderr,"direction %s not recognized.  Must be sense_force, antisense_force, sense_filter, antisense_filter, or auto\n",optarg);
	return 9;
      }
      break;

    case 'j':
      if (user_quality_shift == true) {
	fprintf(stderr,"Cannot specify both -j (--quality-print-shift) and --quality-protocol\n");
	return 9;
      } else {
	quality_shift = atoi(check_valid_int(optarg));
	user_quality_shift = true;
      }
      break;

#endif
    case 'Y': strictp = false; break;
    case 'I': invertmode = atoi(check_valid_int(optarg)); break;
    case 'i': user_ngap = atoi(check_valid_int(optarg)); break;
    case 'l': wraplength = atoi(check_valid_int(optarg)); break;

    case '?': fprintf(stderr,"For usage, run 'gmap --help'\n"); return 9;
    default: return 9;
    }
  }

  if (printtype == SPLICESITES || printtype == INTRONS) {
    if (maxpaths_report > 1 || (sense_try != +1 && sense_filter != +1)) {
      fprintf(stderr,"For splicesites or introns output, you should probably add flags '-n 1' and either '-z sense_force' or '-z sense_filter'.\n");
    }
  }
  
  if (user_ngap >= 0) {
    ngap = user_ngap;
  } else if (printtype == EXONS_CDNA || printtype == EXONS_GENOMIC) {
    /* If user didn't specify, then set to zero */
    ngap = 0;
  } else if (printtype == EXONS_CDNA_WINTRONS || printtype == EXONS_GENOMIC_WINTRONS) {
    /* If user didn't specify, then set to infinity */
    ngap = 10000000;		/* largest possible intron */
  };

  if (maxintronlen > maxtotallen_bound) {
    maxintronlen = maxtotallen_bound;
  }

#ifdef HAVE_PTHREAD
#ifdef USE_DIAGPOOL
  if (diag_debug == true && nworkers > 0) {
    fprintf(stderr,"For diag output, must specify 0 threads\n");
    exit(9);
  }
#endif
#endif

  if (user_cmdline != NULL) {
    part_modulus = 0;
    part_interval = 1;
    input_buffer_size = 0;
    nchromosomes = 1;
    dbroot = (char *) NULL;
  } else if (user_selfalign_p == true) {
    nchromosomes = 1;
    dbroot = (char *) NULL;
  } else if (user_pairalign_p == true) {
    nchromosomes = 1;
    dbroot = (char *) NULL;
  } else if (user_genomicsegs != NULL) {
    /* Ignore -D and -d flags */
    nchromosomes = 1;
    dbroot = (char *) NULL;
  } else if (dbroot == NULL) {
    fprintf(stderr,"Need to specify the -d, -g, -1, -2, or --cmdline flag\n");
    print_program_usage();
    return 9;
  } else if (!strcmp(dbroot,"?")) {
    Datadir_avail_gmap_databases(stdout,user_genomedir);
    return 1;
  }

#ifndef PMAP
  if (printtype == SAM) {
    if (sam_read_group_id == NULL && sam_read_group_name != NULL) {
      sam_read_group_id = sam_read_group_name;
    } else if (sam_read_group_id != NULL && sam_read_group_name == NULL) {
      sam_read_group_name = sam_read_group_id;
    }
  }
#endif

  return 0;
}


static Inbuffer_T
open_input_stream (int *nread, List_T user_genomes, int argc, char **argv) {
  Inbuffer_T inbuffer;
  int nextchar = '\0';
  FILE *input = NULL;
  char **files;
  int nfiles;

  Genome_T genome, genomealt;
  Sequence_T genomeseq;
  char *p;

  if (user_cmdline != NULL) {
    /* Query (and genome) sequences come from the command line */
    p = user_cmdline;
    while (*p != '\0' && *p != ',') {
      p++;
    }
    if (*p == '\0') {
      fprintf(stderr,"--cmdline requires two strings separated by a comma");
      exit(9);
    } else {
      genomeseq = Sequence_genomic_new(user_cmdline,(int) (p - user_cmdline),/*copyp*/true);
      p++;
    }
    genome = genomealt = Genome_from_sequence(genomeseq);
    Sequence_free(&genomeseq);

    /* Read in a single queryseq from command line */
    inbuffer = Inbuffer_cmdline(/*queryseq*/p,/*querylength*/strlen(p),genome,genomealt);
    *nread = 1;

  } else if (argc == 0) {
    /* Query sequence comes from stdin */
    /* Open input stream and peek at first char */
    fprintf(stderr,"Reading from stdin\n");
    input = stdin;
    files = (char **) NULL;
    nfiles = 0;

    /* Read in first batch of sequences */
    inbuffer = Inbuffer_new(nextchar,input,read_files_command,files,nfiles,input_buffer_size,
			    user_pairalign_p,user_genomes,part_modulus,part_interval);
    *nread = Inbuffer_fill_init(inbuffer);

  } else {
    /* Query sequence comes from a file or files */
    input = NULL;
    files = argv;
    nfiles = argc;

    /* Read in first batch of sequences */
    inbuffer = Inbuffer_new(nextchar,input,read_files_command,files,nfiles,input_buffer_size,
			    user_pairalign_p,user_genomes,part_modulus,part_interval);
    *nread = Inbuffer_fill_init(inbuffer);
  }

  return inbuffer;
}


int
main (int argc, char *argv[]) {
  int cmdline_status;
  FILE *input;
  Sequence_T genomeseq;
  List_T user_genomes = NULL, p;
  Genome_T genome;
  int nextchar = '\0';

  char *genomesubdir = NULL, *snpsdir = NULL, *modedir = NULL, *mapdir = NULL, *iitfile = NULL, *fileroot = NULL;
  char *idx_filesuffix1, *idx_filesuffix2;
  int divno;
  Univinterval_T interval;

  int nread;
  double runtime;

  Splicestringpool_T splicestringpool;

#ifdef HAVE_PTHREAD
  int ret, i;
  pthread_attr_t thread_attr_join;
#ifdef WORKER_DETACH
  pthread_attr_t thread_attr_detach;
#endif
#endif

#ifdef HAVE_SIGACTION
  struct sigaction signal_action;
#endif

  extern int optind;

#ifdef MEMUSAGE
  Mem_usage_init();
  Mem_usage_set_threadname("main");
#endif


  cmdline_status = parse_command_line(argc,argv,optind);
  argc -= optind;
  argv += optind;
      
  if (cmdline_status == 0) {
    /* okay to continue */
  } else if (cmdline_status == 1) {
    /* only information needed */
    exit(0);
  } else {
    exit(cmdline_status);
  }

  check_compiler_assumptions();

  if (exception_raise_p == false) {
    fprintf(stderr,"Allowing signals and exceptions to pass through.  If using shared memory, need to remove segments manually.\n");
    Except_inactivate();
  } else {
#ifdef HAVE_SIGACTION
    signal_action.sa_handler = signal_handler;
    signal_action.sa_flags = 0;
    sigfillset(&signal_action.sa_mask); /* After first signal, block all other signals */

    /* Note: SIGKILL and SIGSTOP cannot be caught */

    sigaction(SIGFPE,&signal_action,NULL);
    sigaction(SIGSEGV,&signal_action,NULL);
    sigaction(SIGTRAP,&signal_action,NULL);
    sigaction(SIGUSR1,&signal_action,NULL);
    sigaction(SIGABRT,&signal_action,NULL); /* abnormal termination (abort) */
    sigaction(SIGBUS,&signal_action,NULL);  /* bus error */
    sigaction(SIGFPE,&signal_action,NULL);  /* arithmetic exception */
    sigaction(SIGHUP,&signal_action,NULL);  /* hangup */
    sigaction(SIGILL,&signal_action,NULL);  /* illegal hardware instruction */
    sigaction(SIGINT,&signal_action,NULL);  /* terminal interruption (control-C) */
    sigaction(SIGPIPE,&signal_action,NULL);  /* write to pipe with no readers */
    sigaction(SIGQUIT,&signal_action,NULL);  /* terminal quit (control-backslash) */
    sigaction(SIGSEGV,&signal_action,NULL);  /* invalid memory reference */
    sigaction(SIGSYS,&signal_action,NULL);  /* invalid system call */
    sigaction(SIGTERM,&signal_action,NULL);  /* Unix kill command */
    sigaction(SIGTRAP,&signal_action,NULL);  /* hardware fault */
    sigaction(SIGXCPU,&signal_action,NULL);  /* CPU limit exceeded */
    sigaction(SIGXFSZ,&signal_action,NULL);  /* file size limit exceeded */
#endif
  }

  /* Check query sequences so we know if we can memory map the genome index */
  if (user_genomicsegs != NULL) {
    /* However, with user genomicsegs, we will not be memory mapping
       the genome index, and we need to provide these genomes to the
       inbuffer */
    if ((input = FOPEN_READ_TEXT(user_genomicsegs)) == NULL) {
      fprintf(stderr,"Can't open file %s\n",user_genomicsegs);
      exit(9);
    } else {
      nextchar = '\0';
      while ((genomeseq = Sequence_read_unlimited(&nextchar,input)) != NULL) {
	genome = Genome_from_sequence(genomeseq);
	user_genomes = List_push(user_genomes,(void *) genome);
	Sequence_free(&genomeseq);
      }
      fclose(input);
      user_genomes = List_reverse(user_genomes);
    }
  }

  inbuffer = open_input_stream(&nread,user_genomes,argc,argv);

  if (nread > 1) {
    multiple_sequences_p = true;

  } else {
    /* fprintf(stderr,"Note: only 1 sequence detected.  Ignoring batch (-B) command\n"); */
    multiple_sequences_p = false;
    expand_offsets_p = false;
#ifdef HAVE_MMAP
    offsetsstrm_access = USE_MMAP_ONLY;
    positions_access = USE_MMAP_ONLY;
    genome_access = USE_MMAP_ONLY;
#else
    offsetsstrm_access = USE_ALLOCATE;
    positions_access = USE_ALLOCATE;
    genome_access = USE_ALLOCATE;
#endif
  }

  if (dbroot == NULL) {
    dbversion = (char *) NULL;

    any_circular_p = false;
    circularp = (bool *) MALLOC(1*sizeof(bool));
    circularp[0] = false;
    circular_typeint = -1;

    altlocp = (bool *) MALLOC(sizeof(bool));
    altlocp[0] = false;
    alias_starts = (Univcoord_T *) CALLOC(1,sizeof(Univcoord_T));
    alias_ends = (Univcoord_T *) CALLOC(1,sizeof(Univcoord_T));

  } else {
    /* Prepare genomic data */
    genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);

    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    if ((chromosome_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      fprintf(stderr,"IIT file %s is not valid\n",iitfile);
      exit(9);
#ifdef LARGE_GENOMES
    } else if (Univ_IIT_coord_values_8p(chromosome_iit) == false) {
      fprintf(stderr,"This program gmapl is designed for large genomes.\n");
      fprintf(stderr,"For small genomes of less than 2^32 (4 billion) bp, please run gmap instead.\n");
      exit(9);
#endif
    } else {
      FREE(iitfile);
      nchromosomes = Univ_IIT_total_nintervals(chromosome_iit);
      circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
      circularp = Univ_IIT_circularp(&any_circular_p,chromosome_iit);

      iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				strlen(fileroot)+strlen(".altscaffold.iit")+1,sizeof(char));
      sprintf(iitfile,"%s/%s.altscaffold.iit",genomesubdir,fileroot);
      if ((altscaffold_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
	/* fprintf(stderr,"No altscaffold file found\n"); */
	altlocp = (bool *) CALLOC(nchromosomes+1,sizeof(bool));
	alias_starts = (Univcoord_T *) CALLOC(nchromosomes+1,sizeof(Univcoord_T));
	alias_ends = (Univcoord_T *) CALLOC(nchromosomes+1,sizeof(Univcoord_T));

      } else {
	fprintf(stderr,"Found altscaffold file found\n");
	altlocp = Univ_IIT_altlocp(&alias_starts,&alias_ends,chromosome_iit,altscaffold_iit);
	Univ_IIT_free(&altscaffold_iit);
      }
      FREE(iitfile);
    }

  }

  if (map_iitfile == NULL) {
    /* Skip */
  } else if (!strcmp(map_iitfile,"?")) {
    Datadir_avail_maps(stdout,user_mapdir,genomesubdir,fileroot);
    exit(0);
  } else {
    mapdir = Datadir_find_mapdir(user_mapdir,genomesubdir,fileroot);
    iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
			      strlen(map_iitfile)+strlen(".iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.iit",mapdir,map_iitfile);
    if ((map_iit = IIT_read(iitfile,/*name*/map_iitfile,/*readonlyp*/true,/*divread*/READ_ALL,
			    /*divstring*/NULL,/*add_iit_p*/true)) == NULL) {
      fprintf(stderr,"Map file %s.iit not found in %s.  Available files:\n",map_iitfile,mapdir);
      Datadir_list_directory(stderr,mapdir);
      fprintf(stderr,"Either install file %s.iit or specify a directory for the IIT file\n",iitfile);
      fprintf(stderr,"using the -M flag.\n");
      exit(9);
    } else {
      map_divint_crosstable = Univ_IIT_divint_crosstable(chromosome_iit,map_iit);
    }

    check_map_iit(map_iit,chromosome_iit);

    FREE(iitfile);
    FREE(mapdir);
    FREE(map_iitfile);
  }

  if (splicing_file != NULL) {
    if (user_splicingdir == NULL) {
      if ((splicing_iit = IIT_read(splicing_file,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				   /*divstring*/NULL,/*add_iit_p*/false)) != NULL) {
	fprintf(stderr,"Reading splicing file %s locally...",splicing_file);
      } else {
	iitfile = (char *) CALLOC(strlen(user_splicingdir)+strlen("/")+strlen(splicing_file)+1,sizeof(char));
	sprintf(iitfile,"%s/%s",user_splicingdir,splicing_file);
	if ((splicing_iit = IIT_read(splicing_file,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				     /*divstring*/NULL,/*add_iit_p*/false)) != NULL) {
	  fprintf(stderr,"Reading splicing file %s locally...",splicing_file);
	  FREE(iitfile);
	}
      }
    }

    if (splicing_iit == NULL) {
      mapdir = Datadir_find_mapdir(/*user_mapdir*/NULL,genomesubdir,fileroot);
      iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
				strlen(splicing_file)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",mapdir,splicing_file);
      if ((splicing_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				   /*divstring*/NULL,/*add_iit_p*/true)) != NULL) {
	fprintf(stderr,"Reading splicing file %s...",iitfile);
	FREE(iitfile);
	FREE(mapdir);
      } else {
	fprintf(stderr,"Splicing file %s.iit not found locally or in %s.  Available files:\n",splicing_file,mapdir);
	Datadir_list_directory(stderr,mapdir);
	fprintf(stderr,"Either install file %s or specify a full directory path\n",splicing_file);
	exit(9);
      }
    }
  }

  /* Complement_init(); */
  Dynprog_init(mode);
#ifdef PMAP
  Backtranslation_init();
#endif


  if (user_cmdline != NULL || user_selfalign_p == true || user_pairalign_p == true ||
      user_genomes != NULL) {
    /* Either genome already read from command line, or no genome
       needed, or we will read in each individual genome */

    global_genome = (Genome_T) NULL;
    global_genomealt = (Genome_T) NULL;
    genomelength = 0;

  } else {
    if (snps_root == NULL) {
      global_genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,
				 chromosome_iit,genome_access,sharedp,/*revcompp*/false);
      global_genomealt = global_genome;
      
    } else {
      /* Map against genome with SNPs */
      if (user_snpsdir == NULL) {
	snpsdir = genomesubdir;
      } else {
	snpsdir = user_snpsdir;
      }
      
      global_genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,
				 chromosome_iit,genome_access,sharedp,/*revcompp*/false);
      global_genomealt = Genome_new(snpsdir,fileroot,snps_root,
				    chromosome_iit,genome_access,sharedp,/*revcompp*/false);
    }
    genomelength = Genome_genomelength(global_genome);

    if (user_modedir != NULL) {
      modedir = user_modedir;
    } else {
      modedir = genomesubdir;
    }

    if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
      idx_filesuffix1 = "metct";
      idx_filesuffix2 = "metga";
    } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
      idx_filesuffix1 = "a2iag";
      idx_filesuffix2 = "a2itc";
    } else if (mode == TTOC_STRANDED || mode == TTOC_NONSTRANDED) {
      idx_filesuffix1 = "a2itc";
      idx_filesuffix2 = "a2iag";
    } else {
      idx_filesuffix1 = IDX_FILESUFFIX; /* "ref" */
      idx_filesuffix2 = (char *) NULL;
    }

    if ((indexdb_fwd = Indexdb_new_genome(&index1part,&index1interval,
					  /*genomesubdir*/modedir,snpsdir,
					  fileroot,idx_filesuffix1,snps_root,
					  required_index1part,required_index1interval,
					  offsetsstrm_access,positions_access,
					  sharedp,multiple_sequences_p,/*preload_shared_memory_p*/false,
					  /*unload_shared_memory_p*/false)) == NULL) {

      if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
	fprintf(stderr,"Cannot find %s index file.  Need to run cmetindex first\n",idx_filesuffix1);
      } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED ||
		 mode == TTOC_STRANDED || mode == TTOC_NONSTRANDED) {
	fprintf(stderr,"Cannot find %s index file.  Need to run atoiindex first\n",idx_filesuffix1);
      } else {
	fprintf(stderr,"Cannot find %s index file\n",idx_filesuffix1);
      }
      exit(9);
    }

    if (idx_filesuffix2 == NULL) {
      indexdb_rev = indexdb_fwd;
    } else if ((indexdb_rev = Indexdb_new_genome(&index1part,&index1interval,
						 /*genomesubdir*/modedir,snpsdir,
						 fileroot,idx_filesuffix2,snps_root,
						 required_index1part,required_index1interval,
						 offsetsstrm_access,positions_access,
						 sharedp,multiple_sequences_p,/*preload_shared_memory_p*/false,
						 /*unload_shared_memory_p*/false)) == NULL) {
      if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
	fprintf(stderr,"Cannot find %s index file.  Need to run cmetindex first\n",idx_filesuffix2);
      } else {
	fprintf(stderr,"Cannot find %s index file.  Need to run atoiindex first\n",idx_filesuffix2);
      }
      exit(9);
    }


    if (user_chrsubsetname != NULL) {
      if ((divno = Univ_IIT_find_one(chromosome_iit,user_chrsubsetname)) < 0) {
	fprintf(stderr,"Cannot find chrsubset %s in chromosome IIT file.  Ignoring.\n",user_chrsubsetname);
      } else {
	interval = Univ_IIT_interval(chromosome_iit,divno);
	chrsubset_start = Univinterval_low(interval);
	chrsubset_end = Univinterval_high(interval);
      }
    }
  }

  FREE(genomesubdir);
  FREE(fileroot);
  FREE(dbroot);

  Genome_setup(circular_typeint);
  Request_setup(global_genome,global_genomealt);


  if (splicing_file != NULL && global_genome != NULL) {
    /* TODO: Handle case for observed distances */
    /* min_extra_end no longer used by gregion.c */
    /* min_extra_end = shortsplicedist; */

    splicing_divint_crosstable = Univ_IIT_divint_crosstable(chromosome_iit,splicing_iit);
    if ((donor_typeint = IIT_typeint(splicing_iit,"donor")) >= 0 && 
	(acceptor_typeint = IIT_typeint(splicing_iit,"acceptor")) >= 0) {
      fprintf(stderr,"found donor and acceptor tags, so treating as splicesites file\n");
      splicestringpool = Splicestringpool_new();
      splicesites = Splicetrie_retrieve_via_splicesites(&distances_observed_p,&splicetypes,&splicedists,
							&splicestrings,&splicefrags_ref,&splicefrags_alt,
							&nsplicesites,splicing_iit,splicing_divint_crosstable,
							donor_typeint,acceptor_typeint,chromosome_iit,
							global_genome,global_genomealt,shortsplicedist,splicestringpool);
      if (nsplicesites == 0) {
	fprintf(stderr,"\nWarning: No splicesites observed for genome %s.  Are you sure this splicesite file was built for this genome?  Please compare chromosomes below:\n",
		dbroot);
	fprintf(stderr,"Chromosomes in the genome: ");
	Univ_IIT_dump_labels(stderr,chromosome_iit);
	fprintf(stderr,"Chromosomes in the splicesites IIT file: ");
	IIT_dump_divstrings(stderr,splicing_iit);
	exit(9);

      } else {
	Splicetrie_npartners(&nsplicepartners_skip,&nsplicepartners_obs,&nsplicepartners_max,splicesites,splicetypes,splicedists,
			     splicestrings,nsplicesites,chromosome_iit,shortsplicedist,distances_observed_p);
	Splicetrie_build_via_splicesites(&triecontents_obs,&trieoffsets_obs,&triecontents_max,&trieoffsets_max,
					 nsplicepartners_skip,nsplicepartners_obs,nsplicepartners_max,splicetypes,
					 splicestrings,nsplicesites);
	FREE(nsplicepartners_max);
	FREE(nsplicepartners_obs);
	FREE(nsplicepartners_skip);
	/* Splicestring_gc(splicestrings,nsplicesites); */
	FREE(splicestrings);
      }
      Splicestringpool_free(&splicestringpool);

    } else {
      fprintf(stderr,"no donor or acceptor tags found, so treating as introns file\n");
      splicestringpool = Splicestringpool_new();
      splicesites = Splicetrie_retrieve_via_introns(&splicetypes,&splicedists,
						    &splicestrings,&splicefrags_ref,&splicefrags_alt,
						    &nsplicesites,splicing_iit,splicing_divint_crosstable,
						    chromosome_iit,global_genome,global_genomealt,splicestringpool);
      if (nsplicesites == 0) {
	fprintf(stderr,"\nWarning: No splicesites observed for genome %s.  Are you sure this splicesite file was built for this genome?  Please compare chromosomes below:\n",
		dbroot);
	fprintf(stderr,"Chromosomes in the genome: ");
	Univ_IIT_dump_labels(stderr,chromosome_iit);
	fprintf(stderr,"Chromosomes in the splicesites IIT file: ");
	IIT_dump_divstrings(stderr,splicing_iit);
	exit(9);
      } else {
	Splicetrie_build_via_introns(&triecontents_obs,&trieoffsets_obs,splicesites,splicetypes,
				     splicestrings,nsplicesites,chromosome_iit,splicing_iit,splicing_divint_crosstable);
	triecontents_max = (Triecontent_T *) NULL;
	trieoffsets_max =  (Trieoffset_T *) NULL;
	/* Splicestring_gc(splicestrings,nsplicesites); */
	FREE(splicestrings);
      }
      Splicestringpool_free(&splicestringpool);

    }
    
    fprintf(stderr,"done\n");
  }


  Translation_setup(translation_code,alt_initiation_codons_p);

#ifdef PMAP
  Alphabet_setup(alphabet,alphabet_size,index1part_aa);
  Oligoindex_pmap_setup(genomecomp);
  Oligop_setup(alphabet,alphabet_size,index1part_aa);
  Indexdb_setup(index1part_aa);
  Stage1_setup(index1part_aa,maxextension,maxtotallen_bound,circular_typeint);
#else
  Oligoindex_hr_setup(mode);
  Oligo_setup(mode);
  Indexdb_setup(index1part);
  Stage1_setup(index1part,maxextension,maxtotallen_bound,circular_typeint);
#endif

  Stage2_setup(/*splicingp*/novelsplicingp == true || knownsplicingp == true,cross_species_p,
	       suboptimal_score_start,suboptimal_score_end,sufflookback,nsufflookback,maxintronlen,mode,
	       /*snps_p*/global_genomealt == global_genome ? false : true);
  Dynprog_single_setup(user_open,user_extend,user_dynprog_p,homopolymerp);
  Dynprog_genome_setup(novelsplicingp,splicing_iit,splicing_divint_crosstable,
		       donor_typeint,acceptor_typeint,
		       user_open,user_extend,user_dynprog_p);
  Dynprog_end_setup(splicesites,splicetypes,splicedists,nsplicesites,
		    trieoffsets_obs,triecontents_obs,trieoffsets_max,triecontents_max,
		    user_open,user_extend,user_dynprog_p);
  Pair_setup(novelsplicingp,splicing_iit,trim_indel_score,print_margin_p,
	     gff3_separators_p,sam_insert_0M_p,force_xs_direction_p,
	     md_lowercase_variant_p,/*snps_p*/global_genomealt == global_genome ? false : true,
	     gff3_phase_swap_p,cdstype,sam_cigar_extended_p,cigar_action);

  Stage3_setup(/*splicingp*/novelsplicingp == true || knownsplicingp == true,novelsplicingp,
	       require_splicedir_p,shortsplicedist,splicing_iit,splicing_divint_crosstable,
	       donor_typeint,acceptor_typeint,splicesites,circularp,altlocp,alias_starts,alias_ends,
	       min_intronlength,max_deletionlength,/*min_indel_end_matches*/6,
	       maxpeelback_distalmedial,nullgap,extramaterial_end,extramaterial_paired,
	       extraband_single,extraband_end,extraband_paired,
	       ngap,maxintronlen,maxintronlen_ends,end_trimming_score,minendexon,
	       user_open,user_extend,user_dynprog_p,homopolymerp,gff3_fasta_annotation_type,
	       sam_flippedp,stage3debug,/*genome_totallength*/genomelength);

  Splicetrie_setup(splicesites,splicefrags_ref,splicefrags_alt,
		   trieoffsets_obs,triecontents_obs,trieoffsets_max,triecontents_max,
		   /*snpp*/false,amb_closest_p,/*amb_clip_p*/true,/*min_shortend*/2);
  Output_setup(chromosome_iit,nofailsp,failsonlyp,quiet_if_excessive_p,maxpaths_report,
	       failedinput_root,quality_shift,
	       printtype,invertmode,wraplength,ngap,nointronlenp,sam_paired_p,sam_flippedp,
	       cds_startpos,fulllengthp,truncatep,strictp,checksump,
	       dbversion,user_chrsubsetname,contig_iit,altstrain_iit,
	       /*chimeras_allowed_p*/chimera_margin > 0 ? true : false,
	       map_iit,map_divint_crosstable,map_exons_p,map_bothstrands_p,
	       nflanking,print_comment_p,sam_read_group_id);

  Inbuffer_setup(/*single_cell_p*/false,/*filter_if_both_p*/false);

  Filestring_setup(/*split_simple_p*/false);
  SAM_header_setup(argc,argv,optind,nworkers,orderedp,
		   chromosome_iit,printtype,global_genome,
		   sam_headers_p,sam_read_group_id,sam_read_group_name,
		   sam_read_group_library,sam_read_group_platform);
  Outbuffer_setup(any_circular_p,quiet_if_excessive_p,
		  /*paired_end_p*/false,appendp,/*output_file*/NULL,
		  /*split_simple_p*/false,split_output_root,failedinput_root);
  outbuffer = Outbuffer_new(output_buffer_size,nread);
  Inbuffer_set_outbuffer(inbuffer,outbuffer);

  fprintf(stderr,"Starting alignment\n");
  stopwatch = Stopwatch_new();
  Stopwatch_start(stopwatch);


#if !defined(HAVE_PTHREAD)
  /* Serial version */
  single_thread();

#else
  /* Pthreads version */
  if (nworkers == 0) {
    single_thread();
    
  } else if (multiple_sequences_p == false) {
    single_thread();
    
  } else {
#ifdef WORKER_DETACH
    pthread_attr_init(&thread_attr_detach);
    if ((ret = pthread_attr_setdetachstate(&thread_attr_detach,PTHREAD_CREATE_DETACHED)) != 0) {
      fprintf(stderr,"ERROR: pthread_attr_setdetachstate %d\n",ret);
      exit(1);
    }
#endif
    pthread_attr_init(&thread_attr_join);
    if ((ret = pthread_attr_setdetachstate(&thread_attr_join,PTHREAD_CREATE_JOINABLE)) != 0) {
      fprintf(stderr,"ERROR: pthread_attr_setdetachstate %d\n",ret);
      exit(1);
    }
    
    worker_thread_ids = (pthread_t *) CALLOC(nworkers,sizeof(pthread_t));
    Except_init_pthread();
    pthread_key_create(&global_request_key,NULL);

    if (orderedp == true) {
      pthread_create(&output_thread_id,&thread_attr_join,Outbuffer_thread_ordered,
		     (void *) outbuffer);
    } else {
      pthread_create(&output_thread_id,&thread_attr_join,Outbuffer_thread_anyorder,
		     (void *) outbuffer);
    }

    for (i = 0; i < nworkers; i++) {
#ifdef WORKER_DETACH
      pthread_create(&(worker_thread_ids[i]),&thread_attr_detach,worker_thread,(void *) NULL);
#else
      /* Need to have worker threads finish before we call Inbuffer_free() */
      pthread_create(&(worker_thread_ids[i]),&thread_attr_join,worker_thread,(void *) NULL);
#endif
    }
    
    pthread_join(output_thread_id,NULL);
    for (i = 0; i < nworkers; i++) {
      pthread_join(worker_thread_ids[i],NULL);
    }

    pthread_key_delete(global_request_key);
    /* Do not delete global_except_key, because worker threads might still need it */
    /* Except_term_pthread(); */

    FREE(worker_thread_ids);

  }
#endif /* HAVE_PTHREAD */


  /* Single CPU or Pthreads version */
  runtime = Stopwatch_stop(stopwatch);
  Stopwatch_free(&stopwatch);

  nread = Outbuffer_nread(outbuffer);
  /* nbeyond = Outbuffer_nbeyond(outbuffer); */
  fprintf(stderr,"Processed %u queries in %.2f seconds (%.2f queries/sec)\n",
	  nread,runtime,(double) nread/runtime);

  Outbuffer_free(&outbuffer);
  Inbuffer_free(&inbuffer);	/* Also closes inputs */

  Outbuffer_close_files();

#ifdef PMAP
  Backtranslation_term();
#endif
  Dynprog_term(mode);


  if (nsplicesites > 0) {
    if (splicetrie_precompute_p == true) {
      FREE(triecontents_max);
      FREE(trieoffsets_max);
      FREE(triecontents_obs);
      FREE(trieoffsets_obs);
    } else {
      FREE(nsplicepartners_max);
      FREE(nsplicepartners_obs);
      FREE(nsplicepartners_skip);
      /* Splicestring_gc(splicestrings,nsplicesites); */
      FREE(splicestrings);
    }
    FREE(splicefrags_ref);
    FREE(splicedists);
    FREE(splicetypes);
    FREE(splicesites);
  }

  if (splicing_iit != NULL) {
    FREE(splicing_divint_crosstable);
    IIT_free(&splicing_iit);
  }

#if 0
  /* Oligoindex_localdb_cleanup(); */
  if (localdb != NULL) {
    Localdb_free(&localdb);
  }
#endif

#ifdef PMAP
 if (indexdb_rev != NULL) {
    Indexdb_free(&indexdb_rev);
  }
  if (indexdb_fwd != NULL) {
    Indexdb_free(&indexdb_fwd);
  }
#else
  if (indexdb_rev != indexdb_fwd) {
    Indexdb_free(&indexdb_rev);
  }
  if (indexdb_fwd != NULL) {
    Indexdb_free(&indexdb_fwd);
  }
#endif
  if (dbversion != NULL) {
    FREE(dbversion);
  }
  if (altstrain_iit != NULL) {
    IIT_free(&altstrain_iit);
  }

  if (user_genomes != NULL) {
    for (p = user_genomes; p != NULL; p = List_next(p)) {
      genome = (Genome_T) List_head(p);
      Genome_free(&genome);
    }
    List_free(&user_genomes);
  }

  if (global_genomealt != NULL && global_genomealt != global_genome) {
    Genome_free(&global_genomealt);
  }
  if (global_genome != NULL) {
    Genome_free(&global_genome);
  }

  if (map_iit != NULL) {
    IIT_free(&map_iit);
  }
  if (contig_iit != NULL) {
    Univ_IIT_free(&contig_iit);
  }
  if (altlocp != NULL) {
    FREE(alias_ends);
    FREE(alias_starts);
    FREE(altlocp);
  }
  if (circularp != NULL) {
    FREE(circularp);
  }
  if (chromosome_iit != NULL) {
    Univ_IIT_free(&chromosome_iit);
  }

  Outbuffer_cleanup();

  Access_controlled_cleanup();

  return 0;
}


static void
print_program_usage () {
#ifdef PMAP
    fprintf(stdout,"\
Usage: pmap [OPTIONS...] <FASTA files...>, or\n\
       cat <FASTA files...> | pmap [OPTIONS...]\n\
");
#else
    fprintf(stdout,"\
Usage: gmap [OPTIONS...] <FASTA files...>, or\n\
       cat <FASTA files...> | gmap [OPTIONS...]\n\
");
#endif
    fprintf(stdout,"\n");

    fprintf(stdout,"Input options (must include -d or -g)\n");
    fprintf(stdout,"\
  -D, --dir=directory            Genome directory.  Default (as specified by --with-gmapdb to the configure program) is\n \
                                   %s\n\
",GMAPDB);
    fprintf(stdout,"\
  -d, --db=STRING                Genome database.  If argument is '?' (with\n\
                                   the quotes), this command lists available databases.\n\
");
    fprintf(stdout,"\n");

#ifdef PMAP
    fprintf(stdout,"\
  -a, --alphabet=STRING          Alphabet to use in PMAP genome database\n\
                                   (allowed values in order of preference: 20, 15a, 12a).\n\
                                   If not specified, the program will find the first available\n\
                                   alphabet in the genome database in preference order\n\
");
#endif

#if 0
    /* No longer supported */
    fprintf(stdout,"\
    -G, --genomefull               Use full genome (all ASCII chars allowed;\n \
                                   built explicitly during setup), not\n\
                                   compressed version\n\
");
#endif

    fprintf(stdout,"\
  -k, --kmer=INT                 kmer size to use in genome database (allowed values: 16 or less).\n\
                                   If not specified, the program will find the highest available\n\
                                   kmer size in the genome database\n\
  --sampling=INT                 Sampling to use in genome database.  If not specified, the program\n\
                                   will find the smallest available sampling value in the genome database\n\
                                   within selected k-mer size\n\
  -g, --gseg=filename            User-supplied genomic segments.  If multiple segments are provided, then\n\
                                   every query sequence is aligned against every genomic segment\n\
  -1, --selfalign                Align one sequence against itself in FASTA format via stdin\n\
                                   (Useful for getting protein translation of a nucleotide sequence)\n\
  -2, --pairalign                Align two sequences in FASTA format via stdin, first one being\n\
                                   genomic and second one being cDNA\n\
\n\
  --cmdline=STRING,STRING        Align these two sequences provided on the command line,\n\
                                   first one being genomic and second one being cDNA\n\
  -q, --part=INT/INT             Process only the i-th out of every n sequences\n\
                                   e.g., 0/100 or 99/100 (useful for distributing jobs\n\
                                   to a computer farm).\n\
");
    fprintf(stdout,"\
  --input-buffer-size=INT        Size of input buffer (program reads this many sequences\n\
                                   at a time for efficiency) (default %d)\n\
",input_buffer_size);
    fprintf(stdout,"\n");

    fprintf(stdout,"Computation options\n");
#ifdef HAVE_MMAP
    fprintf(stdout,"\
  -B, --batch=INT                Batch mode (default = 2)\n\
                                 Mode     Positions       Genome\n\
                                   0      mmap            mmap\n\
                                   1      mmap & preload  mmap\n\
                      (default)    2      mmap & preload  mmap & preload\n\
                                   3      allocate        mmap & preload\n\
                                   4      allocate        allocate\n\
                                   5      allocate        allocate     (same as 4)\n\
                           Note: For a single sequence, all data structures use mmap\n\
                           If mmap not available and allocate not chosen, then will use fileio (very slow)\n\
");
#else
    fprintf(stdout,"\
  -B, --batch=INT                Batch mode (default = 4, modes 0-3 disallowed because program configured without mmap)\n\
                                 Mode     Positions       Genome\n\
                      (default)    4      allocate        allocate\n\
                                   5      allocate        allocate     (same as 4)\n\
");
#endif

  fprintf(stdout,"\
  --use-shared-memory=INT        If 1, then allocated memory is shared among all processes on this node\n\
                                   If 0 (default), then each process has private allocated memory\n\
");

    fprintf(stdout,"\
  --nosplicing                   Turns off splicing (useful for aligning genomic sequences\n\
                                   onto a genome)\n\
");
    fprintf(stdout,"\
  --max-deletionlength=INT       Max length for a deletion (default %d).  Above this size,\n\
                                   a genomic gap will be considered an intron rather than a deletion.\n\
                                   If the genomic gap is less than --max-deletionlength and greater\n\
                                   than --min-intronlength, a known splice site or splice site probabilities\n\
                                   of 0.80 on both sides will be reported as an intron.\n\
",max_deletionlength);
    fprintf(stdout,"\
  --min-intronlength=INT         Min length for one internal intron (default %d).  Below this size,\n\
                                   a genomic gap will be considered a deletion rather than an intron.\n\
                                   If the genomic gap is less than --max-deletionlength and greater\n\
                                   than --min-intronlength, a known splice site or splice site probabilities\n\
                                   of 0.80 on both sides will be reported as an intron.\n\
",min_intronlength);
    fprintf(stdout,"\
  --max-intronlength-middle=INT  Max length for one internal intron (default %d).  Note: for backward\n\
                                   compatibility, the -K or --intronlength flag will set both\n\
                                   --max-intronlength-middle and --max-intronlength-ends.\n\
                                   Also see --split-large-introns below.\n\
",maxintronlen);
    fprintf(stdout,"\
  --max-intronlength-ends=INT    Max length for first or last intron (default %d).  Note: for backward\n\
                                   compatibility, the -K or --intronlength flag will set both\n\
                                   --max-intronlength-middle and --max-intronlength-ends.\n\
",maxintronlen_ends);
    fprintf(stdout,"\
  --split-large-introns          Sometimes GMAP will exceed the value for --max-intronlength-middle,\n\
                                   if it finds a good single alignment.  However, you can force GMAP\n\
                                   to split such alignments by using this flag\n\
");
    fprintf(stdout,"\
  --end-trimming-score=INT       Trim ends if the alignment score is below this value\n\
                                   where a match scores +1 and a mismatch scores -3\n\
                                   The value should be 0 (default) or negative.  A negative\n\
                                   allows some mismatches at the ends of the alignment\n\
  --trim-end-exons=INT           Trim end exons with fewer than given number of matches\n\
                                   (in nt, default %d)\n\
",minendexon);
    fprintf(stdout,"\
  -w, --localsplicedist=INT      Max length for known splice sites at ends of sequence\n\
                                   (default %d)\n\
",shortsplicedist);
    fprintf(stdout,"\
  -L, --totallength=INT          Max total intron length (default %d)\n\
",maxtotallen_bound);
    fprintf(stdout,"\
  -x, --chimera-margin=INT       Amount of unaligned sequence that triggers\n\
                                   search for the remaining sequence (default %d).\n\
                                   Enables alignment of chimeric reads, and may help\n\
                                   with some non-chimeric reads.  To turn off, set to\n\
                                   zero.\n\
",chimera_margin);
    fprintf(stdout,"\
  --no-chimeras                  Turns off finding of chimeras.  Same effect as --chimera-margin=0\n\
");

#if 0
    fprintf(stdout,"\
  -w, --reference=filename       Reference cDNA sequence for relative alignment\n\
");
#endif

#ifdef HAVE_PTHREAD
    fprintf(stdout,"\
  -t, --nthreads=INT             Number of worker threads\n\
");
#else
  fprintf(stdout,"\
  -t, --nthreads=INT             Number of worker threads.  Flag is ignored in this version of GMAP, which has pthreads disabled\n\
");
#endif
    fprintf(stdout,"\
  -c, --chrsubset=string         Limit search to given chromosome\n\
  --strand=STRING                Genome strand to try aligning to (plus, minus, or both default)\n\
  -z, --direction=STRING         cDNA direction (sense_force, antisense_force,\n\
                                   sense_filter, antisense_filter,or auto (default))\n\
");
    fprintf(stdout,"\
  --canonical-mode=INT           Reward for canonical and semi-canonical introns\n\
                                   0=low reward, 1=high reward (default), 2=low reward for\n\
                                   high-identity sequences and high reward otherwise\n\
  --cross-species                Use a more sensitive search for canonical splicing, which helps especially\n\
                                   for cross-species alignments and other difficult cases\n\
  --allow-close-indels=INT       Allow an insertion and deletion close to each other\n\
                                   (0=no, 1=yes (default), 2=only for high-quality alignments)\n\
");
    fprintf(stdout,"\
  --microexon-spliceprob=FLOAT   Allow microexons only if one of the splice site probabilities is\n\
                                   greater than this value (default %.2f)\n\
",microexon_spliceprob);

    fprintf(stdout,"\
  --indel-open                   In dynamic programming, opening penalty for indel\n\
  --indel-extend                 In dynamic programming, extension penalty for indel\n\
                                   Values for --indel-open and --indel-extend should be in [-127,-1].\n\
                                   If value is < -127, then will use -127 instead.\n\
                                   If --indel-open and --indel-extend are not specified, values are chosen\n\
                                   adaptively, based on the differences between the query and reference\n\
");

#if 0
    fprintf(stdout,"\
  --homopolymer                  In dynamic programming, favor indels in regions of homopolymers,\n\
                                   e.g., AAAAAA.  Useful for some platforms, such as Pacific Biosciences\n\
");
#endif

#ifndef PMAP
    fprintf(stdout,"\
  --cmetdir=STRING               Directory for methylcytosine index files (created using cmetindex)\n\
                                   (default is location of genome index files specified using -D, -V, and -d)\n\
  --atoidir=STRING               Directory for A-to-I RNA editing index files (created using atoiindex)\n\
                                   (default is location of genome index files specified using -D, -V, and -d)\n\
  --mode=STRING                  Alignment mode: standard (default), cmet-stranded, cmet-nonstranded,\n\
                                    atoi-stranded, atoi-nonstranded, ttoc-stranded, or ttoc-nonstranded.\n\
                                    Non-standard modes requires you to have previously run the cmetindex\n\
                                    or atoiindex programs (which also cover the ttoc modes) on the genome\n\
");
#endif

#if 0
    /* Causes seg faults, so do not advertise */
    fprintf(stdout,"\
  -s, --splicing=STRING          Look for splicing involving known sites\n\
                                   (in <STRING>.iit)\n\
");
#endif

#ifndef PMAP
    fprintf(stdout,"\
  -p, --prunelevel               Pruning level: 0=no pruning (default), 1=poor seqs,\n\
                                   2=repetitive seqs, 3=poor and repetitive\n\
");
#endif
    fprintf(stdout,"\n");

    fprintf(stdout,"\
Output types\n\
  -S, --summary                  Show summary of alignments only\n\
  -A, --align                    Show alignments\n\
  -3, --continuous               Show alignment in three continuous lines\n\
  -4, --continuous-by-exon       Show alignment in three lines per exon\n\
  -E, --exons=STRING             Print exons (\"cdna\" or \"genomic\")\n\
                                   Will also print introns with \"cdna+introns\" or\n\
                                   \"genomic+introns\"\n\
");

#ifdef PMAP    
    fprintf(stdout,"\
  -P, --protein_gen              Print protein sequence (genomic)\n\
  -Q, --nucleotide               Print inferred nucleotide sequence from protein\n\
");
#else
    fprintf(stdout,"\
  -P, --protein_dna              Print protein sequence (cDNA)\n\
  -Q, --protein_gen              Print protein sequence (genomic)\n\
");
#endif

#ifdef PMAP
    fprintf(stdout,"\
  -f, --format=INT               Other format for output (also note the -A and -S options\n\
                                   and other options listed under Output types):\n\
                                   mask_introns,\n\
                                   mask_utr_introns,\n\
                                   psl_pro (or 0) = PSL format in protein coords,\n\
                                   psl_nt (or 1) = PSL format in nucleotide coords,\n\
                                   gff3_gene (or 2) = GFF3 gene format,\n\
                                   gff3_match_cdna (or 3) = GFF3 cDNA_match format,\n\
                                   gff3_match_est (or 4) = GFF3 EST_match format,\n\
                                   map_exons (or 7) = IIT FASTA exon map format,\n\
                                   map_ranges (or 8) = IIT FASTA range map format,\n\
                                   coords (or 9) = coords in table format\n\
");
#else
    fprintf(stdout,"\
  -f, --format=INT               Other format for output (also note the -A and -S options\n\
                                   and other options listed under Output types):\n\
                                   mask_introns,\n\
                                   mask_utr_introns,\n\
                                   psl (or 1) = PSL (BLAT) format,\n\
                                   gff3_gene (or 2) = GFF3 gene format,\n\
                                   gff3_match_cdna (or 3) = GFF3 cDNA_match format,\n\
                                   gff3_match_est (or 4) = GFF3 EST_match format,\n\
                                   splicesites (or 6) = splicesites output (for GSNAP splicing file),\n\
                                   introns = introns output (for GSNAP splicing file),\n\
                                   map_exons (or 7) = IIT FASTA exon map format,\n\
                                   map_ranges (or 8) = IIT FASTA range map format,\n\
                                   coords (or 9) = coords in table format,\n\
                                   sampe = SAM format (setting paired_read bit in flag),\n\
                                   samse = SAM format (without setting paired_read bit),\n\
                                   bedpe = indels and gaps in BEDPE format\n\
");
#endif
    fprintf(stdout,"\n");

    fprintf(stdout,"\
Output options\n\
  -n, --npaths=INT               Maximum number of paths to show (default %d).  If set to 1, GMAP\n\
                                   will not report chimeric alignments, since those imply\n\
                                   two paths.  If you want a single alignment plus chimeric\n\
                                   alignments, then set this to be 0.\n\
",maxpaths_report);
    fprintf(stdout,"\
  --suboptimal-score=FLOAT       Report only paths whose score is within this value of the\n\
                                   best path.\n\
                                 If specified between 0.0 and 1.0, then treated as a fraction\n\
                                   of the score of the best alignment (matches minus penalties for\n\
                                   mismatches and indels).  Otherwise, treated as an integer\n\
                                   number to be subtracted from the score of the best alignment.\n\
                                   Default value is 0.50.\n\
  -O, --ordered                  Print output in same order as input (relevant\n\
                                   only if there is more than one worker thread)\n\
  -5, --md5                      Print MD5 checksum for each query sequence\n\
  -o, --chimera-overlap          Overlap to show, if any, at chimera breakpoint\n\
  --failsonly                    Print only failed alignments, those with no results\n\
  --nofails                      Exclude printing of failed alignments\n\
\n\
  -V, --snpsdir=STRING           Directory for SNPs index files (created using snpindex) (default is\n\
                                   location of genome index files specified using -D and -d)\n \
  -v, --use-snps=STRING          Use database containing known SNPs (in <STRING>.iit, built\n\
                                   previously using snpindex) for tolerance to SNPs\n\
");

  fprintf(stdout,"\
  --split-output=STRING          Basename for multiple-file output, separately for nomapping,\n\
                                   uniq, mult, (and chimera, if --chimera-margin is selected)\n\
  --failed-input=STRING          Print completely failed alignments as input FASTA or FASTQ format\n\
                                   to the given file.  If the --split-output flag is also given, this file\n\
                                   is generated in addition to the output in the .nomapping file.\n\
  --append-output                When --split-output or --failedinput is given, this flag will append output\n\
                                   to the existing files.  Otherwise, the default is to create new files.\n\
");
  fprintf(stdout,"\
  --output-buffer-size=INT       Buffer size, in queries, for output thread (default %d).  When the number\n\
                                   of results to be printed exceeds this size, worker threads wait\n\
                                   until the backlog is cleared\n\
",output_buffer_size);


  fprintf(stdout,"\
  --translation-code=INT         Genetic code used for translating codons to amino acids and computing CDS\n\
                                   Integer value (default=1) corresponds to an available code at\n\
                                   http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi\n\
  --alt-start-codons             Also, use the alternate initiation codons shown in the above Web site\n\
                                   By default, without this option, only ATG is considered an initiation codon\n\
");

#ifdef PMAP    
    fprintf(stdout,"\
  -Y, --tolerant                 Translates genome with corrections for frameshifts\n\
");
#else
    fprintf(stdout,"\
  -F, --fulllength               Assume full-length protein, starting with Met\n\
  -a, --cdsstart=INT             Translate codons from given nucleotide (1-based)\n\
  -T, --truncate                 Truncate alignment around full-length protein, Met to Stop\n\
                                 Implies -F flag.\n\
  -Y, --tolerant                 Translates cDNA with corrections for frameshifts\n\
");
#endif

    fprintf(stdout,"\n");

#ifndef PMAP
  fprintf(stdout,"Options for GFF3 output\n");
  fprintf(stdout,"\
  --gff3-add-separators=INT      Whether to add a ### separator after each query sequence\n\
                                   Values: 0 (no), 1 (yes, default)\n\
  --gff3-swap-phase=INT          Whether to swap phase (0 => 0, 1 => 2, 2 => 1) in gff3_gene format\n\
                                   Needed by some analysis programs, but deviates from GFF3 specification\n\
                                   Values: 0 (no, default), 1 (yes)\n\
  --gff3-fasta-annotation=INT    Whether to include annotation from the FASTA header into the GFF3 output\n\
                                   Values: 0 (default): Do not include\n\
                                           1: Wrap all annotation as Annot=\"<header>\"\n\
                                           2: Include key=value pairs, replacing brackets with quotation marks\n\
                                              and replacing spaces between key=value pairs with semicolons\n\
  --gff3-cds=STRING              Whether to use cDNA or genomic translation for the CDS coordinates\n\
                                   Values: cdna (default), genomic\n\
");
  fprintf(stdout,"\n");

  fprintf(stdout,"Options for SAM output\n");
  fprintf(stdout,"\
  --no-sam-headers               Do not print headers beginning with '@'\n\
  --sam-use-0M                   Insert 0M in CIGAR between adjacent insertions and deletions\n\
                                   Required by Picard, but can cause errors in other tools\n\
  --sam-extended-cigar           Use extended CIGAR format (using X and = symbols instead of M,\n\
                                   to indicate matches and mismatches, respectively\n\
  --sam-flipped                  Flip the query and genomic positions in the SAM output.\n\
                                   Potentially useful with the -g flag when short reads are picked as query\n\
                                   sequences and longer reads as picked as genomic sequences\n\
  --force-xs-dir                 For RNA-Seq alignments, disallows XS:A:? when the sense direction\n\
                                   is unclear, and replaces this value arbitrarily with XS:A:+.\n\
                                   May be useful for some programs, such as Cufflinks, that cannot\n\
                                   handle XS:A:?.  However, if you use this flag, the reported value\n\
                                   of XS:A:+ in these cases will not be meaningful.\n\
  --md-lowercase-snp             In MD string, when known SNPs are given by the -v flag,\n\
                                   prints difference nucleotides as lower-case when they,\n\
                                   differ from reference but match a known alternate allele\n\
  --action-if-cigar-error        Action to take if there is a disagreement between CIGAR length and sequence length\n\
                                   Allowed values: ignore, warning (default), noprint, abort\n\
                                   Note that the noprint option does not print the CIGAR string at all if there\n\
                                   is an error, so it may break a SAM parser\n\
  --read-group-id=STRING         Value to put into read-group id (RG-ID) field\n\
  --read-group-name=STRING       Value to put into read-group name (RG-SM) field\n\
  --read-group-library=STRING    Value to put into read-group library (RG-LB) field\n\
  --read-group-platform=STRING   Value to put into read-group library (RG-PL) field\n\
");
  fprintf(stdout,"\n");

  /* Quality score options */
  fprintf(stdout,"Options for quality scores\n");
  fprintf(stdout,"\
  --quality-protocol=STRING      Protocol for input quality scores.  Allowed values:\n\
                                   illumina (ASCII 64-126) (equivalent to -J 64 -j -31)\n\
                                   sanger   (ASCII 33-126) (equivalent to -J 33 -j 0)\n\
                                 Default is sanger (no quality print shift)\n\
                                 SAM output files should have quality scores in sanger protocol\n\
\n\
                                 Or you can specify the print shift with this flag:\n\
  -j, --quality-print-shift=INT  Shift FASTQ quality scores by this amount in output\n\
                                   (default is 0 for sanger protocol; to change Illumina input\n\
                                   to Sanger output, select -31)\n\
");
#endif

    fprintf(stdout,"\
External map file options\n\
  -M, --mapdir=directory         Map directory\n\
  -m, --map=iitfile              Map file.  If argument is '?' (with the quotes),\n\
                                   this lists available map files.\n\
  -e, --mapexons                 Map each exon separately\n\
  -b, --mapboth                  Report hits from both strands of genome\n\
  -u, --flanking=INT             Show flanking hits (default 0)\n\
  --print-comment                Show comment line for each hit\n\
");
    fprintf(stdout,"\n");

    fprintf(stdout,"\
Alignment output options\n\
  --nolengths                    No intron lengths in alignment\n\
  --nomargin                     No left margin in GMAP standard output (with the -A flag)\n\
  -I, --invertmode=INT           Mode for alignments to genomic (-) strand:\n\
                                   0=Don't invert the cDNA (default)\n\
                                   1=Invert cDNA and print genomic (-) strand\n\
                                   2=Invert cDNA and print genomic (+) strand\n\
");
    fprintf(stdout,"\
  -i, --introngap=INT            Nucleotides to show on each end of intron (default %d)\n\
",ngap);
    fprintf(stdout,"\
  -l, --wraplength=INT           Wrap length for alignment (default %d)\n\
",wraplength);
    fprintf(stdout,"\n");

    fprintf(stdout,"\
Filtering output options\n\
  --min-trimmed-coverage=FLOAT   Do not print alignments with trimmed coverage less\n\
                                   this value (default=0.0, which means no filtering)\n\
                                   Note that chimeric alignments will be output regardless\n\
                                   of this filter\n\
  --min-identity=FLOAT           Do not print alignments with identity less\n\
                                   this value (default=0.0, which means no filtering)\n\
                                   Note that chimeric alignments will be output regardless\n\
                                   of this filter\n\
\n\
Help options\n\
  --check                        Check compiler assumptions\n\
  --version                      Show version\n\
  --help                         Show this help message\n\
");

    return;
}
