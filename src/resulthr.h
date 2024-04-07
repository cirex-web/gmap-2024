/* $Id: 3d7880d00a430fa5411696233577ce7f7acc3fa4 $ */
#ifndef RESULTHR_INCLUDED
#define RESULTHR_INCLUDED

typedef struct Result_T *Result_T;

#include "bool.h"
#include "samflags.h" 		/* for SAM_split_output_type */
#include "intlistpool.h"
#include "univcoord.h"
#include "listpool.h"
#include "pathpool.h"
#include "transcriptpool.h"
#include "hitlistpool.h"


/* PAIRED_UNSPECIFIED assigned only by Stage1hr_paired_read */

typedef enum {EXACT, SUB, SUBSTRINGS,
	      HALFSPLICE_DONOR, HALFSPLICE_ACCEPTOR, SPLICE, SAMECHR_SPLICE, TRANSLOC_SPLICE} Hittype_T;

typedef enum {CONCORDANT, PAIRED_UNSPECIFIED, PAIRED_INVERSION, PAIRED_SCRAMBLE, PAIRED_TOOLONG,
	      CONCORDANT_TRANSLOCATIONS, UNPAIRED, UNSPECIFIED} Pairtype_T;

typedef enum {SINGLEEND_NOMAPPING, PAIREDEND_NOMAPPING,
	      SINGLEEND_UNIQ, SINGLEEND_TRANSLOC, SINGLEEND_MULT,
	      PAIRED_UNIQ_INV, PAIRED_UNIQ_SCR, PAIRED_UNIQ_TOOLONG, PAIRED_MULT,
	      CONCORDANT_UNIQ, CONCORDANT_TRANSLOC, CONCORDANT_MULT,
	      HALFMAPPING_UNIQ, HALFMAPPING_TRANSLOC, HALFMAPPING_MULT,
	      UNPAIRED_UNIQ, UNPAIRED_TRANSLOC, UNPAIRED_MULT} Resulttype_T;


#define T Result_T

extern char *
Pairtype_string (Pairtype_T pairtype);
extern Resulttype_T
Result_resulttype (T this);
extern bool
Result_tr_consistent_p (T this);
extern char *
Resulttype_string (Resulttype_T resulttype);
extern int
Result_id (T this);
extern int
Result_worker_id (T this);
extern void **
Result_array (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq, T this);
extern void **
Result_array2 (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq, T this);
extern T
Result_single_read_new (int id, void **resultarray, int npaths_primary, int npaths_altloc,
			int first_absmq, int second_absmq);
extern T
Result_paired_read_new (int id, void **resultarray, int npaths_primary, int npaths_altloc,
			int first_absmq, int second_absmq,
			Pairtype_T final_pairtype);
extern T
Result_paired_as_singles_new (int id, void **hits5, int npaths5_primary, int npaths5_altloc,
			      int first_absmq5, int second_absmq5,
			      void **hits3, int npaths3_primary, int npaths3_altloc,
			      int first_absmq3, int second_absmq3);
extern void
Result_free (T *old, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	     Hitlistpool_T hitlistpool);

#undef T
#endif

