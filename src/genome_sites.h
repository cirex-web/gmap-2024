/* $Id: genome_sites.h 223511 2020-11-14 15:50:08Z twu $ */
#ifndef GENOME_SITES_INCLUDED
#define GENOME_SITES_INCLUDED

#include "bool.h"
#include "types.h"
#include "univcoord.h"
#include "genomicpos.h"
#include "genome.h"

#define T Genome_T

extern void
Genome_sites_setup (T genome_in, T genomealt_in);

extern int
Genome_donor_sites (int *sites, int *knowni, int *old_knownpos, int *old_knowni,
		    Univcoord_T left, int pos5, int pos3);
extern int
Genome_acceptor_sites (int *sites, int *knowni, int *old_knownpos, int *old_knowni,
		       Univcoord_T left, int pos5, int pos3);

extern int
Genome_antidonor_sites (int *sites, int *knowni, int *old_knownpos, int *old_knowni,
			Univcoord_T left, int pos5, int pos3);

extern int
Genome_antiacceptor_sites (int *sites, int *knowni, int *old_knownpos, int *old_knowni,
			   Univcoord_T left, int pos5, int pos3);

#if 0
extern int
Genome_donor_sites_novel (int *sites, Univcoord_T left, int pos5, int pos3);
extern int
Genome_acceptor_sites_novel (int *sites, Univcoord_T left, int pos5, int pos3);
extern int
Genome_antidonor_sites_novel (int *sites, Univcoord_T left, int pos5, int pos3);
extern int
Genome_antiacceptor_sites_novel (int *sites, Univcoord_T left, int pos5, int pos3);
#endif


#if 0
extern Univcoord_T
Genome_prev_donor_position (Univcoord_T pos, Univcoord_T prevpos);
extern Univcoord_T
Genome_prev_acceptor_position (Univcoord_T pos, Univcoord_T prevpos);
extern Univcoord_T
Genome_prev_antidonor_position (Univcoord_T pos, Univcoord_T prevpos);
extern Univcoord_T
Genome_prev_antiacceptor_position (Univcoord_T pos, Univcoord_T prevpos);
#endif

#undef T
#endif

