/* $Id: univdiag.h 225774 2023-01-27 19:13:10Z twu $ */
#ifndef UNIVDIAG_INCLUDED
#define UNIVDIAG_INCLUDED

typedef struct Univdiag_T *Univdiag_T;

#include "bool.h"
#include "list.h"
#include "genomicpos.h"
#include "types.h"

#define T Univdiag_T

#if 0
/* Now in univdiagpool.c */
extern T
Univdiag_new (int qstart, int qend, Univcoord_T univdiagonal);
#endif

extern void
Univdiag_free (T *old);
extern void
Univdiag_gc (List_T *list);
extern void
Univdiag_list_gc (List_T *paths);
extern int
Univdiag_list_length (List_T path);

extern int
Univdiag_ascending_cmp (const void *a, const void *b);
extern int
Univdiag_descending_cmp (const void *a, const void *b);
extern int
Univdiag_diagonal_cmp (const void *a, const void *b);

#undef T
#endif


