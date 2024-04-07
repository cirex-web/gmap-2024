#ifndef INTERSECT_APPROX_INDICES_UINT8_INCLUDED
#define INTERSECT_APPROX_INDICES_UINT8_INCLUDED

#include "univcoord.h"

extern int *
Intersect_approx_indices_uint8 (int *ndiagpairs,
				const Univcoord_T *set1, const int length1, int diagterm1,
				const Univcoord_T *set2, const int length2, int diagterm2,
				const Univcoord_T below_slop, const Univcoord_T above_slop);

extern void
Intersect_approx_indices_uint8_setup ();

#endif
