/* $Id: uinttableuint.h 224642 2021-08-25 22:03:29Z twu $ */
#ifndef UINTTABLEUINT_INCLUDED
#define UINTTABLEUINT_INCLUDED
#include "bool.h"

#define T Uinttableuint_T
typedef struct T *T;

extern T
Uinttableuint_new (int hint);
extern void 
Uinttableuint_free (T *table);
extern int   
Uinttableuint_length (T table);
extern unsigned int
Uinttableuint_put (T table, const unsigned int key, unsigned int value);
extern unsigned int
Uinttableuint_get (T table, const unsigned int key);
extern unsigned int
Uinttableuint_remove (T table, const unsigned int key);
extern void   
Uinttableuint_map (T table,
	       void (*apply)(const unsigned int key, unsigned int *value, void *cl),
	       void *cl);
extern unsigned int *
Uinttableuint_keys (T table, bool sortp, unsigned int end);
extern unsigned int *
Uinttableuint_keys_by_timeindex (T table);
extern unsigned int *
Uinttableuint_values (T table, unsigned int end);

#undef T
#endif
