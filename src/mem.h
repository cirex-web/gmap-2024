/* $Id: 9bf89de2d61db39e1d1b261a3882a0b34d0e8bab $ */
#ifndef MEM_INCLUDED
#define MEM_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"		/* For HAVE_ALLOCA, HAVE_ALLOCA_H */
#endif

#include <stddef.h>
#include "except.h"

#ifdef HAVE_ALLOCA
#include <string.h>		/* For memset */
#endif

#ifdef HAVE_ALLOCA_H
#include <alloca.h>
#endif

#include "simd.h"
/* Uses SSE4_1 for MALLOC_ALIGN, needed by SIMD merge */


#define MAX_QUERYLENGTH_STACK 10000

extern const Except_T Mem_Failed;

extern void
Mem_trap_start (void **location, const char *file, int line);
extern void
Mem_trap_check (const char *file, int line);

/* MEMUSAGE checks for memory leak in each request */
/* #define MEMUSAGE 1 */
#ifdef MEMUSAGE

typedef enum {MAIN_STORAGE, INPUT_USAGE, OUTPUT_USAGE, WORKER_STORAGE} Memusage_Class_T;

extern void
Mem_usage_init ();
extern void
Mem_usage_set_threadname (char *threadname_in);
extern void
Mem_usage_reset_heap_baseline (long int x);
extern void
Mem_usage_reset_stack_max ();
extern void
Mem_usage_reset_heap_max ();
extern void
Mem_usage_std_stack_add (long int x, const char *file, int line);
extern void
Mem_usage_std_stack_subtract (const char *file, int line);
extern void
Mem_usage_std_heap_add (long int x);
extern long int
Mem_usage_report_std_stack ();
extern long int
Mem_usage_report_std_heap ();
extern long int
Mem_usage_report_std_stack_max ();
extern long int
Mem_usage_report_std_heap_max ();
extern long int
Mem_usage_report_keep ();
extern long int
Mem_usage_report_in ();
extern long int
Mem_usage_report_out ();
#endif


extern void *Mem_alloc (size_t nbytes, const char *file, int line);
extern void *Mem_alloc_keep (size_t nbytes, const char *file, int line);
extern void *Mem_alloc_in (size_t nbytes, const char *file, int line);
extern void *Mem_alloc_out (size_t nbytes, const char *file, int line);
extern void *Mem_alloc_no_exception (size_t nbytes,
				     const char *file, int line);
extern void *Mem_calloc (size_t count, size_t nbytes, const char *file, int line);
extern void *Mem_calloc_keep (size_t count, size_t nbytes, const char *file, int line);
extern void *Mem_calloc_in (size_t count, size_t nbytes, const char *file, int line);
extern void *Mem_calloc_out (size_t count, size_t nbytes, const char *file, int line);
extern void *Mem_calloc_no_exception (size_t count, size_t nbytes,
				      const char *file, int line);
extern void Mem_free (void *ptr, const char *file, int line);
extern void Mem_free_keep (void *ptr, const char *file, int line);
extern void Mem_free_in (void *ptr, const char *file, int line);
extern void Mem_free_out (void *ptr, const char *file, int line);
extern void *Mem_resize (void *ptr, size_t nbytes,
			 const char *file, int line);

#define MTRAP(location) Mem_trap_start((location), __FILE__, __LINE__)
#define MCHECK() Mem_trap_check(__FILE__, __LINE__)


#define MALLOC(nbytes) Mem_alloc((nbytes), __FILE__, __LINE__)
#define CALLOC(count, nbytes) Mem_calloc((count), (nbytes), __FILE__, __LINE__)
#define FREE(ptr) ((void)(Mem_free((ptr),__FILE__, __LINE__), (ptr) = 0))

#ifdef HAVE_ALLOCA
/* Use stack.  Need macros and not a function containing alloca, since the memory would be freed when the function returns */
#ifdef MEMUSAGE
#define ALLOCA(nbytes) alloca((nbytes)); Mem_usage_std_stack_add((nbytes),__FILE__,__LINE__)
#define MALLOCA(nbytes) alloca((nbytes)); Mem_usage_std_stack_add((nbytes),__FILE__,__LINE__)
#define CALLOCA(count, nbytes) memset(alloca((count) * (nbytes)),0,(count) * (nbytes)); Mem_usage_std_stack_add((count) * (nbytes),__FILE__,__LINE__)
#define FREEA(ptr) Mem_usage_std_stack_subtract(__FILE__,__LINE__)
#else
#define ALLOCA(nbytes) alloca((nbytes))
#define MALLOCA(nbytes) alloca((nbytes))
#define CALLOCA(count, nbytes) memset(alloca((count) * (nbytes)),0,(count) * (nbytes))
#define FREEA(ptr)
#endif
#else
/* Use heap: Same as MALLOC, CALLOC, and FREE */
#define MALLOCA(nbytes) Mem_alloc((nbytes), __FILE__, __LINE__)
#define CALLOCA(count, nbytes) Mem_calloc((count), (nbytes), __FILE__, __LINE__)
#define FREEA(ptr) ((void)(Mem_free((ptr),__FILE__, __LINE__), (ptr) = 0))
#endif


/* PAD_UINT8, PAD_UINT4: Pad lengths at end for row-based storage */

/* #define DEBUG_ALIGN 1 */

#ifdef DEBUG_ALIGN
#define debug_align(x) x
#else
#define debug_align(x)
#endif

#ifdef HAVE_ARM
#define MALLOC_ALIGN(new,x) posix_memalign((void **) &(new),16,x)
#ifdef DEBUG_ALIGN
#define FREE_ALIGN(x) printf("Freeing %p in standard pool at %s:%d\n",x,__FILE__,__LINE__); FREE(x)
#else
#define FREE_ALIGN(x) FREE(x)
#endif
#define CHECK_ALIGN(x) assert((unsigned long) x % 16 == 0)
#define PAD_UINT8(x) (x)
#define PAD_UINT4(x) (((x + 3)/4) * 4)

#elif defined(HAVE_AVX512)
#define MALLOC_ALIGN(new,x) new = _mm_malloc(x,64)
#define FREE_ALIGN(x) _mm_free(x)
#define CHECK_ALIGN(x) assert((unsigned long) x % 64 == 0)
#define PAD_UINT8(x) (((x + 7)/8) * 8)
#define PAD_UINT4(x) (((x + 15)/16) * 16)

#elif defined(HAVE_AVX2)
#define MALLOC_ALIGN(new,x) new = _mm_malloc(x,32);
#define FREE_ALIGN(x) _mm_free(x)
#define CHECK_ALIGN(x) assert((unsigned long) x % 32 == 0)
#define PAD_UINT8(x) (((x + 3)/4) * 4)
#define PAD_UINT4(x) (((x + 7)/8) * 8)

#elif defined(HAVE_SSE4_1)
#define MALLOC_ALIGN(new,x) new = _mm_malloc(x,16);
#ifdef DEBUG_ALIGN
#define FREE_ALIGN(x) printf("Freeing %p in standard pool at %s:%d\n",x,__FILE__,__LINE__); _mm_free(x)
#else
#define FREE_ALIGN(x) _mm_free(x)
#endif
#define CHECK_ALIGN(x) assert((unsigned long) x % 16 == 0)
#define PAD_UINT8(x) (x)
#define PAD_UINT4(x) (((x + 3)/4) * 4)

#else
#define MALLOC_ALIGN(new,x) new = MALLOC(x)
#define FREE_ALIGN(x) FREE(x)
#define CHECK_ALIGN(x)
#define PAD_UINT8(x) (x)
#define PAD_UINT4(x) (x)

#endif


#ifdef MEMUSAGE
#define MALLOC_KEEP(nbytes) Mem_alloc_keep((nbytes), __FILE__, __LINE__)
#define MALLOC_IN(nbytes) Mem_alloc_in((nbytes), __FILE__, __LINE__)
#define MALLOC_OUT(nbytes) Mem_alloc_out((nbytes), __FILE__, __LINE__)

#define CALLOC_KEEP(count, nbytes) Mem_calloc_keep((count), (nbytes), __FILE__, __LINE__)
#define CALLOC_IN(count, nbytes) Mem_calloc_in((count), (nbytes), __FILE__, __LINE__)
#define CALLOC_OUT(count, nbytes) Mem_calloc_out((count), (nbytes), __FILE__, __LINE__)

#define FREE_KEEP(ptr) ((void)(Mem_free_keep((ptr),__FILE__, __LINE__), (ptr) = 0))
#define FREE_IN(ptr) ((void)(Mem_free_in((ptr),__FILE__, __LINE__), (ptr) = 0))
#define FREE_OUT(ptr) ((void)(Mem_free_out((ptr),__FILE__, __LINE__), (ptr) = 0))

#else

#define MALLOC_KEEP(nbytes) Mem_alloc((nbytes), __FILE__, __LINE__)
#define MALLOC_IN(nbytes) Mem_alloc((nbytes), __FILE__, __LINE__)
#define MALLOC_OUT(nbytes) Mem_alloc((nbytes), __FILE__, __LINE__)

#define CALLOC_KEEP(count, nbytes) Mem_calloc((count), (nbytes), __FILE__, __LINE__)
#define CALLOC_IN(count, nbytes) Mem_calloc((count), (nbytes), __FILE__, __LINE__)
#define CALLOC_OUT(count, nbytes) Mem_calloc((count), (nbytes), __FILE__, __LINE__)

#define FREE_KEEP(ptr) ((void)(Mem_free((ptr),__FILE__, __LINE__), (ptr) = 0))
#define FREE_IN(ptr) ((void)(Mem_free((ptr),__FILE__, __LINE__), (ptr) = 0))
#define FREE_OUT(ptr) ((void)(Mem_free((ptr),__FILE__, __LINE__), (ptr) = 0))

#endif


#define MALLOC_NO_EXCEPTION(nbytes) \
	Mem_alloc_no_exception((nbytes), __FILE__, __LINE__)
#define CALLOC_NO_EXCEPTION(count, nbytes) \
	Mem_calloc((count), (nbytes), __FILE__, __LINE__)
#define  NEW(p) ((p) = MALLOC(sizeof *(p)))
#define NEW0(p) ((p) = CALLOC(1, sizeof *(p)))
#define RESIZE(ptr, nbytes) 	((ptr) = Mem_resize((ptr), \
	(nbytes), __FILE__, __LINE__))
#endif
