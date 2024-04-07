static char rcsid[] = "$Id: d8367f5685dc09e31593a933835466e84c86ad5e $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "intlist.h"
#include <stdio.h>		/* For sprintf */
#include <stdlib.h>
#include <string.h>		/* For strlen */
#include "mem.h"


#define T Intlist_T

#if !defined(HAVE_INLINE)
T
Intlist_push (T list, int x) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->first = x;
  new->rest = list;
  return new;
}

T
Intlist_pop (T list, int *x) {
  T head;

  if (list) {
    head = list->rest;
    *x = list->first;
    FREE(list);
    return head;
  } else {
    return list;
  }
}
  
int
Intlist_head (T list) {
  return list->first;
}

void
Intlist_head_set (T list, int x) {
  list->first = x;
  return;
}

void
Intlist_head_incr (T list, int x) {
  list->first += x;
  return;
}

int
Intlist_second_value (T list) {
  return list->rest->first;
}

T
Intlist_next (T list) {
  if (list) {
    return list->rest;
  } else {
    return NULL;
  }
}

void
Intlist_free (T *list) {
  T prev;

  while ((prev = *list) != NULL) {
    *list = prev->rest;
    FREE(prev);
  }
}

T
Intlist_reverse (T list) {
  T head = NULL, next;

  for ( ; list; list = next) {
    next = list->rest;
    list->rest = head;
    head = list;
  }
  return head;
}

int
Intlist_length (T list) {
  int n;
  
  for (n = 0; list; list = list->rest) {
    n++;
  }
  return n;
}


#endif	/* inline */


T
Intlist_push_in (T list, int x) {
  T new = (T) MALLOC_IN(sizeof(*new));
  
  new->first = x;
  new->rest = list;
  return new;
}

void
Intlist_delete (T prev, T this) {
  prev->rest = this->rest;
  FREE(this);
  return;
}


void
Intlist_free_in (T *list) {
  T prev;

  while ((prev = *list) != NULL) {
    *list = prev->rest;
    FREE_IN(prev);
  }
}

T
Intlist_keep_one (T list, int i) {
  T head;

  while (--i >= 0) {
    /* Pop */
    head = list->rest;
    FREE(list);
    list = head;
  }

  Intlist_free(&list->rest);
  return list;
}



int
Intlist_max (T list) {
  int m = 0;

  while (list) {
    if (list->first > m) {
      m = list->first;
    }
    list = list->rest;
  }

  return m;
}

int
Intlist_min (T list) {
  int m;

  if (list == NULL) {
    return 0;

  } else {
    m = list->first;
    list = list->rest;
    while (list) {
      if (list->first < m) {
	m = list->first;
      }
      list = list->rest;
    }

    return m;
  }
}


int
Intlist_sum (T list) {
  int sum = 0;

  while (list) {
    sum += list->first;
    list = list->rest;
  }

  return sum;
}


bool
Intlist_vary (T list) {
  int m;

  if (list == NULL) {
    return false;

  } else {
    m = list->first;
    list = list->rest;
    while (list) {
      if (list->first != m) {
	return true;
      }
      list = list->rest;
    }

    return false;
  }
}

bool
Intlist_exists_p (T list, int x) {
  while (list) {
    if (list->first == x) {
      return true;
    }
    list = list->rest;
  }
  return false;
}

bool
Intlist_first_equals_p (T list, int x) {
  if (list == NULL) {
    return false;
  } else if (list->first == x) {
    return true;
  } else {
    return false;
  }
}

int *
Intlist_to_array (int *n, T list) {
  int *array;
  int i;

  *n = Intlist_length(list);
  if (*n == 0) {
    return NULL;
  } else {
    array = (int *) CALLOC(*n,sizeof(int));
    for (i = 0; i < *n; i++) {
      array[i] = list->first;
      list = list->rest;
    }
    return array;
  }
}

int *
Intlist_to_array_plus_extra (int *n, T list) {
  int *array;
  int i;

  *n = Intlist_length(list);
  if (*n == 0) {
    return NULL;
  } else {
    array = (int *) CALLOC((*n)+1,sizeof(int));
    for (i = 0; i < *n; i++) {
      array[i] = list->first;
      list = list->rest;
    }
    return array;
  }
}

void
Intlist_fill_array (int *array, T list) {
  int i = 0;

  while (list) {
    array[i++] = list->first;
    list = list->rest;
  }

  return;
}


void
Intlist_fill_array_and_free (int *array, T *list) {
  T prev;
  int i = 0;

  while ((prev = *list) != NULL) {
    array[i++] = prev->first;
    *list = prev->rest;
    FREE(prev);
  }

  return;
}

int *
Intlist_to_array_out (int *n, T list) {
  int *array;
  int i;

  *n = Intlist_length(list);
  if (*n == 0) {
    return NULL;
  } else {
    array = (int *) CALLOC_OUT(*n,sizeof(int));
    for (i = 0; i < *n; i++) {
      array[i] = list->first;
      list = list->rest;
    }
    return array;
  }
}


char *
Intlist_to_char_array (int *n, T list) {
  char *array;
  int i;

  *n = Intlist_length(list);
  if (*n == 0) {
    return NULL;
  } else {
    array = (char *) MALLOC((*n + 1)*sizeof(char));
    for (i = 0; i < *n; i++) {
      array[i] = (char) list->first;
      list = list->rest;
    }
    array[*n] = '\0';
    return array;
  }
}

char *
Intlist_to_char_array_in (int *n, T list) {
  char *array;
  int i;

  *n = Intlist_length(list);
  if (*n == 0) {
    return NULL;
  } else {
    array = (char *) MALLOC_IN((*n + 1)*sizeof(char));
    for (i = 0; i < *n; i++) {
      array[i] = (char) list->first;
      list = list->rest;
    }
    array[*n] = '\0';
    return array;
  }
}

T
Intlist_from_array (int *array, int n) {
  T list = NULL, p;

  while (--n >= 0) {
    p = (T) MALLOC(sizeof(*p));
    p->first = array[n];
    p->rest = list;
    list = p;
  }

  return list;
}

T
Intlist_copy (T list) {
  T head, *p = &head;

  for ( ; list; list = list->rest) {
    *p = (T) MALLOC(sizeof(**p));
    (*p)->first = list->first;
    p = &(*p)->rest;
  }
  *p = NULL;
  return head;
}

T
Intlist_append (T list, T tail) {
  T *p = &list;

  while (*p) {
    p = &(*p)->rest;
  }
  *p = tail;
  return list;
}

int
Intlist_last_value (T this) {
  T last = NULL, r;

  for (r = this; r != NULL; r = r->rest) {
    last = r;
  }
  return last->first;
}

int
Intlist_penultimate_value (T this) {
  T last = NULL, r;

  for (r = this; r->rest != NULL; r = r->rest) {
    last = r;
  }
  return last->first;
}

int
Intlist_index (T this, int index) {
  while (index-- > 0) {
    this = this->rest;
  }
  return this->first;
}


#if 0
/* A specialized version is included in gsnap-counts.c */
T
Intlist_from_string (char *string) {
  T this = NULL;
  char *p = string;
  int x;

  while (sscanf(p,"%d",&x) > 0) {
    this = Intlist_push(this,x);
    while (*p != '\0' && *p != ',') {
      p++;
    }
    if (*p == ',') {
      p++;
    }
  }
  return Intlist_reverse(this);
}
#endif


char *
Intlist_to_string (T this) {
  char *string, Buffer[256];
  T p;
  int n, i, strlength;

  if ((n = Intlist_length(this)) == 0) {
    string = (char *) CALLOC(1,sizeof(char));
    string[0] = '\0';
  } else {
    strlength = 0;
    for (i = 0, p = this; i < n-1; i++, p = Intlist_next(p)) {
      sprintf(Buffer,"%d,",Intlist_head(p));
      strlength += strlen(Buffer);
    }
    sprintf(Buffer,"%d",Intlist_head(p));
    strlength += strlen(Buffer);

    string = (char *) CALLOC(strlength + 1,sizeof(char));
    string[0] = '\0';
    for (i = 0, p = this; i < n-1; i++, p = Intlist_next(p)) {
      sprintf(Buffer,"%d,",Intlist_head(p));
      strcat(string,Buffer);
    }
    sprintf(Buffer,"%d",Intlist_head(p));
    strcat(string,Buffer);
  }

  return string;
}


struct Cell_T {
  int elt;
  int keyvalue;
};

static int
cell_ascending (const void *a, const void *b) {
  struct Cell_T *x = (struct Cell_T *) a;
  struct Cell_T *y = (struct Cell_T *) b;

  if (x->keyvalue < y->keyvalue) {
    return -1;
  } else if (y->keyvalue < x->keyvalue) {
    return +1;
  } else {
    return 0;
  }
}

static int
cell_ascending_dual (const void *a, const void *b) {
  struct Cell_T *x = (struct Cell_T *) a;
  struct Cell_T *y = (struct Cell_T *) b;

  if (x->keyvalue < y->keyvalue) {
    return -1;
  } else if (y->keyvalue < x->keyvalue) {
    return +1;
  } else if (x->elt < y->elt) {
    return -1;
  } else if (y->elt < x->elt) {
    return +1;
  } else {
    return 0;
  }
}

static int
cell_descending (const void *a, const void *b) {
  struct Cell_T *x = (struct Cell_T *) a;
  struct Cell_T *y = (struct Cell_T *) b;

  if (x->keyvalue > y->keyvalue) {
    return -1;
  } else if (y->keyvalue > x->keyvalue) {
    return +1;
  } else {
    return 0;
  }
}

int *
Intlist_array_ascending_by_key (int *n, T this, T key) {
  int *sorted, i;
  struct Cell_T *cells;
  T p, q;

  /* find length */
  for (*n = 0, p = this; p; p = p->rest) {
    (*n)++;
  }

  cells = (struct Cell_T *) CALLOC(*n,sizeof(struct Cell_T));
  for (p = this, q = key, i = 0; p != NULL; p = p->rest, q = q->rest, i++) {
    cells[i].elt = p->first;
    cells[i].keyvalue = q->first;
  }
  qsort(cells,*n,sizeof(struct Cell_T),cell_ascending);

  sorted = (int *) CALLOC(*n,sizeof(int));
  for (i = 0; i < *n; i++) {
    sorted[i] = cells[i].elt;
  }

  FREE(cells);
  return sorted;
}


void
Intlist_array_dual_ascending_by_key (int *sorted, int *keyarray, int n, T this, T key) {
  int i;
  struct Cell_T *cells;
  T p, q;

  cells = (struct Cell_T *) MALLOCA(n * sizeof(struct Cell_T));
  for (p = this, q = key, i = 0; p != NULL; p = p->rest, q = q->rest, i++) {
    cells[i].elt = p->first;
    cells[i].keyvalue = q->first;
  }
  qsort(cells,n,sizeof(struct Cell_T),cell_ascending_dual);

  for (i = 0; i < n; i++) {
    sorted[i] = cells[i].elt;
    keyarray[i] = cells[i].keyvalue;
  }

  FREEA(cells);
  return;
}


T
Intlist_list_ascending_by_key (T this, T key) {
  T sorted = NULL, p, q;
  int n, i;
  struct Cell_T *cells;

  /* find length */
  for (n = 0, p = this; p; p = p->rest) {
    n++;
  }

  cells = (struct Cell_T *) CALLOC(n,sizeof(struct Cell_T));
  for (p = this, q = key, i = 0; p != NULL; p = p->rest, q = q->rest, i++) {
    cells[i].elt = p->first;
    cells[i].keyvalue = q->first;
  }
  qsort(cells,n,sizeof(struct Cell_T),cell_descending);

  for (i = 0; i < n; i++) {
    sorted = Intlist_push(sorted,cells[i].elt);
  }
  /* No need to reverse list because we sorted by descending key */

  FREE(cells);

  return sorted;
}


T
Intlist_list_descending_by_key (T this, T key) {
  T sorted = NULL, p, q;
  int n, i;
  struct Cell_T *cells;

  /* find length */
  for (n = 0, p = this; p; p = p->rest) {
    n++;
  }

  cells = (struct Cell_T *) CALLOC(n,sizeof(struct Cell_T));
  for (p = this, q = key, i = 0; p != NULL; p = p->rest, q = q->rest, i++) {
    cells[i].elt = p->first;
    cells[i].keyvalue = q->first;
  }
  qsort(cells,n,sizeof(struct Cell_T),cell_ascending);

  for (i = 0; i < n; i++) {
    sorted = Intlist_push(sorted,cells[i].elt);
  }
  /* No need to reverse list because we sorted by ascending key */

  FREE(cells);

  return sorted;
}


T
Intlist_sort_ascending (T this) {
  T sorted = NULL, p;
  int n, i;
  struct Cell_T *cells;

  /* find length */
  for (n = 0, p = this; p; p = p->rest) {
    n++;
  }

  cells = (struct Cell_T *) CALLOC(n,sizeof(struct Cell_T));
  for (p = this, i = 0; p != NULL; p = p->rest, i++) {
    cells[i].keyvalue = p->first;
  }
  qsort(cells,n,sizeof(struct Cell_T),cell_descending);

  for (i = 0; i < n; i++) {
    sorted = Intlist_push(sorted,cells[i].keyvalue);
  }
  /* No need to reverse list because we sorted by descending key */

  FREE(cells);

  return sorted;
}

bool
Intlist_equal (T x, T y) {

  while (x && y && x->first == y->first) {
    x = x->rest;
    y = y->rest;
  }
  if (x || y) {
    return false;
  } else {
    return true;
  }
}

bool
Intlist_intersect_p (T x, T y) {
  T q;
  int value;

  while (x) {
    value = x->first;
    for (q = y; q != NULL; q = q->rest) {
      if (q->first == value) {
	return true;
      }
    }
    x = x->rest;
  }

  return false;
}

#if 0
/* Reverses elements in x */
T
Intlist_intersect (T x, T y) {
  T list = NULL, new, q;
  int value;

  while (x) {
    value = x->first;
    for (q = y; q != NULL; q = q->rest) {
      if (q->first == value) {
	/* list = Intlist_push(list,value) */
	new = (T) MALLOC(sizeof(*new));
	new->first = value;
	new->rest = list;
	list = new;
      }
    }
    x = x->rest;
  }

  return list;
}
#endif


T
Intlist_filter (T list, bool *deletep) {
  T newlist = NULL, next;
  int i = 0;

  while (list) {
    next = list->rest;
    if (deletep[i++] == false) {
      /* Push cell onto newlist */
      list->rest = newlist;
      newlist = list;
    } else {
      FREE(list);
    }
    list = next;
  }

  /* return Intlist_reverse(newlist); */
  return newlist;
}
