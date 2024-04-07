static char rcsid[] = "$Id: list.c 225763 2023-01-26 18:05:13Z twu $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "list.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mem.h"

#define T List_T

#if !defined(HAVE_INLINE)
T
List_push (T list, void *x) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->first = x;
  new->rest = list;
  return new;
}


T
List_push_keep (T list, void *x) {
  T new = (T) MALLOC_KEEP(sizeof(*new));
  
  new->first = x;
  new->rest = list;
  return new;
}

T
List_push_out (T list, void *x) {
  T new = (T) MALLOC_OUT(sizeof(*new));
  
  new->first = x;
  new->rest = list;
  return new;
}


/* Returns new tail, which must be reassigned to tail */
T
List_unshift_in (T *head, T tail, void *x) {
  T new = (T) MALLOC_IN(sizeof(*new));

  new->first = x;
  new->rest = (T) NULL;

  if (*head == NULL) {
    *head = new;
  } else {
    tail->rest = new;
  }

  return new;
}


/* Returns new tail, which must be reassigned to tail */
T
List_unshift_out (T *head, T tail, void *x) {
  T new = (T) MALLOC_OUT(sizeof(*new));

  new->first = x;
  new->rest = (T) NULL;

  if (*head == NULL) {
    *head = new;
  } else {
    tail->rest = new;
  }

  return new;
}


T
List_pop (T list, void **x) {
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


void *
List_head (T list) {
  return list->first;
}

T
List_next (T list) {
  if (list) {
    return list->rest;
  } else {
    return NULL;
  }
}

void
List_free (T *list) {
  T prev;

  while ((prev = *list) != NULL) {
    *list = prev->rest;
    FREE(prev);
  }

  return;
}

int
List_length (T list) {
  int n;
  
  for (n = 0; list; list = list->rest) {
    n++;
  }

  return n;
}

T
List_reverse (T list) {
  T head = NULL, next;

  for ( ; list; list = next) {
    next = list->rest;
    list->rest = head;
    head = list;
  }
  return head;
}

T
List_transfer_one (T dest, T *source) {
  T next;

  next = (*source)->rest;
  (*source)->rest = dest;
  dest = *source;
  *source = next;
  return dest;
}
#endif


T
List_pop_in (T list, void **x) {
  T head;

  if (list) {
    head = list->rest;
    *x = list->first;
    FREE_IN(list);
    return head;
  } else {
    return list;
  }
}

T
List_pop_out (T list, void **x) {
  T head;

  if (list) {
    head = list->rest;
    *x = list->first;
    FREE_OUT(list);
    return head;
  } else {
    return list;
  }
}

void
List_head_set (T this, void *x) {
  this->first = x;
  return;
}

void
List_tail_set (T this, T rest) {
  this->rest = rest;
  return;
}

void
List_free_keep (T *list) {
  T prev;

  while ((prev = *list) != NULL) {
    *list = (*list)->rest;
    FREE_KEEP(prev);
  }
}

void
List_free_out (T *list) {
  T prev;

  while ((prev = *list) != NULL) {
    *list = (*list)->rest;
    FREE_OUT(prev);
  }
}

T
List_truncate (T list, int n) {
  T head = list;

  while (--n > 0) {
    list = list->rest;
  }
  if (list) {
    list->rest = (T) NULL;
  }
  return head;
}

void **
List_to_array (T list, void *end) {
  void **array;
  int i, n = List_length(list);

#if 0
  if (n == 0) {
    return (void *) NULL;
  } else {
#endif
    array = (void **) MALLOC((n+1)*sizeof(*array));
    for (i = 0; i < n; i++) {
      array[i] = list->first;
      list = list->rest;
    }
    array[i] = end;
    return array;
#if 0
  }
#endif
}

void
List_fill_array (void **array, T list) {
  int i = 0;

  while (list) {
    array[i++] = list->first;
    list = list->rest;
  }
  return;
}

void
List_fill_array_and_free (void **array, T *list) {
  T prev;
  int i = 0;

  while ((prev = *list) != NULL) {
    array[i++] = prev->first;
    *list = prev->rest;
    FREE(prev);
  }

  return;
}

T
List_fill_array_with_handle (struct T *new, void **array, int nelts) {
  T list = NULL;
  int i;

  for (i = nelts; i > 0; i--) {
    new[i].first = array[i-1];
    new[i].rest = list;
    list = &(new[i]);
  }

  /* Add initial list element as a handle */
  new[0].first = (void *) NULL;
  new[0].rest = list;

  return &(new[0]);
}
    

void **
List_to_array_out (T list, void *end) {
  void **array;
  int i, n = List_length(list);

#if 0
  if (n == 0) {
    return (void *) NULL;
  } else {
#endif
    array = (void **) MALLOC_OUT((n+1)*sizeof(*array));
    for (i = 0; i < n; i++) {
      array[i] = list->first;
      list = list->rest;
    }
    array[i] = end;
    return array;
#if 0
  }
#endif
}

void **
List_to_array_n (int *n, T list) {
  void **array;
  int i;

  *n = List_length(list);
  if (*n == 0) {
    return NULL;
  } else {
    array = (void **) MALLOC((*n)*sizeof(*array));
    for (i = 0; i < *n; i++) {
      array[i] = list->first;
      list = list->rest;
    }
    return array;
  }
}

void **
List_to_array_out_n (int *n, T list) {
  void **array;
  int i;

  *n = List_length(list);
  if (*n == 0) {
    return NULL;
  } else {
    array = (void **) MALLOC_OUT((*n)*sizeof(*array));
    for (i = 0; i < *n; i++) {
      array[i] = list->first;
      list = list->rest;
    }
    return array;
  }
}

T
List_copy (T list) {
  T head, *p = &head;

  for ( ; list; list = list->rest) {
    *p = (T) MALLOC(sizeof(**p));
    (*p)->first = list->first;
    p = &(*p)->rest;
  }
  *p = NULL;
  return head;
}

void
List_dump (T list) {
  while (list) {
    printf("%p\n",list);
    list = list->rest;
  }
  return;
}

T
List_append (T list, T tail) {
  T *p = &list;

  while (*p) {
    p = &(*p)->rest;
  }
  *p = tail;
  return list;
}

T
List_insert_end (T list, void *x) {
  T tail = (T) MALLOC(sizeof(*tail));
  T *p = &list;

  tail->first = x;
  tail->rest = NULL;
  
  while (*p) {
    p = &(*p)->rest;
  }
  *p = tail;

  return list;
}


T
List_drop_last (T list, void **x) {
  T *p = &list;

  while ((*p)->rest) {
    p = &(*p)->rest;
  }
  *x = (*p)->first;
  /* FREE(*p); */
  *p = NULL;

  return list;
}

void
List_last_set (T list, void *x) {
  T *p = &list;

  while ((*p)->rest) {
    p = &(*p)->rest;
  }
  (*p)->first = x;

  return;
}


void *
List_last_value (T this, T end) {
  T last = NULL, r;

  for (r = this; r != end; r = r->rest) {
    last = r;
  }
  return last->first;
}

T
List_last_item (T this, T end) {
  T last = NULL, r;

  for (r = this; r != end; r = r->rest) {
    last = r;
  }
  return last;
}

void *
List_index (T this, int index) {
  while (index-- > 0) {
    this = this->rest;
  }
  return this->first;
}


#if 1
void
List_insert (T *listptr, void *x) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->first = x;
  new->rest = *listptr;
  *listptr = new;

  return;
}
#else
/* Doesn't work if p is NULL */
T
List_insert (T p, void *x) {
  T new = (T) MALLOC(sizeof(*new));
  
  new->first = x;
  new->rest = p->rest;
  p->rest = new;
  return new;
}
#endif

void
List_reinsert (T *listptr, T cell) {
  cell->rest = *listptr;
  *listptr = cell;

  return;
}


void
List_delete (T this, T prev) {
  prev->rest = this->rest;
  return;
}


/* Use by doing

    newsource = List_next(source);
    dest = List_push_existing(dest,source);
    source = newsource;

    or dest = List_transfer_one(dest,&source);
*/
T
List_push_existing (T dest, T source) {
  source->rest = dest;
  return source;
}

T
List_from_string (char *string) {
  T this = NULL;
  char *p = string, *scout = string, *substring;
  int substringlen;

  while (*++scout != '\0') {
    if (*scout == ',') {
      substringlen = (scout-p)/sizeof(char);
      substring = (char *) MALLOC((substringlen+1)*sizeof(char));
      strncpy(substring,p,substringlen);
      substring[substringlen] = '\0';
      this = List_push(this,substring);
      scout++;
      p = scout;
    }
  }

  substringlen = (scout-p)/sizeof(char);
  substring = (char *) MALLOC((substringlen+1)*sizeof(char));
  strncpy(substring,p,substringlen);
  substring[substringlen] = '\0';
  this = List_push(this,substring);

  return List_reverse(this);
}
