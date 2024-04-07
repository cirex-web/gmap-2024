static char rcsid[] = "$Id$";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "auxinfo.h"

#include <stdio.h>
#include "assert.h"
#include "mem.h"
#include "path-eval.h"		/* For Path_local_cmp */

#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* Auxinfo_assign_chrinfo */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

static EF64_T chromosome_ef64;

#define T Auxinfo_T


/* Modified from List_free */
void
Auxinfo_free_wpaths (T *old, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool,
		     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		     Listpool_T listpool, Pathpool_T pathpool,
		     Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  T prev;

  while ((prev = *old) != NULL) {
    *old = prev->rest;

    Intlistpool_free_list(&prev->best_sense_partners,intlistpool
			  intlistpool_trace(__FILE__,__LINE__));
    Intlistpool_free_list(&prev->best_antisense_partners,intlistpool
			  intlistpool_trace(__FILE__,__LINE__));

    Hitlistpool_free_list(&prev->best_sense_paths,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));
    Hitlistpool_free_list(&prev->best_antisense_paths,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));

    Path_gc(&prev->unextended_sense_paths,intlistpool,univcoordlistpool,
	    listpool,pathpool,transcriptpool,hitlistpool);
    Path_gc(&prev->unextended_antisense_paths,intlistpool,univcoordlistpool,
	    listpool,pathpool,transcriptpool,hitlistpool);

    Path_gc(&prev->complete_sense_paths,intlistpool,univcoordlistpool,
	    listpool,pathpool,transcriptpool,hitlistpool);
    Path_gc(&prev->complete_antisense_paths,intlistpool,univcoordlistpool,
	    listpool,pathpool,transcriptpool,hitlistpool);

    if (prev->right_univdiags != NULL) {
      Univdiagpool_gc(&prev->right_univdiags,univdiagpool
		      univdiagpool_trace(__FILE__,__LINE__));
    }
    if (prev->left_univdiags != NULL) {
      Univdiagpool_gc(&prev->left_univdiags,univdiagpool
		      univdiagpool_trace(__FILE__,__LINE__));
    }
      
    Auxinfopool_free_auxinfo(&prev,auxinfopool
			     auxinfopool_trace(__FILE__,__LINE__)); /* Allocated by Pathpool_new_path */
  }

  return;
}


void
Auxinfo_gc_wpaths (T *array, int n, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool,
		   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		   Listpool_T listpool, Pathpool_T pathpool,
		   Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  int i;
  T this;

  for (i = 0; i < n; i++) {
    this = array[i];
    Auxinfo_free_wpaths(&this,univdiagpool,auxinfopool,
			intlistpool,univcoordlistpool,
			listpool,pathpool,transcriptpool,hitlistpool);
  }

  FREE(array);
  return;
}  


void
Auxinfo_free (T *old, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool,
	      Intlistpool_T intlistpool, Hitlistpool_T hitlistpool) {
  T prev;

  while ((prev = *old) != NULL) {
    *old = prev->rest;

    Intlistpool_free_list(&prev->best_sense_partners,intlistpool
			  intlistpool_trace(__FILE__,__LINE__));
    Intlistpool_free_list(&prev->best_antisense_partners,intlistpool
			  intlistpool_trace(__FILE__,__LINE__));

    Hitlistpool_free_list(&prev->best_sense_paths,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));
    Hitlistpool_free_list(&prev->best_antisense_paths,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));

    Hitlistpool_free_list(&prev->unextended_sense_paths,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));
    Hitlistpool_free_list(&prev->unextended_antisense_paths,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));

    Hitlistpool_free_list(&prev->complete_sense_paths,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));
    Hitlistpool_free_list(&prev->complete_antisense_paths,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));

    if (prev->right_univdiags != NULL) {
      Univdiagpool_gc(&prev->right_univdiags,univdiagpool
		      univdiagpool_trace(__FILE__,__LINE__));
    }
    if (prev->left_univdiags != NULL) {
      Univdiagpool_gc(&prev->left_univdiags,univdiagpool
		      univdiagpool_trace(__FILE__,__LINE__));
    }
      
    Auxinfopool_free_auxinfo(&prev,auxinfopool
			     auxinfopool_trace(__FILE__,__LINE__)); /* Allocated by Pathpool_new_path */
  }

  return;
}


void
Auxinfo_gc (T *array, int n, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool,
	    Intlistpool_T intlistpool, Hitlistpool_T hitlistpool) {
  int i;
  T this;

  for (i = 0; i < n; i++) {
    this = array[i];
    Auxinfo_free(&this,univdiagpool,auxinfopool,intlistpool,hitlistpool);
  }

  FREE(array);
  return;
}  


static List_T
best_paths (List_T paths, Hitlistpool_T hitlistpool) {
  List_T best;

  int npaths, i, j, k;
  Path_T *patharray, path0, path;

  if ((npaths = List_length(paths)) == 0) {
    return (List_T) NULL;

  } else if (npaths == 1) {
    path = (Path_T) List_head(paths);
    return Hitlist_push(NULL,hitlistpool,(void *) path
			hitlistpool_trace(__FILE__,__LINE__));

  } else {
    patharray = (Path_T *) List_to_array(paths,NULL);

    /* Eliminate duplicates */
    qsort(patharray,npaths,sizeof(Path_T),Path_structure_cmp);

    k = 0;
    i = 0;
    while (i < npaths) {
      j = i + 1;
      while (j < npaths && Path_structure_cmp(&(patharray[j]),&(patharray[i])) == 0) {
	j++;
      }

      /* Found an identical group by structure.  Keep just the first one */
      patharray[k++] = patharray[i];
	
      i = j;
    }
    npaths = k;

    /* Find best ones */
    qsort(patharray,npaths,sizeof(Path_T),Path_cmp);

    path0 = patharray[0];
    debug8(printf("Best: ")); debug8(Path_print(path0));
    best = Hitlist_push(NULL,hitlistpool,(void *) path0
			hitlistpool_trace(__FILE__,__LINE__));
    for (i = 1; i < npaths; i++) {
      path = patharray[i];
      if (Path_cmp(&path,&path0) <= 0) {
	debug8(printf("Tie: ")); debug8(Path_print(path));
	best = Hitlist_push(best,hitlistpool,(void *) path
			    hitlistpool_trace(__FILE__,__LINE__));
      } else {
	debug8(printf("Not best: ")); debug8(Path_print(path));
      }
    }
    debug8(printf("\n"));

    FREE(patharray);
    return List_reverse(best);
  }
}


void
Auxinfo_set_best_paths (T this, Hitlistpool_T hitlistpool) {
  if (this->complete_sense_paths != NULL) {
    this->complete_sense_p = true;
    this->best_sense_paths = best_paths(this->complete_sense_paths,hitlistpool);
  } else {
    this->best_sense_paths = best_paths(this->unextended_sense_paths,hitlistpool);
  }
  if (this->complete_antisense_paths != NULL) {
    this->complete_antisense_p = true;
    this->best_antisense_paths = best_paths(this->complete_antisense_paths,hitlistpool);
  } else {
    this->best_antisense_paths = best_paths(this->unextended_antisense_paths,hitlistpool);
  }
  return;
}


void
Auxinfo_set_best_sense_paths (T this, Hitlistpool_T hitlistpool, bool only_complete_p) {

  if (only_complete_p == true) {
    this->complete_sense_p = true;
    this->best_sense_paths = best_paths(this->complete_sense_paths,hitlistpool);
  } else if (this->complete_sense_paths != NULL) {
    this->complete_sense_p = true;
    this->best_sense_paths = best_paths(this->complete_sense_paths,hitlistpool);
  } else {
    this->best_sense_paths = best_paths(this->unextended_sense_paths,hitlistpool);
  }

  return;
}

void
Auxinfo_set_best_antisense_paths (T this, Hitlistpool_T hitlistpool, bool only_complete_p) {

  if (only_complete_p == true) {
    this->complete_antisense_p = true;
    this->best_antisense_paths = best_paths(this->complete_antisense_paths,hitlistpool);
  } else if (this->complete_antisense_paths != NULL) {
    this->complete_antisense_p = true;
    this->best_antisense_paths = best_paths(this->complete_antisense_paths,hitlistpool);
  } else {
    this->best_antisense_paths = best_paths(this->unextended_antisense_paths,hitlistpool);
  }

  return;
}



/* Standard version */
T
Auxinfo_new (Method_T method, int qstart, int qend, Auxinfopool_T auxinfopool) {
  T new = Auxinfopool_new_auxinfo(auxinfopool
				  auxinfopool_trace(__FILE__,__LINE__));
  new->chrnum = 0;

  new->solvedp = false;
  new->complete_sense_p = false;
  new->complete_antisense_p = false;

  new->best_sense_partners = (Intlist_T) NULL;
  new->best_antisense_partners = (Intlist_T) NULL;

  new->best_sense_paths = (List_T) NULL;
  new->best_antisense_paths = (List_T) NULL;

  new->unextended_sense_paths = (List_T) NULL;
  new->unextended_antisense_paths = (List_T) NULL;
  new->complete_sense_paths = (List_T) NULL;
  new->complete_antisense_paths = (List_T) NULL;

  new->method = method;
  new->qstart = qstart;
  new->qend = qend;
  new->nmismatches = -1;

  new->right_univdiags = (List_T) NULL;
  new->left_univdiags = (List_T) NULL;

  new->rest = (T) NULL;

  return new;
}


T
Auxinfo_new_tr (Auxinfopool_T auxinfopool, Path_T path) {
  T new = Auxinfopool_new_auxinfo(auxinfopool
				  auxinfopool_trace(__FILE__,__LINE__));

  new->chrnum = path->chrnum;
  new->chroffset = path->chroffset;
  new->chrhigh = path->chrhigh;

  new->solvedp = false;
  new->complete_sense_p = false;
  new->complete_antisense_p = false;

  new->best_sense_partners = (Intlist_T) NULL;
  new->best_antisense_partners = (Intlist_T) NULL;

  new->best_sense_paths = (List_T) NULL;
  new->best_antisense_paths = (List_T) NULL;

  new->unextended_sense_paths = (List_T) NULL;
  new->unextended_antisense_paths = (List_T) NULL;
  new->complete_sense_paths = (List_T) NULL;
  new->complete_antisense_paths = (List_T) NULL;

  new->method = TR_METHOD;	/* Some tr method, e.g., TR_EXACT1 */
  new->qstart = -1;
  new->qend = -1;
  new->nmismatches = -1;

  new->right_univdiags = (List_T) NULL;
  new->left_univdiags = (List_T) NULL;

  new->rest = (T) NULL;

  return new;
}

T
Auxinfo_new_univdiags (Method_T method, int qstart, int qend, int nmismatches,
		       List_T right_univdiags, List_T left_univdiags, Auxinfopool_T auxinfopool) {
  T new = Auxinfopool_new_auxinfo(auxinfopool
				  auxinfopool_trace(__FILE__,__LINE__));
  new->chrnum = 0;

  new->solvedp = false;
  new->complete_sense_p = false;
  new->complete_antisense_p = false;

  new->best_sense_partners = (Intlist_T) NULL;
  new->best_antisense_partners = (Intlist_T) NULL;

  new->best_sense_paths = (List_T) NULL;
  new->best_antisense_paths = (List_T) NULL;

  new->unextended_sense_paths = (List_T) NULL;
  new->unextended_antisense_paths = (List_T) NULL;
  new->complete_sense_paths = (List_T) NULL;
  new->complete_antisense_paths = (List_T) NULL;

  new->method = method;
  new->qstart = qstart;
  new->qend = qend;
  new->nmismatches = nmismatches;

  new->right_univdiags = right_univdiags;
  new->left_univdiags = left_univdiags;

  new->rest = (T) NULL;

  return new;
}


T
Auxinfo_new_sense_path (Path_T path, Hitlistpool_T hitlistpool, Auxinfopool_T auxinfopool) {
  T new = Auxinfopool_new_auxinfo(auxinfopool
				  auxinfopool_trace(__FILE__,__LINE__));
  new->chrnum = path->chrnum;
  new->chroffset = path->chroffset;
  new->chrhigh = path->chrhigh;

  new->solvedp = false;
  new->complete_sense_p = false;
  new->complete_antisense_p = false;

  new->best_sense_partners = (Intlist_T) NULL;
  new->best_antisense_partners = (Intlist_T) NULL;

  new->best_sense_paths = Hitlist_push(NULL,hitlistpool,(void *) path
				       hitlistpool_trace(__FILE__,__LINE__));
  new->best_antisense_paths = (List_T) NULL;

  new->unextended_sense_paths = (List_T) NULL;
  new->unextended_antisense_paths = (List_T) NULL;
  new->complete_sense_paths = (List_T) NULL;
  new->complete_antisense_paths = Hitlist_push(NULL,hitlistpool,(void *) path
				      hitlistpool_trace(__FILE__,__LINE__));

  new->method = path->method;
  new->qstart = 0;
  new->qend = path->querylength;
  new->nmismatches = -1;

  new->right_univdiags = (List_T) NULL;
  new->left_univdiags = (List_T) NULL;

  new->rest = (T) NULL;

  return new;
}


T
Auxinfo_new_antisense_path (Path_T path, Hitlistpool_T hitlistpool, Auxinfopool_T auxinfopool) {
  T new = Auxinfopool_new_auxinfo(auxinfopool
				  auxinfopool_trace(__FILE__,__LINE__));
  new->chrnum = path->chrnum;
  new->chroffset = path->chroffset;
  new->chrhigh = path->chrhigh;

  new->solvedp = false;
  new->complete_sense_p = false;
  new->complete_antisense_p = false;

  new->best_sense_partners = (Intlist_T) NULL;
  new->best_antisense_partners = (Intlist_T) NULL;

  new->best_sense_paths = (List_T) NULL;
  new->best_antisense_paths = Hitlist_push(NULL,hitlistpool,(void *) path
					   hitlistpool_trace(__FILE__,__LINE__));

  new->unextended_sense_paths = (List_T) NULL;
  new->unextended_antisense_paths = (List_T) NULL;
  new->complete_sense_paths = (List_T) NULL;
  new->complete_antisense_paths = Hitlist_push(NULL,hitlistpool,(void *) path
				      hitlistpool_trace(__FILE__,__LINE__));

  new->method = path->method;
  new->qstart = 0;
  new->qend = path->querylength;
  new->nmismatches = -1;

  new->right_univdiags = (List_T) NULL;
  new->left_univdiags = (List_T) NULL;

  new->rest = (T) NULL;

  return new;
}


/* Modified from List_append */
T
Auxinfo_append (T this1, T this2) {
  T first, second, *p;

  if ((this1->complete_sense_paths != NULL || this1->complete_antisense_paths != NULL) &&
      (this2->complete_sense_paths != NULL || this2->complete_antisense_paths != NULL)) {
    first = this1; second = this2;
  } else if ((this2->complete_sense_paths != NULL || this2->complete_antisense_paths != NULL) &&
	     (this1->complete_sense_paths != NULL || this1->complete_antisense_paths != NULL)) {
    first = this2; second = this1;
  } else if (this1->method < this2->method) {
    first = this1; second = this2;
  } else if (this2->method < this1->method) {
    first = this2; second = this1;
  } else {
    fprintf(stderr,"Unexpected repeated methods on the same univdiagonal\n");
    abort();
  }

  p = &first;

  while (*p) {
    p = &(*p)->rest;
  }
  *p = second;

  return first;
}


void
Auxinfo_assign_chrinfo (Univcoord_T * univdiagonals, T *auxinfo, int n, int querylength) {
  int i, j;
  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh;
  T auxinfo_i, auxinfo_j;

  debug9(printf("Entering Auxinfo_assign_chrinfo with %d univdiagonals/auxinfo\n",n));
  i = 0;
  while (i < n && univdiagonals[i] < (Univcoord_T) querylength) {
    auxinfo_i = auxinfo[i];
    chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			 /*univdiagonals[i] - querylength*/0,univdiagonals[i]);
    auxinfo_i->chrnum = chrnum;
    auxinfo_i->chroffset = chroffset;
    auxinfo_i->chrhigh = chrhigh;

    debug9(printf("i* %d %u %d %u..%u\n",i,univdiagonals[i],chrnum,chroffset,chrhigh));
    
    j = i + 1;
    while (j < n && univdiagonals[j] < chrhigh) {
      auxinfo_j = auxinfo[j];
      auxinfo_j->chrnum = chrnum;
      auxinfo_j->chroffset = chroffset;
      auxinfo_j->chrhigh = chrhigh;
      debug9(printf("j* %d %u %d %u..%u\n",j,univdiagonals[j],chrnum,chroffset,chrhigh));

      j++;
    }
    
    i = j;
  }

  while (i < n) {
    auxinfo_i = auxinfo[i];
    chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			 univdiagonals[i] - querylength,univdiagonals[i]);
    auxinfo_i->chrnum = chrnum;
    auxinfo_i->chroffset = chroffset;
    auxinfo_i->chrhigh = chrhigh;

    debug9(printf("i %d %u %d %u..%u\n",i,univdiagonals[i],chrnum,chroffset,chrhigh));
    
    j = i + 1;
    while (j < n && univdiagonals[j] < chrhigh) {
      auxinfo_j = auxinfo[j];
      auxinfo_j->chrnum = chrnum;
      auxinfo_j->chroffset = chroffset;
      auxinfo_j->chrhigh = chrhigh;
      debug9(printf("j %d %u %d %u..%u\n",j,univdiagonals[j],chrnum,chroffset,chrhigh));

      j++;
    }
    
    i = j;
  }

  return;
}



#if 0
void
Auxinfo_collect_unextended_paths (List_T *unextended_sense_paths,
				  List_T *unextended_antisense_paths,
				  T *auxinfo_array, int nunivdiagonals,
				  Hitlistpool_T hitlistpool) {
  T auxinfo;
  Path_T path;
  List_T p;
  int i;

  *unextended_sense_paths = *unextended_antisense_paths = (List_T) NULL;

  for (i = 0; i < nunivdiagonals; i++) {
    auxinfo = auxinfo_array[i];
    for (p = auxinfo->unextended_sense_paths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      *unextended_sense_paths =
	Hitlist_push(*unextended_sense_paths,hitlistpool,(void *) path
		     hitlistpool_trace(__FILE__,__LINE__));
    }

    for (p = auxinfo->unextended_antisense_paths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      *unextended_antisense_paths =
	Hitlist_push(*unextended_antisense_paths,hitlistpool,(void *) path
		     hitlistpool_trace(__FILE__,__LINE__));
    }
  }

  *unextended_sense_paths = List_reverse(*unextended_sense_paths);
  *unextended_antisense_paths = List_reverse(*unextended_antisense_paths);

  return;
}
#endif


void
Auxinfo_collect_paths (bool *foundp, List_T *sense_paths, List_T *antisense_paths,
		       T *auxinfo_array, int nunivdiagonals,
		       Hitlistpool_T hitlistpool) {

  T this;
  Path_T path;
  List_T p;
  int i;

  *sense_paths = *antisense_paths = (List_T) NULL;

  for (i = 0; i < nunivdiagonals; i++) {
    this = auxinfo_array[i];

    /* Collect the best ones */
    for (p = this->best_sense_paths; p != NULL; p = List_next(p)) {
      *foundp = true;
      path = (Path_T) List_head(p);
      *sense_paths = Hitlist_push(*sense_paths,hitlistpool,(void *) path
				  hitlistpool_trace(__FILE__,__LINE__));
    }
    for (p = this->best_antisense_paths; p != NULL; p = List_next(p)) {
      *foundp = true;
      path = (Path_T) List_head(p);
      *antisense_paths = Hitlist_push(*antisense_paths,hitlistpool,(void *) path
				  hitlistpool_trace(__FILE__,__LINE__));
    }
  }

  *sense_paths = List_reverse(*sense_paths);
  *antisense_paths = List_reverse(*antisense_paths);

  return;
}


void
Auxinfo_setup (EF64_T chromosome_ef64_in) {

  chromosome_ef64 = chromosome_ef64_in;

  return;
}
