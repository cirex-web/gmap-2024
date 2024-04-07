static char rcsid[] = "$Id: get-genome.c 226315 2023-02-28 18:22:20Z twu $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>		/* For rindex */
#include <ctype.h>

#include "bool.h"
#include "mem.h"
#include "access.h"
#include "getline.h"

#include "genomicpos.h"
#include "sequence.h"
#include "chrom.h"
#include "genome.h"
#include "iit-read.h"
#include "datadir.h"
#include "parserange.h"
#include "separator.h"
#include "complement.h"
#include "getopt.h"
#include "intlist.h"


#define INFINITY -1


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* Program Options */

#ifdef HAVE_MMAP
static Access_mode_T genome_access = USE_MMAP_ONLY;
#else
static Access_mode_T genome_access = USE_ALLOCATE;
#endif


#if 0
static bool replacex = false;
#endif
static int print_snps_mode = 3;	/* 1 = alt version only, 2 = snps marked, 3 = both */
static char *snps_root = NULL;

static char *user_genomedir = NULL;
static char *user_snpsdir = NULL;

static bool uppercasep = false;
static bool revcomp = false;
static bool coordp = false;
static char *genomesubdir = NULL;
static char *dbroot = NULL;
static char *dbversion = NULL;
static char *transcriptome_dbroot = NULL;

static char *user_typestring = NULL;
static int wraplength = 60;
static char *header = NULL;
static bool levelsp = false;
static bool rawp = false;

static bool genes_file_p = false;
static char *user_mapdir = NULL;
static char *map_iitfile = NULL;
static int nflanking = 0;
static bool exonsp = false;
static bool sequencep = false;
static bool genome_sequence_p = false; /* introns plus exons */
static bool uniquep = false;
static bool force_label_p = false;

static bool exactp = false;
static bool sortp = false;
static bool signedp = false;

static bool vareffect_p = false;

static bool codingp = false;

/* Dump options */
static bool dumpallp = false;
static bool add_circular_p = false;
static bool stream_chars_p = false;
static bool stream_ints_p = false;
static bool dumpchrp = false;
static bool dumpchr_forsam_p = false;
static bool dumpsegsp = false;

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'}, /* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* dbroot */
  {"transcriptdb", required_argument, 0, 'c'}, /* transcriptome_dbroot */
  {"genes", no_argument, 0, 0},	     /* genes_file_p */

  /* Output options */
  {"coords", no_argument, 0, 'C'}, /* coordp */
  {"uppercase", no_argument, 0, 'U'}, /* uppercasep */
  {"wraplength", required_argument, 0, 'l'}, /* wraplength */
  {"header", required_argument, 0, 'h'}, /* header */
  {"snpsdir", required_argument, 0, 'V'}, /* user_snpsdir */
  {"usesnps", required_argument, 0, 'v'}, /* snps_root */
  {"snpformat", required_argument, 0, 'f'}, /* print_snps */

  /* External map options */
  {"mapdir", required_argument, 0, 'M'}, /* user_mapdir */
  {"map", required_argument, 0, 'm'},	/* map_iitfile */
  {"ranks", no_argument, 0, 'k'}, /* levelsp */
  {"raw", no_argument, 0, 'r'}, /* rawp */
  {"flanking", required_argument, 0, 'u'}, /* nflanking */
  {"exons", no_argument, 0, 'E'},	   /* exonsp */
  {"sequence", no_argument, 0, 'S'},	   /* sequencep */
  {"gsequence", no_argument, 0, 0},	   /* genome_sequence_p */
  {"nunique", no_argument, 0, 0},	   /* uniquep */
  {"exact", no_argument, 0, 0},		/* exactp */
  {"signed", no_argument, 0, 's'},	   /* signedp */
  {"aslabel", no_argument, 0, 0},	   /* force_label_p */
  {"coding", no_argument, 0, 0},	   /* codingp */

  /* Dump options */
  {"dump", no_argument, 0, 'A'},	/* dumpallp */
  {"add-circular", no_argument, 0, 0}, /* add_circular_p */
  {"stream-chars", no_argument, 0, 0},	/* stream_chars_p */
  {"chromosomes", no_argument, 0, 'L'},	/* dumpchrp */
  {"forsam", no_argument, 0, 0},	/* dumpchr_forsam_p */
  {"contigs", no_argument, 0, 'I'}, /* dumpsegsp */
  
  /* Special (hidden) options */
  {"vareffect", no_argument, 0, 0}, /* vareffect_p */

  /* Help options */
  {"version", no_argument, 0, '^'}, /* print_program_version */
  {"help", no_argument, 0, '?'}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"get-genome: retrieval utility for genomic segments\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Thomas D. Wu and Colin K. Watanabe, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

static void
print_program_usage () {
  fprintf(stdout,"\
Usage: get-genome [OPTIONS...] -d genome [genome:]range, or\n\
       get-genome [OPTIONS...] -d genome chromosome:range, or\n\
       get-genome [OPTIONS...] -d genome contig[:range]\n\
where\n\
   range is startposition..endposition (endpos < startpos means - strand)\n\
         or startposition+length (+ strand)\n\
         or startposition+-length (- strand)\n\
\n\
Input options\n\
  -D, --dir=STRING        Genome directory\n\
  -d, --db=STRING         Genome database.  If argument is '?' (with\n\
                            the quotes), this command lists available databases.\n\
  -c, --transcriptdb=STRING\n\
                          Retrieve the gene structure of the given contig argument\n\
                            (Assumes transcriptome was built with gmap_build -c)\n\
  --genes                 Argument is a genes file\n\
\n\
Output options\n\
  -C, --coords            Show coordinates only\n\
  -U, --uppercase         Convert sequence to uppercase\n\
  -l, --wraplength=INT    Wrap length for sequence (default=60)\n\
  -h, --header=STRING     Desired header line\n\
  -V, --snpsdir=STRING    Directory for SNPs index files (created using snpindex)\n\
  -v, --usesnps=STRING    Use snp version (built by snpindex)\n\
  -f, --snpformat=INT     Print snp information from database built previously\n\
                            using snpindex (0=none, 1=alternate version only\n\
                            2=both versions merged (using N), 3=both versions separate (default)\n\
\n\
External map file options\n\
  -M, --mapdir=directory  Map directory\n\
  -m, --map=iitfile       Map file.  If argument is '?' (with the quotes),\n\
                            this lists available map files.\n\
  -S, --sequence          For a gene map file, prints the coding sequence\n\
  --gsequence             For a gene map file, prints the gene sequence (exons plus introns), one per line\n\
  -E, --exons             For a gene map file, prints the sequence, one exon per line\n\
  --nunique               For a gene map file, also prints the number of unique positions\n\
  -k, --ranks             Prints levels for non-overlapping printing of map hits\n\
  -r, --raw               Prints sequence as ASCII numeric codes\n\
  -u, --flanking=INT      Show flanking hits (default 0)\n\
  -t, --maptype=STRING    Show only intervals with given type\n\
  -s, --signed            Show only intervals with same direction as query.  If flanking hits\n\
                            are also requested, show only flanking hits downstream in direction of\n\
                            query.\n\
  --aslabel               Consider all queries to be labels, even if numeric\n\
  --coding                Print entry only if position overlaps a coding exon\n\
\n\
Dump options\n\
  -A, --dump              Dump entire genome in FASTA format\n\
  --stream                Dump entire genome as a single stream of ACGTX bytes\n\
  -L, --chromosomes       List all chromosomes with universal coordinates\n\
  --forsam                List all chromosomes for use in a SAM file\n\
  -I, --contigs           List all contigs with universal coordinates\n\
\n\
Help options\n\
  -^, --version           Show version\n\
  -?, --help              Show this help message\n\
");
  return;
}


/* Printing functions */

static void
print_two_coords (Univcoord_T left, Chrpos_T length, Univ_IIT_T chromosome_iit) {
  char *chromosome;
  Chrpos_T chrpos;

  printf("%llu%s%llu\t",(unsigned long long) left+1,SEPARATOR,(unsigned long long) left+length);
  chromosome = Univ_IIT_string_from_position(&chrpos,left,chromosome_iit);
  printf("%s:%u\t",chromosome,chrpos+1U);
  FREE(chromosome);

  chromosome = Univ_IIT_string_from_position(&chrpos,left+length-1,chromosome_iit);
  printf("%s:%u\n",chromosome,chrpos+1U);
  FREE(chromosome);
  return;
}



static void
print_sequence (Genome_T genome, Genome_T genomealt, Univcoord_T genomicstart, Chrpos_T genomiclength,
		Univ_IIT_T chromosome_iit, bool whole_chromosome_p, bool circular_part_p) {
  char *chromosome1, *chromosome2, *ptr, c;
  Sequence_T genomicseg, genomicseg_alt, genomicseg_snp;
  Chrpos_T chrpos;

  /* Handle reference strain */
  if (stream_chars_p == true || stream_ints_p == true) {
    /* Don't print a header */
  } else if (vareffect_p == true) {
    /* Don't print a header */
  } else if (user_typestring != NULL) {
    /* Don't print a header */
  } else if (header != NULL) {
    printf(">%s\n",header);
  } else if (circular_part_p == true) {
    /* Don't print a header */
  } else if (revcomp == true) {
    printf(">");
    chromosome1 = Univ_IIT_string_from_position(&chrpos,genomicstart+genomiclength-1,chromosome_iit);
    if (whole_chromosome_p == true) {
      printf("%s\n",chromosome1);
    } else {
      printf("%s:%u",chromosome1,chrpos+1U);
      printf("%s",SEPARATOR);
      chromosome2 = Univ_IIT_string_from_position(&chrpos,genomicstart+1,chromosome_iit);
      if (strcmp(chromosome2,chromosome1) == 0) {
	printf("%u",chrpos);
      } else {
	printf("%s:%u",chromosome2,chrpos);
      }
      printf(" %s:%llu%s%llu\n",dbversion,(unsigned long long) genomicstart+genomiclength,SEPARATOR,(unsigned long long) genomicstart+1U);
      FREE(chromosome2);
    }
    FREE(chromosome1);
  } else {
    printf(">");
    chromosome1 = Univ_IIT_string_from_position(&chrpos,genomicstart,chromosome_iit);
    if (whole_chromosome_p == true) {
      printf("%s\n",chromosome1);
    } else {
      printf("%s:%u",chromosome1,chrpos+1U);
      printf("%s",SEPARATOR);
      chromosome2 = Univ_IIT_string_from_position(&chrpos,genomicstart+genomiclength-1,chromosome_iit);
      if (strcmp(chromosome2,chromosome1) == 0) {
	printf("%u",chrpos+1U);
      } else {
	printf("%s:%u",chromosome2,chrpos+1U);
      }
      printf(" %s:%llu%s%llu\n",dbversion,(unsigned long long) genomicstart+1U,SEPARATOR,(unsigned long long) genomicstart+genomiclength);
      FREE(chromosome2);
    }
    FREE(chromosome1);
  }

  if (vareffect_p == true) {
    chromosome1 = Univ_IIT_string_from_position(&chrpos,genomicstart,chromosome_iit);
    genomicseg = Genome_get_segment(genome,genomicstart,genomiclength,chromosome_iit,revcomp);
    ptr = Sequence_fullpointer(genomicseg);
    while ((c = *ptr++) != '\0') {
      chrpos++;			/* Converts to a 1-based coordinate */
      switch (c) {
      case 'A':
	printf("%s\t%u\t%u\tA/C\n",chromosome1,chrpos,chrpos);
	printf("%s\t%u\t%u\tA/G\n",chromosome1,chrpos,chrpos);
	printf("%s\t%u\t%u\tA/T\n",chromosome1,chrpos,chrpos);
	break;
      case 'C':
	printf("%s\t%u\t%u\tC/A\n",chromosome1,chrpos,chrpos);
	printf("%s\t%u\t%u\tC/G\n",chromosome1,chrpos,chrpos);
	printf("%s\t%u\t%u\tC/T\n",chromosome1,chrpos,chrpos);
	break;
      case 'G':
	printf("%s\t%u\t%u\tG/A\n",chromosome1,chrpos,chrpos);
	printf("%s\t%u\t%u\tG/C\n",chromosome1,chrpos,chrpos);
	printf("%s\t%u\t%u\tG/T\n",chromosome1,chrpos,chrpos);
	break;
      case 'T':
	printf("%s\t%u\t%u\tT/A\n",chromosome1,chrpos,chrpos);
	printf("%s\t%u\t%u\tT/C\n",chromosome1,chrpos,chrpos);
	printf("%s\t%u\t%u\tT/G\n",chromosome1,chrpos,chrpos);
	break;
      default:
	printf("%s\t%u\t%u\tN/A\n",chromosome1,chrpos,chrpos);
	printf("%s\t%u\t%u\tN/C\n",chromosome1,chrpos,chrpos);
	printf("%s\t%u\t%u\tN/G\n",chromosome1,chrpos,chrpos);
	printf("%s\t%u\t%u\tN/T\n",chromosome1,chrpos,chrpos);
	break;
      }
    }
    Sequence_free(&genomicseg);
    FREE(chromosome1);

  } else if (stream_chars_p == true) {
    genomicseg = Genome_get_segment(genome,genomicstart,genomiclength,chromosome_iit,revcomp);
    Sequence_stdout_stream_chars(genomicseg);
    Sequence_free(&genomicseg);

  } else if (stream_ints_p == true) {
    genomicseg = Genome_get_segment(genome,genomicstart,genomiclength,chromosome_iit,revcomp);
    Sequence_stdout_stream_ints(genomicseg);
    Sequence_free(&genomicseg);

  } else if (snps_root == NULL || print_snps_mode == 0 || print_snps_mode == 2) {
    if (user_typestring == NULL) {
      /* Cannot use Genome_dump here because it does not allow for revcomp */
      genomicseg = Genome_get_segment(genome,genomicstart,genomiclength,chromosome_iit,revcomp);
      Sequence_stdout(genomicseg,uppercasep,wraplength,/*trimmedp*/false);
      Sequence_free(&genomicseg);
    }

  } else if (print_snps_mode == 1) {
    /* Handle both reference and alternate versions */
    genomicseg = Genome_get_segment(genome,genomicstart,genomiclength,chromosome_iit,revcomp);
    genomicseg_alt = Genome_get_segment_alt(genomealt,genomicstart,genomiclength,chromosome_iit,revcomp);
    genomicseg_snp = Genome_get_segment_snp(genomealt,genomicstart,genomiclength,chromosome_iit,revcomp);
    if (user_typestring == NULL) {
      if (rawp == true) {
	Sequence_stdout_raw(genomicseg);
      } else {
	Sequence_stdout_alt(genomicseg,genomicseg_alt,genomicseg_snp,uppercasep,wraplength);
      }
      Sequence_free(&genomicseg_snp);
      Sequence_free(&genomicseg_alt);
      Sequence_free(&genomicseg);
    }

  } else {
    /* Handle both reference and alternate versions */
    genomicseg = Genome_get_segment(genome,genomicstart,genomiclength,chromosome_iit,revcomp);
    genomicseg_snp = Genome_get_segment_snp(genomealt,genomicstart,genomiclength,chromosome_iit,revcomp);
    if (user_typestring == NULL) {
      if (rawp == true) {
	Sequence_stdout_raw(genomicseg);
      } else {
	Sequence_stdout_two(genomicseg,genomicseg_snp,uppercasep,wraplength);
      }
      Sequence_free(&genomicseg_snp);
      Sequence_free(&genomicseg);
    }
  }


#if 0
  /* Handle alternate strains */
  if (altstrainp == true || user_typestring != NULL) {
    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".altstrain.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.altstrain.iit",genomesubdir,fileroot);
    global_altstrain_iit = altstrain_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,
						    /*divread*/READ_ALL,/*divstring*/NULL,/*add_iit_p*/false);
    FREE(iitfile);
  } else {
    global_altstrain_iit = altstrain_iit = (IIT_T) NULL;
  }

  /* We rely upon the fact that gbuffer1 still holds the genomic segment.  This code is duplicated in gmap.c. */
  if (altstrain_iit != NULL) {
    if (user_typestring == NULL) {
      /* No user-specified strain.  Get all indices. */
      indexarray = IIT_get(&nindices,altstrain_iit,/*divstring*/NULL,
			   genomicstart+1,genomicstart+genomiclength-1,/*sortp*/false);
    } else {
      user_type = IIT_typeint(altstrain_iit,user_typestring);
      if (user_type < 0) {
	/* Invalid user-specified strain.  Print nothing. */
	fprintf(stderr,"No such type as %s.  Allowed strains are:\n",user_typestring);
	IIT_dump_typestrings(stderr,altstrain_iit);
	indexarray = NULL;
	nindices = 0;
	Sequence_free(&genomicseg);
      } else {
	/* Valid user-specified strain.  Get subset of indices. */
	indexarray = IIT_get_typed(&nindices,altstrain_iit,/*divstring*/NULL,
				   genomicstart+1,genomicstart+genomiclength-1,user_type,/*sortp*/false);
	if (nindices == 0) {
	  /* Print reference strain */
	  if (rawp == true) {
	    Sequence_stdout_raw(genomicseg);
	  } else {
	    Sequence_stdout(genomicseg,uppercasep,wraplength,/*trimmedp*/false);
	  }
	  Sequence_free(&genomicseg);
	}
      }
    }

    if (nindices > 0) {
      /* Sort according to type and genome position*/
      qsort(indexarray,nindices,sizeof(int),index_compare);

      for (j = 0; j < nindices; ) {
	i = j++;
	type = Interval_type(interval = IIT_interval(altstrain_iit,indexarray[i]));
	strain = IIT_typestring(altstrain_iit,type);
	sourcelength = IIT_annotation_strlen(altstrain_iit,indexarray[i]);
	if (sourcelength > Interval_length(interval)) {
	  extra = sourcelength - Interval_length(interval);
	} else {
	  extra = 0;
	}
	while (j < nindices && Interval_type(interval = IIT_interval(altstrain_iit,indexarray[j])) == type) {
	  sourcelength = IIT_annotation_strlen(altstrain_iit,indexarray[j]);
	  if (sourcelength > Interval_length(interval)) {
	    extra += sourcelength - Interval_length(interval);
	  }
	  j++;
	}
	/* Patch from i to j */
	gbuffer3 = (char *) CALLOC(genomiclength+extra+1,sizeof(char));
	genomicseg = Genome_patch_strain(&(indexarray[i]),j-i,altstrain_iit,
					 genomicstart,genomiclength,
					 revcomp,gbuffer1,gbuffer2,gbuffer3,
					 genomiclength+extra);

	if (header != NULL) {
	  printf(">%s\n",header);
	} else if (revcomp == true) {
	  printf(">%s:%u%s%u_rc variant:%s\n",
		 dbversion,genomicstart+1,SEPARATOR,genomicstart+genomiclength,strain);
	} else {
	  printf(">%s:%u%s%u variant:%s\n",
		 dbversion,genomicstart+1,SEPARATOR,genomicstart+genomiclength,strain);
	}
	if (rawp == true) {
	  Sequence_stdout_raw(genomicseg);
	} else {
	  Sequence_stdout(genomicseg,uppercasep,wraplength,/*trimmedp*/false);
	}
	Sequence_free(&genomicseg);
	FREE(gbuffer3);
      }
      FREE(indexarray);
    }
    IIT_free(&altstrain_iit);
  }
#endif

  return;
}


static int *
get_matches (int *nmatches, int *sign, char *divstring, Chrpos_T coordstart, Chrpos_T coordend, bool revcomp,
	     int **leftflanks, int *nleftflanks, int **rightflanks, int *nrightflanks,
	     char *typestring, IIT_T *iit, char *filename) {
  int *matches;
  int typeint;

  if (signedp == false) {
    *sign = 0;
  } else if (revcomp == true) {
    *sign = -1;
  } else if (coordend == coordstart) {
    *sign = 0;
  } else {
    *sign = +1;
  }

  if (*iit == NULL) {
    /* This call can give a warning if there are no entries in the map IIT file for the given div */
    if ((*iit = IIT_read(filename,/*name*/NULL,true,/*divread*/READ_ONE,divstring,/*add_iit_p*/true)) == NULL) {
      fprintf(stderr,"File %s appears to be an invalid IIT file\n",filename);
      exit(9);
    }
  }

  if (typestring == NULL) {
    /* Treat query as coordinates, without a typestring */
    if (exactp == true) {
      matches = IIT_get_exact_multiple(&(*nmatches),*iit,divstring,coordstart,coordend,/*type*/0);
    } else {
      matches = IIT_get_signed(&(*nmatches),*iit,divstring,coordstart,coordend,*sign,sortp);
    }
    if (nflanking > 0) {
      IIT_get_flanking(&(*leftflanks),&(*nleftflanks),&(*rightflanks),&(*nrightflanks),*iit,divstring,
		       coordstart,coordend,nflanking,*sign);
    }

  } else if ((typeint = IIT_typeint(*iit,typestring)) < 0) {
    fprintf(stderr,"No such type as %s.  Ignoring the type.\n",typestring);
    if (exactp == true) {
      matches = IIT_get_exact_multiple(&(*nmatches),*iit,divstring,coordstart,coordend,/*type*/0);
    } else {
      matches = IIT_get_signed(&(*nmatches),*iit,divstring,coordstart,coordend,*sign,sortp);
    }
    if (nflanking > 0) {
      IIT_get_flanking(&(*leftflanks),&(*nleftflanks),&(*rightflanks),&(*nrightflanks),*iit,divstring,
		       coordstart,coordend,nflanking,*sign);
    }

  } else {
    /* Treat query as coordinates, with a typestring */
    if (exactp == true) {
      debug(printf("Calling IIT_get_exact_multiple\n"));
      matches = IIT_get_exact_multiple(&(*nmatches),*iit,divstring,coordstart,coordend,typeint);
    } else {
      debug(printf("Calling IIT_get_typed_signed with sign %d\n",*sign));
      matches = IIT_get_typed_signed(&(*nmatches),*iit,divstring,coordstart,coordend,typeint,*sign,sortp);
    }
    if (nflanking > 0) {
      IIT_get_flanking_typed(&(*leftflanks),&(*nleftflanks),&(*rightflanks),&(*nrightflanks),*iit,divstring,
			     coordstart,coordend,nflanking,typeint,*sign);
    }
  }

  return matches;
}


static char complCode[128] = COMPLEMENT_LC;

#if 0
static void
make_complement_buffered (char *complement, char *sequence, Chrpos_T length) {
  int i, j;

  for (i = length-1, j = 0; i >= 0; i--, j++) {
    complement[j] = complCode[(int) sequence[i]];
  }
  complement[length] = '\0';
  return;
}
#endif

static void
make_complement_inplace (char *sequence, Chrpos_T length) {
  char temp;
  Chrpos_T i, j;

  for (i = 0, j = length-1; i < length/2; i++, j--) {
    temp = complCode[(int) sequence[i]];
    sequence[i] = complCode[(int) sequence[j]];
    sequence[j] = temp;
  }
  if (i == j) {
    sequence[i] = complCode[(int) sequence[i]];
  }

  return;
}



static void
genemap_print_exons (char *annot, Univcoord_T chroffset, Genome_T genome) {
  char *p;
  char *gbuffer;
  Chrpos_T exonstart, exonend, exonlow, exonhigh, exonlength;

  /* Skip header */
  p = annot;
  while (*p != '\0' && *p != '\n') {
    putchar(*p);
    p++;
  }
  if (*p == '\n') p++;
  printf("\n");

  while (*p != '\0') {
    if (sscanf(p,"%u %u",&exonstart,&exonend) != 2) {
      fprintf(stderr,"Can't parse exon coordinates in %s\n",p);
      abort();
    } else {
      if (exonstart <= exonend) {
	exonlow = exonstart;
	exonhigh = exonend;
	exonlength = exonhigh - exonlow + 1;
	gbuffer = (char *) CALLOC(exonlength+1,sizeof(char));
	Genome_fill_buffer(genome,/*left*/chroffset-1U+exonlow,exonlength,gbuffer);
      } else {
	exonlow = exonend;
	exonhigh = exonstart;
	exonlength = exonhigh - exonlow + 1;
	gbuffer = (char *) CALLOC(exonlength+1,sizeof(char));
	Genome_fill_buffer(genome,/*left*/chroffset-1U+exonlow,exonlength,gbuffer);
	make_complement_inplace(gbuffer,exonlength);
      }
      printf("%s\n",gbuffer);
      FREE(gbuffer);
    }

    while (*p != '\0' && *p != '\n') p++;
    if (*p == '\n') p++;
  }

  return;
}


static void
genemap_print_sequence (char *annot, Univcoord_T chroffset, Genome_T genome) {
  char *p;
  char *gbuffer;
  Chrpos_T exonstart, exonend, exonlow, exonhigh, exonlength;
  int col, i;

  /* Print header */
  /* printf(">"); -- accession already includes this */
  p = annot;
  while (*p != '\0' && *p != '\n') {
    putchar(*p);		/* Gene information */
    p++;
  }
  if (*p == '\n') p++;
  printf("\n");

  col = 0;
  while (*p != '\0') {
    if (sscanf(p,"%u %u",&exonstart,&exonend) != 2) {
      fprintf(stderr,"Can't parse exon coordinates in %s\n",p);
      abort();
    } else {
      if (exonstart <= exonend) {
	exonlow = exonstart;
	exonhigh = exonend;
	exonlength = exonhigh - exonlow + 1;
	gbuffer = (char *) CALLOC(exonlength+1,sizeof(char));
	Genome_fill_buffer(genome,/*left*/chroffset-1U+exonlow,exonlength,gbuffer);
      } else {
	exonlow = exonend;
	exonhigh = exonstart;
	exonlength = exonhigh - exonlow + 1;
	gbuffer = (char *) CALLOC(exonlength+1,sizeof(char));
	Genome_fill_buffer(genome,/*left*/chroffset-1U+exonlow,exonlength,gbuffer);
	make_complement_inplace(gbuffer,exonlength);
      }
      
      for (i = 0; i < (int) exonlength; i++) {
	putchar(gbuffer[i]);
	if (++col >= wraplength) {
	  printf("\n");
	  col = 0;
	}
      }
      FREE(gbuffer);
    }

    while (*p != '\0' && *p != '\n') p++;
    if (*p == '\n') p++;
  }

  if (col > 0) {
    printf("\n");
  }

  return;
}


static void
genemap_print_gsequence (char *annot, Univcoord_T chroffset, Genome_T genome) {
  char *p;
  char *gbuffer;
  Chrpos_T exonstart, exonend, exonlow, exonhigh, exonlength, intronlength,
    prev_exonlow = 0, prev_exonhigh = 0;

  /* Print header */
  /* printf(">"); -- accession already includes this */
  p = annot;
  while (*p != '\0' && *p != '\n') {
    putchar(*p);
    p++;
  }
  if (*p == '\n') p++;
  printf("\n");

  while (*p != '\0') {
    if (sscanf(p,"%u %u",&exonstart,&exonend) != 2) {
      fprintf(stderr,"Can't parse exon coordinates in %s\n",p);
      abort();
    } else {
      if (exonstart <= exonend) {
	exonlow = exonstart;
	exonhigh = exonend;
	if (prev_exonhigh > 0) {
	  intronlength = exonlow - prev_exonhigh - 1;
	  gbuffer = (char *) CALLOC(intronlength+1,sizeof(char));
	  Genome_fill_buffer(genome,/*left*/chroffset-1U+(prev_exonhigh+1),
			     intronlength,gbuffer);
	  printf("%s\n",gbuffer);
	  FREE(gbuffer);
	}
	exonlength = exonhigh - exonlow + 1;
	gbuffer = (char *) CALLOC(exonlength+1,sizeof(char));
	Genome_fill_buffer(genome,/*left*/chroffset-1U+exonlow,exonlength,gbuffer);
	printf("%s\n",gbuffer);
	FREE(gbuffer);

      } else {
	exonlow = exonend;
	exonhigh = exonstart;
	if (prev_exonlow > 0) {
	  intronlength = prev_exonlow - 1 - exonhigh;
	  gbuffer = (char *) CALLOC(intronlength+1,sizeof(char));
	  Genome_fill_buffer(genome,/*left*/chroffset-1U+(exonhigh+1),
			     intronlength,gbuffer);
	  make_complement_inplace(gbuffer,intronlength);
	  printf("%s\n",gbuffer);
	  FREE(gbuffer);
	}
	exonlength = exonhigh - exonlow + 1;
	gbuffer = (char *) CALLOC(exonlength+1,sizeof(char));
	Genome_fill_buffer(genome,/*left*/chroffset-1U+exonlow,exonlength,gbuffer);
	make_complement_inplace(gbuffer,exonlength);
	printf("%s\n",gbuffer);
	FREE(gbuffer);
      }

      prev_exonlow = exonlow;
      prev_exonhigh = exonhigh;
    }

    while (*p != '\0' && *p != '\n') p++;
    if (*p == '\n') p++;
  }

  return;
}


static void
genemap_print_unique (char *annot, Intlist_T unique_positions, Intlist_T unique_splicep, bool print_geneline_p) {
  char *p;
  Chrpos_T exonstart, exonend;
  int nunique, uniquep;

  /* Skip header */
  p = annot;
  if (print_geneline_p == true) {
    while (*p != '\0' && *p != '\n') {
      putchar(*p);
      p++;
    }
    if (*p == '\n') p++;
    printf("\n");
  } else {
    while (*p != '\0' && *p != '\n') p++;
    if (*p == '\n') p++;
  }

  while (*p != '\0') {
    if (sscanf(p,"%u %u",&exonstart,&exonend) != 2) {
      fprintf(stderr,"Can't parse exon coordinates in %s\n",p);
      abort();
    } else {
      unique_positions = Intlist_pop(unique_positions,&nunique);
      if (unique_splicep == NULL) {
	printf("%u %u %d NA\n",exonstart,exonend,nunique);
      } else {
	unique_splicep = Intlist_pop(unique_splicep,&uniquep);
	if (uniquep == 0) {
	  printf("%u %u %d F\n",exonstart,exonend,nunique);
	} else {
	  printf("%u %u %d T\n",exonstart,exonend,nunique);
	}
      }
    }

    while (*p != '\0' && *p != '\n') p++;
    if (*p == '\n') p++;
  }

  return;
}


static void
print_interval (char *divstring, int index, IIT_T iit, int ndivs, Univ_IIT_T chromosome_iit,
		Genome_T genome, int fieldint) {
  Interval_T interval;
  Univinterval_T chr_interval;
  char *label, *annotation, *restofheader;
  bool allocp;
  bool annotationonlyp = false, signed_output_p = true;
  int divno;
  Univcoord_T chroffset;
  Intlist_T unique_positions, unique_splicep;

  if (exonsp == true || sequencep == true || genome_sequence_p == true || uniquep == true) {
    label = IIT_label(iit,index,&allocp);
    printf(">%s ",label);
    if (allocp == true) {
      FREE(label);
    }

    annotation = IIT_annotation(&restofheader,iit,index,&allocp);
    printf("%s",restofheader);
    if (IIT_version(iit) < 5) {
      printf("\n");
    }

    if (ndivs > 1) {
      if (divstring == NULL) {
	/* For example, if interval was retrieved by label */
	divstring = IIT_divstring_from_index(iit,index);
      }
    }

    if ((divno = Univ_IIT_find_one(chromosome_iit,divstring)) < 0) {
      fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file\n",divstring);
      /* exit(9); */
    } else {
      chr_interval = Univ_IIT_interval(chromosome_iit,divno);
      chroffset = Univinterval_low(chr_interval);
      if (exonsp == true) {
	genemap_print_exons(annotation,chroffset,genome);
      } else if (sequencep == true) {
	genemap_print_sequence(annotation,chroffset,genome);
      } else if (genome_sequence_p == true) {
	genemap_print_gsequence(annotation,chroffset,genome);
      } else if (uniquep == true) {
	unique_positions = IIT_unique_positions(iit,index,divno);
	unique_splicep = IIT_unique_splicep(iit,index,divno);
	genemap_print_unique(annotation,unique_positions,unique_splicep,/*print_geneline_p*/true);
	/* Intlist_free(&unique_splicep); -- popped by above procedure */
	/* Intlist_free(&unique_positions); -- popped by above procedure */
      }
    }
    if (allocp == true) {
      FREE(restofheader);
    }

  } else {
    if (annotationonlyp == false) {
      label = IIT_label(iit,index,&allocp);
      printf(">%s ",label);
      if (allocp == true) {
	FREE(label);
      }
      
      if (ndivs > 1) {
	if (divstring == NULL) {
	  /* For example, if interval was retrieved by label */
	  divstring = IIT_divstring_from_index(iit,index);
	}
	printf("%s:",divstring);
      }

      interval = IIT_interval(iit,index);
      if (signed_output_p == false) {
	printf("%u..%u",Interval_low(interval),Interval_high(interval));
      } else if (Interval_sign(interval) < 0) {
	printf("%u..%u",Interval_high(interval),Interval_low(interval));
      } else {
	printf("%u..%u",Interval_low(interval),Interval_high(interval));
      }
      if (Interval_type(interval) > 0) {
	printf(" %s",IIT_typestring(iit,Interval_type(interval)));
      }
#if 0
      /* Unnecessary because of "\n" after restofheader below */
      if (IIT_version(iit) < 5) {
	printf("\n");
      }
#endif

    }

    if (fieldint < 0) {
      annotation = IIT_annotation(&restofheader,iit,index,&allocp);
      printf("%s\n",restofheader);
      printf("%s",annotation);
      if (allocp == true) {
	FREE(restofheader);
      }

    } else {
      annotation = IIT_fieldvalue(iit,index,fieldint);
      printf("%s",annotation);
      FREE(annotation);
    }
  }

  return;
}


int
main (int argc, char *argv[]) {
  char *snpsdir = NULL;
  char *iitfile;
  Genome_T genome = NULL, genomealt = NULL;
  Univcoord_T genomicstart, chroffset;
  Univinterval_T interval;
  Chrpos_T genomiclength, chrlength, chrstart, chrend;
  char *mapdir = NULL, *fileroot = NULL;
  char *typestring_ptr = NULL;
  char *divstring, *divstring2;
  Univ_IIT_T chromosome_iit, contig_iit;
  IIT_T map_iit = NULL;
  Divread_T divread;
  char *segment, *coords, *typestring, *line,
    *annot, *accession, *space, *string, *p, *q;
  FILE *fp;
  int annot_length, line_length;
  List_T lines, l;
  bool whole_chromosome_p, chr_valid_p;
  int divno;
  Univinterval_T chr_interval;

  int index;
  Intlist_T unique_positions, unique_splicep;
  char *label, *annotation, *restofheader;

  int fieldint = -1;
  int *matches, nmatches, ndivs, i, *leftflanks, *rightflanks, nleftflanks = 0, nrightflanks = 0;
  int sign;

  int circular_typeint;
  bool *circularp = NULL;
  bool any_circular_p;

  char *chr;
  int indx;
  bool allocp;

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"D:d:c:CUl:h:V:v:f:M:m:kru:ESALI^?",
			    long_options,&long_option_index)) != -1) {
    switch (opt) {
    case 0:
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"version")) {
	print_program_version();
	exit(0);
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	exit(0);

      } else if (!strcmp(long_name,"genes")) {
	genes_file_p = true;

      } else if (!strcmp(long_name,"gsequence")) {
	genome_sequence_p = true;

      } else if (!strcmp(long_name,"nunique")) {
	uniquep = true;

      } else if (!strcmp(long_name,"exact")) {
	exactp = true;

      } else if (!strcmp(long_name,"aslabel")) {
	force_label_p = true;

      } else if (!strcmp(long_name,"forsam")) {
	dumpchrp = true;
	dumpchr_forsam_p = true;

      } else if (!strcmp(long_name,"vareffect")) {
	vareffect_p = true;

      } else if (!strcmp(long_name,"stream-chars")) {
	stream_chars_p = true;

      } else if (!strcmp(long_name,"stream-ints")) {
	stream_ints_p = true;

      } else if (!strcmp(long_name,"add-circular")) {
	add_circular_p = true;

      } else if (!strcmp(long_name,"coding")) {
	codingp = true;

      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'get-genome --help'\n",long_name);
	exit(9);
      }
      break;

    case 'D': user_genomedir = optarg; break;
    case 'd': 
      dbroot = (char *) CALLOC(strlen(optarg)+1,sizeof(char));
      strcpy(dbroot,optarg);
      break;
    case 'c': transcriptome_dbroot = optarg; break;

    case 'C': coordp = true; break;
    case 'U': uppercasep = true; break;
    case 'l': wraplength = atoi(optarg); break;
    case 'h': header = optarg; break;
    case 'V': user_snpsdir = optarg; break;
    case 'v': snps_root = optarg; break;
    case 'f': print_snps_mode = atoi(optarg); break;

    case 'M': user_mapdir = optarg; break;
    case 'm': map_iitfile = optarg; break;
    case 'k': levelsp = true; break;
    case 'r': rawp = true; break;
    case 'u': nflanking = atoi(optarg); break;
    case 'E': exonsp = true; break;
    case 'S': sequencep = true; break;
    case 's': signedp = true; break;

    case 'A': dumpallp = true; break;
    case 'L': dumpchrp = true; dumpchr_forsam_p = false; break;
    case 'I': dumpsegsp = true; break;

    case '^': print_program_version(); exit(0);
    case '?': print_program_usage(); exit(0);
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;
  
  if (dbroot == NULL) {
    print_program_usage();
    exit(9);
  } else if (!strcmp(dbroot,"?")) {
    Datadir_avail_gmap_databases(stdout,user_genomedir);
    exit(0);
  } else {
    genomesubdir = Datadir_find_genomesubdir(&fileroot,&dbversion,user_genomedir,dbroot);
  }

  if (transcriptome_dbroot != NULL) {
    iitfile = (char *) MALLOC((strlen(genomesubdir)+strlen("/")+
			       strlen(fileroot)+strlen(".transcripts/")+
			       strlen(transcriptome_dbroot)+strlen(".genes.iit")+1)*sizeof(char));
    sprintf(iitfile,"%s/%s.transcripts/%s.genes.iit",genomesubdir,fileroot,
	    transcriptome_dbroot);
    if ((map_iit = IIT_read(iitfile,/*name*/NULL,true,/*divread*/READ_ALL,/*divstring*/NULL,
			    /*add_iit_p*/true)) == NULL) {
      fprintf(stderr,"Cannot open IIT file %s\n",iitfile);
      exit(9);
    } else if ((index = IIT_find_one(map_iit,argv[0])) < 0) {
      fprintf(stderr,"Cannot find %s in map file %s\n",argv[0],iitfile);
    } else {
      annotation = IIT_annotation(&restofheader,map_iit,index,&allocp);
      printf("%s",annotation);
      if (allocp == true) {
	FREE(restofheader);
      }
    }

    IIT_free(&map_iit);
    FREE(iitfile);

    return 0;
  }


  if (user_snpsdir == NULL) {
    snpsdir = genomesubdir;
  } else {
    snpsdir = user_snpsdir;
  }

  if (stream_chars_p == true || stream_ints_p == true) {
    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    chromosome_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false);
    circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
    circularp = Univ_IIT_circularp(&any_circular_p,chromosome_iit);
    FREE(iitfile);

    genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,
			chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
    Genome_setup(circular_typeint);

    for (indx = 1; indx <= Univ_IIT_total_nintervals(chromosome_iit); indx++) {
      chr = Univ_IIT_label(chromosome_iit,indx,&allocp);
      interval = Univ_IIT_interval(chromosome_iit,indx);
      genomicstart = Univinterval_low(interval);
      genomiclength = Univinterval_length(interval);

      /* Do not print header */
      if (add_circular_p == true && circularp[indx] == true) {
	Genome_dump(genome,genomicstart,genomicstart + 2*genomiclength,wraplength);
      } else {
	Genome_dump(genome,genomicstart,genomicstart + genomiclength,wraplength);
      }

      if (allocp == true) {
	FREE(chr);
      }
    }

    Genome_free(&genome);

    Univ_IIT_free(&chromosome_iit);

    return 0;


  } else if (dumpallp == true) {
    if (map_iitfile != NULL) {
      if (Access_file_exists_p(map_iitfile) == true) {
	iitfile = (char *) CALLOC(strlen(map_iitfile)+1,sizeof(char));
	strcpy(iitfile,map_iitfile);

      } else {
	mapdir = Datadir_find_mapdir(user_mapdir,genomesubdir,fileroot);
	iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
				  strlen(map_iitfile)+strlen(".iit")+1,sizeof(char));
	sprintf(iitfile,"%s/%s.iit",mapdir,map_iitfile);

	if (Access_file_exists_p(iitfile) == false) {
	  /* Try main level of genome */
	  FREE(iitfile);
	  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				    strlen(map_iitfile)+strlen(".iit")+1,sizeof(char));
	  sprintf(iitfile,"%s/%s.iit",genomesubdir,map_iitfile);
	  
	  if (Access_file_exists_p(iitfile) == false) {
	    fprintf(stderr,"Map file %s.iit not found in %s.  Available files:\n",map_iitfile,mapdir);
	    Datadir_list_directory(stderr,mapdir);
	    fprintf(stderr,"Either install file %s.iit or specify a full directory path\n",map_iitfile);
	    fprintf(stderr,"using the -M flag to gmap.\n");
	    exit(9);
	  }
	}
	FREE(mapdir);
      }
    
      if (sequencep == true || genome_sequence_p == true || exonsp == true) {
	/* User is requesting all sequences or exons in the map file */
	if ((map_iit = IIT_read(iitfile,/*name*/NULL,true,/*divread*/READ_ALL,/*divstring*/NULL,
				/*add_iit_p*/true)) == NULL) {
	  fprintf(stderr,"Cannot open IIT file %s\n",iitfile);
	  exit(9);

	} else {
	  FREE(iitfile);

	  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
	  sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
	  chromosome_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false);
	  circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
	  FREE(iitfile);

	  genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,
			      chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
	  Genome_setup(circular_typeint);

	  ndivs = IIT_ndivs(map_iit);
	  for (indx = 1; indx <= Univ_IIT_total_nintervals(chromosome_iit); indx++) {
	    chr = Univ_IIT_label(chromosome_iit,indx,&allocp);
	    chrlength = Univ_IIT_length(chromosome_iit,indx);
	    matches = IIT_get(&nmatches,map_iit,chr,/*coordstart*/0,/*coordend*/chrlength,/*sortp*/false);
	    for (i = 0; i < nmatches; i++) {
	      print_interval(chr,matches[i],map_iit,ndivs,chromosome_iit,genome,/*fieldint*/-1);
	    }
	    FREE(matches);
	    if (allocp) {
	      FREE(chr);
	    }
	  }

	  Genome_free(&genome);
	  Univ_IIT_free(&chromosome_iit);
	  IIT_free(&map_iit);
	  FREE(dbversion);
	  FREE(genomesubdir);
	  FREE(fileroot);
	  FREE(dbroot);

	  return 0;
	}

      } else {
	/* User is requesting a dump of the map file contents */
	/* Don't need to find competing labels, since they are given, so divread can be READ_NONE */
	if ((map_iit = IIT_read(iitfile,/*name*/NULL,true,/*divread*/READ_NONE,/*divstring*/NULL,
				/*add_iit_p*/true)) == NULL) {
	  fprintf(stderr,"Cannot open IIT file %s\n",iitfile);
	  exit(9);

	} else {
	  FREE(iitfile);

	  IIT_dump(map_iit,/*sortp*/false);
	  IIT_free(&map_iit);
	  FREE(dbversion);
	  FREE(genomesubdir);
	  FREE(fileroot);
	  FREE(dbroot);
	  return 0;
	}
      }
      
    } else {
      /* User is requesting a dump of the genome */

      iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
      sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
      chromosome_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false);
      circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
      circularp = Univ_IIT_circularp(&any_circular_p,chromosome_iit); /* Needed if --with-circular is specified */
      FREE(iitfile);

      if (snps_root == NULL || print_snps_mode == 0) {
	genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,
			    chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
      } else if (print_snps_mode == 2) {
	genome = Genome_new(snpsdir,fileroot,snps_root,
			    chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
      } else if (print_snps_mode == 1 || print_snps_mode == 3) {
	genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,
			    chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
	genomealt = Genome_new(snpsdir,fileroot,snps_root,
			       chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
      }
      Genome_setup(circular_typeint);

      for (indx = 1; indx <= Univ_IIT_total_nintervals(chromosome_iit); indx++) {
	chr = Univ_IIT_label(chromosome_iit,indx,&allocp);

#if 0
	/* Previously method was too slow */
	Parserange_universal(&segment,&revcomp,&whole_chromosome_p,&genomicstart,&genomiclength,&chrstart,&chrend,
			     &chroffset,&chrlength,with_colon,genomesubdir,fileroot);
	print_sequence(genome,genomealt,genomicstart,genomiclength,chromosome_iit,
		       /*whole_chromosome_p*/true,/*circular_part_p*/false);
#endif

	printf(">%s\n",chr);
	interval = Univ_IIT_interval(chromosome_iit,indx);
	genomicstart = Univinterval_low(interval);
	genomiclength = Univinterval_length(interval);
	if (add_circular_p == true && circularp[indx] == true) {
	  Genome_dump(genome,genomicstart,genomicstart + 2*genomiclength,wraplength);
	} else {
	  Genome_dump(genome,genomicstart,genomicstart + genomiclength,wraplength);
	}
	if (allocp == true) {
	  FREE(chr);
	}
      }

      FREE(circularp);
      if (genomealt != NULL) {
	Genome_free(&genomealt);
      }
      Genome_free(&genome);

      Univ_IIT_free(&chromosome_iit);
      FREE(dbversion);
      FREE(genomesubdir);
      FREE(fileroot);
      FREE(dbroot);

      return 0;
    }

  } else if (dumpchrp == true) {
    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    chromosome_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false);
    /* circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular"); */
    FREE(iitfile);

    if (dumpchr_forsam_p == true) {
      Univ_IIT_dump_sam(/*fp*/stdout,chromosome_iit,/*sam_read_group_id*/NULL,/*sam_read_group_name*/NULL,
			/*sam_read_group_library*/NULL,/*sam_read_group_platform*/NULL);
    } else {
      Univ_IIT_dump_table(chromosome_iit,/*zerobasedp*/false);
    }
    Univ_IIT_free(&chromosome_iit);
    return 0;

  } else if (dumpsegsp == true) {
    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
    chromosome_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false);
    /* circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular"); */
    FREE(iitfile);

    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".contig.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.contig.iit",genomesubdir,fileroot);
    contig_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false);
    FREE(iitfile);

    Univ_IIT_dump_contigs(contig_iit,chromosome_iit,/*directionalp*/true);
    Univ_IIT_free(&contig_iit);
    FREE(dbversion);
    FREE(genomesubdir);
    FREE(fileroot);
    FREE(dbroot);
    return 0;

  }

#if 0
  if (replacex == true) {
    Genome_replace_x();
  }
#endif

  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
  chromosome_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false);
  circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
  FREE(iitfile);


  if (genes_file_p == true) {
    if (snps_root == NULL || print_snps_mode == 0) {
      genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,
			  chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
    } else if (print_snps_mode == 2) {
      genome = Genome_new(snpsdir,fileroot,snps_root,
			  chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
    } else if (print_snps_mode == 1 || print_snps_mode == 3) {
      genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,
			  chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
      genomealt = Genome_new(snpsdir,fileroot,snps_root,
			     chromosome_iit,genome_access,/*sharedp*/false,/*recompp*/false);
    }
    Genome_setup(circular_typeint);

    for (i = 0; i < argc; i++) {
      if ((fp = fopen(argv[i],"r")) == NULL) {
	fprintf(stderr,"Cannot open file %s\n",argv[i]);
	exit(9);
      } else {
	lines = (List_T) NULL;
	annot_length = 0;
	while ((line = Getline_wlinefeed(&line_length,fp)) != NULL) {
	  if (line[0] != '>') {
	    if (chr_valid_p == true) {
	      lines = List_push(lines,(void *) line);
	      annot_length += line_length;
	    }

	  } else {
	    /* Process previous gene */
	    if (lines != NULL) {
	      /* Concatenate lines */
	      lines = List_reverse(lines);
	      p = annot = (char *) MALLOC((annot_length+1)*sizeof(char));
	      for (l = lines; l != NULL; l = List_next(l)) {
		string = (char *) List_head(l);
		strcpy(p,string);
		p += strlen(string);
		FREE(string);
	      }

	      genemap_print_sequence(annot,chroffset,genome);
	      FREE(annot);
	      List_free(&lines);
	    }

	    /* Parse header to get accession and chroffset */
	    p = line;
	    while (*p != '\0' && *p != ' ') {
	      p++;
	    }
	    if (*p == '\0') {
	      fprintf(stderr,"Could not find an accession in %s",line);
	      exit(9);
	    } else {
	      accession = (char *) MALLOC((p - line + 1) * sizeof(char));
	      strncpy(accession,line,p - line);
	      accession[p - line] = '\0';

	      q = ++p;		/* Skip space */
	      while (*q != '\0' && *q != ':') {
		q++;
	      }
	      if (*q == '\0') {
		fprintf(stderr,"Could not find a divstring in %s",line);
		exit(9);
	      } else {
		divstring = (char *) MALLOC((q - p + 1) * sizeof(char));
		strncpy(divstring,p,q - p);
		divstring[q - p] = '\0';
		if ((divno = Univ_IIT_find_one(chromosome_iit,divstring)) < 0) {
		  fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file, so skipping %s",
			  divstring,line);
		  chr_valid_p = false;
		} else {
		  chr_interval = Univ_IIT_interval(chromosome_iit,divno);
		  chroffset = Univinterval_low(chr_interval);
		  chr_valid_p = true;
		}
		FREE(divstring);
	      }
	    }
	    FREE(line);

	    lines = (List_T) NULL;
	    annot_length = 0;
	    if (chr_valid_p == false) {
	      FREE(accession);
	    } else {
	      lines = List_push(lines,(void *) accession);
	      annot_length += strlen(accession);
	      space = (char *) MALLOC(2*sizeof(char));
	      space[0] = ' ';
	      space[1] = '\0';
	      lines = List_push(lines,(void *) space);
	      annot_length += strlen(space);
	    }
	  }
	}

	if (lines != NULL) {
	  /* Concatenate lines */
	  lines = List_reverse(lines);
	  p = annot = (char *) MALLOC((annot_length+1)*sizeof(char));
	  for (l = lines; l != NULL; l = List_next(l)) {
	    string = (char *) List_head(l);
	    strcpy(p,string);
	    p += strlen(string);
	    FREE(string);
	  }

	  genemap_print_sequence(annot,chroffset,genome);
	  FREE(annot);
	  List_free(&lines);
	}
      }
      fclose(fp);
    }

  } else if (argc > 1 && map_iitfile != NULL) {
    /* User is requesting multiple labels */

    mapdir = Datadir_find_mapdir(user_mapdir,genomesubdir,fileroot);
    iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
			      strlen(map_iitfile)+strlen(".iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.iit",mapdir,map_iitfile);

    if (Access_file_exists_p(iitfile) == false) {
      /* Try main level of genome */
      FREE(iitfile);
      iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				strlen(map_iitfile)+strlen(".iit")+1,sizeof(char));
      sprintf(iitfile,"%s/%s.iit",genomesubdir,map_iitfile);

      if (Access_file_exists_p(iitfile) == false) {
	fprintf(stderr,"Map file %s.iit not found in %s or %s.  Available files:\n",map_iitfile,mapdir,genomesubdir);
	Datadir_list_directory(stderr,mapdir);
	fprintf(stderr,"Either install file %s.iit or specify a full directory path\n",map_iitfile);
	fprintf(stderr,"using the -M flag to gmap.\n");
	exit(9);
      }
    }
    
    /* Don't need to find competing labels, since they are given, so divread can be READ_NONE */
    if ((map_iit = IIT_read(iitfile,/*name*/NULL,true,/*divread*/READ_NONE,/*divstring*/NULL,
			    /*add_iit_p*/true)) == NULL) {
      fprintf(stderr,"Cannot open IIT file %s\n",iitfile);
      exit(9);
    } else {
      ndivs = IIT_ndivs(map_iit);
    }

    if (uniquep == false) {
      for (i = 0; i < argc; i++) {
	if ((index = IIT_find_one(map_iit,argv[i])) < 0) {
	  fprintf(stderr,"Cannot find %s in map file %s\n",argv[i],iitfile);
	} else {
	  print_interval(/*divstring*/NULL,index,map_iit,ndivs,chromosome_iit,/*genome*/NULL,/*fieldint*/-1);
	}
      }

    } else {
      /* User is requesting a comparison of labels for identifying uniqueness */
      nmatches = argc;
      matches = (int *) MALLOC(nmatches * sizeof(int));
      for (i = 0; i < nmatches; i++) {
	if ((matches[i] = IIT_find_one(map_iit,argv[i])) < 0) {
	  fprintf(stderr,"Cannot find %s in map file %s\n",argv[i],iitfile);
	  return 9;
	}
      }
      
      for (i = 0; i < nmatches; i++) {
	unique_positions = IIT_unique_positions_given_others(map_iit,matches[i],matches,nmatches);
	unique_splicep = IIT_unique_splicep_given_others(map_iit,matches[i],matches,nmatches);

	label = IIT_label(map_iit,matches[i],&allocp);
	printf(">%s",label);		/* Separator needed for multiple labels */
	/* TODO: Print coordinates */
	printf("\n");
	if (allocp == true) {
	  FREE(label);
	}

	annotation = IIT_annotation(&restofheader,map_iit,matches[i],&allocp);
	genemap_print_unique(annotation,unique_positions,unique_splicep,/*print_geneline_p*/false);
	if (allocp == true) {
	  FREE(annotation);
	}
      }
    }

    FREE(iitfile);
    FREE(mapdir);

  } else if (argc >= 1) {
    /* User is providing a single query */
    if (coordp == true) {
      debug(printf("coordp is true\n"));
      if (Parserange_universal(&segment,&revcomp,&whole_chromosome_p,&genomicstart,&genomiclength,&chrstart,&chrend,
			       &chroffset,&chrlength,argv[0],genomesubdir,fileroot) == true) {
	debug(printf("Query %s parsed as: genomicstart = %llu, genomiclength = %u, revcomp = %d\n",
		     argv[0],(unsigned long long) genomicstart,genomiclength,revcomp));
	print_two_coords(genomicstart,genomiclength,chromosome_iit);
      }
      FREE(segment);
      
    } else if (map_iitfile == NULL) {
      debug(printf("No map file\n"));

      if (snps_root == NULL || print_snps_mode == 0) {
	genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,
			    chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
      } else if (print_snps_mode == 2) {
	genome = Genome_new(snpsdir,fileroot,snps_root,
			    chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
      } else if (print_snps_mode == 1 || print_snps_mode == 3) {
	genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,
			    chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
	genomealt = Genome_new(snpsdir,fileroot,snps_root,
			       chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
      }
      Genome_setup(circular_typeint);

      if (Parserange_universal(&segment,&revcomp,&whole_chromosome_p,&genomicstart,&genomiclength,&chrstart,&chrend,
			       &chroffset,&chrlength,argv[0],genomesubdir,fileroot) == true) {
	debug(printf("Query %s parsed as: genomicstart = %llu, genomiclength = %u, revcomp = %d\n",
		     argv[0],(unsigned long long) genomicstart,genomiclength,revcomp));
	print_sequence(genome,genomealt,genomicstart,genomiclength,chromosome_iit,
		       whole_chromosome_p,/*circular_part_p*/false);
      }
      FREE(segment);

    } else if (!strcmp(map_iitfile,"?")) {
      Datadir_avail_maps(stdout,user_mapdir,genomesubdir,fileroot);
      exit(0);

    } else {
      debug(printf("Map file\n"));

      mapdir = Datadir_find_mapdir(user_mapdir,genomesubdir,fileroot);
      iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
				strlen(map_iitfile)+strlen(".iit")+1,sizeof(char));
      sprintf(iitfile,"%s/%s.iit",mapdir,map_iitfile);

      if (Access_file_exists_p(iitfile) == false) {
	/* Try main level of genome */
	FREE(iitfile);
	iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				  strlen(map_iitfile)+strlen(".iit")+1,sizeof(char));
	sprintf(iitfile,"%s/%s.iit",genomesubdir,map_iitfile);

	if (Access_file_exists_p(iitfile) == false) {
	  fprintf(stderr,"Map file %s.iit not found in %s or %s.  Available files:\n",map_iitfile,mapdir,genomesubdir);
	  Datadir_list_directory(stderr,mapdir);
	  fprintf(stderr,"Either install file %s.iit or specify a full directory path\n",map_iitfile);
	  fprintf(stderr,"using the -M flag to gmap.\n");
	  exit(9);
	}
      }

      if (exonsp == true || sequencep == true || genome_sequence_p == true) {
	if (snps_root == NULL || print_snps_mode == 0) {
	  genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,
			      chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
	} else if (print_snps_mode == 2) {
	  genome = Genome_new(snpsdir,fileroot,snps_root,
			      chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
	} else if (print_snps_mode == 1 || print_snps_mode == 3) {
	  genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,
			      chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
#if 0
	  genomealt = Genome_new(snpsdir,fileroot,snps_root,
				 chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
#endif
	}
	Genome_setup(circular_typeint);
      }

      segment = (char *) NULL;
      if (force_label_p == false &&
	  Parserange_universal(&segment,&revcomp,&whole_chromosome_p,&genomicstart,&genomiclength,&chrstart,&chrend,
			       &chroffset,&chrlength,argv[0],genomesubdir,fileroot) == true) {
	debug(printf("Query %s parsed as: genomicstart = %llu, genomiclength = %u, revcomp = %d\n",
		     argv[0],(unsigned long long) genomicstart,genomiclength,revcomp));
	divstring = Univ_IIT_string_from_position(&chrstart,genomicstart,chromosome_iit);
	divstring2 = Univ_IIT_string_from_position(&chrend,genomicstart+genomiclength-1U,chromosome_iit);
	if (strcmp(divstring,divstring2)) {
	  fprintf(stderr,"Coordinates cross chromosomal boundary\n");
	  exit(9);
	} else {
	  chrstart += 1U;
	  chrend += 1U;
	  debug(printf("Query translated to %s:%u..%u\n",divstring,chrstart,chrend));
	}
	if (argc <= 1) {
	  debug(printf("No typestring\n"));
	  typestring_ptr = (char *) NULL;
	} else {
	  debug(printf("Typestring is %s\n",argv[1]));
	  typestring_ptr = argv[1];
	}
	matches = get_matches(&nmatches,&sign,divstring,chrstart,chrend,revcomp,
			      &leftflanks,&nleftflanks,&rightflanks,&nrightflanks,
			      typestring_ptr,&map_iit,/*filename*/iitfile);
	ndivs = IIT_ndivs(map_iit);
	if (nflanking > 0) {
	  if (sign != +1) {
	    for (i = nleftflanks-1; i >= 0; i--) {
	      if (codingp == false || IIT_gene_overlapp(map_iit,leftflanks[i],chrstart,chrend) == true) {
		print_interval(divstring,leftflanks[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	      }
	    }
	  }
	  printf("====================\n");
	  FREE(leftflanks);
	}
	
	for (i = 0; i < nmatches; i++) {
	  if (codingp == false || IIT_gene_overlapp(map_iit,matches[i],chrstart,chrend) == true) {
	    print_interval(divstring,matches[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	  }
	}
	
	if (nflanking > 0) {
	  printf("====================\n");
	  if (sign != -1) {
	    for (i = 0; i < nrightflanks; i++) {
	      if (codingp == false || IIT_gene_overlapp(map_iit,rightflanks[i],chrstart,chrend) == true) {
		print_interval(divstring,rightflanks[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	      }
	    }
	  }
	  FREE(rightflanks);
	}
	
	FREE(divstring2);
	FREE(divstring);

      } else {
	/* Must have been a label */
#if 0
	if ((*iit = IIT_read(filename,/*name*/NULL,true,/*divread*/READ_NONE,/*divstring*/NULL,/*add_iit_p*/true)) == NULL) {
	}
#endif
	
	if (uniquep == true) {
	  divread = READ_ALL;	/* For subsequent search of overlapping genes */
	} else {
	  divread = READ_NONE;
	}
	if ((map_iit = IIT_read(iitfile,/*name*/NULL,true,divread,/*divstring*/NULL,
				/*add_iit_p*/true)) == NULL) {
	  fprintf(stderr,"Cannot open IIT file %s\n",iitfile);
	  exit(9);
	}
	matches = IIT_find(&nmatches,map_iit,argv[0]);
	debug(printf("Looking up label %s => %d matches\n",argv[0],nmatches));
	ndivs = IIT_ndivs(map_iit);

#if 0
	/* Not doing flanking */
	if (nflanking > 0) {
	  for (i = nleftflanks-1; i >= 0; i--) {
	    print_interval(/*divstring*/NULL,leftflanks[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	  }
	  printf("====================\n");
	  FREE(leftflanks);
	}
#endif

	for (i = 0; i < nmatches; i++) {
	  print_interval(/*divstring*/NULL,matches[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	}
	FREE(matches);

#if 0
	/* Not doing flanking */
	if (nflanking > 0) {
	  printf("====================\n");
	  for (i = 0; i < nrightflanks; i++) {
	    print_interval(/*divstring*/NULL,rightflanks[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	  }
	  FREE(rightflanks);
	}
#endif

      }
      FREE(segment);
      FREE(iitfile);
      FREE(mapdir);
    }

  } else {
    /* argv == 0.  Read from stdin */

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

      if ((map_iit = IIT_read(iitfile,/*name*/map_iitfile,/*readonlyp*/true,
			      /*divread*/READ_ALL,/*divstring*/NULL,/*add_iit_p*/false)) == NULL) {
	FREE(iitfile);
	iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				  strlen(map_iitfile)+strlen(".iit")+1,sizeof(char));
	sprintf(iitfile,"%s/%s.iit",genomesubdir,map_iitfile);
	if ((map_iit = IIT_read(iitfile,/*name*/map_iitfile,/*readonlyp*/true,
				/*divread*/READ_ALL,/*divstring*/NULL,/*add_iit_p*/false)) == NULL) {

	  fprintf(stderr,"Map file %s.iit not found in %s or %s.  Available files:\n",map_iitfile,mapdir,genomesubdir);
	  Datadir_list_directory(stderr,mapdir);
	  fprintf(stderr,"Either install file %s.iit or specify a full directory path\n",map_iitfile);
	  fprintf(stderr,"using the -M flag to gmap.\n");
	  exit(9);
	}

#if 0
      } else if (map_typestring != NULL) {
	map_type = IIT_typeint(map_iit,map_typestring);
	if (map_type < 0) {
	  /* Invalid user-specified strain.  Print nothing. */
	  fprintf(stderr,"No such type as %s.  Allowed strains are:\n",map_typestring);
	  IIT_dump_typestrings(stderr,map_iit);
	  exit(9);
	}
#endif
      } else {
	ndivs = IIT_ndivs(map_iit);
      }
    }
    FREE(iitfile);

    if (snps_root == NULL || print_snps_mode == 0) {
      genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,
			  chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
    } else if (print_snps_mode == 2) {
      genome = Genome_new(snpsdir,fileroot,snps_root,
			  chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
    } else if (print_snps_mode == 1 || print_snps_mode == 3) {
      genome = Genome_new(genomesubdir,fileroot,/*snps_root*/NULL,
			  chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
      genomealt = Genome_new(snpsdir,fileroot,snps_root,
			     chromosome_iit,genome_access,/*sharedp*/false,/*revcompp*/false);
    }
    Genome_setup(circular_typeint);


    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
			      strlen(".contig.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.contig.iit",genomesubdir,fileroot);
    contig_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false);
    FREE(iitfile);

    while ((line = Getline_wlength(&line_length,stdin)) != NULL) {
#if 0
      if ((p = rindex(Buffer,'\n')) != NULL) {
	*p = '\0';
      }
#endif

      coords = (char *) MALLOC((line_length+1)*sizeof(char));
      typestring = (char *) MALLOC((line_length+1)*sizeof(char));

      if (sscanf(line,"%s %s",coords,typestring) < 2) {
	typestring_ptr = (char *) NULL;
      } else {
	typestring_ptr = &(typestring[0]);
      }
      
      fprintf(stdout,"# Query: %s\n",coords);
      segment = (char *) NULL;
      if (force_label_p == false &&
	  Parserange_universal_iit(&segment,&revcomp,&genomicstart,&genomiclength,&chrstart,&chrend,
				   &chroffset,&chrlength,coords,chromosome_iit,contig_iit) == true) {
	debug(printf("Query %s parsed as: genomicstart = %llu, genomiclength = %u, revcomp = %d\n",
		     coords,(unsigned long long) genomicstart,genomiclength,revcomp));
	divstring = Univ_IIT_string_from_position(&chrstart,genomicstart,chromosome_iit);
	divstring2 = Univ_IIT_string_from_position(&chrend,genomicstart+genomiclength-1U,chromosome_iit);
	if (strcmp(divstring,divstring2)) {
	  fprintf(stderr,"Coordinates cross chromosomal boundary\n");
	  exit(9);
	} else {
	  chrstart += 1U;
	  chrend += 1U;
	  debug(printf("Query translated to %s:%u..%u\n",divstring,chrstart,chrend));
	}
	if (map_iit != NULL) {
	  matches = get_matches(&nmatches,&sign,divstring,chrstart,chrend,revcomp,
				&leftflanks,&nleftflanks,&rightflanks,&nrightflanks,
				typestring_ptr,&map_iit,/*filename*/NULL);
	  if (nflanking > 0) {
	    if (sign != +1) {
	      for (i = nleftflanks-1; i >= 0; i--) {
		print_interval(divstring,leftflanks[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	      }
	    }
	    printf("====================\n");
	    FREE(leftflanks);
	  }

	  for (i = 0; i < nmatches; i++) {
	    print_interval(divstring,matches[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	  }

	  if (nflanking > 0) {
	    printf("====================\n");
	    if (sign != -1) {
	      for (i = 0; i < nrightflanks; i++) {
		print_interval(divstring,rightflanks[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	      }
	    }
	    FREE(rightflanks);
	  }

	  fflush(stdout);
	  FREE(matches);
	  
	} else if (coordp == true) {
	  print_two_coords(genomicstart,genomiclength,chromosome_iit);
	  
	} else {
	  print_sequence(genome,genomealt,genomicstart,genomiclength,chromosome_iit,
			 /*whole_chromosome_p*/false,/*circular_part_p*/false);
	}

	FREE(divstring2);
	FREE(divstring);

      } else {
	/* Must have been a label */

	if (map_iit != NULL) {
	  matches = IIT_find(&nmatches,map_iit,line);
	  debug(printf("Looking up label %s => %d matches\n",line,nmatches));

	  /* Not doing flanking */

	  for (i = 0; i < nmatches; i++) {
	    print_interval(/*divstring*/NULL,matches[i],map_iit,ndivs,chromosome_iit,genome,fieldint);
	  }

	  /* Not doing flanking */

	  fflush(stdout);
	  FREE(matches);
	}
      }
      fprintf(stdout,"# End\n");

      FREE(segment);

      FREE(typestring);
      FREE(coords);
      FREE(line);
    }

    Univ_IIT_free(&contig_iit);

    if (map_iitfile != NULL) {
      FREE(iitfile);
      FREE(mapdir);
    }

  }
    
  if (genomealt != NULL) {
    Genome_free(&genomealt);
  }
  Genome_free(&genome);

  Univ_IIT_free(&chromosome_iit);

  if (map_iitfile != NULL) {
    IIT_free(&map_iit);
  }

  FREE(dbversion);
  FREE(genomesubdir);
  FREE(fileroot);
  FREE(dbroot);

  return 0;
}

