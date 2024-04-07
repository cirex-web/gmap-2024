/* $Id: littleendian.h 223349 2020-10-28 02:49:25Z twu $ */
#ifndef LITTLEENDIAN_INCLUDED
#define LITTLEENDIAN_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"		/* For HAVE_64_BIT */
#endif

#include <stdio.h>
#include "types.h"

extern void
Littleendian_write_char (unsigned char value, int fd);
extern void
Littleendian_write_uint (UINT4 value, int fd);
#ifdef HAVE_64_BIT
extern void
Littleendian_write_uint8 (UINT8 value, int fd);
#endif


#define FREAD_CHAR(p,fp) fread(p,sizeof(unsigned char),1,fp)
#define FREAD_USHORTS(a,n,fp) fread(a,sizeof(UINT2),n,fp)
#define FREAD_INT(p,fp) fread(p,sizeof(int),1,fp)
#define FREAD_UINT(p,fp) fread(p,sizeof(UINT4),1,fp)
#define FREAD_INTS(a,n,fp) fread(a,sizeof(int),n,fp)
#define FREAD_UINTS(a,n,fp) fread(a,sizeof(UINT4),n,fp)
#define FREAD_FLOATS(a,n,fp) fread(a,sizeof(float),n,fp)
#ifdef HAVE_64_BIT
#define FREAD_UINT8(p,fp) fread(p,sizeof(UINT8),1,fp)
#define FREAD_UINT8S(a,n,fp) fread(a,sizeof(UINT8),n,fp)
#endif

#define FWRITE_CHAR(x,fp) fwrite(&(x),sizeof(unsigned char),1,fp)
#define FWRITE_USHORT(x,fp) fwrite(&(x),sizeof(UINT2),1,fp)
#define FWRITE_INT(x,fp) fwrite(&(x),sizeof(int),1,fp)
#define FWRITE_UINT(x,fp) fwrite(&(x),sizeof(UINT4),1,fp)
#define FWRITE_DOUBLE(x,fp) fwrite(&(x),sizeof(double),1,fp)
#define WRITE_CHAR(x,fd) Littleendian_write_char(x,fd)
#define WRITE_UINT(x,fd) Littleendian_write_uint(x,fd)
#ifdef HAVE_64_BIT
#define WRITE_UINT8(x,fd) Littleendian_write_uint8(x,fd)
#define WRITE_UINT8_AS_UINT(x,fd) Littleendian_write_uint8_as_uint(x,fd)
#endif

#define FWRITE_USHORTS(a,n,fp) fwrite(a,sizeof(UINT2),n,fp)
#define FWRITE_INTS(a,n,fp) fwrite(a,sizeof(int),n,fp)
#define FWRITE_UINTS(a,n,fp) fwrite(a,sizeof(UINT4),n,fp)
#define FWRITE_FLOATS(a,n,fp) fwrite(a,sizeof(float),n,fp)

#ifdef HAVE_64_BIT
#define FWRITE_UINT8(x,fp) fwrite(&(x),sizeof(UINT8),1,fp)
#define FWRITE_UINT8S(a,n,fp) fwrite(a,sizeof(UINT8),n,fp)
#endif

#define FREAD_CHARS(a,n,fp) fread(a,sizeof(char),n,fp)
#define FWRITE_CHARS(a,n,fp) fwrite(a,sizeof(char),n,fp)

#endif

