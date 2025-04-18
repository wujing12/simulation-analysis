/* -*- mode: c; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- 
 *
 * $Id: xdrfile_xtc.c,v 1.5 2009/05/18 09:06:38 spoel Exp $
 *
 * Copyright (c) Erik Lindahl, David van der Spoel 2003,2004.
 * Coordinate compression (c) by Frans van Hoesel. 
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 */
 
#include <stdlib.h>
#include "xdrfile.h"
#include "xdrfile_xtc.h"
	
#define MAGIC 1995

enum { FALSE, TRUE };

static int xtc_header(XDRFILE *xd,int *N_atom,int *step,float *time,mybool bRead)
{
	int result,magic,n=1;
	
	/* Note: read is same as write. He he he */
	magic  = MAGIC;
	if ((result = xdrfile_write_int(&magic,n,xd)) != n)
		{
			if (bRead)
				return exdrENDOFFILE;
			else
				return exdrINT;
		}
	if (magic != MAGIC)
		return exdrMAGIC;
	if ((result = xdrfile_write_int(N_atom,n,xd)) != n)
		return exdrINT;
	if ((result = xdrfile_write_int(step,n,xd)) != n)
		return exdrINT;
	if ((result = xdrfile_write_float(time,n,xd)) != n)
		return exdrFLOAT;
	
	return exdrOK;
}

static int xtc_coord(XDRFILE *xd,int *N_atom,matrix box,rvec *x,float *prec,
					 mybool bRead)
{
	int i,j,result;
    
	/* box */
	result = xdrfile_read_float(box[0],DIM*DIM,xd);
	if (DIM*DIM != result)
		return exdrFLOAT;
	else 
		{
			if (bRead)
				{
					result = xdrfile_decompress_coord_float(x[0],N_atom,prec,xd); 
					if (result != *N_atom)
						return exdr3DX;
				}
			else
				{
					result = xdrfile_compress_coord_float(x[0],*N_atom,*prec,xd); 
					if (result != *N_atom)
						return exdr3DX;
				}
		}
	return exdrOK;
}

int read_xtc_natoms(char *fn,int *N_atom)
{
	XDRFILE *xd;
	int step,result;
	float time;
	
	xd = xdrfile_open(fn,"r");
	if (NULL == xd)
		return exdrFILENOTFOUND;
	result = xtc_header(xd,N_atom,&step,&time,TRUE);
	xdrfile_close(xd);
	
	return result;
}

int read_xtc(XDRFILE *xd,
			 int N_atom,int *step,float *time,
			 matrix box,rvec *x,float *prec)
/* Read subsequent frames */
{
	int result;
  
	if ((result = xtc_header(xd,&N_atom,step,time,TRUE)) != exdrOK)
		return result;
	  
	if ((result = xtc_coord(xd,&N_atom,box,x,prec,1)) != exdrOK)
		return result;
  
	return exdrOK;
}

int write_xtc(XDRFILE *xd,
			  int N_atom,int step,float time,
			  matrix box,rvec *x,float prec)
/* Write a frame to xtc file */
{
	int result;
  
	if ((result = xtc_header(xd,&N_atom,&step,&time,FALSE)) != exdrOK)
		return result;

	if ((result = xtc_coord(xd,&N_atom,box,x,&prec,0)) != exdrOK)
		return result;
  
	return exdrOK;
}
