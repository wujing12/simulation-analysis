/* -*- mode: c; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- 
 *
 * $Id: xdrfile_xtc.h,v 1.5 2009/05/18 09:06:38 spoel Exp $
 *
 * Copyright (c) Erik Lindahl, David van der Spoel 2003,2004.
 * Coordinate compression (c) by Frans van Hoesel. 
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 */

#ifndef _xdrfile_xtc_h
#define _xdrfile_xtc_h

#ifdef CPLUSPLUS
extern "C" {
#endif

#include "xdrfile.h"
  
  /* All physics return exdrOK if succesfull. 
   * (error codes defined in xdrfile.h).
   */  
   
  /* This physics returns the number of atoms in the xtc file in *N_atom */
  extern int read_xtc_natoms(char *fn,int *N_atom);
  
  /* Read one frame of an open xtc file */
  extern int read_xtc(XDRFILE *xd,int N_atom,int *step,float *time,
		      matrix box,rvec *x,float *prec);
  
  /* Write a frame to xtc file */
  extern int write_xtc(XDRFILE *xd,
		       int N_atom,int step,float time,
		       matrix box,rvec *x,float prec);
  
#ifdef CPLUSPLUS
}
#endif

#endif
