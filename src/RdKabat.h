/*************************************************************************

   Program:    KabatMan
   File:       RdKabat.h
   
   Version:    V2.18
   Date:       10.09.97
   Function:   Include file for using ReadNextKabatEntry()
   
   Copyright:  (c) Andrew C. R. Martin, UCL 1994-7
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure and Modelling Unit,
               Department of Biochemistry and Molecular Biology,
               University College,
               Gower Street,
               London.
   EMail:      martin@biochem.ucl.ac.uk
               
**************************************************************************

   This program is copright. 

   Any copying without the express permission of the author is illegal.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V0.1  12.04.94 Development version
   V0.2  18.04.94 Added insertsnuc to structure
   V1.0  27.04.94 Original
   V1.1  11.05.94 Changes for routines moved into bioplib
   V2.0  30.06.94 Skipped
   V2.1  11.07.94 Skipped
   V2.2  21.07.94 Skipped
   V2.3  23.01.95 Skipped
   V2.4  10.02.95 Skipped
   V2.5  07.03.95 Skipped
   V2.6  16.03.95 Skipped
   V2.7  16.05.95 Skipped
   V2.8  22.06.95 Skipped
   V2.9  23.06.95 Skipped
   V2.10 27.06.95 Skipped
   V2.11 15.12.95 Skipped
   V2.12 02.04.96 Added kadbid to KABATENTRY
   V2.13 11.04.96 Skipped
   V2.14 18.04.96 Skipped
   V2.15 22.04.96 Skipped
   V2.16 08.05.96 Skipped
   V2.17 29.05.96 Skipped
   V2.18 10.09.97 Skipped

*************************************************************************/
#ifndef _RDKABAT_H
#define _RDKABAT_H

/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "bioplib/SysDefs.h"
#include "bioplib/macros.h"
#include "bioplib/seq.h"

#ifdef NOBIOPLIB
#include "libroutines.p"
#endif

/************************************************************************/
/* Defines
*/
#define SMALLBUFF     40
#define LARGEBUFF    320
#define MAXSEQ       800
#define SEQBUFF     1600
#define MAXKABATSEQ  200

/************************************************************************/
/* Structure definitions
*/
typedef struct
{
   char aatable[SMALLBUFF],
        aaname[SMALLBUFF],
        codname[SMALLBUFF],
        reference[LARGEBUFF],
        antigen[LARGEBUFF],
        species[SMALLBUFF],
        class[SMALLBUFF],
        strain[SMALLBUFF],
        source[LARGEBUFF],
        insertsaa[LARGEBUFF],
        insertsnuc[LARGEBUFF],
        notesaa[LARGEBUFF],
        kabatnum[SEQBUFF],
        kabatseq[SEQBUFF],
        sequence[MAXSEQ],
        kadbid[SMALLBUFF];
}  KABATENTRY;

#endif /* _RDKABAT_H                                                    */