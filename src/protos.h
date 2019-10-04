/*************************************************************************

   Program:    KabatMan
   File:       protos.h
   
   Version:    V2.19
   Date:       14.10.98
   Function:   Include all prototype files
   
   Copyright:  (c) UCL / Andrew C. R. Martin, UCL 1994-8
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure and Modelling Unit,
               Department of Biochemistry and Molecular Biology,
               University College,
               Gower Street,
               London.
   EMail:      martin@biochem.ucl.ac.uk
               andrew@stagleys.demon.co.uk
               
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
   V0.1  13.04.94 Development version
   V1.0  27.04.94 Original
   V1.1  11.05.94 libroutine.p only included if bioplib not available
   V1.2  11.05.94 Skipped
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
   V2.12 02.04.96 Skipped
   V2.13 11.04.96 Skipped
   V2.14 18.04.96 Skipped
   V2.15 22.04.96 Skipped
   V2.16 07.05.96 Added KabCho.p
   V2.17 29.05.96 Skipped
   V2.18 10.09.97 Added subgroup.p
   V2.19 14.10.98 Skipped

*************************************************************************/
/* Includes
*/
#include "kabatman.p"
#include "RdKabat.p"
#include "BuildSelect.p"
#include "BuildWhere.p"
#include "ExecSearch.p"
#include "KabCho.p"
#include "subgroup.p"

#ifdef NOBIOPLIB
#include "libroutines.p"
#endif
