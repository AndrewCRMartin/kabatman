/*************************************************************************

   Program:    KabatMan
   File:       kabatman.h
   
   Version:    V2.22
   Date:       31.07.00
   Function:   Database program for reading Kabat sequence files
   
   Copyright:  (c) UCL / Andrew C. R. Martin 1994-2000
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
   V0.1  12.04.94 Development version
   V1.0  27.04.94 Original
   V1.1  11.05.94 Skipped
   V1.2  11.05.94 Skipped
   V2.0  30.06.94 Skipped
   V2.1  11.07.94 Skipped
   V2.2  21.07.94 Added gOldFormat
   V2.3  23.01.95 Skipped
   V2.4  10.02.95 Skipped
   V2.5  07.03.95 Skipped
   V2.6  16.03.95 Skipped
   V2.7  16.05.95 Skipped
   V2.8  22.06.95 Increased sizes of buffers in SELECTION and WHERE 
                  structures
   V2.9  23.06.95 Skipped
   V2.10 27.06.95 Skipped
   V2.11 15.12.95 Skipped
   V2.12 02.04.96 Added idlight and idheavy to DATA.
                  Added ID***** and URL***** field types
                  Added URLFORMAT define
                  Added gHTML support
   V2.13 11.04.96 Added gFileDate
   V2.14 18.04.96 Added gURLFormat
   V2.15 22.04.96 Increased MAXCHOTHRES from 20
   V2.16 07.05.96 Added gCanonChothNum
   V2.17 29.05.96 Skipped
   V2.18 09.09.97 Added FIELD_SUBGROUP
   V2.19 14.10.98 Added gDelim
   V2.20 xx.xx.xx Skipped
   V2.21 13.07.00 Skipped
   V2.22 31.07.00 Added LOOP definitions for Contact CDR definitions

*************************************************************************/
#ifndef _KABATMAN_H
#define _KABATMAN_H

/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <time.h>

#include "bioplib/macros.h"
#include "bioplib/SysDefs.h"
#include "bioplib/general.h"

#include "RdKabat.h"

/************************************************************************/
/* Defines and macros
*/
#define URLFORMAT \
"<a href=http://immuno.bme.nwu.edu/scripts/noninter.tcl?qt=%s>%s</a>"  
                                 /* Format for URL to get entry         */
#define MAXBUFF      160         /* General pupose buffer size          */
#define DEF_FOF      "kabat.fof" /* Default file of files               */
#define DEF_KABAT    "kabat.dat" /* Stored Kabat data                   */
#define DEF_CHOTHIA  "chothia.dat"  /* Default Chothia data file        */
#define DEF_INFO     1           /* Default info level                  */
#define DEF_VARIABILITY (REAL)0.0   /* Default variability              */
#define MAXLFILES    3           /* Max number of LC files per HC       */
#define MINSEQ       75          /* Min sequence size to bother keeping */
#define STACKDEPTH   10          /* Max set operation stack depth       */
#define ENV_KABATDIR "KABATDIR"  /* Environment variable for Kabat      */
                                 /* directory                           */
#define MAXCHOTHRES  40          /* Max number of key residues per class*/

#define FIELD_NAME       1        /* These types define the items which */
#define FIELD_ANTIGEN    2        /* may be specified in SELECT and     */
#define FIELD_L1         3        /* WHERE clauses.                     */
#define FIELD_L2         4
#define FIELD_L3         5
#define FIELD_H1         6
#define FIELD_H2         7
#define FIELD_H3         8
#define FIELD_CLASS      9
#define FIELD_SOURCE    10
#define FIELD_REF       11
#define FIELD_LENGTH    12
#define FIELD_RES       13
#define FIELD_PIR       14
#define FIELD_COMPLETE  15
#define FIELD_VAR       16
#define FIELD_LIGHT     17
#define FIELD_HEAVY     18
#define FIELD_CANONICAL 19
#define FIELD_IDLIGHT   20
#define FIELD_IDHEAVY   21
#define FIELD_URLLIGHT  22
#define FIELD_URLHEAVY  23
#define FIELD_SUBGROUP  24

#define OPER_AND        1        /* Types for logical set operators     */
#define OPER_OR         2
#define OPER_NOT        3

#define COMP_EQ         1        /* These define the comparisons to be  */
#define COMP_NE         2        /* made within each part of the WHERE  */
#define COMP_LT         3        /* clause                              */
#define COMP_GT         4
#define COMP_LE         5
#define COMP_GE         6
#define COMP_SIM        7

#define LOOP_KABAT      1        /* Loop definitions                    */
#define LOOP_ABM        2
#define LOOP_CHOTHIA    3
#define LOOP_CONTACT    4

#define CLASS_LAMBDA    1        /* Light chain classes                 */
#define CLASS_KAPPA     2

/* A linked list of DATA structures is used to store the actual Kabat
   data
*/
typedef struct _data
{
   struct _data *next;
   BOOL         active[STACKDEPTH];
   char         **LNumbers,
                **HNumbers;
   char         antigen[LARGEBUFF],
                class[SMALLBUFF],
                name[SMALLBUFF],
                source[SMALLBUFF],
                fsource[LARGEBUFF],
                reference[LARGEBUFF],
                light[LARGEBUFF],
                heavy[LARGEBUFF],
                idlight[SMALLBUFF],
                idheavy[SMALLBUFF];
}  DATA;
   
/* An array of FIELD structures links each of the field strings to a 
   FIELD_xxxx definition
*/
typedef struct
{
   int  type,       /* Field type number                                */
        length;     /* Number of characters to compare                  */
   char *name;      /* Text for comparison                              */
}  FIELD;

/* A linked list of SELECTION structures is used to store the elements
   of the SELECT statement
*/
typedef struct _selection
{
   struct _selection *next;
   int               type;
   char              param[MAXBUFF];
}  SELECTION;

/* A linked list of WHERE structures is used to store the elements of
   the WHERE statement
*/
typedef struct _where
{
   struct _where *next;
   int           type,
                 logic,
                 comparison;
   BOOL          SetOper;
   char          param[MAXBUFF],
                 data[MAXBUFF*2];
}  WHERE;

/* An array of LOOP structures is used to store the alternative loop
   defintions
*/
typedef struct
{
   char *name,                 /* The loop name (L1,...)                */
        *KabatS,   *KabatE,    /* Kabat defintition (e.g. L24 L34)      */
        *AbMS,     *AbME,      /* AbM definition                        */
        *ChothiaS, *ChothiaE,  /* Chothia definition                    */
        *ContactS, *ContactE;  /* Contact definitions                   */
}  LOOP;

typedef struct _chothia
{
   struct _chothia *next;
   int             length;
   char            LoopID[8],
                   class[8],
                   resnum[MAXCHOTHRES][8],
                   restype[MAXCHOTHRES][24];
}  CHOTHIA;

/************************************************************************/
/* Globals
*/
#ifdef MAIN        /*-------------- Global definitions -----------------*/
char  **gFlagList = NULL,
      gFOF[MAXBUFF],
      gKabatFile[MAXBUFF],
      gChothiaFile[MAXBUFF],
      gFileDate[MAXBUFF],
      gURLFormat[MAXBUFF],
      gDelim      = ',';
DATA  *gData      = NULL;                   /* Kabat data linked list   */
int   gInfoLevel  = DEF_INFO,               /* Information level        */
      gLoopMode   = LOOP_KABAT;             /* Loop definition mode     */
BOOL  gShowInserts= FALSE,                  /* Show inserts in loops?   */
      gCanonChothNum = FALSE;               /* Use Chothia numbering in
                                               canonical files?         */
FIELD gField[]    =                         /* Link field names/numbers */
{  {  FIELD_NAME,      4, "NAME"},
   {  FIELD_ANTIGEN,   7, "ANTIGEN"},
   {  FIELD_L1,        2, "L1"},
   {  FIELD_L2,        2, "L2"},
   {  FIELD_L3,        2, "L3"},
   {  FIELD_H1,        2, "H1"},
   {  FIELD_H2,        2, "H2"},
   {  FIELD_H3,        2, "H3"},
   {  FIELD_CLASS,     5, "CLASS"},
   {  FIELD_SOURCE,    6, "SOURCE"},
   {  FIELD_REF,       3, "REFERENCE"},
   {  FIELD_LENGTH,    3, "LENGTH"},
   {  FIELD_RES,       3, "RESIDUE"},
   {  FIELD_PIR,       3, "PIRFILE"},
   {  FIELD_COMPLETE,  8, "COMPLETE"},
   {  FIELD_VAR,       3, "VARIABILITY"},
   {  FIELD_LIGHT,     5, "LIGHT"},
   {  FIELD_HEAVY,     5, "HEAVY"},
   {  FIELD_CANONICAL, 3, "CANONICAL"},
   {  FIELD_IDLIGHT,   3, "IDLIGHT"},
   {  FIELD_IDHEAVY,   3, "IDHEAVY"},
   {  FIELD_URLLIGHT,  4, "URLLIGHT"},
   {  FIELD_URLHEAVY,  4, "URLHEAVY"},
   {  FIELD_SUBGROUP,  4, "SUBGROUP"},
   {  0,               0, NULL}
}  ;
FIELD gSetOper[]  =                         /* Link set oper names/nums */
{  {  OPER_AND,        3, "AND"},
   {  OPER_OR,         2, "OR"},
   {  OPER_NOT,        3, "NOT"},
   {  0,               0, NULL}
}  ;
LOOP gLoopDefs[]   = 
{  {  "L1", "L24", "L34",  "L24", "L34",  "L24", "L34",  "L30", "L36"},
   {  "L2", "L50", "L56",  "L50", "L56",  "L50", "L56",  "L46", "L55"},
   {  "L3", "L89", "L97",  "L89", "L97",  "L89", "L97",  "L89", "L96"},
   {  "H1", "H31", "H35B", "H26", "H35B", "H26", "H32",  "H30", "H35B"},
   {  "H2", "H50", "H65",  "H50", "H58",  "H52", "H56",  "H47", "H58"},
   {  "H3", "H95", "H102", "H95", "H102", "H95", "H102", "H93", "H101"},
   {  NULL, NULL , NULL  , NULL , NULL  , NULL , NULL  }
}  ;
   
SELECTION *gSelectClause  = NULL,           /* SELECT statement list    */
          *gCurrentSelect = NULL;
WHERE     *gWhereClause   = NULL,           /* WHERE statement list     */
          *gCurrentWhere  = NULL;
CHOTHIA   *gChothia       = NULL;           /* CHOTHIA class types      */
BOOL      gOldFormat      = FALSE;          /* Old Kabat format         */
REAL      gVariability    = 0.0;            /* Variability              */
BOOL      gHTML           = FALSE;          /* HTML Output format       */

#else              /*------------- External  references ----------------*/
extern char      **gFlagList,
                 gFOF[MAXBUFF],
                 gKabatFile[MAXBUFF],
                 gChothiaFile[MAXBUFF],
                 gFileDate[MAXBUFF],
                 gURLFormat[MAXBUFF],
                 gDelim;
extern DATA      *gData;
extern int       gInfoLevel,
                 gLoopMode;
extern BOOL      gShowInserts,
                 gCanonChothNum;
extern FIELD     gField[],
                 gSetOper[];
extern LOOP      gLoopDefs[];
extern SELECTION *gSelectClause,
                 *gCurrentSelect;
extern WHERE     *gWhereClause,
                 *gCurrentWhere;
extern CHOTHIA   *gChothia;
extern BOOL      gOldFormat;
extern REAL      gVariability;
extern BOOL      gHTML;

#endif             /*-------------- End of global data -----------------*/


#endif /* _KABATMAN_H                                                   */


