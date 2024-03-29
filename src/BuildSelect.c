/*************************************************************************

   Program:    KabatMan
   File:       BuildSelect.c
   
   Version:    V2.26
   Date:       04.10.19
   Function:   Database program for reading Kabat sequence files
   
   Copyright:  (c) UCL / Andrew C. R. Martin, UCL 1994-2019
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure and Modelling Unit,
               Department of Biochemistry and Molecular Biology,
               University College,
               Gower Street,
               London.
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is copyright. 

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
   V0.1  20.04.94 Development version
   V1.0  27.04.94 Original
   V1.1           Skipped
   V1.2  13.05.94 Skipped
   V2.0  30.06.94 Skipped
   V2.1  11.07.94 Skipped
   V2.2  21.07.94 Skipped
   V2.3  23.01.95 Skipped
   V2.4  10.02.95 Skipped
   V2.5  07.03.95 Skipped
   V2.6  16.03.95 Skipped
   V2.7  16.05.95 Skipped
   V2.8  22.06.95 Skipped
   V2.9  23.06.95 Clean compile under gcc -Wall
   V2.10 27.06.95 Skipped
   V2.11 15.12.95 Skipped
   V2.12 02.04.96 Skipped
   V2.13 11.04.96 Skipped
   V2.14 18.04.96 Skipped
   V2.15 22.04.96 Skipped
   V2.16 08.05.96 Skipped
   V2.17 29.05.96 Skipped
   V2.18 10.09.97 Skipped
   V2.19 14.10.98 Skipped
   V2.20 xx.xx.xx Skipped
   V2.21 13.07.00 Skipped
   V2.22 31.07.00 Skipped
   V2.23 03.04.02 Skipped
   V2.24 28.02.05 blGetWord() takes extra parameter
   V2.25 24.08.06 Skipped
   V2.26 04.10.19 Changed all bioplib calls to blXXX()

*************************************************************************/
/* Includes
*/
#include "kabatman.h"

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
#include "protos.h"
/************************************************************************/
/*>BOOL BuildSelect(char *buffer)
   ------------------------------
   Input:   char      *buffer         String containing select clause
   Returns: BOOL                      Success (FALSE indicates syntax 
                                      error or memory allocation error)
   Globals: SELECTION *gSelectClause  The selection linked list

   Parses a select clause and creates a linked list of selection actions.
   The word `SELECT' within the input buffer is ignored.

   20.04.94 Original   By: ACRM
   23.06.95 Added missing return value
   28.02.05 Added word length parameter for blGetWord()
*/
BOOL BuildSelect(char *buffer)
{
   char      *pch,
             word[MAXBUFF];
   int       i;

   pch=buffer;

   /* Step through the buffer pulling a word at a time out of the buffer*/
   do
   {
      pch=blGetWord(pch,word,MAXBUFF);

      if(blUpstrcmp(word,"SELECT")) /* If the word is not `SELECT'        */
      {
         /* Run through the field list to see if we get a match         */
         for(i=0; gField[i].type; i++)
         {
            if(!blUpstrncmp(word,gField[i].name,gField[i].length))
            {
               /* We've got a match in our array of fields, so allocate 
                  the next entry in our linked list of selections.
               */
               if(gSelectClause == NULL)
               {
                  INIT(gSelectClause,SELECTION);
                  gCurrentSelect = gSelectClause;
               }
               else
               {
                  ALLOCNEXT(gCurrentSelect,SELECTION);
               }

               /* Check allocation was OK                               */
               if(gCurrentSelect==NULL)
               {
                  fprintf(stderr,"Error: No memory to build selection\n");
                  return(FALSE);
               }

               /* Fill in the field type for this selection sub-clause  */
               gCurrentSelect->type = gField[i].type;

               /* Now see if there is a parameter specified             */
               FillParameter(gCurrentSelect->param, word);

               /* Break out of FIELD array                              */
               break;
            }  /* Found field in FIELD array                            */
         }  /* Step through FIELD array                                 */

         /* Check to see if we failed to find this item; if so, return
            FALSE
         */
         if(!gField[i].type)
            return(FALSE);
      }  /* This was not `SELECT'                                       */
   }  while(pch);  /* Step through words in buffer                      */

   return(TRUE);
}

/************************************************************************/
/*>void FillParameter(char *param, char *word)
   -------------------------------------------
   Input:   char *word    The word which may contain a parameter to copy
   Output:  char *param   The copied parameter or a blank string

   Copies a parameter (if one exists) out of word into parameter. 
   Thus if word = TEST(PARAM), then param will be PARAM. If no parameter
   is specified, param will contain a blank string.

   20.04.94 Original    By: ACRM
*/
void FillParameter(char *param, char *word)
{
   char *pw1,
        *pw2;

   param[0] = '\0';
   
   if((pw1 = strchr(word,'('))!=NULL)             /* See if there's a ( */
   {
      if(*(++pw1))                                /* Character after (  */
      {
         if((pw2 = strchr(pw1,')'))!=NULL)
            *pw2 = '\0';                          /* Terminate at )     */
         
         strcpy(param,pw1);
      }
   }
}
   
/************************************************************************/
/*>void ClearSelect(void)
   ----------------------
   Clears the current selection

   21.04.94 Original    By: ACRM
*/
void ClearSelect(void)
{
   if(gSelectClause != NULL)
   {
      FREELIST(gSelectClause,SELECTION);
      gSelectClause = gCurrentSelect = NULL;
   }
}
