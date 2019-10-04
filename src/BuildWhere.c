/*************************************************************************

   Program:    KabatMan
   File:       BuildWhere.c
   
   Version:    V2.21
   Date:       13.07.00
   Function:   Database program for reading Kabat sequence files
   
   Copyright:  (c) UCL / Andrew C. R. Martin, UCL 1994-2000
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure and Modelling Unit,
               Department of Biochemistry and Molecular Biology,
               University College,
               Gower Street,
               London.
   Phone:      +44 (0) 1372 275775 (Home)
   EMail:      martin@biochem.ucl.ac.uk
               andrew@stagleys.demon.co.uk
               
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
   V2.6  16.03.95 Allows FORTRAN style comparison operators
   V2.7  16.05.95 Skipped
   V2.8  22.06.95 Increased word buffer size in BuildWhere()
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
   V2.21 13.07.00 Fixed long-standing bug in SetWhereData() which meant
                  you couldn't do a string comparison with anything 
                  containing a ' or "

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
/*>BOOL BuildWhere(char *buffer)
   -----------------------------
   Input:   char   *buffer        String containing where clause
   Returns: BOOL                  Success (FALSE indicates syntax 
                                  error or memory allocation error)
   Globals: WHERE  *gWhereClause  The where clause linked list

   Parses a where clause and creates a linked list of where conditions.
   The word `WHERE' within the input buffer is ignored.

   Currently the logic of the WHERE clause must be expressed in RPN.

   20.04.94 Original   By: ACRM
   22.06.95 Doubled length of word buffer
   23.06.95 Added missing return value
*/
BOOL BuildWhere(char *buffer)
{
   char  *pch = buffer,
         word[2*MAXBUFF];
   BOOL  error = FALSE;

   /* Step through the buffer pulling a word at a time out of the buffer*/
   do
   {
      pch=GetWord(pch,word);

      if(upstrcmp(word,"WHERE"))  /* If the word is not `WHERE'         */
      {
         if(!CheckForSetOper(word, &error))
         {
            if(error) return(FALSE);
            
            pch = HandleWhereSubClause(pch,word,&error);
            if(error) return(FALSE);
         }
      }  /* This was not `WHERE'                                        */
   }  while(pch);  /* Step through words in buffer                      */

   return(TRUE);
}

/************************************************************************/
/*>BOOL CheckForSetOper(char *word, BOOL *error)
   ---------------------------------------------
   Input:   char  *word     The word which may be a set operator
   Output:  BOOL  *error    TRUE if an error occurred
   Returns: BOOL            TRUE if the word was a set operator

   Checks to see if a word was a set operator. If so, adds an item to the
   WHERE clause linked list containing this set operator.

   20.04.94 Original    By: ACRM
*/
BOOL CheckForSetOper(char *word, BOOL *error)
{
   int   i;

   *error = FALSE;

   /* Run through the list of logical set operators                     */
   for(i=0; gSetOper[i].type; i++)
   {
      if(!upstrncmp(word,gSetOper[i].name,gSetOper[i].length))
      {
         /* We've got a match in our array of fields, so allocate 
            the next entry in our linked list of selections.
         */
         if(gWhereClause == NULL)
         {
            INIT(gWhereClause,WHERE);
            gCurrentWhere = gWhereClause;
         }
         else
         {
            ALLOCNEXT(gCurrentWhere,WHERE);
         }
         
         /* Check allocation was OK                                     */
         if(gCurrentWhere==NULL)
         {
            fprintf(stderr,"Error: No memory to build where clause\n");
            *error = TRUE;
            return(FALSE);
         }
         
         /* Fill in the field type for this selection sub-clause        */
         gCurrentWhere->SetOper = TRUE;
         gCurrentWhere->type    = gSetOper[i].type;

         return(TRUE);
      }  /* Got a match                                                 */
   }

   return(FALSE);
}

/************************************************************************/
/*>char *HandleWhereSubClause(char *buffer, char *word, BOOL *error)
   -----------------------------------------------------------------
   Input:   char  *buffer      Current start of next word in buffer
            char  *word        The current word to be tested (modified on
                               exit)
   Output:  BOOL  *error       TRUE if an error occurred
   Returns: char  *            Pointer to start of next word in buffer 
                               after processing this sub-clause
   Globals: WHERE gWhereClause Linked list of where sub-clauses and 
                               actions

   Checks to see if a word is a valid field type. If so, reads the 
   comparison operator and data to be compared from the input buffer and
   places the data into the WHERE linked list

   20.04.94 Original    By: ACRM
*/
char *HandleWhereSubClause(char *buffer, char *word, BOOL *error)
{
   int   i;
   char  *pch = buffer;

   *error = FALSE;
   
   /* Run through the field list to see if we get a match               */
   for(i=0; gField[i].type; i++)
   {
      if(!upstrncmp(word,gField[i].name,gField[i].length))
      {
         /* We've got a match in our array of fields, so allocate 
            the next entry in our linked list of selections.
         */
         if(gWhereClause == NULL)
         {
            INIT(gWhereClause,WHERE);
            gCurrentWhere = gWhereClause;
         }
         else
         {
            ALLOCNEXT(gCurrentWhere,WHERE);
         }
         
         /* Check allocation was OK                                     */
         if(gCurrentWhere==NULL)
         {
            fprintf(stderr,"Error: No memory to build where clause\n");
            *error = TRUE;
            return(NULL);
         }
         
         /* Fill in the field type for this selection sub-clause        */
         gCurrentWhere->SetOper = FALSE;
         gCurrentWhere->type    = gField[i].type;

         /* See if there's a parameter to fill in                       */
         FillParameter(gCurrentWhere->param, word);
         
         /* Now get the comparison to be performed                      */
         pch = GetWord(pch,word);
         if(!word[0] || !SetComparison(gCurrentWhere,word))
         {
            *error = TRUE;
            return(FALSE);
         }
                
         /* Now get the data for comparison                             */
         pch = GetWord(pch,word);

         if(!word[0])
         {
            *error = TRUE;
            return(NULL);
         }

         SetWhereData(gCurrentWhere,word);

         /* Return the new pointer to the next word                     */
         return(pch);
      }  /* Found field in FIELD array                                  */
   }  /* Step through FIELD array                                       */

   /* If we get here, we failed to find this item so return error       */
   *error = TRUE;

   return(NULL);
}
   
/************************************************************************/
/*>BOOL SetComparison(WHERE *p, char *word)
   ----------------------------------------
   Input:   char  *word     Data comparison symbol
   Output:  WHERE *p        p->comparison set to appropriate type
   Returns: BOOL            Success?

   Sets the comparison field of the structure according to the given
   symbol. = or == are accepted for equal and != or <> are accepted for
   not equal.

   20.04.94 Original    By: ACRM
   21.04.94 Added INC to substring options
   16.03.95 Added FORTRAN style (EQ, NE, etc.) comparisons
*/
BOOL SetComparison(WHERE *p, char *word)
{
   UPPER(word);
   
   if(!strcmp(word, "=")  || 
      !strcmp(word, "==") || 
      !strcmp(word, "EQ"))
   {
      p->comparison = COMP_EQ;
      return(TRUE);
   }
   else if(!strcmp(word, "!=") || 
           !strcmp(word, "<>") || 
           !strcmp(word, "NE"))
   {
      p->comparison = COMP_NE;
      return(TRUE);
   }
   else if(!strcmp(word, "<") ||
           !strcmp(word, "LT"))
   {
      p->comparison = COMP_LT;
      return(TRUE);
   }
   else if(!strcmp(word, "<=") ||
           !strcmp(word, "LE"))
   {
      p->comparison = COMP_LE;
      return(TRUE);
   }
   else if(!strcmp(word, ">") ||
           !strcmp(word, "GT"))
   {
      p->comparison = COMP_GT;
      return(TRUE);
   }
   else if(!strcmp(word, ">=") ||
           !strcmp(word, "GE"))
   {
      p->comparison = COMP_GE;
      return(TRUE);
   }
   else if(!strncmp(word, "CONT",4) || 
           !strncmp(word, "LIKE",4) || 
           !strncmp(word, "SIM",3)  || 
           !strncmp(word, "INC",3)  || 
           !strncmp(word, "SUB",3))
   {
      p->comparison = COMP_SIM;
      return(TRUE);
   }

   return(FALSE);
}

/************************************************************************/
/*>void SetWhereData(WHERE *wh, char *word)
   ----------------------------------------
   Input:   char  *word    The data to copy into the WHERE item
   Output:  WHERE *wh      wh->data will contain the word

   Copies data into the where item data field, removing any leading or
   trailing ' or "

   20.04.94 Original    By: ACRM
   13.07.00 The check for trailing ' or " was actually stopping the string
            at the first occurrence rather than the last. Fixed this!
*/
void SetWhereData(WHERE *wh, char *word)
{  
   char *p;
   int  len;
   
   /* If the first character is a ' or " then skip it and end the string at
      the last ' or "
   */
   if(word[0] == '"')
   {
      strcpy(wh->data, word+1);
      len = strlen(wh->data);
      for(p=wh->data + len-1; p>=wh->data; p--)
      {
         if(*p=='"')
         {
            *p = '\0';
            break;
         }
      }
   }
   else if(word[0] == '\'')
   {
      strcpy(wh->data, word+1);
      len = strlen(wh->data);
      for(p=wh->data + len-1; p>=wh->data; p--)
      {
         if(*p=='\'')
         {
            *p = '\0';
            break;
         }
      }
   }
   else
   {
      strcpy(wh->data, word);
   }
}

/************************************************************************/
/*>void ClearWhere(void)
   ---------------------
   Clears the current where statement

   21.04.94 Original    By: ACRM
*/
void ClearWhere(void)
{
   DATA *d;
   
   if(gWhereClause != NULL)
   {
      FREELIST(gWhereClause,WHERE);
      gWhereClause = gCurrentWhere = NULL;

      for(d=gData; d!=NULL; NEXT(d))
         d->active[0] = FALSE;
   }
}
