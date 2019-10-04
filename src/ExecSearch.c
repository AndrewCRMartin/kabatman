/*************************************************************************

   Program:    KabatMan
   File:       ExecSearch.c
   
   Version:    V2.18
   Date:       09.09.97
   Function:   Database program for reading Kabat sequence files
   
   Copyright:  (c) Andrew C. R. Martin 1994-7
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure and Modelling Unit,
               Department of Biochemistry and Molecular Biology,
               University College,
               Gower Street,
               London.
   Phone:      +44 (0) 1372 275775 (Home)
   EMail:      martin@biochem.ucl.ac.uk
               
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
   V0.1  12.04.94 Development version
   V1.0  27.04.94 Original
   V1.1           Skipped
   V1.2  13.05.94 Added Chothia canonical function
   V2.0  30.06.94 Skipped
   V2.1  11.07.94 Skipped
   V2.2  21.07.94 Skipped
   V2.3  23.01.95 Added variability support
   V2.4  10.02.95 Skipped
   V2.5  07.03.95 Skipped
   V2.6  16.03.95 Corrected handling of PIR option when file name is not
                  specified.
   V2.7  16.05.95 Skipped
   V2.8  22.06.95 Added optimisations to fuzzystrstr when given blank
                  string(s) or a substring longer than the actual string
   V2.9  23.06.95 Clean compile under gcc -Wall
   V2.10 27.06.95 Bug fix in fuzzystrstr()
   V2.11 15.12.95 Changed to use QueryStrStr() rather than strstr() to
                  allow single char wildcard matches
   V2.12 02.04.96 Added ID and URL handling
   V2.13 11.04.96 Writes dataset date with number of hits
   V2.14 18.04.96 Takes the URL from a variable 
   V2.15 22.04.96 Skipped
   V2.16 08.05.96 FindCanonical() modified to allow Chothia numbering in
                  the Chothia data file
   V2.17 29.05.96 Skipped
   V2.18 09.09.97 Added subgroup handling

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
/*>BOOL ExecuteSearch(char *filename)
   ----------------------------------
   Input:   char *filename   Filename for output or NULL for stdout
   Returns: BOOL             Success

   Actually execute the search. Steps through the WHERE clause linked list
   calling HandleMatch() or HandleLogical() as appropriate to add data to
   the stack or to perform a logical set operation on the current stack.

   20.04.94 Original    By: ACRM
   21.04.94 Added filename parameter
   26.04.94 Prints error message if stack depth wrong.
   23.06.95 Added missing return value
*/
BOOL ExecuteSearch(char *filename)
{
   int   StackDepth = 0;
   WHERE *wh;
   FILE  *fp = stdout;

   for(wh=gWhereClause; wh!=NULL; NEXT(wh))
   {
      if(wh->SetOper)        /* This is a logical operator              */
      {
         if(!HandleLogical(wh,&StackDepth))
            return(FALSE);
      }
      else                   /* This is a standard comparison           */
      {
         if(!HandleMatch(wh,&StackDepth))
            return(FALSE);
      }
   }

   if(StackDepth != 1)
   {
      fprintf(stderr,"Error: Stack depth (%d) should be 1\n",StackDepth);
      return(FALSE);
   }

   if(filename != NULL)
   {
      if((fp=fopen(filename,"w"))==NULL)
      {
         fprintf(stderr,"Unable to open redirection file: %s\n",filename);
         fp = stdout;
      }
   }
   
   DisplaySearch(fp,StackDepth);
   
   if(fp != stdout)
      fclose(fp);

   return(TRUE);
}
   
      
/************************************************************************/
/*>BOOL HandleLogical(WHERE *wh, int *StackDepth)
   ----------------------------------------------
   Input:   WHERE *wh          An entry from the WHERE clause linked list
   I/O:     int   *StackDepth  The stack depth before and after the 
                               logical operation
   Returns: BOOL               Success (fails if stack too shallow)

   Handles a logical set operation in the WHERE linked list. Takes one
   or two sets from the stack and replaces them with a single set 
   containing the result of the logical operation.

   20.04.94 Original    By: ACRM
*/
BOOL HandleLogical(WHERE *wh, int *StackDepth)
{
   DATA *d;

   switch(wh->type)
   {
   case OPER_NOT:
      if(*StackDepth < 1) return(FALSE);
      for(d=gData; d!=NULL; NEXT(d))
         TOGGLE(d->active[(*StackDepth)-1]);
      break;
   case OPER_AND:
      if(*StackDepth < 2) return(FALSE);
      for(d=gData; d!=NULL; NEXT(d))
         d->active[(*StackDepth)-2] = (d->active[(*StackDepth)-2] &&
                                       d->active[(*StackDepth)-1]);
      (*StackDepth)--;
      break;
   case OPER_OR:
      if(*StackDepth < 2) return(FALSE);
      for(d=gData; d!=NULL; NEXT(d))
         d->active[(*StackDepth)-2] = (d->active[(*StackDepth)-2] ||
                                       d->active[(*StackDepth)-1]);
      (*StackDepth)--;
      break;
   default:
      break;
   }

   return(TRUE);
}


/************************************************************************/
/*>BOOL HandleMatch(WHERE *wh, int *StackDepth)
   --------------------------------------------
   Input:   WHERE *wh         An item from the WHERE clause linked list
   I/O:     int   *StackDepth The stack depth before and after adding
                              this set.
   Returns: BOOL              Success?

   Adds a new set to the stack by comparing each item in the data list
   with the specification in this WHERE clause.

   20.04.94 Original    By: ACRM
   26.04.94 Added light and heavy.
            Modified for new DoStrTest()
   27.04.94 Added L1...H3
   11.05.94 Added Canonical class handling
   02.04.96 Added ID (Kadbid) handling
   10.09.97 Added subgroup handling
*/
BOOL HandleMatch(WHERE *wh, int *StackDepth)
{
   DATA *d;
   char loop[40],
        res,
        class[8];
   int  len,
        idata;

   if(++(*StackDepth) >= STACKDEPTH)
   {
      fprintf(stderr,"Error: Search stack depth (%d) exceeded.\n",
              STACKDEPTH);
      return(FALSE);
   }

   for(d=gData; d!=NULL; NEXT(d))
   {
      switch(wh->type)
      {
      case FIELD_NAME:
         d->active[(*StackDepth)-1] = DoStrTest(d->name,
                                                wh->comparison, 
                                                wh->data, FALSE);
         break;
      case FIELD_ANTIGEN:
         d->active[(*StackDepth)-1] = DoStrTest(d->antigen,
                                                wh->comparison, 
                                                wh->data, FALSE);
         break;
      case FIELD_CLASS:
         d->active[(*StackDepth)-1] = DoStrTest(d->class,
                                                wh->comparison, 
                                                wh->data, FALSE);
         break;
      case FIELD_SOURCE:
         d->active[(*StackDepth)-1] = DoStrTest(d->source,
                                                wh->comparison, 
                                                wh->data, FALSE);
         break;
      case FIELD_REF:
         d->active[(*StackDepth)-1] = DoStrTest(d->reference,
                                                wh->comparison, 
                                                wh->data, FALSE);
         break;
      case FIELD_LENGTH:
         FillLoop(wh->param,d,loop);                /* Get loop         */
         len = TrueSeqLen(loop);                    /* Get length       */
         if(!sscanf(wh->data,"%d",&idata)) idata=0; /* Get required len */
         d->active[(*StackDepth)-1] = DoIntTest(len,
                                                wh->comparison, idata);
         break;
      case FIELD_RES:
         res = GetResidue(d, wh->param);
         d->active[(*StackDepth)-1] = DoCharTest(res,
                                                 wh->comparison, 
                                                 wh->data[0]);
         break;
      case FIELD_COMPLETE:
         d->active[(*StackDepth)-1] = DoBoolTest(IsComplete(d),
                                                 wh->comparison,
                                                 wh->data);
         break;
      case FIELD_LIGHT:
         d->active[(*StackDepth)-1] = DoStrTest(d->light,
                                                wh->comparison,
                                                wh->data, TRUE);
         break;
      case FIELD_HEAVY:
         d->active[(*StackDepth)-1] = DoStrTest(d->heavy,
                                                wh->comparison, 
                                                wh->data, TRUE);
         break;
      case FIELD_L1:
         FillLoop("L1",d,loop);
         d->active[(*StackDepth)-1] = DoStrTest(loop,
                                                wh->comparison, 
                                                wh->data, TRUE);
         break;
      case FIELD_L2:
         FillLoop("L2",d,loop);
         d->active[(*StackDepth)-1] = DoStrTest(loop,
                                                wh->comparison, 
                                                wh->data, TRUE);
         break;
      case FIELD_L3:
         FillLoop("L3",d,loop);
         d->active[(*StackDepth)-1] = DoStrTest(loop,
                                                wh->comparison, 
                                                wh->data, TRUE);
         break;
      case FIELD_H1:
         FillLoop("H1",d,loop);
         d->active[(*StackDepth)-1] = DoStrTest(loop,
                                                wh->comparison, 
                                                wh->data, TRUE);
         break;
      case FIELD_H2:
         FillLoop("H2",d,loop);
         d->active[(*StackDepth)-1] = DoStrTest(loop,
                                                wh->comparison, 
                                                wh->data, TRUE);
         break;
      case FIELD_H3:
         FillLoop("H3",d,loop);
         d->active[(*StackDepth)-1] = DoStrTest(loop,
                                                wh->comparison, 
                                                wh->data, TRUE);
         break;
      case FIELD_VAR:
         break;
      case FIELD_CANONICAL:
         if(FindCanonical(d,wh->param,class))
         {
            d->active[(*StackDepth)-1] = DoStrTest(class,
                                                   wh->comparison, 
                                                   wh->data, FALSE);
         }
         else
         {
            fprintf(stderr,"Error: Unable to get canonical \
information\n");
            return(FALSE);
         }
         break;
      case FIELD_IDLIGHT:
         d->active[(*StackDepth)-1] = DoStrTest(d->idlight,
                                                wh->comparison, 
                                                wh->data, FALSE);
         break;
      case FIELD_IDHEAVY:
         d->active[(*StackDepth)-1] = DoStrTest(d->idheavy,
                                                wh->comparison, 
                                                wh->data, FALSE);
         break;
      case FIELD_SUBGROUP:
         GetSubgroup(d,wh->param,class);
         d->active[(*StackDepth)-1] = DoStrTest(class,
                                                wh->comparison, 
                                                wh->data, FALSE);
         break;
      default:
         return(FALSE);
         break;
      }
   }

   return(TRUE);
}


/************************************************************************/
/*>BOOL DoStrTest(char *text, int comparison, char *subtext, BOOL fuzzy)
   ---------------------------------------------------------------------
   Input:   char *text         The data text
            int  comparison    The comparison type
            char *subtext      The text we are searching for
            BOOL fuzzy         Ignore -'s in either string
   Returns: BOOL               Match?

   Compares 2 strings using the mode specified by comparison. Not case
   sensitive.

   20.04.94 Original    By: ACRM
   26.04.94 Added fuzzy parameter
   15.12.95 Changed call to strstr() to use QueryStrStr()
*/
BOOL DoStrTest(char *text, int comparison, char *subtext, BOOL fuzzy)
{
   UPPER(text);
   UPPER(subtext);
   
   switch(comparison)
   {
   case COMP_EQ:
      if(!strcmp(text,subtext))          return(TRUE);
      break;
   case COMP_NE:
      if(strcmp(text,subtext))           return(TRUE);
      break;
   case COMP_SIM:
      if(fuzzy)
      {
         if(fuzzystrstr(text,subtext))   return(TRUE);
      }
      else
      {
         if(QueryStrStr(text,subtext)!=NULL)  return(TRUE);
      }
      break;
   default:
      break;
   }
   return(FALSE);
}


/************************************************************************/
/*>BOOL fuzzystrstr(char *text, char *subtext)
   -------------------------------------------
   Input:   char  *text       String in which to search
            char  *subtext    Substring for which to look
   Returns: BOOL              Found?

   Like strstr(), but ignores -'s in either string. Also returns BOOL
   (found?) rather than a pointer into the string.

   26.04.94 Original    By: ACRM
   22.06.95 Added checks on blank strings and string lengths
   23.06.95 Initialise SubText to NULL
   27.06.95 malloc() for SubText was being done wrong!
   15.12.95 Changed call to strstr() to use QueryStrStr()
*/
BOOL fuzzystrstr(char *text, char *subtext)
{
   char *Text,
        *SubText = NULL,
        *ptr = NULL;
   int  i,
        j;

   /* Both are blank strings                                            */
   if(!text[0] && !subtext[0])
      return(TRUE);
   
   /* One is a blank string                                             */
   if(!text[0] || !subtext[0])
      return(FALSE);

   /* Test if substring is longer than main string                      */
   if(strlen(text) < strlen(subtext))
      return(FALSE);
   
   if((Text = malloc((strlen(text)+1)*sizeof(char)))!=NULL)
   {
      if((SubText = malloc((strlen(subtext)+1)*sizeof(char)))!=NULL)
      {
         for(i=0, j=0; text[i]; i++)
            if(text[i] != '-')    Text[j++]    = text[i];
         Text[j]    = '\0';
         
         for(i=0, j=0; subtext[i]; i++)
            if(subtext[i] != '-') SubText[j++] = subtext[i];
         SubText[j] = '\0';

         ptr = QueryStrStr(Text,SubText);
      }         
   }

   if(Text    != NULL) free(Text);
   if(SubText != NULL) free(SubText);

   if(ptr     != NULL) return(TRUE);

   return(FALSE);
}


/************************************************************************/
/*>BOOL DoIntTest(int number, int comparison, int testnum)
   -------------------------------------------------------
   Input:   int  number        The data value
            int  comparison    The comparison type
            int  testnum       The value we are searching for
   Returns: BOOL               Match?

   Compares 2 integers using the mode specified by comparison.

   20.04.94 Original    By: ACRM
*/
BOOL DoIntTest(int number, int comparison, int testnum)
{
   switch(comparison)
   {
   case COMP_EQ:
      if(number==testnum) return(TRUE);
      break;
   case COMP_NE:
      if(number!=testnum) return(TRUE);
      break;
   case COMP_LT:
      if(number<testnum) return(TRUE);
      break;
   case COMP_LE:
      if(number<=testnum) return(TRUE);
      break;
   case COMP_GT:
      if(number>testnum) return(TRUE);
      break;
   case COMP_GE:
      if(number>=testnum) return(TRUE);
      break;
   default:
      break;
   }

   return(FALSE);
}


/************************************************************************/
/*>BOOL DoCharTest(char ch, int comparison, char testch)
   -----------------------------------------------------
   Input:   char ch            The data character
            int  comparison    The comparison type
            char testch        The character we are searching for
   Returns: BOOL               Match?

   Compares 2 characters using the mode specified by comparison. Not case
   sensitive.

   20.04.94 Original    By: ACRM
   26.04.94 Was ignoring the comparison mode
*/
BOOL DoCharTest(char ch, int comparison, char testch)
{
   char ch1, ch2;
   ch1 = (islower(ch)     ? toupper(ch)     : ch);
   ch2 = (islower(testch) ? toupper(testch) : testch);
   
   switch(comparison)
   {
   case COMP_EQ:
      if(ch1==ch2) return(TRUE);
      break;
   case COMP_NE:
      if(ch1!=ch2) return(TRUE);
      break;
   default:
      break;
   }
   
   return(FALSE);
}


/************************************************************************/
/*>BOOL DoBoolTest(BOOL condition, int comparison, char *test)
   -----------------------------------------------------------
   Input:   BOOL condition     The data logical value
            int  comparison    The comparison type
            char *test         The condition we are searching for
                               ("TRUE" or "FALSE")
   Returns: BOOL               Match?

   Compares a logical flag with a string representation using the
   specified comparison

   21.04.94 Original    By: ACRM
*/
BOOL DoBoolTest(BOOL condition, int comparison, char *test)
{
   BOOL TVal;
   
   if(test[0] == 'T' || test[0] == 't') TVal = TRUE;
   else                                 TVal = FALSE;

   switch(comparison)
   {
   case COMP_EQ:
      if(TVal == condition) return(TRUE);
      break;
   case COMP_NE:
      if(TVal != condition) return(TRUE);
      break;
   default:
      break;
   }

   return(FALSE);
}


/************************************************************************/
/*>BOOL IsComplete(DATA *d)
   ------------------------
   Input:   DATA    *d      A data entry
   Returns: BOOL            Complete?

   Tests whether an entry is complete (i.e. has both light and heavy 
   chains)

   20.04.94 Original    By: ACRM
*/
BOOL IsComplete(DATA *d)
{
   if(strlen(d->light) && strlen(d->heavy)) return(TRUE);
   
   return(FALSE);
}


/************************************************************************/
/*>void DisplaySearch(FILE *fp, int StackDepth)
   --------------------------------------------
   Input:   FILE *fp           Output file pointer
            int  StackDepth    Current stack depth

   Displays the selection specified by the top item on the stack.
   
   20.04.94 Original    By: ACRM
   21.04.94 Added commas between fields. Added output file pointer.
            Only ends line if something has been printed to screen.
            Added output of light and heavy
   25.04.94 Fixed bug in display of length
   11.05.94 Added canonical class
   12.05.94 Added hit count
   25.01.95 Added call to RemoveDupes()
   16.03.95 Corrected handling of PIR option when file name is not
            specified.
   02.04.96 Added ID***** and URL***** handling
            Added gHTML varibale testing
   11.04.96 Writes dataset date with number of hits
   10.09.97 Added subgroup
*/
void DisplaySearch(FILE *fp, int StackDepth)
{
   DATA      *d;
   SELECTION *p;
   FILE      *fpPIR     = NULL;
   char      loop[40],
             class[8],
             res;
   int       len,
             NHits      = 0;
   BOOL      first,
             GotPrint   = FALSE;

   if(StackDepth < 1 || StackDepth >= STACKDEPTH)
      return;
   
   if(gVariability > 0.0)
      RemoveDupes(StackDepth);
   
   for(d=gData; d!=NULL; NEXT(d))
   {
      if(d->active[StackDepth-1])
      {
         first = TRUE;
         NHits++;
         
         for(p=gSelectClause; p!=NULL; NEXT(p))
         {
            if(first)
               first = FALSE;
            else
               fprintf(fp,", ");
               
            switch(p->type)
            {
            case FIELD_NAME:
               fprintf(fp,"%s",d->name);
               GotPrint = TRUE;
               break;
            case FIELD_ANTIGEN:
               fprintf(fp,"%s",d->antigen);
               GotPrint = TRUE;
               break;
            case FIELD_L1:
               FillLoop("L1", d, loop);
               fprintf(fp,"%s",loop);
               GotPrint = TRUE;
               break;
            case FIELD_L2:
               FillLoop("L2", d, loop);
               fprintf(fp,"%s",loop);
               GotPrint = TRUE;
               break;
            case FIELD_L3:
               FillLoop("L3", d, loop);
               fprintf(fp,"%s",loop);
               GotPrint = TRUE;
               break;
            case FIELD_H1:
               FillLoop("H1", d, loop);
               fprintf(fp,"%s",loop);
               GotPrint = TRUE;
               break;
            case FIELD_H2:
               FillLoop("H2", d, loop);
               fprintf(fp,"%s",loop);
               GotPrint = TRUE;
               break;
            case FIELD_H3:
               FillLoop("H3", d, loop);
               fprintf(fp,"%s",loop);
               GotPrint = TRUE;
               break;
            case FIELD_CLASS:
               fprintf(fp,"%s",d->class);
               GotPrint = TRUE;
               break;
            case FIELD_SOURCE:
               fprintf(fp,"%s",d->source);
               GotPrint = TRUE;
               break;
            case FIELD_REF:
               fprintf(fp,"%s",d->reference);
               GotPrint = TRUE;
               break;
            case FIELD_LENGTH:
               FillLoop(p->param, d, loop);
               len = TrueSeqLen(loop);
               fprintf(fp,"%d",len);
               GotPrint = TRUE;
               break;
            case FIELD_RES:
               res = GetResidue(d, p->param);
               fprintf(fp,"%c",res);
               GotPrint = TRUE;
               break;
            case FIELD_PIR:
               fpPIR = fp;
               if(p->param[0])
               {
                  if((fpPIR=fopen(p->param,"w"))==NULL)
                     fpPIR = fp;
               }
               WriteAsPIR(fpPIR,d);
               break;
            case FIELD_LIGHT:
               fprintf(fp,"%s",d->light);
               GotPrint = TRUE;
               break;
            case FIELD_HEAVY:
               fprintf(fp,"%s",d->heavy);
               GotPrint = TRUE;
               break;
            case FIELD_CANONICAL:
               if(FindCanonical(d,p->param,class))
               {
                  fprintf(fp,"%s",class);
                  GotPrint = TRUE;
               }
               break;
            case FIELD_IDLIGHT:
               fprintf(fp,"%s",d->idlight);
               GotPrint = TRUE;
               break;
            case FIELD_IDHEAVY:
               fprintf(fp,"%s",d->idheavy);
               GotPrint = TRUE;
               break;
            case FIELD_URLLIGHT:
               if(d->idlight[0])
                  fprintf(fp,gURLFormat,d->idlight,d->idlight);
               else
                  fprintf(fp,"??????");
               GotPrint = TRUE;
               break;
            case FIELD_URLHEAVY:
               if(d->idheavy[0])
                  fprintf(fp,gURLFormat,d->idheavy,d->idheavy);
               else
                  fprintf(fp,"??????");
               GotPrint = TRUE;
               break;
            case FIELD_SUBGROUP:
               GetSubgroup(d,p->param,class);
               fprintf(fp,"%s",class);
               GotPrint = TRUE;
               break;
            default:
               break;
            }
         }
         if(GotPrint) fprintf(fp,"\n");
      }
   }

   if(gHTML) fprintf(fp,"<p><i>");
   fprintf(fp,"\nNumber of hits = %d (Dataset created %s)\n",
           NHits,gFileDate);
   if(gHTML) fprintf(fp,"</i><p>\n");

   if(fpPIR!=NULL && fpPIR!=stdout)
      fclose(fpPIR);
   fpPIR = NULL;
}


/************************************************************************/
/*>void FillLoop(char *loopname, DATA *d, char *loop)
   --------------------------------------------------
   Input:   char  *loopname     The loop name (L1...H3)
            DATA  *d            A data structure
   Output:  char  *loop         The sequence of the specified loop

   Extracts the sequence for a specified loop from a data structure

   25.04.94 Original    By: ACRM
*/
void FillLoop(char *loopname, DATA *d, char *loop)
{
   int  i,
        j,
        count,
        StartOff,
        EndOff;
   char start[8],
        end[8],
        *chain,
        **KabatIndex;
   
   /* Default (error) situation, output a blank string                  */
   loop[0] = '\0';

   /* Search the loop definition table for this loop                    */
   for(i=0; gLoopDefs[i].name != NULL; i++)
   {
      if(!upstrncmp(loopname,gLoopDefs[i].name,2))
      {
         /* Found the required loop, get the start and end              */
         switch(gLoopMode)
         {
         case LOOP_KABAT:
            strcpy(start, gLoopDefs[i].KabatS);
            strcpy(end,   gLoopDefs[i].KabatE);
            break;
         case LOOP_ABM:
            strcpy(start, gLoopDefs[i].AbMS);
            strcpy(end,   gLoopDefs[i].AbME);
            break;
         case LOOP_CHOTHIA:
            strcpy(start, gLoopDefs[i].ChothiaS);
            strcpy(end,   gLoopDefs[i].ChothiaE);
            break;
         default:
            return;
         }

         /* Store pointers depending on chain                           */
         if(start[0] == 'L')
         {
            chain      = d->light;
            KabatIndex = d->LNumbers;
         }
         else
         {
            chain      = d->heavy;
            KabatIndex = d->HNumbers;
         }
         

         /* Build residues into the output loop array                   */
         StartOff = GetKabatOffset(KabatIndex, start);
         EndOff   = GetKabatOffset(KabatIndex, end);
         
         for(count=0,j=StartOff; j<=EndOff; j++)
         {
            if(gShowInserts || chain[j] != '-')
               loop[count++] = chain[j];
         }
         loop[count] = '\0';

         /* Break out of the search                                     */
         break;
      }
   }
}


/************************************************************************/
/*>char GetResidue(DATA *d, char *resid)
   -------------------------------------
   Input:   DATA  *d      Item in Kabat data linked list
            char  *resid  Residue label (<chain><label>)
   Returns: char          Amino acid code for this residue

   Finds the 1-letter amino acid code for a residue specified as a chain
   and number (e.g. L24, H100D, etc). Returns X if no data available
   or the residue id is specified incorrectly.

   25.04.94 Original    By: ACRM
*/
char GetResidue(DATA *d, char *resid)
{
   char RetVal = 'X',
        ResID[8];

   strncpy(ResID, resid, 8);
   ResID[7] = '\0';

   UPPER(ResID);
   
   if(ResID[0] == 'L' && strlen(d->light))
      RetVal = d->light[GetKabatOffset(d->LNumbers,ResID)];
   else if(ResID[0] == 'H' && strlen(d->heavy))
      RetVal = d->heavy[GetKabatOffset(d->HNumbers,ResID)];
      
   return(RetVal);
}


/************************************************************************/
/*>BOOL FindCanonical(DATA *d, char *LoopID, char *class)
   ----------------------------------------------------
   Input:   DATA *d      A pointer to a data structure
            char *LoopID The loop name (L1...H3)
   Output:  char *class  The canonical class for this loop
   Returns: BOOL         Was the canonical class data available?

   Finds the Chothia canonical class for a spaecified loop.

   11.05.94 Original    By: ACRM
   12.05.94 Changed loop mode to AbM rather than Chothia
   07.05.96 The Canonical data are now stored with Chothia numbering
            if the keyword CHOTHIANUMBERS appears at the top of the
            data file. If this is so, the flag gCanonChothiaNum is
            set and the ChoKab() routine is called to convert each
            residue number before testing
*/
BOOL FindCanonical(DATA *d, char *LoopID, char *class)
{
   CHOTHIA *p;
   char    res,
           chain,
           ResID[8],
           LoopSeq[40];
   BOOL    Matched;
   int     i,
           LoopLen,
           LoopMode;

   /* Initialise class to unknown                                       */
   strcpy(class,"?");

   if(gChothia == NULL)
      return(FALSE);

   /* Get the chain id from the loop id                                 */
   chain = LoopID[0];
   if(islower(chain))
      chain = toupper(chain);
   
   /* Store the current loop mode and switch to AbM                     */
   LoopMode  = gLoopMode;
   gLoopMode = LOOP_ABM;    
   
   /* Get the loop length                                               */
   FillLoop(LoopID,d,LoopSeq);
   LoopLen = TrueSeqLen(LoopSeq);
   
   for(p=gChothia; p!=NULL; NEXT(p))      /* Go through Chothia data    */
   {
      if(LoopLen == p->length)            /* Check loop length          */
      {
         if(!upstrcmp(LoopID,p->LoopID))  /* Correct loop ID            */
         {
            Matched = TRUE;               /* Assume all residues match  */
            
            /* Step through the residue numbers for this canonical      */
            for(i=0; strncmp(p->resnum[i],"-1",2); i++)
            {
               /* Get the residue type for this data entry              */
               if(isalpha(p->resnum[i][0]))
                  strcpy(ResID,p->resnum[i]);
               else
                  sprintf(ResID,"%c%s",chain,p->resnum[i]);

               if(gCanonChothNum)
               {
                  /* If the Chothia data file uses Chothia numbering
                     rather than Kabat numbering, call the ChoKab() 
                     routine to convert this Chothia number to a Kabat
                     number (based on the length of L1/H1) before looking
                     up the contents of this residue
                  */
                  if(ResID[0] == 'L' || ResID[0] == 'l')
                  {
                     FillLoop("L1",d,LoopSeq);
                     res = GetResidue(d, ChoKab("L1", 
                                                TrueSeqLen(LoopSeq),
                                                ResID));
                  }
                  else
                  {
                     FillLoop("H1",d,LoopSeq);
                     res = GetResidue(d, ChoKab("H1", 
                                                TrueSeqLen(LoopSeq),
                                                ResID));
                  }
               }
               else
               {
                  res = GetResidue(d, ResID);
               }
               
               /* See if this type features in the allowed types        */
               if(strchr(p->restype[i],res)==NULL)
               {
                  Matched = FALSE;
                  break;
               }
            }
            
            /* If we still have a match, then copy in the class data, 
               restore the loop mode and return
               */
            if(Matched)
            {
               strcpy(class,p->class);
               gLoopMode  = LoopMode;
               return(TRUE);
            }
         }  /* Correct loop id                                          */
      }  /* Correct loop length                                         */
   }  /* Step through Chothia data                                      */

   /* We didn't find it, but we still return TRUE (after restoring the
      loop mode
   */
   gLoopMode  = LoopMode;
   return(TRUE);
}


/************************************************************************/
/*>void RemoveDupes(int StackDepth)
   --------------------------------
   Input:    int   StackDepth        Depth of stack. Must be 1

   Runs through the top of the stack and compares the sequences to
   remove near-duplicates.

   25.01.95 Original   By: ACRM
   23.06.95 Removed redundant variables
*/
void RemoveDupes(int StackDepth)
{
   DATA *d, *e;

   /* Check that there is only one item in the stack                    */
   if(StackDepth != 1)
      return;

   for(d=gData; d!=NULL; NEXT(d))
   {
      if(d->active[0])
      {
         for(e=d->next; e!=NULL; NEXT(e))
         {
            if(e->active[0])
            {
               if(TooSimilar(d,e,gVariability))
               {
                  e->active[0] = FALSE;
               }
            }
         }
      }
   }
}


/************************************************************************/
/*>BOOL TooSimilar(DATA *d, DATA *e, REAL Cutoff)
   ----------------------------------------------
   Input:    DATA     *d       Pointer to a data entry
             DATA     *e       Pointer to a data entry
             REAL     Cutoff   Percentage identity cutoff
   Returns:  BOOL              True of identity > Cutoff 

   Calculates a percentage identity between two antibodies. A missing
   chain in one of the sequences is not penalised. Each mismatch
   scores a penalty of 100/MeanLength % and deletions in one sequence 
   wrt the other score a double penalty.

   Returns immediately the identity falls below the specified Cutoff

   25.01.95 Original    By: ACRM
   26.01.95 Modified from SimilarityScore()
*/
BOOL TooSimilar(DATA *d, DATA *e, REAL Cutoff)
{
   int  i,
        Length1,
        Length2;
   REAL MeanLength,
        LScore = (REAL)100.0,
        HScore = (REAL)100.0,
        Score,
        Penalty,
        TwoPenalty,
        TwoCutoff,
        TwoCutMinus100;
   BOOL LDone = FALSE,
        HDone = FALSE;

   TwoCutoff      = (REAL)2.0*Cutoff;
   TwoCutMinus100 = TwoCutoff - (REAL)100.0;

   LDone = (d->light[0] && e->light[0]);
   HDone = (d->heavy[0] && e->heavy[0]);

   /* Compare the light chains                                          */
   if(LDone)
   {
      /* Calculate the mean sequence length and mismatch penalty        */
      Length1    = TrueSeqLen(d->light);
      Length2    = TrueSeqLen(e->light);
      MeanLength = (REAL)(Length1+Length2)/(REAL)2.0;
      Penalty    = (REAL)100.0/MeanLength;
      TwoPenalty = (REAL)2.0 * Penalty;

      /* For each residue in the light chain, decrement the score if there
         is a mismatch
      */
      for(i=0; i<strlen(d->light); i++)
      {
         if(d->light[i] != e->light[i])   /* A mismatch                 */
         {
            /* If one was a deletion wrt to the other, double penalty   */
            if((d->light[i] == '-') || (e->light[i] == '-'))
               LScore -= TwoPenalty;
            else
               LScore -= Penalty;

            /* This is an optimisation of 
               ((LScore+100.0)/2.0 < Cutoff)
            */
            if(LScore < TwoCutMinus100)
               return(FALSE);
         }
      }
   }
   
   /* Compare the heavy chains                                          */
   if(HDone)
   {
      /* Calculate the mean sequence length and mismatch penalty        */
      Length1    = TrueSeqLen(d->heavy);
      Length2    = TrueSeqLen(e->heavy);
      MeanLength = (REAL)(Length1+Length2)/(REAL)2.0;
      Penalty    = (REAL)100.0/MeanLength;
      TwoPenalty = (REAL)2.0 * Penalty;

      /* For each residue in the heavy chain, decrement the score if there
         is a mismatch
      */
      for(i=0; i<strlen(d->heavy); i++)
      {
         if(d->heavy[i] != e->heavy[i])   /* A mismatch                 */
         {
            /* If one was a deletion wrt to the other, double penalty   */
            if((d->heavy[i] == '-') || (e->heavy[i] == '-'))
               HScore -= TwoPenalty;
            else
               HScore -= Penalty;

            /* This is an optimisation of 
               ((HScore+LScore)/2.0 < Cutoff)
            */
            if((HScore+LScore) < TwoCutoff)
               return(FALSE);
         }
      }
   }

   /* Calculate the score based on which comparisons were performed     */
   if(LDone && HDone)
   {
      Score = (LScore + HScore) / (REAL)2.0;
   }
   else if(LDone)
   {
      Score = LScore;
   }
   else if(HDone)
   {
      Score = HScore;
   }
   else
   {
      Score = (REAL)100.0;
   }

   if(Score < Cutoff)
      return(FALSE);

   return(TRUE);
}


/************************************************************************/
/*>void GetSubgroup(DATA *d, char *chain, char *subgroup)
   ------------------------------------------------------
   This is an interface to Sophie Deret's routine for finding the 
   subgroup of a human sequence.


*/
void GetSubgroup(DATA *d, char *chain, char *subgroup)
{
   char class[16];

   *chain = (islower(*chain) ? toupper(*chain) : *chain);
   
   if(!strncmp(d->source,"HUMAN",5))
   {
      if(*chain == 'L')
      {
         DoGetSubgroup(d->light, class, subgroup);
         if(!strstr(d->class, class))
            strcpy(subgroup,"?");
      }
      else if(*chain == 'H')
      {
         DoGetSubgroup(d->heavy, class, subgroup);
      }
   }
   else
   {
      strcpy(subgroup,"?");
   }
}


/************************************************************************/
void DoGetSubgroup(char *sequence, char *class, char *subgroup)
{
   long classnum, sgpenum;
   char localseq[LARGEBUFF];
   int  i, j;
   
   /* Copy the sequence, skipping - characters                          */
   for(i=0,j=0; sequence[i]; i++)
   {
      if(sequence[i] == '?')
         localseq[j++] = 'X';
      else if(sequence[i] != '-')
         localseq[j++] = sequence[i];
   }
   localseq[j] = '\0';

   /* Call Sophie Deret's code to find the class & subgroup             */
   det_sgpe(localseq, &classnum, &sgpenum);

   /* Change the class/subgroup number to text                          */
   switch(classnum)
   {
   case 0:
      strcpy(class,"HEAVY");
      break;
   case 1:
      strcpy(class,"KAPPA");
      break;
   case 2:
      strcpy(class,"LAMBDA");
      break;
   default:
      strcpy(class,"?");
   }

   switch(sgpenum)
   {
   case 1:
      strcpy(subgroup,"I");
      break;
   case 2:
      strcpy(subgroup,"II");
      break;
   case 3:
      strcpy(subgroup,"III");
      break;
   case 4:
      strcpy(subgroup,"IV");
      break;
   case 5:
      strcpy(subgroup,"V");
      break;
   case 6:
      strcpy(subgroup,"VI");
      break;
   default:
      strcpy(subgroup,"?");
   }      
}


/************************************************************************/
/*>void WriteAsPIR(FILE *fp, DATA *d)
   ----------------------------------
   Input:   FILE  *fp      File pointer for output file
            DATA  *d       POinter to this data item

   Writes the data from the entry in simple PIR format

   21.04.94 Original    By: ACRM
   25.04.94 Puts chain name in title if only one chain
            Code truncated at 6 chars
*/
void WriteAsPIR(FILE *fp, DATA *d)
{
   char buffer[40],
        *chp;
   int  i,
        len,
        count;
   
   /* Remove any trailing 'xxx from the name                            */
   strncpy(buffer,d->name,40);
   if((chp=strchr(buffer,'\''))!=NULL) *chp = '\0';

   /* Do the header lines                                               */
   fprintf(fp,">P1;");
   for(i=0; i<6; i++)
   {
      if(!buffer[i]) break;
      fprintf(fp,"%c",  buffer[i]);
   }
   fprintf(fp,"\n");

   fprintf(fp,"%s - (%s) %s", buffer, d->source, d->fsource);
   if(strlen(d->light)==0)
      fprintf(fp," (HEAVY CHAIN)");
   else if(strlen(d->heavy)==0)
      fprintf(fp," (LIGHT CHAIN)");
   fprintf(fp,"\n");

   /* Do the light chain                                                */
   if((len = strlen(d->light))!=0)
   {
      for(i=0, count=0; i<len; i++)
      {
         if(d->light[i] != '-')
         {
            fprintf(fp,"%c",d->light[i]);
            if(!((++count)%40)) fprintf(fp,"\n");
         }
      }
      fprintf(fp,"*\n");
   }
   
   /* Do the heavy chain                                                */
   if((len = strlen(d->heavy))!=0)
   {
      for(i=0, count=0; i<len; i++)
      {
         if(d->heavy[i] != '-')
         {
            fprintf(fp,"%c",d->heavy[i]);
            if(!(++count%40)) fprintf(fp,"\n");
         }
      }
      fprintf(fp,"*\n");
   }
}


