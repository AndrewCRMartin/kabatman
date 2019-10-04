/*************************************************************************

   Program:    KabatMan
   File:       RdKabat.c
   
   Version:    V2.26
   Date:       04.10.19
   Function:   Read data from a new format Kabat sequence file
   
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
   V0.1  13.04.94 Development version
   V1.0  27.04.94 Original
   V1.1  05.04.94 Fixed bug in InsertSeq()
   V1.2  13.05.94 Now reads blank sequence entries as ? rather than -
   V2.0  30.06.94 Reads the new or old format files
   V2.1  11.07.94 Skipped
   V2.2  21.07.94 Read strain from ANNOTA STRN records.
                  Reads Multiple references
   V2.3  23.01.95 Skipped
   V2.4  10.02.95 Modified for the 1995 Kabat release format
   V2.5  07.03.95 Modified sequence storage to replace ? at end of 
                  sequence with -
   V2.6  16.03.95 Skipped
   V2.7  16.05.95 Skipped
   V2.8  22.06.95 Skipped
   V2.9  23.06.95 Clean compile under gcc -Wall
   V2.10 27.06.95 Skipped
   V2.11 15.12.95 Skipped
   V2.12 02.04.96 Added kadbid stuff
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
   V2.23 03.04.02 Added refdate stuff
   V2.24 28.02.05 Skipped
   V2.25 24.08.06 Modified GetKabatOffset() so it can return the label
                  from the offset as well as the offset from the label
   V2.26 04.10.19 Changed all bioplib calls to blXXX()

*************************************************************************/
/* Includes
*/
#include "RdKabat.h"
#include <stdlib.h>

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
#include "RdKabat.p"
static int FixSequence(KABATENTRY *KabatEntry, BOOL *insert, BOOL
                       OldFormat);
static int InsertSeq(KABATENTRY *KabatEntry, char *resnum, int *i,
                     BOOL OldFormat);

/************************************************************************/
/*>int ReadNextKabatEntry(FILE *fp, KABATENTRY *Kabat, BOOL *insert,
                          BOOL OldFormat)
   -----------------------------------------------------------------
   Input:   FILE        *fp          Kabat file pointer
            BOOL        OldFormat    Read old format files
   Output:  KABATENTRY  *Kabat       Completed KABATENTRY data structure
            BOOL        *insert      Has an insertion occurred?
   Returns: int                      Length of sequence 
                                      0: End of file
                                     -1: Unable to handle insertion

   Read the next entry out of a Kabat data file

   12.04.94 Original    By: ACRM
   13.04.94 Calls ClearKabatEntry(). Changed to return int rather 
            than BOOL
   22.04.94 Added insert parameter
   29.06.94 Checks KADBID rather than A TABLE
   30.06.94 Added OldFormat option
   23.06.95 Removed redundant variables
*/
int ReadNextKabatEntry(FILE *fp, KABATENTRY *Kabat, BOOL *insert,
                       BOOL OldFormat)
{
   char buffer[SEQBUFF],    /* A BIG buffer to read sequence lines      */
        *bp;

   ClearKabatEntry(Kabat);

   while(fgets(buffer,SEQBUFF,fp))
   {
      TERMINATE(buffer);
      KILLLEADSPACES(bp,buffer);

      if(strlen(bp) != 0)
      {
         if(OldFormat)
         {
            if(!strncmp(buffer,"AA TABLE",8))
            {
               return(ReadOldKabatEntry(fp,buffer,SEQBUFF,Kabat,insert));
            }
         }
         else
         {
            if(!strncmp(buffer,"KADBID",6))
            {
               return(ReadKabatEntry(fp,buffer,SEQBUFF,Kabat,insert));
            }
         }
      }
   }

   return(0);
}

/************************************************************************/
/*>static int FixSequence(KABATENTRY *KabatEntry, BOOL *insert, 
                          BOOL OldFormat)
   ------------------------------------------------------------
   I/O:     KABATENTRY  *KabatEntry   Input:  Raw Kabat entry from file
                                      Output: Kabat entry with 1-letter
                                              sequence data complete
   Output:  BOOL        *insert       Flag to indicate inserts c.f.
                                      Kabat standard numbering
   Input:   BOOL        OldFomat      Flag for old format
   Returns: int                       Number of residues in sequence
                                      -1: Unable to handle insertion

   Process the 3 letter code sequence from the raw Kabat data into a
   1-letter code sequence handling insertions.
   Used by ReadKabatEntry()
   
   16.03.94 Original    By: ACRM
   18.03.94 Added insert parameter
   13.04.94 No longer increments sequence count for blank entries
   14.04.94 Now handles blanks as -, but for both counts the length
            separately from the sequence array. Thus returns the
            true length rather than the sequence buffer usage.
   18.04.94 Added -1 return
   22.04.94 Sets *insert to FALSE on entry
   25.04.94 Replaces the |s when stepping on to next entry
            Uses own copy of the sequence and residue number buffers
   27.04.94 Allows unseparted bars for deletions
   13.05.94 Treats a blank entry as ? rather than -. Returns
            KnownSeqLen() rather than TrueSeqLen()
   21.07.94 Added OldFormat flag
   07.03.95 Runs from the end of the sequence replacing all ?s with -s
            until a real residue is hit. This stops all sequences which
            are truncated before or within CDR-3 from having maximum
            length loops when the LEN() query is used
*/
static int FixSequence(KABATENTRY *KabatEntry, BOOL *insert, BOOL
                       OldFormat)
{
   char *sequence,
        *seqstart,
        *resnum,
        *resstart,
        *sp,
        *np;
   int  i      = 0;
   BOOL ok     = TRUE;


   if((seqstart = 
       (char *)malloc((1+strlen(KabatEntry->kabatseq))*sizeof(char)))
      ==NULL)
   {
      fprintf(stderr,"Error: No memory for Kabat numbering index\n");
      return(0);
   }
   sequence = seqstart;
   strcpy(sequence,KabatEntry->kabatseq);

   if((resstart = 
       (char *)malloc((1+strlen(KabatEntry->kabatnum))*sizeof(char)))
      ==NULL)
   {
      fprintf(stderr,"Error: No memory for Kabat numbering index\n");
      free(seqstart);
      return(0);
   }
   resnum   = resstart;
   strcpy(resnum  ,KabatEntry->kabatnum);

   *insert  = FALSE;

   while((sp = strchr(sequence,'|')) != NULL)
   {
      /* Terminate string at next bar                                   */
      *sp = '\0';

      /* Get the residue number                                         */
      if((np = strchr(resnum,'|')) != NULL)
         *np = '\0';
      else
         ok = FALSE;

      /* Do sequence conversion                                         */
      if(*sequence == ' ')
      {
         KabatEntry->sequence[i++] = '?';
      }
      else if(*sequence == '-' || *sequence == '\0')
      {
         KabatEntry->sequence[i++] = '-';
      }
      else if(*sequence == '#')
      {
         *insert = TRUE;
         if(!InsertSeq(KabatEntry,resnum,&i,OldFormat))
         {
            free(seqstart);
            free(resstart);
            return(-1);
         }
      }
      else
      {
         KabatEntry->sequence[i++] = blThrone(sequence);
      }

      /* Step on to next entry                                          */
      *sp = '|';
      *np = '|';
      sequence = sp+1;
      if(ok) resnum = np+1;
   }

   /* Do the last residue                                               */
   if(*sequence == '-' || *sequence == ' ')
   {
      KabatEntry->sequence[i++] = '-';
   }
   else if(*sequence == '#')
   {
      *insert = TRUE;
      if(!InsertSeq(KabatEntry,resnum,&i,OldFormat))
      {
         free(seqstart);
         free(resstart);
         return(-1);
      }
   }
   else
   {
      KabatEntry->sequence[i++] = blThrone(sequence);
   }

   /* Terminate the sequence                                            */
   KabatEntry->sequence[i] = '\0';
   
   free(seqstart);
   free(resstart);

   /* Now step back from the end of the sequence replacing ? with - until
      a real residue is hit
   */
   sp = KabatEntry->sequence + i - 1;
   while(sp >= KabatEntry->sequence && (*sp == '?' || *sp == '-'))
   {
      *sp = '-';
      sp--;
   }

   return(blKnownSeqLen(KabatEntry->sequence));
}

/************************************************************************/
/*>static int InsertSeq(KABATENTRY *KabatEntry, char *resnum, int *i,
                        BOOL OldFormat)
   ------------------------------------------------------------------
   I/O:     KABATENTRY  *KabatEntry   Kabat entry in which an insertion
                                      is to be processed
            int         *i            Current sequential residue number
   Input:   char        *resnum       Residue id at which the insertion
                                      is being handled
            BOOL        OldFormat     Flag to read old format
   Returns: int                       Number of residues inserted
                                      (0 = error)

   Handle an insertion in a Kabat entry.
   Used by ReadKabatEntry()   

   16.03.94 Original    By: ACRM
   18.03.94 Returns flag to indicate insertion found
   18.04.94 Now only reads insertions from the INSERTSAA or INSERTSNUC
            lines. NOTESAA is unsafe.
   22.04.94 Now returns the number of residues inserted (0 if error)
   27.04.94 Fixed bug in finding next number where more than one insert
   05.05.94 Fixed bug when resnum has a letter
   21.07.94 For new format only scans the insertsaa data.
            Added parameter for oldformat
*/
static int InsertSeq(KABATENTRY *KabatEntry, char *resnum, int *i,
                     BOOL OldFormat)
{
   char *cp,
        *end,
        *start,
        buffer[LARGEBUFF];
   int  len,
        count = 0;
   
   /* First scan the INSERTSAA record to see if the insertion is
      specified here
   */
   strncpy(buffer,KabatEntry->insertsaa,LARGEBUFF);
   len = strlen(buffer);
   
   if((cp=strstr(buffer,resnum)) != NULL)
   {
      /* Found it, so extract the sequence.
         First terminate at the next number (if found)
      */
      cp += strlen(resnum);
      for(end=cp; *end && !isdigit(*end); end++);
      *end = '\0';

      start = cp;
      len   = strlen(start);
      
      /* Now step over any non-alpha characters; what's left should be
         the sequence
      */
      while(!isalpha(*cp)) cp++;
      
      end = cp+3;
      while(end-start <= len)
      {
         /* Terminate this entry                                        */
         *end = '\0';

         /* Store 1-letter code                                         */
         KabatEntry->sequence[(*i)++] = blThrone(cp);
         count++;

         /* Step on to the next entry                                   */
         cp  = end+1;
         end = cp+3;
      }

      return(count);
   }

   if(OldFormat)
   {
      /* We didn't find it in INSERTSAA, so look in INSERTSNUC          */
      strncpy(buffer,KabatEntry->insertsnuc,LARGEBUFF);
      len = strlen(buffer);

      if((cp=strstr(buffer,resnum)) != NULL)
      {
         /* Found it, so extract the sequence.
            First terminate at the next number (if found)
            */
         for(end=cp; *end && !isdigit(*end); end++);
         end = '\0';
         
         /* Now step over any non-alpha characters; what's left should be
            the sequence
            */
         while(!isalpha(*cp)) cp++;
         
         end = cp+3;
         while(end-buffer <= len)
         {
            /* Terminate this entry                                     */
            *end = '\0';
            
            /* Store 1-letter code                                      */
            KabatEntry->sequence[(*i)++] = blDNAtoAA(cp);
            count++;
            
            /* Step on to the next entry                                */
            cp  = end+1;
            end = cp+3;
         }
         
         return(count);
      }
   }
   
   
   /* Didn't find it in INSERTSNUC either. If it's the first position
      in the sequence, put in a -; otherwise return an error
   */
   if(!(*i))
   {
      KabatEntry->sequence[(*i)++] = '-';
      return(1);
   }
   
   return(0);
}


/************************************************************************/
/*>void ClearKabatEntry(KABATENTRY *KabatEntry)
   --------------------------------------------
   I/O:     KABATENTRY  *KabatEntry    Kabat entry to be cleared

   Clears a Kabat entry structure

   17.03.94 Original    By: ACRM
   13.04.94 Changed name from ClearEntry()
   02.04.96 Added kadbid
   03.04.02 Added refdate
*/
void ClearKabatEntry(KABATENTRY *KabatEntry)
{
   KabatEntry->aatable[0]   = '\0';
   KabatEntry->aaname[0]    = '\0';
   KabatEntry->codname[0]   = '\0';
   KabatEntry->reference[0] = '\0';
   KabatEntry->antigen[0]   = '\0';
   KabatEntry->species[0]   = '\0';
   KabatEntry->class[0]     = '\0';
   KabatEntry->strain[0]    = '\0';
   KabatEntry->source[0]    = '\0';
   KabatEntry->insertsaa[0] = '\0';
   KabatEntry->notesaa[0]   = '\0';
   KabatEntry->kabatnum[0]  = '\0';
   KabatEntry->kabatseq[0]  = '\0';
   KabatEntry->sequence[0]  = '\0';
   KabatEntry->kadbid[0]    = '\0';
   KabatEntry->refdate      = 9999;
}



/************************************************************************/
/*>int GetKabatOffset(char **table, char *label, int count)
   --------------------------------------------------------
   Input:   char  *table    Lookup table or NULL (use internal tables)
            int   count     <0 or offset for which to find the label
   I/O      char  *label    If (count<0) this is the label for which to 
                            search (input only)
                            If (count>=0) the chain is taken from this
                            label and the full label for the offset 
                            specified in count is returned
   Returns: int             Offset into table (-1 if error)

   If (count < 0):

   Finds the offset of a Kabat residue number into the sequence array.
   If the table is specified as NULL, the standard kabat tables will be
   used, otherwise the table given will be used (this must be an array
   of pointers, terminated with a NULL pointer).

   The label must be specified as <chain><number> (e.g. L27A), although
   the <chain> will be ignored if a table is specified.

   If (count >= 0):

   Takes the chain name from 'label'. Then works out the label for the
   offset specified in 'count' and replaces 'label' with the full label.
   Returns the offset into the table (same as the input 'count') or
   (-1) if the count is outside the range of labels.

   22.04.94 Original    By: ACRM
   24.08.06 Now takes an additional parameter, count. If count<0 then
            the routine behaves as before. If count>=0 then this is the
            offset and the routine takes the chain identier from 'label'
            and replaces label with the full label for the specified
            offset
*/
int GetKabatOffset(char **table, char *label, int count)
{
   static char *KabatLightNumbers[] = 
   {
      "0","1","2","3","4","5","6","7","8","9","10","11","12","13","14",
      "15","16","17","18","19","20","21","22","23","24","25","26","27",
      "27A","27B","27C","27D","27E","27F","28","29","30","31","32","33",
      "34","35","36","37","38","39","40","41","42","43","44","45","46",
      "47","48","49","50","51","52","53","54","55","56","57","58","59",
      "60","61","62","63","64","65","66","67","68","69","70","71","72",
      "73","74","75","76","77","78","79","80","81","82","83","84","85",
      "86","87","88","89","90","91","92","93","94","95","95A","95B","95C",
      "95D","95E","95F","96","97","98","99","100","101","102","103","104",
      "105","106","106A","107","108","109",NULL
   }  ;

   static char *KabatHeavyNumbers[] = 
   {
      "0","1","2","3","4","5","6","7","8","9","10","11","12","13","14",
      "15","16","17","18","19","20","21","22","23","24","25","26","27",
      "28","29","30","31","32","33","34","35","35A","35B","36","37","38",
      "39","40","41","42","43","44","45","46","47","48","49","50","51",
      "52","52A","52B","52C","53","54","55","56","57","58","59","60","61",
      "62","63","64","65","66","67","68","69","70","71","72","73","74",
      "75","76","77","78","79","80","81","82","82A","82B","82C","83","84",
      "85","86","87","88","89","90","91","92","93","94","95","96","97",
      "98","99","100","100A","100B","100C","100D","100E","100F","100G",
      "100H","100I","100J","100K","101","102","103","104","105","106",
      "107","108","109","110","111","112","113",NULL
   }  ;

   char chain,
        buffer[16],
        *pch;
   int  i;

   if(count < 0)
   {
      /* Upcase the input label                                         */
      strncpy(buffer,label,16);
      UPPER(buffer);
      
      /* Get the chain label out of the label                           */
      chain = buffer[0];
      pch   = buffer + 1;
      
      if(table == NULL)  /* No table, use the internal one              */
      {
         table = (chain=='L' ? KabatLightNumbers : 
                  (chain=='H' ? KabatHeavyNumbers : NULL));
      }
      
      if(table == NULL)  /* Illegal chain specifier                     */
         return(-1);
      
      /* Search for this string and return offset into the table        */
      for(i=0; table[i]!=NULL; i++)
         if(!strcmp(table[i],pch))
            return(i);
   }
   else
   {
      /* Upcase the input label                                         */
      strncpy(buffer,label,16);
      UPPER(buffer);
      
      /* Get the chain label out of the label                           */
      chain = buffer[0];
      pch   = buffer + 1;
      
      if(table == NULL)  /* No table, use the internal one              */
      {
         table = (chain=='L' ? KabatLightNumbers : 
                  (chain=='H' ? KabatHeavyNumbers : NULL));
      }
      
      if(table == NULL)  /* Illegal chain specifier                     */
         return(-1);
      
      /* Check the size of the table. If the count is less than the table
         size, find the residue name and output it. Return the count
      */
      for(i=0; table[i]!=NULL; i++);
      if(count < i)
      {
         strcpy(label, table[count]);
         return(count);
      }
   }
      
   /* Not found, return (-1)                                            */
   return(-1);
}


/************************************************************************/
/*>char **BuildKabatNumbering(KABATENTRY Kabat, BOOL OldFormat)
   ------------------------------------------------------------
   Input:   KABATENTRY    Kabat      Kabat entry structure
            BOOL          OldFormat  Flag for old format insertions
   Returns: char          *          Array of Kabat numbering

   Builds a sensible numbering scheme for Kabat entries which have
   insertions

   22.04.94 Original    By: ACRM
   25.04.94 Uses own copy of sequence and resnum buffer.
            Correctly initialise i and j counters
            Call to InsertSeq() stores offset in junk
            Added call to CheckKabatNumbering()
   21.07.94 Added OldFormat flag
*/
char **BuildKabatNumbering(KABATENTRY Kabat, BOOL OldFormat)
{
   KABATENTRY TempKabat;
   char       **KabatIndex = NULL,
              *sp,
              *np,
              *sequence,
              *seqstart,
              *resnum,
              *resstart;
   BOOL       ok           = TRUE;
   int        i            = 0,
              j            = 0,
              junk,
              NInsert;
   
   if((KabatIndex = (char **)malloc(MAXKABATSEQ * sizeof(char *)))==NULL)
   {
      fprintf(stderr,"Error: No memory for Kabat numbering index\n");
      return(NULL);
   }

   TempKabat = Kabat;
   if((seqstart = 
       (char *)malloc((1+strlen(TempKabat.kabatseq))*sizeof(char)))==NULL)
   {
      fprintf(stderr,"Error: No memory for Kabat numbering index\n");
      return(NULL);
   }
   sequence = seqstart;
   strcpy(sequence,TempKabat.kabatseq);

   if((resstart = 
       (char *)malloc((1+strlen(TempKabat.kabatnum))*sizeof(char)))==NULL)
   {
      fprintf(stderr,"Error: No memory for Kabat numbering index\n");
      free(seqstart);
      return(NULL);
   }
   resnum = resstart;
   strcpy(resnum,TempKabat.kabatnum);

   while((sp = strchr(sequence,'|')) != NULL)
   {
      /* Terminate string at next bar                                   */
      *sp = '\0';

      /* Get the residue number                                         */
      if((np = strchr(resnum,'|')) != NULL)
         *np = '\0';
      else
         ok = FALSE;

      if((KabatIndex[i] = (char *)malloc(5*sizeof(char)))==NULL)
      {
         fprintf(stderr,"Error: No memory for Kabat numbering index\n");
         free(resstart);
         free(seqstart);
         return(NULL);
      }
      strcpy(KabatIndex[i++], resnum);

      /* Use InsertSeq() to find how many insertions required           */
      if(*sequence == '#')
      {
         junk = i;
         if((NInsert = InsertSeq(&TempKabat,resnum,&junk,OldFormat))==0)
         {
            free(resstart);
            free(seqstart);
            return(NULL);
         }
         
         for(j=1; j<NInsert; j++)
         {
            if((KabatIndex[i++] = CreateKabatNumber(resnum,j))==NULL)
            {
               fprintf(stderr,"Error: No memory for Kabat numbering \
index\n");
               free(resstart);
               free(seqstart);
               return(NULL);
            }
         }
      }

      /* Step on to next entry                                          */
      *sp = '|';
      *np = '|';
      sequence = sp+1;
      if(ok) resnum = np+1;
   }

   /* Do the last residue                                               */
   if((KabatIndex[i] = (char *)malloc(5*sizeof(char)))==NULL)
   {
      fprintf(stderr,"Error: No memory for Kabat numbering index\n");
      free(resstart);
      free(seqstart);
      return(NULL);
   }
   strcpy(KabatIndex[i++], resnum);
   
   /* Use InsertSeq() to find how many insertions required              */
   if(*sequence == '#')
   {
      junk = i;
      if((NInsert = InsertSeq(&TempKabat,resnum,&junk,OldFormat))==0)
      {
         free(resstart);
         free(seqstart);
         return(NULL);
      }
      
      for(j=1; j<NInsert; j++)
      {
         if((KabatIndex[i++] = CreateKabatNumber(resnum,j))==NULL)
         {
            fprintf(stderr,"Error: No memory for Kabat numbering \
index\n");
            free(resstart);
            free(seqstart);
            return(NULL);
         }
      }
   }

   /* Terminate the sequence                                            */
   KabatIndex[i] = NULL;

   /* Finally check for silly positions of insertions                   */
   CheckKabatNumbering(KabatIndex);
   
   free(resstart);
   free(seqstart);
   return(KabatIndex);
}

/************************************************************************/
/*>void CheckKabatNumbering(char **KabatIndex)
   -------------------------------------------
   I/O:     char  **KabatIndex     Array of Kabat residue numbers

   Checks that insertion numbering makes sense. If an insertion has
   been specified within the region of Kabat lettered inserts rather
   than at the end, a letter could repeat.
   e.g. Normal Kabat heavy chain numbering has 100, 100A...100K. If an
   extra insertion of 5 residues is specified at position 100J, the
   inserted letters will be 100, 100A...100O, 100K i.e. 100K occurs
   twice (This is seen in Chicken heavy chain HC86).

   This routine sorts out insertion lettering so it is always
   correct.

   25.04.94 Original    By: ACRM
   23.06.95 Removed redundant variables
*/
void CheckKabatNumbering(char **KabatIndex)
{
   int  i,
        offset = 0,
        len    = 0;
   char label[8];
   BOOL Insert = FALSE;
   
   for(i=0; KabatIndex[i] != NULL; i++)
   {
      if(strchr(KabatIndex[i],'A'))
      {
         /* We've got the start of a new inserted region. Copy the
            residue number into label, terminating it at the insert
            character
         */
         offset = 0;
         strcpy(label,KabatIndex[i]);
         *strchr(label,'A') = '\0';
         len = strlen(label);
         Insert = TRUE;
      }

      /* Compare the current Kabat label with the first in an inserted
         region.
      */
      if(Insert)
      {
         if(!strncmp(KabatIndex[i],label,len))
         {
            /* We got a match, set the insert character to the next in 
               sequence
               */
            KabatIndex[i][len] = 'A'+offset;
            offset++;
         }
         else
         {
            Insert = FALSE;
         }
      }
   }
}

/************************************************************************/
/*>char *CreateKabatNumber(char *resnum, int offset)
   -------------------------------------------------
   Input:   char  *resnum     Residue number from which to derive new one
            int   offset      Offset from this base residue
   Returns: char  *           Pointer to malloc'd char array containing
                              new label

   Allocates space and creates a new Kabat label based on a known label
   and an offset from that label. Thus given "24" and 2 as parameters
   the string "24B" will be returned. Also "27F",2 will generate "27H"

   N.B. Assumes sequential (ASCII) coding of characters
*/
char *CreateKabatNumber(char *resnum, int offset)
{
   char insert = ('A'-1),     /* One before 'A'                         */
        *buffer;
   int  i;

   /* Allocate some space for our output                                */
   if((buffer = (char *)malloc(5 * sizeof(char)))==NULL)
      return(NULL);

   /* Copy in the current label                                         */
   strncpy(buffer,resnum,5);
   UPPER(buffer);

   /* First see if there's a letter in the current label                */
   for(i=0; i<strlen(buffer); i++)
   {
      if(isalpha(buffer[i]))
      {
         insert = buffer[i];
         break;
      }
   }

   /* Add the offset to the insert character                            */
   insert += offset;

   /* Patch this character into the buffer for return                   */
   for(i=0; i<strlen(buffer); i++)
   {
      if(isalpha(buffer[i]) || buffer[i] == ' ')
      {
         buffer[i] = insert;
         return(buffer);
      }
   }
   buffer[i] = insert;
   buffer[i+1] = '\0';
   
   return(buffer);
}

/************************************************************************/
/*>int ReadKabatEntry(FILE *fp, char *buffer, int bufflen, 
                      KABATENTRY *KabatEntry, BOOL *insert)
   --------------------------------------------------------
   Input:   FILE       *fp          Kabat file pointer
            int        bufflen      Size of buffer
   I/O:     char       *buffer      Input:  First line of Kabat entry
                                    Output: Last line of Kabat entry
   Output:  KABATENTRY *KabatEntry  Completed kabat data entry
   Returns: int                     Number of residues read
                                    -1: Unable to handle insertion

   Read a complete entry from a Kabat data file.
   
   16.03.94 Original    By: ACRM
   18.03.94 Added insert parameter. FixSequence() called here. Returns
            number of residues read
   23.03.94 Ensures strings are terminated.
   13.04.94 Changed name from ReadEntry() to ReadKabatEntry()
   18.04.94 Added -1 return; added reading of INSERTSNUC
   29.06.94 Re-written for new format Kabat files
   04.07.94 A few small fixes to new format reading
   18.07.94 Reads strain from ANNOTA STRN record
   21.07.94 Reads all references separating each with a bar
            Added FALSE OldFormat flag in call to FixSequence()
   10.02.95 Reads SEQTPA rather than SEQRES records for Jan 1995 release
            of Kabat database.
   02.04.96 Added KADBID and changed while to do/while since the
            KADBID stuff is in the buffer when we enter this routine
   03.04.02 Added reading of earliest reference date into refdate
*/
int ReadKabatEntry(FILE *fp, char *buffer, int bufflen, 
                   KABATENTRY *KabatEntry, BOOL *insert)
{
   BOOL GotJournal = FALSE;

   do
   {
      TERMINATE(buffer);

      if(!strncmp(buffer,"RECEND|",7))
         break;
      else if(!strncmp(buffer,"AANAME",6))
      {
         strncpy(KabatEntry->aaname,buffer+12,SMALLBUFF);
         KabatEntry->aaname[SMALLBUFF-1] = '\0';
      }
      else if(!strncmp(buffer,"AAREFA    ",10))
      {
         int len = strlen(KabatEntry->reference);

         if(GotJournal)
         {
            strncat(KabatEntry->reference,"|",LARGEBUFF-len);
            len++;
         }
         
         strncat(KabatEntry->reference,buffer+12,LARGEBUFF-len);
         KabatEntry->reference[LARGEBUFF-1] = '\0';

         GotJournal = FALSE;
      }
      else if(!strncmp(buffer,"AAREFJ    ",10))
      {
         int  len     = strlen(KabatEntry->reference);
         char *year_p = buffer+12;
         int  pubdate;
         BOOL GotDate = FALSE;

         strncat(KabatEntry->reference,buffer+12,LARGEBUFF-len);
         KabatEntry->reference[LARGEBUFF-1] = '\0';

         GotJournal = TRUE;

         /* 03.04.02 Grab the year out of this reference                */
         while(!GotDate && (year_p!=NULL) && *year_p)
         {
            /* Find the first occurrence of (                           */
            if((year_p = strchr(year_p, '('))!=NULL)
            {
               /* See if there is a ) 5 characters on                   */
               if((year_p+5) < (buffer+bufflen))
               {
                  if(*(year_p+5) == ')')
                  {
                     /* Set this to a null and read the intervening 
                        characters as the date
                     */
                     *(year_p+5) = '\0';
                     if(sscanf(year_p+1, "%d", &pubdate))
                     {
                        if(pubdate < KabatEntry->refdate)
                        {
                           KabatEntry->refdate = pubdate;
                        }
                        GotDate=TRUE;
                     }
                     else
                     {
                        year_p++;
                     }
                  }
                  else
                  {
                     year_p++;
                  }
               }
            }
         }
      }
      else if(!strncmp(buffer,"ANNOTA SPEC",11))
      {
         int len = strlen(KabatEntry->antigen);
         
         strncat(KabatEntry->antigen,buffer+12,LARGEBUFF-len);
         KabatEntry->antigen[LARGEBUFF-1] = '\0';
      }
/*
      else if(!strncmp(buffer,"SPECIE",6))
      {
         strncpy(KabatEntry->species,buffer+12,SMALLBUFF);
         KabatEntry->species[SMALLBUFF-1] = '\0';
      }
*/
      else if(!strncmp(buffer,"ANNOTA CLAS",11))
      {
         strncpy(KabatEntry->class,buffer+12,SMALLBUFF);
         KabatEntry->class[SMALLBUFF-1] = '\0';
      }
      else if(!strncmp(buffer,"ANNOTA STRN",11))
      {
         strncpy(KabatEntry->strain,buffer+12,SMALLBUFF);
         KabatEntry->strain[SMALLBUFF-1] = '\0';
      }
      else if(!strncmp(buffer,"SPECIE",6))
      {
         strncpy(KabatEntry->source,buffer+12,LARGEBUFF);
         KabatEntry->source[LARGEBUFF-1] = '\0';
      }
      else if(!strncmp(buffer,"ANNOTA AAIN",11))
      {
         int len = strlen(KabatEntry->insertsaa);
         
         strncat(KabatEntry->insertsaa,buffer+12,LARGEBUFF - len);
         KabatEntry->insertsaa[LARGEBUFF-1] = '\0';
      }
      else if(!strncmp(buffer,"ANNOTA NNIN",11))
      {
         int len = strlen(KabatEntry->insertsnuc);

         strncat(KabatEntry->insertsnuc,buffer+12,LARGEBUFF - len);
         KabatEntry->insertsnuc[LARGEBUFF-1] = '\0';
      }
      else if(!strncmp(buffer,"NOTES AA",8))
      {
         int len = strlen(KabatEntry->notesaa);

         strncat(KabatEntry->notesaa,buffer+13,LARGEBUFF - len);
         KabatEntry->notesaa[LARGEBUFF-1] = '\0';
      }
/*
      else if(!strncmp(buffer,"KABAT NUM",9))
      {
         strncpy(KabatEntry->kabatnum,buffer+13,SEQBUFF);
         KabatEntry->kabatnum[SEQBUFF-1] = '\0';
      }
*/
      else if(!strncmp(buffer,"SEQTPA",6))
      {
         int len;
         char *b;

         /* First the sequence                                          */
         len = strlen(KabatEntry->kabatseq);
#ifdef KABAT_1994
         buffer[40] = '\0';
#else
         buffer[31] = '\0';
#endif
         if(len != 0)
         {
            strncat(KabatEntry->kabatseq,"|",SEQBUFF - len);
            len++;
         }

#ifdef KABAT_1994
         strncat(KabatEntry->kabatseq,buffer+37,SEQBUFF - len);
#else
         strncat(KabatEntry->kabatseq,buffer+28,SEQBUFF - len);
#endif
         KabatEntry->kabatseq[SEQBUFF-1] = '\0';

         /* Now the numbering                                           */
         len = strlen(KabatEntry->kabatnum);
         buffer[23] = '\0';
         if(len != 0)
         {
            strncat(KabatEntry->kabatnum,"|",SEQBUFF - len);
            len++;
         }

         KILLLEADSPACES(b,buffer+18);
         strncat(KabatEntry->kabatnum,b,SEQBUFF - len);
         KabatEntry->kabatnum[SEQBUFF-1] = '\0';
      }
      else if(!strncmp(buffer,"KADBID",6))  /* Added 02.04.96           */
      {
         strncpy(KabatEntry->kadbid,buffer+12,6);
         KabatEntry->kadbid[6] = '\0';
      }
   }  while(fgets(buffer,bufflen,fp));

   return(FixSequence(KabatEntry, insert, FALSE));
}

/************************************************************************/
/*>int ReadOldKabatEntry(FILE *fp, char *buffer, int bufflen, 
                        KABATENTRY *KabatEntry, BOOL *insert)
   ----------------------------------------------------------
   Input:   FILE       *fp          Kabat file pointer
            int        bufflen      Size of buffer
   I/O:     char       *buffer      Input:  First line of Kabat entry
                                    Output: Last line of Kabat entry
   Output:  KABATENTRY *KabatEntry  Completed kabat data entry
   Returns: int                     Number of residues read
                                    -1: Unable to handle insertion

   Read a complete entry from a Kabat data file.
   
   16.03.94 Original    By: ACRM
   18.03.94 Added insert parameter. FixSequence() called here. Returns
            number of residues read
   23.03.94 Ensures strings are terminated.
   13.04.94 Changed name from ReadEntry() to ReadKabatEntry()
   18.04.94 Added -1 return; added reading of INSERTSNUC
   30.06.94 Changed name from ReadKabatEntry()
   21.07.94 Added TRUE OldFormat flag in call to FixSequence()
*/
int ReadOldKabatEntry(FILE *fp, char *buffer, int bufflen, 
                      KABATENTRY *KabatEntry, BOOL *insert)
{
   strncpy(KabatEntry->aatable,buffer+13,SMALLBUFF);
   
   while(fgets(buffer,bufflen,fp))
   {
      TERMINATE(buffer);
      
      if(!strncmp(buffer,"//",2))
         break;
      else if(!strncmp(buffer,"AMINO NAME",10))
      {
         strncpy(KabatEntry->aaname,buffer+13,SMALLBUFF);
         KabatEntry->aaname[SMALLBUFF-1] = '\0';
      }
      else if(!strncmp(buffer,"REFERENCE",9))
      {
         strncpy(KabatEntry->reference,buffer+13,LARGEBUFF);
         KabatEntry->reference[LARGEBUFF-1] = '\0';
      }
      else if(!strncmp(buffer,"AB SPECIFI",10))
      {
         strncpy(KabatEntry->antigen,buffer+13,LARGEBUFF);
         KabatEntry->antigen[LARGEBUFF-1] = '\0';
      }
      else if(!strncmp(buffer,"SPECIES",7))
      {
         strncpy(KabatEntry->species,buffer+13,SMALLBUFF);
         KabatEntry->species[SMALLBUFF-1] = '\0';
      }
      else if(!strncmp(buffer,"CLASS",5))
      {
         strncpy(KabatEntry->class,buffer+13,SMALLBUFF);
         KabatEntry->class[SMALLBUFF-1] = '\0';
      }
      else if(!strncmp(buffer,"STRAIN",6))
      {
         strncpy(KabatEntry->strain,buffer+13,SMALLBUFF);
         KabatEntry->strain[SMALLBUFF-1] = '\0';
      }
      else if(!strncmp(buffer,"SOURCE",6))
      {
         strncpy(KabatEntry->source,buffer+13,LARGEBUFF);
         KabatEntry->source[LARGEBUFF-1] = '\0';
      }
      else if(!strncmp(buffer,"INSERTSAA",9))
      {
         strncpy(KabatEntry->insertsaa,buffer+13,LARGEBUFF);
         KabatEntry->insertsaa[LARGEBUFF-1] = '\0';
      }
      else if(!strncmp(buffer,"INSERTSNUC",10))
      {
         strncpy(KabatEntry->insertsnuc,buffer+13,LARGEBUFF);
         KabatEntry->insertsnuc[LARGEBUFF-1] = '\0';
      }
      else if(!strncmp(buffer,"NOTES AA",8))
      {
         strncpy(KabatEntry->notesaa,buffer+13,LARGEBUFF);
         KabatEntry->notesaa[LARGEBUFF-1] = '\0';
      }
      else if(!strncmp(buffer,"KABAT NUM",9))
      {
         strncpy(KabatEntry->kabatnum,buffer+13,SEQBUFF);
         KabatEntry->kabatnum[SEQBUFF-1] = '\0';
      }
      else if(!strncmp(buffer,"AA SEQUEN",9))
      {
         strncpy(KabatEntry->kabatseq,buffer+13,SEQBUFF);
         KabatEntry->kabatseq[SEQBUFF-1] = '\0';
      }
   }
   
   return(FixSequence(KabatEntry, insert, TRUE));
}

/************************************************************************/
/*
   #define TEST_RDKABAT
*/

#ifdef TEST_RDKABAT
main()
{
   FILE *fp;
   int len;
   BOOL insert;
   KABATENTRY Kabat;
   
   if((fp=fopen("test.kabat","r"))!=NULL)
   {
      while((len = ReadNextKabatEntry(fp, &Kabat, &insert, FALSE))!=0)
      {
         if(len == (-1))
         {
            printf("Unable to read entry\n");
         }
         else
         {
            printf("Entry:    %s\n",   Kabat.aaname);
            printf("Date:     %d\n",   Kabat.refdate);
            printf("Sequence: %s\n\n", Kabat.sequence);
         }
      }
   }
}
#endif


