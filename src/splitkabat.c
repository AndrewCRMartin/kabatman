/*************************************************************************

   Program:    splitkabat
   File:       splitkabat.c
   
   Version:    V1.2
   Date:       08.02.95
   Function:   Split Kabat new format database files into types
   
   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 1994-5
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0)1372 275775
   EMail:      martin@biochem.ucl.ac.uk
               andrew@stagleys.demon.co.uk
               
**************************************************************************


   This program is copyright. 

   Any copying without the express permission of the author is illegal.

**************************************************************************

   Description:
   ============
   The old format Kabat dump files were maintained in a split form.
   i.e. All human antibody heavy chains were in a file called
   human.ig.hc, kappa chains in the file human.ig.kappa, etc. This
   allowed one to select only the antibodies and increased one's chance
   of matching up light and heavy chains.

   The new format files no longer maintain this separation; they are
   simply split 1000 entries at a time; currently there does appear to
   be some organisation in the first files, but this, undoubtedly, will
   go as the database increases in size.

   This program splits out the antibodies (by rejecting any entries
   with keywords from the gSkipKeys[] array in their DEFINI records)
   and writes out files of the appropriate species and chain type.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  30.06.94 Original
   V1.1  18.07.94 Now checks SPECIES records as well as DEFINI for the
                  type
   V1.2  08.02.95 Changes for new 1995 release format.

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "bioplib/SysDefs.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/
#define NHEADER    16
#define MAXBUFF    160
#define MAXDEFLINE 8

typedef struct
{
   char *source,
        *chain,
        *file;
}  FILEINFO;

/************************************************************************/
/* Globals
*/
FILEINFO gFileInfo[] =              /* Map keywords to filename         */
{  { "CAT ", "HEAVY", "cat.ig.hc"},
   { "CAT ", "LAMBDA", "cat.ig.lambda"},
   { "CAT ", "KAPPA", "cat.ig.kappa"},
   { "CAT\n", "HEAVY", "cat.ig.hc"},
   { "CAT\n", "LAMBDA", "cat.ig.lambda"},
   { "CAT\n", "KAPPA", "cat.ig.kappa"},
   { "CHICKEN", "HEAVY", "chicken.ig.hc"},
   { "CHICKEN", "LAMBDA", "chicken.ig.lambda"},
   { "CHICKEN", "KAPPA", "chicken.ig.kappa"},
   { "DOG ", "HEAVY", "dog.ig.hc"},
   { "DOG ", "LAMBDA", "dog.ig.lambda"},
   { "DOG ", "KAPPA", "dog.ig.kappa"},
   { "DOG\n", "HEAVY", "dog.ig.hc"},
   { "DOG\n", "LAMBDA", "dog.ig.lambda"},
   { "DOG\n", "KAPPA", "dog.ig.kappa"},
   { "FROG", "HEAVY", "frog.ig.hc"},
   { "FROG", "LAMBDA", "frog.ig.lambda"},
   { "FROG", "KAPPA", "frog.ig.kappa"},
   { "GOPHER", "HEAVY", "gopher.ig.hc"},
   { "GOPHER", "LAMBDA", "gopher.ig.lambda"},
   { "GOPHER", "KAPPA", "gopher.ig.kappa"},
   { "HORSE", "HEAVY", "horse.ig.hc"},
   { "HORSE", "LAMBDA", "horse.ig.lambda"},
   { "HORSE", "KAPPA", "horse.ig.kappa"},
   { "HUMAN", "HEAVY", "human.ig.hc"},
   { "HUMAN", "LAMBDA", "human.ig.lambda"},
   { "HUMAN", "KAPPA", "human.ig.kappa"},
   { "MOUSE", "HEAVY", "mouse.ig.hc"},
   { "MOUSE", "LAMBDA", "mouse.ig.lambda"},
   { "MOUSE", "KAPPA", "mouse.ig.kappa"},
   { "RABBIT", "HEAVY", "rabbit.ig.hc"},
   { "RABBIT", "LAMBDA", "rabbit.ig.lambda"},
   { "RABBIT", "KAPPA", "rabbit.ig.kappa"},
   { "RAT", "HEAVY", "rat.ig.hc"},
   { "RAT", "LAMBDA", "rat.ig.lambda"},
   { "RAT", "KAPPA", "rat.ig.kappa"},
   { "SHARK", "HEAVY", "shark.ig.hc"},
   { "SHARK", "LAMBDA", "shark.ig.lambda"},
   { "SHARK", "KAPPA", "shark.ig.kappa"},
   { "SHEEP", "HEAVY", "sheep.ig.hc"},
   { "SHEEP", "LAMBDA", "sheep.ig.lambda"},
   { "SHEEP", "KAPPA", "sheep.ig.kappa"},
   { "VARIOUS", "HEAVY", "various.ig.hc"},
   { "VARIOUS", "LAMBDA", "various.ig.lambda"},
   { "VARIOUS", "KAPPA", "various.ig.kappa"},
   { NULL, NULL, NULL}  
}  ;

char *gSkipKeys[] =               /* Any DEFINI having these is skipped */
{  "MINIGENE",
   "RECEPTOR",
   "PSEUDOGENE",
   "HISTOCOMPATIBILITY",
   "CONSTANT",
   "ADHESION",
   "MICROGLOBULIN",
   "COMPLEMENT",
   "T-CELL",
   "SIGNAL",
   NULL
}  ;
     
/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
FILE *SetOutputFile(char *FileName);
char *FindFileType(FILE *in, char header[NHEADER][MAXBUFF], 
                   int *HeaderLine);
void WriteEntry(FILE *in, FILE *out, char header[NHEADER][MAXBUFF], 
                int HeaderLine);
void Usage(void);
void SkipEntry(FILE *in);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for splitting Kabat files

   29.06.94 Original    By:ACRM
   30.06.94 Added call to SkipEntry()
*/
int main(int argc, char **argv)
{
   FILE *in  = NULL,
        *out = NULL;
   char header[NHEADER][MAXBUFF],
        *FileName        = NULL;
   int  HeaderLine       = 0;
   
   if(argc < 2)
   {
      Usage();
      return(1);
   }
   
   if((in=fopen(argv[1],"r"))==NULL)
   {
      fprintf(stderr,"Unable to open input file: %s\n",argv[1]);
      return(1);
   }
   
   while(!feof(in))
   {
      if((FileName = FindFileType(in, header, &HeaderLine))==((char *)-1))
         break;

      if(FileName)
      {
         if((out = SetOutputFile(FileName))==NULL)
         {
            fprintf(stderr,"Warning: Unable to open output file %s\n",
                    FileName);
            SkipEntry(in);
         }
         else
         {
            WriteEntry(in, out, header, HeaderLine);
         }
      }
      else
      {
         D("No appropriate output file\n");
         SkipEntry(in);
      }
   }
   
   return(0);
}

/************************************************************************/
/*>FILE *SetOutputFile(char *FileName)
   -----------------------------------
   Opens the specified file for append access if it's not already open

   29.06.94 Original    By:ACRM
*/
FILE *SetOutputFile(char *FileName)
{
   static char *CurrentFileName = NULL;
   static FILE *out             = NULL;

   if(FileName == NULL)
   {
      D("NULL filename given\n");
      
      if(out != NULL)
      {
         D("Closing old output file\n");
         fclose(out);
      }
      
      out = NULL;
   }
   else
   {
      if(FileName != CurrentFileName)
      {
         if(out != NULL)
            fclose(out);
         
         out=fopen(FileName,"a");
         CurrentFileName = FileName;
      }
   }
   
   return(out);
}
      
/************************************************************************/
/*>char *FindFileType(FILE *in, char header[NHEADER][MAXBUFF], 
                      int *HeaderLine)
   -----------------------------------------------------------
   Reads in the first lines from a Kabat entry and selects the
   appropriate output file from the DEFINI line

   29.06.94 Original    By:ACRM
   30.06.94 Handles multiple DEFINI lines
            Handles gSkipKeys[]
   18.07.94 Checks specifes lines as well as DEFINI lines
   08.02.95 Just looks for SEQ to end a header rather than SEQRES
            for latest release of the d/b. Upcases the header buffer
            before testing for keywords.
*/
char *FindFileType(FILE *in, char header[NHEADER][MAXBUFF], 
                   int *HeaderLine)
{
   int  i,
        DefLine[MAXDEFLINE],
        NDefLine = 0;
   char buffer[MAXDEFLINE * MAXBUFF];

   *HeaderLine = 0;
   buffer[0] = '\0';

   while(fgets(header[*HeaderLine],MAXBUFF,in))
   {
      /* N.B. TERMINATE() must *not* be called                          */
#ifdef DEBUG
      printf("%d: %s",*HeaderLine, header[*HeaderLine]);
#endif
      
      if(!strncmp(header[*HeaderLine],"DEFINI",6) ||
         !strncmp(header[*HeaderLine],"SPECIE",6))
         DefLine[NDefLine++] = *HeaderLine;
      if(!strncmp(header[*HeaderLine],"AANAME",6))
      {
         (*HeaderLine)++;
         break;
      }
      if(!strncmp(header[*HeaderLine],"SEQ",3))
      {
         NDefLine = 0;
         break;
      }
      
      if(++(*HeaderLine) >= NHEADER)
      {
         NDefLine = 0;
         break;
      }
   }

   if(!(*HeaderLine))
      return((char *)(-1));

   if(NDefLine == 0)
   {
      fprintf(stderr,"File missing DEFINI or AANAME line in entry:\n");
      for(i=0; i<(*HeaderLine); i++)
         fputs(header[i],stderr);
      return(NULL);
   }

   /* Build all DEFINI lines into a single buffer                       */
   for(i=0; i<NDefLine; i++)
      strcat(buffer,header[DefLine[i]]);

   /* Upcase the buffer                                                 */
   for(i=0; i<strlen(buffer); i++)
   {
      if(islower(buffer[i]))
         buffer[i] = toupper(buffer[i]);
   }
   
   /* See if any of the skip keys are present, if so return NULL        */
   for(i=0; gSkipKeys[i] != NULL; i++)
   {
      if(strstr(buffer,gSkipKeys[i]))
         return(NULL);
   }

   /* Search for specific types                                         */
   for(i=0; gFileInfo[i].file != NULL; i++)
   {
      if(strstr(buffer,gFileInfo[i].source) &&
         strstr(buffer,gFileInfo[i].chain))
         return(gFileInfo[i].file);
   }
   
   /* Now do `various' types                                            */
   for(i=0; gFileInfo[i].file != NULL; i++)
   {
      if(!strncmp(gFileInfo[i].source,"VARIOUS",7) &&
         strstr(buffer,gFileInfo[i].chain))
         return(gFileInfo[i].file);
   }

   return(NULL);
}
      
/************************************************************************/
/*>void WriteEntry(FILE *in, FILE *out, char header[NHEADER][MAXBUFF], 
                   int HeaderLine)
   -------------------------------------------------------------------
   First writes the header information, then copies the rest of the
   entry.

   29.06.94 Original    By:ACRM
*/
void WriteEntry(FILE *in, FILE *out, char header[NHEADER][MAXBUFF], 
                int HeaderLine)
{
   int i;
   char buffer[MAXBUFF];
   
   for(i=0; i<HeaderLine; i++)
      fputs(header[i], out);

   while(fgets(buffer,MAXBUFF,in))
   {
      fputs(buffer,out);
      if(!strncmp(buffer,"RECEND",6))
         break;
   }
}
   
/************************************************************************/
/*>void SkipEntry(FILE *in)
   ------------------------
   Skips the rest of an entry for which we haven't got an output file
   30.06.94 Original    By: ACRM
*/
void SkipEntry(FILE *in)
{
   char buffer[MAXBUFF];
   
   while(fgets(buffer,MAXBUFF,in))
   {
      if(!strncmp(buffer,"RECEND",6))
         break;
   }
}
  
/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   29.06.94 Original    By:ACRM
   08.02.95 V1.2
*/
void Usage(void)
{
   fprintf(stderr,"\nSplitKabat V1.2 (c) 1994-5 Dr. Andrew C.R. Martin, \
UCL\n");
   fprintf(stderr,"Usage: splitkabat <kabatfile>\n");
   fprintf(stderr,"Splits new format Kabat files into separate files\n");
   fprintf(stderr,"based on source and chain type.\n\n");
}
