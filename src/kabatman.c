/*************************************************************************

   Program:    KabatMan
   File:       kabatman.c
   
   Version:    V2.18
   Date:       10.09.97
   Function:   Database program for reading Kabat sequence files
   
   Copyright:  (c) Andrew C. R. Martin, UCL 1994-7
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
   KabatMan is a database program for reading the Kabat antibody sequence
   data. Raw input is the Kabat `dump' files which are processed into
   a more easily handled form.

   See kabatman.doc for more information.

**************************************************************************

   Usage:
   ======
   kabatman [-f] [-q] [-o] [-v[v...]]
            -f        Force reading of the raw data files
            -q        Read data quietly
            -o        Read old format files
            -v        Increase verbosity level
            -version  Just print version info

**************************************************************************

   Revision History:
   =================
   V0.1  12.04.94 Development version
   V1.0  27.04.94 Original
   V1.1  11.05.94 Slight changes for installation of general purpose
                  routines in bioplib. Fixed bug in command line parsing.
   V1.2  11.05.94 Added new function for finding Chothia canonical classes
                  Change to RdKabat to read blanks as ? not -
   V2.0  30.06.94 Reads the new Kabat file format by default. Read the
                  old format using the -o command line flag.
   V2.1  11.07.94 Additional check made when matching light & heavy 
                  chains; at least one author's name must match.
   V2.2  21.07.94 Changed to compare multiple references.
   V2.3  23.01.95 Added support for a variability value
   V2.4  10.02.95 Skipped (New 1995 format in RdKabat.c)
   V2.5  07.03.95 Skipped (RdKabat.c puts - instead of ? at end of seq)
   V2.6  16.03.95 -version flag prints header and exits
                  BuildWhere.c - allows FORTRAN style comparison operators
                  ExecSearch.c - corrected handling of PIR option when 
                  file name is not specified.
   V2.7  16.05.95 Allows ignored keyword SOURCE in Chothia data file
   V2.8  22.06.95 Increased buffer sizes for SELECT statement
   V2.9  23.06.95 Clean compile under gcc -Wall
   V2.10 27.06.95 Skipped
   V2.11 15.12.95 Skipped
   V2.12 02.04.96 Added handling for IDs and URLs
   V2.13 11.04.96 Added header lines to data file and datafile date
                  Prints accession codes of skipped entries
   V2.14 18.04.96 Added SET URL
   V2.15 22.04.96 Fixed potential access violation in ReadChothiaData()
                  Fixed initialisation of URL string
   V2.16 07.05.96 ReadChothiaData() sets gCanonChothNum
   V2.17 29.05.96 ReadChothiaData() now frees any pre-existing data
                  Added SET CANONICAL
   V2.18 10.09.97 New SUBGROUP field

*************************************************************************/
/* Includes
*/
#define MAIN
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
static BOOL doStoreData(DATA **pData, KABATENTRY *Kabat, DATA *extra,
                        char chain, BOOL allocate, char *source,
                        BOOL GotInsert);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for the Kabat manager

   12.04.94 Original   By: ACRM
   14.04.94 Rewritten
   26.04.94 Added command line handling
   27.04.94 Added call to DisplayCopyright()
   11.05.94 Added initialisation call to GetWord() to switch on comma mode
   30.06.94 Added OldFormat flag
   21.07.94 Made OldFormat global
   23.06.95 Removed redundant variables
   11.04.96 Initialise gFileDate
*/
int main(int argc, char **argv)
{
   BOOL ForceRead = FALSE;

   strcpy(gFOF,         DEF_FOF);
   strcpy(gKabatFile,   DEF_KABAT);
   strcpy(gChothiaFile, DEF_CHOTHIA);
   strcpy(gURLFormat,   URLFORMAT);
   gInfoLevel         = DEF_INFO;
   gOldFormat         = FALSE;
   gFileDate[0]       = '\0';

   /* This causes GetWord() to return inverted commas as part of the
      words read out of the buffer. (Default mode is to strip them.)
   */
   GetWord(NULL, NULL);

   if(ParseCmdLine(argc, argv, &ForceRead))
   {
      if(ForceRead || !ReadStoredData(gKabatFile))
      {
         if(!ReadKabatData(gFOF))
         {
            fprintf(stderr,"Error: Unable to read Kabat data\n");
            return(1);
         }
         else
         {
            if(!StoreKabatData(gKabatFile))
            {
               fprintf(stderr,"Warning: Unable to store Kabat data\n");
            }
         }
      }
   }

   DisplayCopyright();
   
   if(!ReadChothiaData(gChothiaFile))
   {
      fprintf(stderr,"Warning: Unable to read Chothia data; canonicals \
not available.\n\n");
   }
   
   CommandLoop();
   
   return(0);
}


/************************************************************************/
/*>void DisplayCopyright(void)
   ---------------------------
   Displays a copyright message.

   27.04.94 Original    By: ACRM
   11.05.94 V1.1
   11.05.94 V1.2
   30.06.94 V2.0
   11.07.94 V2.1
   18.07.94 V2.2
   23.01.95 V2.3
   10.02.95 V2.4 
   07.03.95 V2.5
   16.03.95 V2.6
   16.05.95 V2.7
   22.06.95 V2.8
   23.06.95 V2.9
   27.06.95 V2.10   
   15.12.95 V2.11
   02.04.96 V2.12
   11.04.96 V2.13
   18.04.96 V2.14
   22.04.96 V2.15
   08.05.96 V2.16
   29.05.96 V2.17
   10.09.97 V2.18
*/
void DisplayCopyright(void)
{
   printf("\nKabatMan V2.18\n");
   printf("==============\n");
   printf("Copyright (c) 1994-7, Dr. Andrew C.R. Martin, University \
College London.\n\n");
   printf("This program is copyright. Any copying without the permission \
of the\n");
   printf("author is prohibited.\n\n");
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, BOOL *ForceRead)
   ---------------------------------------------------------
   Input:   int   argc         Number of arguments
            char  **argv       Argument list
   Output:  BOOL  *ForceRead   Force reading of Kabat files? (-f)
   Globals: int   gInfoLevel   Information level (-q, -v)
            BOOL  gOldFormat   Old Kabat dump format
   Returns: BOOL               Success?

   Parses the command line.

   26.04.94 Original    By: ACRM
   11.05.94 Correctly returns TRUE on successful completion (!)
   30.06.94 Added -o
   21.07.94 Changed OldFormat to a global variable.
   16.03.95 Added -version handling
   11.04.96 Also allow --version (Posix standard for long flags)
*/
BOOL ParseCmdLine(int argc, char **argv, BOOL *ForceRead)
{
   int i;
   
   argc--;       /* Skip the program name                               */
   argv++;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         if(!strcmp(argv[0],"-version") || !strcmp(argv[0],"--version"))
         {
            DisplayCopyright();
            exit(0);
         }
         
         switch(argv[0][1])
         {
         case 'f':              /* Force reading of Kabat files         */
         case 'F':
            *ForceRead = TRUE;
            break;
         case 'o':              /* Read old format files                */
         case 'O':
            gOldFormat = TRUE;
            break;
         case 'q':              /* Quiet                                */
         case 'Q':
            gInfoLevel = 0;
            break;
         case 'v':              /* Verbose: -vvv, etc to increase level */
         case 'V':
            if(gInfoLevel < 1) gInfoLevel = 1;
            for(i=1; argv[0][i]=='v' || argv[0][i]=='V'; i++)
               gInfoLevel++;
            break;
         default:
            return(FALSE);
         }
      }
      else
      {
         return(FALSE);
      }
      argc--;
      argv++;
   }

   return(TRUE);
}


/************************************************************************/
/*>BOOL ReadStoredData(char *filename)
   -----------------------------------
   Input:   char  *filename     Name of datafile
   Returns: BOOL                Success
   Globals: DATA  *gData        Linked list containing read data (output)

   Read data stored in our own internal format.

   14.04.94 Original    By: ACRM
   21.04.94 Added Reference field
   25.04.94 Added light and heavy chain special numbering
            Added initialisation of LNumbers & HNumbers
            Ensures light and heavy chains set to NULL
   26.04.94 If we fail to find the file in the current directory, look
            in the directory specified by the environment variable named
            by the constant ENV_KABATDIR
   23.06.95 Initialise p and line
   02.04.96 Added idlight & idheavy
   11.04.96 Skip comment lines from start of file, but reads date from
            them if found
*/
BOOL ReadStoredData(char *filename)
{
   DATA *p = NULL;
   FILE *fp;
   char buffer[SEQBUFF],
        FileBuff[MAXBUFF],
        *kabatdir,
        *ptr;
   int  line = 0;
   BOOL StartOfFile = TRUE;
   
   /* Open the data file                                                */
   if((fp=fopen(filename,"r"))==NULL)
   {
      /* Unable to open file so get environment variable                */
      if((kabatdir = getenv(ENV_KABATDIR))==NULL)
         return(FALSE);

      /* Construct filename and try again                               */
      sprintf(FileBuff,"%s/%s",kabatdir,filename);
      if((fp=fopen(FileBuff,"r"))==NULL)
         return(FALSE);
   }

   /* Free any currently stored data                                    */
   if(gData != NULL)
      FREELIST(gData,DATA);
   gData = NULL;
   
   /* Now read the file into a new data linked list                     */
   while(fgets(buffer,SEQBUFF,fp))
   {
      TERMINATE(buffer);

      /* Read comment lines (at start of file only)                     */
      if(StartOfFile)
      {
         if(buffer[0] == '!')
         {
            if((ptr=strstr(buffer,"DATE: "))!=NULL)
            {
               ptr += 6;
               strncpy(gFileDate,ptr,MAXBUFF);
            }
            
            continue;
         }
      }

      /* Check for start of new record and allocate space if found      */
      if(buffer[0] == '>')
      {
         StartOfFile = FALSE;
         
         line = 0;      /* Field counter within record                  */

         if(gData == NULL)
         {
            INIT(gData, DATA);
            p = gData;
         }
         else
         {
            ALLOCNEXT(p,DATA);
         }
         
         if(p==NULL)
         {
            fprintf(stderr,"Error: Unable to allocate memory for stored \
data.\n");
            return(FALSE);
         }

         p->LNumbers = p->HNumbers = NULL;
         p->light[0] = p->heavy[0] = '\0';
      }
      else
      {
         line++;
      }
      
      switch(line)
      {
      case 1:
         strcpy(p->name,       buffer);
         break;
      case 2:
         strcpy(p->fsource,    buffer);
         break;
      case 3:
         strcpy(p->source,     buffer);
         break;
      case 4:
         strcpy(p->class,      buffer);
         break;
      case 5:
         strcpy(p->antigen,    buffer);
         break;
      case 6:
         strcpy(p->reference,  buffer);
         break;
      case 7:
         if(buffer[0] != '*')
         {
            /* We've got a special numbering scheme                     */
            if((p->LNumbers = ReadSpecialNumbering(buffer))==NULL)
            {
               fprintf(stderr,"Error: Unable to allocate memory for \
Kabat numbering of stored data\n");
               return(FALSE);
            }
         }
         break;
      case 8:
         if(buffer[0] != '*')
         {
            /* We've got a special numbering scheme                     */
            if((p->HNumbers = ReadSpecialNumbering(buffer))==NULL)
            {
               fprintf(stderr,"Error: Unable to allocate memory for \
Kabat numbering of stored data\n");
               return(FALSE);
            }
         }
         break;
      case 9:
         strcpy(p->light,      buffer);
         break;
      case 10:
         strcpy(p->heavy,      buffer);
         break;
      case 11:
         strcpy(p->idlight,    buffer);
         break;
      case 12:
         strcpy(p->idheavy,    buffer);
         break;
      default:
         break;
      }
   }

   fclose(fp);
   return(TRUE);
}


/************************************************************************/
/*>char **ReadSpecialNumbering(char *buffer)
   -----------------------------------------
   Input:   char  *buffer      Buffer read from data file
   Returns: char  **           Array of pointers to char containing
                               the residue numbers. NULL if error.

   Allocates an array of strings and copies in the residue numbers for
   a Kabat special numbering scheme from a buffer containing the numbers
   separated by spaces and terminated by a *. This occurs when inserts
   have occurred over the standard Kabat scheme.

   25.04.94 Original    By: ACRM
*/
char **ReadSpecialNumbering(char *buffer)
{
   char **NumBuff = NULL,
        *np       = NULL,
        *numbers  = buffer;
   int  i         = 0;
   
   if((NumBuff = (char **)malloc(MAXKABATSEQ * sizeof(char *)))==NULL)
      return(NULL);
   
   while((np = strchr(numbers,' ')) != NULL)
   {
      /* Terminate string at next space                                 */
      *np = '\0';

      /* Break out if we hit a *                                        */
      if(numbers[0] == '*')
         break;

      /* Allocate space to store number                                 */
      if((NumBuff[i] = (char *)malloc(5 * sizeof(char)))==NULL)
      {
         /* Error in allocation. Free up memory and return              */
         int j;
         
         for(j=0; j<i; j++)
            free(NumBuff[j]);
         free(NumBuff);
         return(NULL);
      }
      
      /* Copy in the number                                             */
      strcpy(NumBuff[i],numbers);

      /* Step on to next entry                                          */
      *np = ' ';
      numbers = np+1;
      i++;
   }

   /* Put a NULL in the last position in the array                      */
   NumBuff[i] = NULL;

   return(NumBuff);
}


/************************************************************************/
/*>BOOL StoreKabatData(char *filename)
   -----------------------------------
   Input:   char   *filename    Filename to write
   Returns: BOOL                Success in opening file?
   Globals: DATA   gData        The Kabat data linked list (input)

   Writes the data contained in the global data linked list to a file

   14.04.94 Original   By: ACRM
   21.04.94 Added reference field
   25.04.94 Writes light chain and heavy chain numbers.
   02.04.96 Added idlight & idheavy
   11.04.96 Writes comment lines at the top which include the date of
            writing and the file format version
*/
BOOL StoreKabatData(char *filename)
{
   DATA   *p;
   FILE   *fp;
   int    i;
   time_t TheTime;

   if(gData != NULL)
   {
      if((fp=fopen(filename,"w"))==NULL)
      {
         fprintf(stderr,"Warning: Unable to open file %s\n",filename);
         return(FALSE);
      }

      /* Determine the current date                                     */
      time(&TheTime);
      strcpy(gFileDate, ctime(&TheTime));
      TERMINATE(gFileDate);

      /* Put a header into the file including a date stamp              */
      fprintf(fp,"! KabatMan data file.\n");
      fprintf(fp,"! KabatMan is Copyright 1994-6, \
Dr. Andrew C.R. Martin, UCL\n");
      fprintf(fp,"! FILE VERSION: 3.0\n");
      fprintf(fp,"! CREATION DATE: %s\n",gFileDate);
      fprintf(fp,"!\n");
      
      for(p=gData; p!=NULL; NEXT(p))
      {
         fprintf(fp,">\n");
         fprintf(fp,"%s\n",p->name);
         fprintf(fp,"%s\n",p->fsource);
         fprintf(fp,"%s\n",p->source);
         fprintf(fp,"%s\n",p->class);
         fprintf(fp,"%s\n",p->antigen);
         fprintf(fp,"%s\n",p->reference);

         if(p->LNumbers)
         {
            for(i=0; p->LNumbers[i]!=NULL; i++)
               fprintf(fp,"%s ",p->LNumbers[i]);
         }
         fprintf(fp,"*\n");
         if(p->HNumbers)
         {
            for(i=0; p->HNumbers[i]!=NULL; i++)
               fprintf(fp,"%s ",p->HNumbers[i]);
         }
         fprintf(fp,"*\n");
         
         fprintf(fp,"%s\n",p->light);
         fprintf(fp,"%s\n",p->heavy);
         fprintf(fp,"%s\n",p->idlight);
         fprintf(fp,"%s\n",p->idheavy);
      }
      
      fclose(fp);
   }
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL ReadKabatData(char *FoF)
   -----------------------------
   Input:   char *FoF        Kabat file of files
   Globals: BOOL gOldFormat  Read old format Kabat files
   Returns: BOOL             Success?
   Globals: DATA *gData      Data linked list

   Read kabat datafiles into the global gData linked list

   12.04.94 Original    By: ACRM
   13.04.94 Modified to work with L files read into memory
            Checks minimum length of sequence
   14.04.94 Modified to take file name as parameter rather than pointer
            Gets source from filename
   18.04.94 Handles entries skipped by ReadNextKabatEntry()
   21.04.94 ReadLFiles() has source parameter
            Gets LCClass array from filename and passes it into 
            ReadLFiles()
   22.04.94 Added Insert handling
   26.04.94 If we fail to find a file in the current directory, look
            in the directory specified by the environment variable named
            by the constant ENV_KABATDIR
            Changed gInfoLevel setting
   30.06.94 Added OldFormat flag
   21.07.94 Made OldFormat global
   11.04.96 On skipped chains, also print accession code
*/
BOOL ReadKabatData(char *FoF)
{
   char       buffer[MAXBUFF],
              FileBuff[MAXBUFF],
              *KabatDir,
              HFile[160],
              LFile[MAXLFILES][160],
              *p = buffer,
              source[MAXBUFF];
   FILE       *fp = NULL,
              *fpH = NULL,
              *fpL[MAXLFILES];
   int        NLFile = 0,
              i,
              LCClass[MAXLFILES],
              nseq;
   KABATENTRY KabatH;
   DATA       *KabatLData[MAXLFILES];
   BOOL       GotInsert;
   
   /* Get the Kabat directory environment for future use                */
   KabatDir = getenv(ENV_KABATDIR);

   /* Open Kabat file of files for reading                              */
   if((fp=fopen(FoF,"r"))==NULL)
   {
      if(KabatDir == NULL)
      {
         printf("Error: Unable to open Kabat file of files: %s\n",FoF);
         return(FALSE);
      }
         
      sprintf(FileBuff,"%s/%s",KabatDir,FoF);
      if((fp=fopen(FileBuff,"r"))==NULL)
      {
         printf("Error: Unable to open Kabat file of files: %s\n",FoF);
         return(FALSE);
      }
   }

   while(fgets(buffer,MAXBUFF-1,fp))
   {
      TERMINATE(buffer);

      /* Process if it's not a blank line                               */
      KILLLEADSPACES(p,buffer);
      if(strlen(p))
      {
         if(gInfoLevel >= 1) printf("Processing files: %s\n",p);
         
         /* Got a line specifying a group of files                      */
         p = buffer;
         p = GetWord(p, HFile);
         for(i=0; i<MAXLFILES; i++)
            p = GetWord(p, LFile[i]);

         /* Open these files                                            */
         if(HFile[0] && HFile[0] != '-')
         {
            if((fpH = fopen(HFile,"r"))==NULL) 
            {
               if(KabatDir == NULL)
               {
                  fprintf(stderr,"Error: Unable to open heavy chain \
file: `%s'\n", HFile);
                  fclose(fp);
                  return(FALSE);
               }
               
               sprintf(FileBuff,"%s/%s",KabatDir,HFile);
               if((fpH=fopen(FileBuff,"r"))==NULL)
               {
                  fprintf(stderr,"Error: Unable to open heavy chain \
file: `%s'\n", HFile);
                  fclose(fp);
                  return(FALSE);
               }
            }
         }
      
         for(NLFile=0; NLFile<MAXLFILES; NLFile++)
         {
            if(LFile[NLFile][0])
            {
               if((fpL[NLFile]=fopen(LFile[NLFile],"r"))==NULL)
               {
                  if(KabatDir == NULL)
                  {
                     fprintf(stderr,"Error: Unable to open light chain \
file: `%s'\n", LFile[NLFile]);
                     fclose(fp);
                     return(FALSE);
                  }
                  
                  sprintf(FileBuff,"%s/%s",KabatDir,LFile[NLFile]);
                  if((fpL[NLFile]=fopen(FileBuff,"r"))==NULL)
                  {
                     fprintf(stderr,"Error: Unable to open light chain \
file: `%s'\n", LFile[NLFile]);
                     fclose(fp);
                     return(FALSE);
                  }
               }
               
               /* Set the class flag based on the filename              */
               LCClass[NLFile] = 0;
               UPPER(LFile[NLFile]);
               if(strstr(LFile[NLFile],"KAPPA")!=NULL) 
                  LCClass[NLFile] = CLASS_KAPPA;
               else if(strstr(LFile[NLFile],"LAMBDA")!=NULL) 
                  LCClass[NLFile] = CLASS_LAMBDA;
            }
            else
            {
               break;
            }
         }

         /* Find out the source from the filename                       */
         if(HFile[0] != '-')
            GetSource(HFile, source);
         else
            GetSource(LFile[0], source);

         /* Read data in from the light chain files                     */
         if(!ReadLFiles(fpL, NLFile, KabatLData, source, LCClass))
         {
            fprintf(stderr,"Error: Failed to make temporary store for \
L-chain data\n");
            fclose(fp);
            return(FALSE);
         }
      
         if(fpH)            /* If we have a heavy chain file            */
         {
            /* Read entries from the heavy chain file                   */
            while((nseq=ReadNextKabatEntry(fpH,&KabatH,&GotInsert,
                                           gOldFormat)) != 0)
            {
               if(nseq == (-1))
               {
                  if(gInfoLevel >= 1)
                     printf("Skipped H chain for %s (%s)\n", 
                            KabatH.aaname, KabatH.kadbid);
               }
               else if(nseq > MINSEQ)
               {
                  /* Store entry and search for matching light chain    */
                  if((gData = StoreHAndMatchL(gData, KabatH, KabatLData,
                                              NLFile, source, GotInsert))
                     ==(DATA *)(-1))
                  {
                     fprintf(stderr,"Error: Failed to store H-chain \
data\n");
                     fclose(fp);
                     return(FALSE);
                  }
               }
            }
         }
      
         /* Store any unmatched light chain entries                     */
         if((gData = StoreUnmatchedL(gData, KabatLData,
                                     NLFile, source))==(DATA *)(-1))
         {
            fprintf(stderr,"Error: Failed to store L-chain data\n");
            fclose(fp);
            return(FALSE);
         }
         
         /* Close files                                                 */
         fclose(fpH);
         
         for(i=0; i<NLFile; i++)
         {
            fclose(fpL[i]);
            fpL[i] = NULL;
         }

         /* Free memory for the L files                                 */
         for(i=0; i<NLFile; i++)
            FREELIST(KabatLData[i], DATA);
         
      }  /* Not a blank line in the file                                */
   }  /* while() line in File of files                                  */

   fclose(fp);
   return(TRUE);
}


/************************************************************************/
/*>BOOL StoreKabatInData(DATA **pData, KABATENTRY Kabat, char chain,
                         char *source, BOOL GotInsert)
   -----------------------------------------------------------------
   I/O:     DATA       **pData   Address of data linked list
   Input:   KABATENTRY Kabat     Item to be appended to Kabat data list
            char       chain     Chain for data (l or h)
            char       *source   Source derived from filename
            BOOL       GotInsert Kabat has an insertion
   Returns: BOOL                 Success

   Interface to the doStoreData() routine which adds an item to the 
   data linked list allocating space for the entry

   12.04.94 Original    By: ACRM
   13.04.94 Added NULL parameter to doStoreData()
   14.04.94 Added source parameter
*/
BOOL StoreKabatInData(DATA **pData, KABATENTRY Kabat, char chain,
                      char *source, BOOL GotInsert)
{
   return(doStoreData(pData, &Kabat, NULL, chain, TRUE, source, 
                      GotInsert));
}


/************************************************************************/
/*>BOOL StoreDataInData(DATA **pData, DATA *p, char chain, char *source)
   ---------------------------------------------------------------------
   I/O:     DATA       **pData Address of data linked list
   Input:   DATA       *p      Item to be appended to Kabat data list
            char       chain   Chain for data (l or h)
            char       *source Source derived from filename
   Returns: BOOL               Success

   Interface to the doStoreData() routine which adds an item to the 
   data linked list allocating space for the entry

   12.04.94 Original    By: ACRM
   13.04.94 Added NULL parameter to doStoreData()
   14.04.94 Added source parameter
   22.04.94 Added insert (FALSE) parameter to doStoreData()
*/
BOOL StoreDataInData(DATA **pData, DATA *p, char chain, char *source)
{
   BOOL ret;
   
   ret = doStoreData(pData, NULL, p, chain, TRUE, source, FALSE);
   
   return(ret);
}


/************************************************************************/
/*>void AddKabatToData(DATA *Data, KABATENTRY Kabat, char chain,
                       BOOL DoInsert)
   -------------------------------------------------------------
   I/O:     DATA       *Data    Data linked list
   Input:   KABATENTRY Kabat    Item to be added to Kabat data list
            char       chain    Chain for data (l or h)
            BOOL       DoInsert Kabat has insertion
   Returns: BOOL                Success

   Interface to the doStoreData() routine which inserts data into the
   last item of the data linked list.

   12.04.94 Original    By: ACRM
   13.04.94 Added NULL parameter to doStoreData()
   14.04.94 Added source parameter to doStoreData()
   22.04.94 Added DoInsert parameter
*/
void AddKabatToData(DATA *Data, KABATENTRY Kabat, char chain, 
                    BOOL DoInsert)
{
   doStoreData(&Data, &Kabat, NULL, chain, FALSE, NULL, DoInsert);
}


/************************************************************************/
/*>void AddDataToData(DATA *Data, DATA *Extra, char chain)
   -------------------------------------------------------
   I/O:     DATA       *Data   Data linked list
   Input:   DATA       *Extra  Item to be added to Kabat data list
            char       chain   Chain for data (l or h)
   Returns: BOOL               Success

   Interface to the doStoreData() routine which inserts data into the
   last item of the data linked list.

   13.04.94 Original    By: ACRM
   14.04.94 Added source parameter to doStoreData()
   22.04.94 Added insert (FALSE) parameter to doStoreData()
*/
void AddDataToData(DATA *Data, DATA *Extra, char chain)
{
   doStoreData(&Data, NULL, Extra, chain, FALSE, NULL, FALSE);
}


/************************************************************************/
/*>static BOOL doStoreData(DATA **pData, KABATENTRY *Kabat, DATA *extra,
                           char chain, BOOL allocate, char *source,
                           BOOL DoInsert)
   -------------------------------------------------------------------
   I/O:     DATA       **pData    Address of data linked list
   Input:   KABATENTRY *Kabat     Item to be appended (or added) to Kabat 
                                  data list or NULL
            DATA       *extra     Item to be appended (or added) to Kabat
                                  data list of NULL
            char       chain      Chain for data (l or h)
            BOOL       allocate   Should memory be allocated
            char       *source    Source derived from filename
            BOOL       GotInsert  Kabat has insert?
   Globals: BOOL       gOldFormat Old format Kabat files
   Returns: BOOL                  Success

   Adds an item to the data linked list allocating space for the entry
   or inserts data into the last entry in the current list

   12.04.94 Original    By: ACRM
   14.04.94 Added source parameter
   22.04.94 Added GotInsert handling
   25.04.94 Added initialisation of LNumbers & HNumbers
            Ensures light and heavy chains set to NULL
   18.07.94 Only copies source name in if not VARIOUS
   02.04.96 Initialises new id strings to NULL
*/
static BOOL doStoreData(DATA **pData, KABATENTRY *Kabat, DATA *extra,
                        char chain, BOOL allocate, char *source,
                        BOOL GotInsert)
{
   static DATA *p=NULL;

   if(Kabat !=NULL || extra != NULL)
   {
      if(allocate)
      {
         if(*pData == NULL)
         {
            INIT((*pData), DATA);
            p = *pData;
         }
         else
         {
            ALLOCNEXT(p,DATA);
         }

         if(p!=NULL)
         {
            p->LNumbers   = p->HNumbers   = NULL;
            p->light[0]   = p->heavy[0]   = '\0';
            p->idlight[0] = p->idheavy[0] = '\0';
         }
      }
   
      if(p==NULL) return(FALSE);
   }

   if(Kabat != NULL) CopyKabatToData(p, *Kabat, chain, GotInsert);
   if(extra != NULL) CopyDataToData(p, extra, chain);

   if(Kabat != NULL || extra != NULL)
   {
      if(source != NULL && upstrncmp(source,"VARIOUS",7))
         strcpy(p->source, source);
      
      UPPER(p->source);
   }
   
   return(TRUE);
}


/************************************************************************/
/*>DATA *StoreHAndMatchL(DATA *Data, KABATENTRY KabatH, DATA *KabL[], 
                        int NLFile, char *source, BOOL GotInsert)
   ----------------------------------------------------------------
   Input:   DATA       *Data     Kabat data linked list
            KABATENTRY KabatH    Kabat heavy chain entry to store
            DATA       *KabL[]   Array of linked lists containing LCs
            int        NLFile    Number of light chain file pointers
            char       *source   Source derived from filename
            BOOL       GotInsert The Kabat entry has an insert
   Returns: DATA       *         Pointer to new linked list (-1 if error)

   Stores a heavy chain entry into the data linked list then searches the
   light chain files for matching entries and adds their sequence data.
   Any stored light chains are flagged.

   12.04.94 Original    By: ACRM
   13.04.94 Modified to use L/C data in memory
   14.04.94 Added source parameter
   21.04.94 DATA.active is now an array
   22.04.94 Added GotInsert parameter and handling
   26.04.94 Changed gInfoLevel setting
   11.07.94 Now checks that an author name matches as well as the 
            sequence name. If authors are stupid enough to give the same
            name to two different antibodies which they have sequenced,
            then there's not a lot we can do!
   21.07.94 Now uses RefCheck() rather than NameCheck(). This allows
            multiple references
   23.06.95 Removed redundant variables
*/
DATA *StoreHAndMatchL(DATA *Data, KABATENTRY KabatH, DATA *KabL[], 
                      int NLFile, char *source, BOOL GotInsert)
{
   int        i;
   DATA       *p;
   
   if(!StoreKabatInData(&Data, KabatH, 'H', source, GotInsert))
      return((DATA *)(-1));

   if(gInfoLevel >= 2)
      printf("Stored H chain for %s\n", KabatH.aaname);
   
   /* For each of the light chain data structures                       */
   for(i=0; i<NLFile; i++)
   {
      /* For each entry in this linked list                             */
      for(p=KabL[i]; p!=NULL; NEXT(p))
      {
         if(!strcmp(KabatH.aaname, p->name) && 
            RefCheck(KabatH.reference, p->reference))
         {
            /* This light chain entry matches the heavy chain entry     */
            AddDataToData(Data, p, 'L');
   
            if(gInfoLevel >= 2)
               printf("Stored L chain for %s\n", p->name);

            /* Flag this one as used                                    */
            p->active[0] = TRUE;
            return(Data);
         }
      }
   }
   
   return(Data);
}


/************************************************************************/
/*>DATA *StoreUnmatchedL(DATA *Data, DATA *KabL[], int NLFile, 
                         char *source)
   -----------------------------------------------------------
   Input:   DATA       *Data    Kabat data linked list
            DATA       *KabL[]  Array of linked lists containing LCs
            int        NLFile   Number of light chain file pointers
            char       *source  Source derived from filename
   Returns: DATA       *        Pointer to new linked list (-1 if error)

   Stores any light chain entries into the data linked list which
   haven't previously bee stored as H chain partners

   12.04.94 Original    By: ACRM
   13.04.94 Modified to use L/C data in memory
   25.04.94 DATA.active is an array!
   26.04.94 Changed gInfoLevel setting
   23.06.95 Removed redundant variables
*/
DATA *StoreUnmatchedL(DATA *Data, DATA *KabL[], int NLFile, char *source)
{
   int        i;
   DATA       *p;

   /* For each of the light chain files                                 */
   for(i=0; i<NLFile; i++)
   {
      /* For each entry in this light chain data structure              */
      for(p=KabL[i]; p!=NULL; NEXT(p))
      {
         if(!p->active[0])
         {
            if(!StoreDataInData(&Data, p, 'L', source))
               return((DATA *)(-1));

            if(gInfoLevel >= 2)
               printf("Stored L chain for %s\n", p->name);
         }
      }
   }
   
   return(Data);
}


/************************************************************************/
/*>void CopyKabatToData(DATA *p, KABATENTRY Kabat, char chain,
                        BOOL GotInsert)
   -----------------------------------------------------------
   Input:   KABATENTRY  Kabat      Full Kabat data entry
            char        chain      Chain indicator (h or l)
            BOOL        GotInsert  Kabat entry has an insert
   Output:  DATA        *p         Completed Kabat data structure

   Copies required data from a complete KABATENTRY structure into a
   DATA structure used by this program.

   12.04.94 Original    By: ACRM
   13.04.94 Clears active flag
   21.04.94 DATA.active is now an array
            DATA.fsource now has the source as described in the file
            Added reference
   22.04.94 Added GotInsert parameter
   18.07.94 Copies into source as well as fsource
   21.07.94 Added gOldFormat flag to calls to BuildKabatNumbering()
   02.04.96 Added kadbid
*/
void CopyKabatToData(DATA *p, KABATENTRY Kabat, char chain, 
                     BOOL GotInsert)
{
   strcpy(p->antigen,    Kabat.antigen);
   strcpy(p->class,      Kabat.class);
   strcpy(p->name,       Kabat.aaname);
   strcpy(p->fsource,    Kabat.source);
   strcpy(p->source,     Kabat.source);
   strcpy(p->reference,  Kabat.reference);

   p->active[0] = FALSE;

   if(chain=='l' || chain=='L')
   {
      strcpy(p->light, Kabat.sequence);
      if(GotInsert) p->LNumbers = BuildKabatNumbering(Kabat,gOldFormat);
      strcpy(p->idlight, Kabat.kadbid);
   }
   else
   {
      strcpy(p->heavy, Kabat.sequence);
      if(GotInsert) p->HNumbers = BuildKabatNumbering(Kabat,gOldFormat);
      strcpy(p->idheavy, Kabat.kadbid);
   }
}


/************************************************************************/
/*>void CopyDataToData(DATA *p, DATA *extra, char chain)
   -----------------------------------------------------
   Input:   DATA        *extra     Data structure entry to copy
            char        chain      Chain indicator (h or l)
   Output:  DATA        *p         Completed Kabat data structure

   Copies required data from a one local DATA structure into another.

   13.04.94 Original    By: ACRM
   21.04.94 DATA.active is now an array
            Added source & reference
   22.04.94 Added LNumbers & HNumbers
   21.07.94 Fixed LNumbers and HNumbers only to be copied if needed (!)
   02.04.96 Added kadbid
*/
void CopyDataToData(DATA *p, DATA *extra, char chain)
{
   strcpy(p->antigen,    extra->antigen);
   strcpy(p->class,      extra->class);
   strcpy(p->name,       extra->name);
   strcpy(p->source,     extra->source);
   strcpy(p->fsource,    extra->fsource);
   strcpy(p->reference,  extra->reference);

   p->active[0] = FALSE;

   if(chain=='l' || chain=='L')
   {
      p->LNumbers =         extra->LNumbers;
      strcpy(p->light,      extra->light);
      strcpy(p->idlight,    extra->idlight);
   }
   else
   {
      p->HNumbers =         extra->HNumbers;
      strcpy(p->heavy,      extra->heavy);
      strcpy(p->idheavy,    extra->idheavy);
   }
}


/************************************************************************/
/*>BOOL ReadLFiles(FILE *fpL[], int NLFile, DATA *KabatLData[], 
                   char *source, int *LCClass)
   ------------------------------------------------------------
   Input:   FILE   *fpL[]         Array of light chain file pointers
            int    NLFile         Number of light chains
            char   *source        Source name from filename
            int    *LCClass       Array of filename based classes
   Globals: BOOL   gOldFormat     Old format files
   Output:  DATA   *KabatLData[]  Array of light chain data linked lists

   Read in the light chain files into DATA linked lists

   13.04.94 Original    By: ACRM
   18.04.94 Handles entries skipped by ReadNextKabatEntry()
   21.04.94 Added source parameter and setting of source value
            Added LCClass parameter and setting of class
   22.04.94 Added Insert handling
   25.04.94 Added initialisation of LNumbers & HNumbers
            Ensures light and heavy chains set to NULL
   26.04.94 Changed gInfoLevel setting
   30.06.94 Handles OldFormat flag
   18.07.94 Only takes source from filename only if not VARIOUS
   21.07.94 Made OldFormat flag global
   23.06.95 Initialised p
   11.04.96 Also prints accession code for skipped entries
*/
BOOL ReadLFiles(FILE *fpL[], int NLFile, DATA *KabatLData[], char *source,
                int *LCClass)
{
   int        i,
              nseq;
   KABATENTRY KabatEntry;
   DATA       *p = NULL;
   BOOL       GotInsert;

   for(i=0; i<NLFile; i++)
      KabatLData[i] = NULL;
   
   for(i=0; i<NLFile; i++)
   {
      while((nseq=ReadNextKabatEntry(fpL[i],&KabatEntry,&GotInsert,
                                     gOldFormat)) != 0)
      {
         if(nseq == (-1))
         {
            if(gInfoLevel >= 1)
               printf("Skipped L chain for %s (%s)\n", 
                      KabatEntry.aaname, KabatEntry.kadbid);
         }
         else if(nseq > MINSEQ)
         {
            if(KabatLData[i] == NULL)
            {
               INIT((KabatLData[i]), DATA);
               p = KabatLData[i];
            }
            else
            {
               ALLOCNEXT(p, DATA);
            }

            if(p==NULL) return(FALSE);

            p->LNumbers = p->HNumbers = NULL;
            p->light[0] = p->heavy[0] = '\0';
            
            CopyKabatToData(p, KabatEntry, 'L', GotInsert);

            /* Only copy source from filename if not VARIOUS            */
            if(upstrncmp(source,"VARIOUS",7))
               strcpy(p->source, source);
            
            UPPER(p->source);

            /* Set class from filename if not specified in file         */
            if(!(p->class[0]))
            {
               if(LCClass[i] == CLASS_LAMBDA)
                  strcpy(p->class,"LAMBDA");
               else if(LCClass[i] == CLASS_KAPPA)
                  strcpy(p->class,"KAPPA");
            }
         }
      }
   }

   return(TRUE);
}


/************************************************************************/
/*>void GetSource(char *filename, char *source)
   --------------------------------------------
   Input:   char  *filename      Complete filename
   Output:  char  *source        The source derived from the filename

   Takes the first part of the true filename (after any /, : or ])
   and truncates at the next . or _   Outputs this as the source.

   14.04.94 Original    By: ACRM
*/
void GetSource(char *filename, char *source)
{
   char *p;
   
   GetFilestem(filename, source);
   if((p=strchr(source,'.')) != NULL)
      *p = '\0';
   if((p=strchr(source,'_')) != NULL)
      *p = '\0';
}


/************************************************************************/
/*>void CommandLoop(void)
   ----------------------
   The main command loop which does the work for the Kabat manager.

   18.04.94 Original    By: ACRM
   21.04.94 Corrected nesting after ; or .
            Added prompt
            Uses upstr(n)cmp() to do comparisons
            Added redirection of output
            Added SET variable handling
   25.04.94 No longer reports blank lines as syntax error
   22.06.95 Doubled buffer size to stop problems with people entering
            whole chains
*/
void CommandLoop(void)
{
   char buffer[2*MAXBUFF],
        *p;
   int Mode = 0;
   
   printf("KABATMAN> ");
   
   while(fgets(buffer,MAXBUFF,stdin))
   {
      /* Tidy up the buffer                                             */
      TERMINATE(buffer);
      KILLLEADSPACES(p,buffer);


      if(p[0] == ';' || p[0] == '.')  /* Cause the search to be run     */
      {
         Mode = 0;
         ExecuteSearch(NULL);
      }
      else if(p[0] == '>')            /* Run search and redirect        */
      {
         Mode = 0;
         KILLLEADSPACES(p,p+1);
         ExecuteSearch(p);
      }
      else if(p[0] == '\0')           /* Blank line                     */
      {
         while(FALSE) ;
      }
      else
      {
         /* See if we're entering a new mode                            */
         if(!upstrncmp(p,"SET",3))
         {
            if(Mode != 0)
               fprintf(stderr,"SET only allowed at main prompt\n");
            else
               HandleSetCommand(p);
         }
         else if(!upstrncmp(p,"SELECT",6))
         {
            Mode = 1;
            ClearSelect();
            BuildSelect(p);
         }
         else if(!upstrncmp(p,"WHERE",5))
         {
            Mode = 2;
            ClearWhere();
            BuildWhere(p);
         }
         else if(!upstrncmp(p,"FROM",4))
         {
            Mode = 3;
            BuildFrom(p);
         }
         else if(!upstrncmp(p,"QUIT",4) || !upstrncmp(p,"EXIT",4))
         {
            ClearSelect();
            ClearWhere();
            return;
         }
         else     /* Add info to the current mode                       */
         {
            switch(Mode)
            {
            case 1:
               BuildSelect(p);
               break;
            case 2:
               BuildWhere(p);
               break;
            case 3:
               BuildFrom(p);
               break;
            default:
               fprintf(stderr,"Error: (Syntax) %s\n",p);
               break;
            }
         }
      }

      switch(Mode)
      {
      case 0:
         printf("KABATMAN> ");
         break;
      case 1:
         printf("SELECT> ");
         break;
      case 2:
         printf("WHERE> ");
         break;
      case 3:
         printf("FROM> ");
         break;
      }
   }

   exit(0);
}


/************************************************************************/
/*>BOOL RefCheck(char *ref1, char *ref2)
   -------------------------------------
   Input:   char  *ref1          First reference list
            char  *ref2          Second reference list
   Returns: BOOL                 TRUE: a match found   
                                 FALSE: no matches

   Runs through two reference lists and returns TRUE if one author
   matches in the 2 lists. This is used to compare the author lists
   in the reference field when a light and heavy chain have been found
   with the same name. If at least one author's name matches, it is
   fairly safe to assume that the chains are co-expressed.

   Each reference within each list is delineated by a |

   21.07.94 Original    By: ACRM
*/
BOOL RefCheck(char *inref1, char *inref2)
{
   char *r1, 
        *r2,
        *r1e = (char *)1,
        *r2e = (char *)1,
        ref1[LARGEBUFF],
        ref2[LARGEBUFF];

   /* Copy the strings into new buffers so we don't mess up the original*/
   strcpy(ref1,inref1);
   strcpy(ref2,inref2);
   
   for(r1 = ref1; r1e != NULL && *r1; r1 = r1e+1)
   {
      r1e = strchr(r1,'|');
      if(r1e != NULL) *r1e = '\0';
      
      for(r2 = ref2; r2e != NULL && *r2; r2 = r2e+1)
      {
         r2e = strchr(r2,'|');
         if(r2e != NULL) *r2e = '\0';

         if(NameCheck(r1, r2))
            return(TRUE);
      }
   }
   
   return(FALSE);
}
  

/************************************************************************/
/*>BOOL NameCheck(char *instr1, char *instr2)
   ------------------------------------------
   Input:   char  *instr1        First string
            char  *instr2        Second string
   Returns: BOOL                 TRUE: one word matched
                                 FALSE: no matches

   Runs through two strings and returns TRUE if at least one word
   matches in the 2 strings. This is used to compare the author lists
   in the reference field when a light and heavy chain have been found
   with the same name. If at least one author's name matches, it is
   pretty safe to assume that the chains are co-expressed.

   11.07.94 Original    By: ACRM
   23.06.95 Removed redundant variables
*/
BOOL NameCheck(char *instr1, char *instr2)
{
   char *p1, *p2,
        string1[LARGEBUFF], 
        string2[LARGEBUFF],
        word1[160], word2[160];
   
   /* Copy the strings into new buffers so we don't mess up the original*/
   strcpy(string1,instr1);
   strcpy(string2,instr2);

   /* First terminate the strings at the first ( . This marks the end of
      the author list and the start of the reference itself (i.e. the
      date)
   */
   if((p1 = strchr(string1,'('))!=NULL) *p1 = '\0';
   if((p2 = strchr(string2,'('))!=NULL) *p2 = '\0';

   /* Step through the first string extracting a word at a time         */
   p1 = string1;
   do
   {
      p1 = GetFullStopWord(p1, word1);

      /* Now step through the second string                             */
      p2 = string2;
      do
      {
         p2 = GetFullStopWord(p2, word2);

         /* Compare the 2 words; if they match, return TRUE             */
         if(!strcmp(word1, word2))
            return(TRUE);

         p2 = EatInitials(p2);
      }  while(p2 != NULL);

      p1 = EatInitials(p1);
   }  while(p1 != NULL);

   return(FALSE);
}


/************************************************************************/
/*>char *GetFullStopWord(char *buffer, char *word)
   -----------------------------------------------
   Input:   char  *buffer     The buffer from which to extract a word
   Output:  char  *word       Word extracted from buffer
   Returns: char  *           Pointer into buffer after word extracted

   Gets a word from a string delimited by a full stop. This is used to
   extract a surname and first initial from a reference string.

   11.07.94 Original    By: ACRM
*/
char *GetFullStopWord(char *buffer, char *word)
{
   char        *p;
   int         i         = 0,
               j         = 0;
   
   /* Return a blank string if the input buffer is NULL                 */
   word[0] = '\0';
   if(buffer==NULL) return(NULL);

   /* Remove leading spaces                                             */
   KILLLEADSPACES(p, buffer);

   /* Copy up to next comma or white space                              */
   for(i=0; p[i]; i++)
   {
      /* If we get a . then our word has ended.                         */
      if(p[i]=='.')
         break;

      word[j++] = p[i];
   }

   /* Terminate output string                                           */
   word[j] = '\0';

   /* Move p on to the next word                                        */
   p += i;                  /* Move p onto the next character           */
   KILLLEADSPACES(p,p);     /* Strip any leading spaces                 */

   if(*p == '\0') p = NULL;
   
   return(p);
}


/************************************************************************/
/*>char *EatInitials(char *buffer)
   -------------------------------
   Input:   char  *buffer     Input string
   Returns: char  *           Pointer into string after initials removed

   GetFullStopWord() removes a surname and first initial from the
   reference string. This routine is used to eat characters up to the
   next white space or comma. It then removes any more white spaces,
   commas or ampersands used to delimit the last author.

   11.07.94 Original    By: ACRM
*/
char *EatInitials(char *buffer)
{
   char *p;
   
   /* Step along until we hit a space or a comma                        */
   for(p=buffer; (*p!='\0' && *p!=' ' && *p!='\t' && *p!=','); p++) ;

   /* Step over any space, comma or ampersand                           */
   while(*p==' ' || *p==',' || *p=='&') p++;
   
   /* If we've hit the end of the string return NULL, else return 
      pointer to current position
   */
   return((*p=='\0') ? (NULL) : p);
}


/************************************************************************/
void BuildFrom(char *buffer)
{
   
}

/************************************************************************/
/*>void HandleSetCommand(char *buffer)
   -----------------------------------
   Input:   char  *buffer       Buffer containing SET command
   Globals: int   gLoopMode     Set by LOOP  {KABAT|ABM|CHOTHIA}
            int   gInfoLevel    Set by LEVEL {value}
            BOOL  gShowInserts  Set by INSERTS {ON|OFF}
            REAL  gVariability  Set by VARiability {value}
            BOOL  gHTML         Set by HTML {ON|OFF}
            char  *gURLFormat   Set by URL {format}
            CHOTHIA *gChothia   The Chothia data; refreshed by 
                                   SET CANONICAL {type}

   Handles commands which set variables

   21.04.94 Original    By: ACRM
   25.04.94 Added INSERTS {ON|OFF}
   26.04.94 Changed default gInfoLevel
   23.01.95 Added handling of SET VARiability {value}
   02.04.96 Added SET HTML {ON|OFF}
   18.04.96 Added SET URL {format}
   29.05.96 Added SET CANONICAL {type}
*/
void HandleSetCommand(char *buffer)
{
   char word[MAXBUFF],
        value[MAXBUFF],
        *p = buffer;
   
   p = GetWord(p,word);
   do
   {
      if(upstrcmp(word,"SET"))
      {
         /* It's not `SET', so it's a variable name, so get the value   */
         p = GetWord(p,value);
         
         /* Carry out required action for each variable                 */
         if(!upstrncmp(word,"LEVEL",5))
         {
            if(sscanf(value,"%d",&gInfoLevel)==0) gInfoLevel = DEF_INFO;
         }
         else if(!upstrncmp(word,"LOOP",4))
         {
            if(!upstrncmp(value,"KAB",3))
               gLoopMode = LOOP_KABAT;
            else if(!upstrncmp(value,"ABM",3))
               gLoopMode = LOOP_ABM;
            else if(!upstrncmp(value,"CHOTH",5))
               gLoopMode = LOOP_CHOTHIA;
            else
               fprintf(stderr,"Error: Unknown loop mode (%s)\n",value);
         }
         else if(!upstrncmp(word,"INSERT",6))
         {
            if(!upstrncmp(value,"ON",2))
               gShowInserts = TRUE;
            else
               gShowInserts = FALSE;
         }
         else if(!upstrncmp(word,"VAR",3))
         {
            if(sscanf(value,"%lf",&gVariability)==0) 
               gVariability = DEF_VARIABILITY;
         }
         else if(!upstrncmp(word,"HTML",4))
         {
            if(!upstrncmp(value,"ON",2))
               gHTML = TRUE;
            else
               gHTML = FALSE;
         }
         else if(!upstrncmp(word,"URL",3))
         {
            int len;

            if(!upstrncmp(value,"DEF",3))
            {
               strcpy(gURLFormat, URLFORMAT);
            }
            else if(strncmp(value,"href=",5))
            {
               fprintf(stderr,"Error: Invalid URL. Default restored\n");
               strcpy(gURLFormat, URLFORMAT);
            }
            else
            {
               strcpy(gURLFormat,"\"<a ");
               strncpy(gURLFormat+4, value, MAXBUFF);
               gURLFormat[MAXBUFF-1] = '\0';
               
               len = strlen(gURLFormat);
               if(len<MAXBUFF-11)
               {
                  strcpy(gURLFormat+len, "%s>%s</a>\"");
               }
fprintf(stderr,"%s\n",gURLFormat);
            }
         }
         else if(!upstrncmp(word,"CAN",3))
         {
            if(!upstrncmp(value,"DEF",3))
            {
               ReadChothiaData(gChothiaFile);
            }
            else
            {
               char *filename;
               int  length = strlen(value) + strlen(gChothiaFile) + 8;
               
               if((filename = (char *)malloc(length * sizeof(char)))
                  ==NULL)
               {
                  fprintf(stderr,"No memory for filename buffer when \
reading key residue file\n");
                  return;
               }
               
               LOWER(value);
               sprintf(filename,"%s.%s",gChothiaFile, value);
               
               ReadChothiaData(filename);

               free(filename);
            }
         }
         else
         {
            fprintf(stderr,"Error: Unknown variable (%s)\n",word);
         }
      }

      /* Get the next word out of the string                            */
      p = GetWord(p,word);
   }  while(p!=NULL);
}


/************************************************************************/
/*>BOOL ReadChothiaData(char *filename)
   ------------------------------------
   Input:   char *filename   The Chothia data filename
   Returns: BOOL             Success?

   Reads a Chothia canonical definition file. This file has the format:
   LOOP loopid class length
   resid types
   resid types
   ...

   11.05.94 Original
   16.05.95 Now allows keyword SOURCE in a LOOP definition which is 
            ignored. Uses upstrncmp() rather than strncmp(). Allows
            a # to start comments.
   22.04.96 Changed check on MAXCHOTHRES from > to >=
   07.05.96 Checks for CHOTHIANUM keyword in the file and, if found,
            sets the gCanonChothNum flag
   08.05.96 Was doing the last else clause after reading the CHOTHIANUM
            keyword.
   29.05.96 Now frees any pre-existing data so this can be called
            multiple times.
*/
BOOL ReadChothiaData(char *filename)
{
   FILE    *fp;
   char    buffer[MAXBUFF],
           FileBuff[MAXBUFF],
           word[40],
           *chp,
           *kabatdir;
   CHOTHIA *p = NULL;
   int     count = 0;
   
   /* Open the data file                                                */
   if((fp=fopen(filename,"r"))==NULL)
   {
      /* Unable to open file so get environment variable                */
      if((kabatdir = getenv(ENV_KABATDIR))==NULL)
         return(FALSE);

      /* Construct filename and try again                               */
      sprintf(FileBuff,"%s/%s",kabatdir,filename);
      if((fp=fopen(FileBuff,"r"))==NULL)
         return(FALSE);
   }

   /* Free the current Chothia data if there is any                     */
   if(gChothia != NULL)
   {
      FREELIST(gChothia, CHOTHIA);
      gChothia = NULL;
   }

   /* This flag indicatess whether the file contains Chothia or Kabat
      numbering
   */
   gCanonChothNum = FALSE;

   while(fgets(buffer,160,fp))
   {
      TERMINATE(buffer);
      
      if(strlen(buffer) && buffer[0] != '!' && buffer[0] != '#')
      {
         /* Ignore the keyword SOURCE                                   */
         if(!upstrncmp(buffer,"SOURCE",6))
         {
            continue;
         }
         else if(!upstrncmp(buffer,"CHOTHIANUM",10))
         {
            gCanonChothNum = TRUE;
         }
         else if(!upstrncmp(buffer,"LOOP",4))   /* Start of entry       */
         {
            /* Terminate the previous list of resnums                   */
            if(p!=NULL)
               strcpy(p->resnum[count],"-1");
            
            /* Allocate space in linked list                            */
            if(gChothia == NULL)
            {
               INIT(gChothia,CHOTHIA);
               p = gChothia;
            }
            else
            {
               ALLOCNEXT(p,CHOTHIA);
            }
            if(p==NULL) return(FALSE);
            
            /* Strip out the word LOOP                                  */
            chp = GetWord(buffer,word);
            /* Get the loop id                                          */
            chp = GetWord(chp,p->LoopID);
            /* Get the class name                                       */
            chp = GetWord(chp,p->class);
            /* Get the loop length                                      */
            chp = GetWord(chp,word);
            sscanf(word,"%d",&(p->length));
            
            /* Set the resnum counter to zero                           */
            count = 0;
         }
         else
         {
            /* Not the start of an entry, so must be a resid/type pair  */
            if(p!=NULL)
            {
               chp = GetWord(buffer,p->resnum[count]);
               chp = GetWord(chp,p->restype[count]);
               if(++count >= MAXCHOTHRES)
               {
                  fprintf(stderr,"Error: (Reading Chothia file) Too many \
Chothia key residues.\n");
                  return(FALSE);
               }
            }
         }
      }
   }
   
   /* Terminate the previous list of resnums                            */
   if(p!=NULL)
      strcpy(p->resnum[count],"-1");
   
   return(TRUE);
}

