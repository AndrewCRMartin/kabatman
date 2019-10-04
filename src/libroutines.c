/*************************************************************************

   Program:    KabatMan
   File:       libroutines.c
   
   Version:    V1.2
   Date:       22.03.96
   Function:   Library routines
   
   Copyright:  (c) SciTech Software 1993 and Andrew C. R. Martin 1994
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure and Modelling Unit,
               Department of Biochemistry and Molecular Biology,
               University College,
               Gower Street,
               London.
   EMail:      martin@biochem.ucl.ac.uk
               
**************************************************************************

   The routines in this file are part of the Bioplib package. They may
   only be used for the purposes of recompiling KabatMan.

   You may obtain a licence form to use these files outside KabatMan via
   the WWW address:
   http://www.biochem.ucl.ac.uk/~martin/text/BioplibLicence.ps

**************************************************************************

   Description:
   ============

   These routines form part of the bioplib C protein handling library.
   For more information on this library of routines, please contact
   the author.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.1 22.03.96 Added QueryStrStr()

*************************************************************************/
/* Includes
*/
#include "bioplib/macros.h"  /* DNAtoAA() */
#include <string.h>          /* DNAtoAA() */
#include <ctype.h>           /* DNAtoAA() */
#include "bioplib/SysDefs.h" /* GetWord() */

#include "libroutines.h"

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

static char *sAACode[4][4] =
{  { "FFLL", "SSSS", "YYXX", "CCXW" },   /* TTX, TCX, TAX, TGX          */
   { "LLLL", "PPPP", "HHQQ", "RRRR" },   /* CTX, CCX, CAX, CGX          */
   { "IIIM", "TTTT", "NNKK", "SSRR" },   /* ATX, ACX, AAX, AGX          */
   { "VVVV", "AAAA", "DDEE", "GGGG" }    /* GTX, GCX, GAX, GGX          */
} ;


/************************************************************************/
/*>char DNAtoAA(char *dna)
   -----------------------
   Input:   char  *dna        DNA/RNA codon
   Returns: char              1-letter amino acid code (X=termination)

   Converts a nucleic acid codon to the 1-letter amino acid equivalent.
   Termination codons are returned as X. No special action is taken
   for initiation codons.

   18.04.94 Original    By: ACRM
*/
char DNAtoAA(char *dna)
{
   char buffer[8], *p;
   int idx1, idx2, idx3;
   
   KILLLEADSPACES(p,dna);
   
   strncpy(buffer,p,8);
   buffer[7] = '\0';
   UPPER(buffer);

   idx1 = NUCINDEX(buffer[0]);
   idx2 = NUCINDEX(buffer[1]);
   idx3 = NUCINDEX(buffer[2]);
   
   return(sAACode[idx1][idx2][idx3]);
}

/************************************************************************/
/*>int TrueSeqLen(char *sequence)
   ------------------------------
   Input:   char  *sequence    A sequence containing deletions
   Returns: int                Length without deletions

   Scans a 1-letter code sequence and calculate the length without
   `-' or ` ' residues

   14.04.94 Original    By: ACRM
*/
int TrueSeqLen(char *sequence)
{
   int length = 0,
       i = 0;
   
   for(i=0; sequence[i]; i++)
   {
      if(sequence[i] != '-' && sequence[i] != ' ')
         length++;
   }
   
   return(length);
}

/************************************************************************/
/*>int KnownSeqLen(char *sequence)
   -------------------------------
   Input:   char  *sequence    A sequence containing deletions
   Returns: int                Length without deletions

   Scans a 1-letter code sequence and calculate the length without
   `-', ` ' or '?' residues

   13.05.94 Original    By: ACRM
*/
int KnownSeqLen(char *sequence)
{
   int length = 0,
       i = 0;
   
   for(i=0; sequence[i]; i++)
   {
      if(sequence[i] != '-' && sequence[i] != ' ' && sequence[i] != '?')
         length++;
   }
   
   return(length);
}

/************************************************************************/
/*>void GetFilestem(char *filename, char *stem)
   --------------------------------------------
   Input:   char  *filename      Complete filename
   Output:  char  *stem          The filestem

   Extracts the filestem from a complete filename. Should work under
   Unix, VMS, MS-DOS, AmigaDOS, etc.

   14.04.94 Original    By: ACRM
*/
void GetFilestem(char *filename, char *stem)
{
   char *p, 
        *q;

   q = filename;

   /* First step past a ] if found (VMS)                                */
   if((p=strchr(q,']'))!=NULL)      q=p+1;

   /* Step past any colons/double colons (VMS, AmigaDOS, Unix, MS-DOS)  */
   while((p=strchr(q,':'))!=NULL)   q=p+1;

   /* Step past any / (Unix, AmigaDOS)                                  */
   while((p=strchr(q,'/'))!=NULL)   q=p+1;
   
   /* Step past any \ (MS-DOS)                                          */
   while((p=strchr(q,'\\'))!=NULL)  q=p+1;
   
   /* We should now have the actual filename, with the path removed.
      Copy it into our output array
   */
   strcpy(stem, q);
   
   /* Terminate at the last . to remove the extension.                  */
   q = stem;
   for(p = q + strlen(q) - 1; p!=q; p--)
   {
      if(*p == '.')
      {
         *p = '\0';
         break;
      }
   }
}

/************************************************************************/
/*>int upstrcmp(char *word1, char *word2)
   --------------------------------------
   Input:   char *word1     First word
            char *word2     Second word
   Returns: int             0 if strings match or offset of first 
                            mismatched character

   Like strcmp(), but upcases each character before comparison

   20.04.94 Original   By: ACRM
*/
int upstrcmp(char *word1, char *word2)
{
   int i;

   for(i=0; word1[i] && word2[i]; i++)
      if((islower(word1[i])?toupper(word1[i]):word1[i]) != 
         (islower(word2[i])?toupper(word1[i]):word2[i])) return(i+1);

   if(word1[i] || word2[i]) return(i+1);
   
   return(0);
}

/************************************************************************/
/*>int upstrncmp(char *word1, char *word2, int ncomp)
   --------------------------------------------------
   Input:   char *word1     First word
            char *word2     Second word
            int  ncomp      Number of characters to compare
   Returns: int             0 if strings match or offset of first 
                            mismatched character

   Like strncmp(), but upcases each character before comparison

   20.04.94 Original   By: ACRM
*/
int upstrncmp(char *word1, char *word2, int ncomp)
{
   int i;

   for(i=0; i<ncomp; i++)
   {
      if(!word1[i] || !word2[i]) return(i+1);
      
      if((islower(word1[i])?toupper(word1[i]):word1[i]) != 
         (islower(word2[i])?toupper(word1[i]):word2[i])) return(i+1);
   }
   
   return(0);
}

/************************************************************************/
/*>char *GetWord(char *buffer, char *word)
   ---------------------------------------
   Input:   char  *buffer        String buffer
   Output:  char  *word          Word extracted from buffer
   Returns: char  *              Pointer into buffer after word

   Extracts a comma or white space delimited word from a buffer and
   returns a pointer into the buffer after the word has been removed.
   ' or " may be used to define a word containing white space.

   By default, the bounding ' or " are not returned as part of the word.
   Calling the routine with 2 NULL pointers will switch the mode such
   that the bounding inverted commas are returned.
   Calling again will switch back to the default mode.

   12.04.94 Original    By: ACRM
   26.04.94 Now handles groups of characters in inverted commas as a
            single word.
   11.05.94 Added special call; if called with 2 NULL pointers, switches
            the treatment of ' and " to be included or not in the output
*/
char *GetWord(char *buffer, char *word)
{
   char        *p;
   int         i         = 0,
               j         = 0,
               Commas    = 0;     /* 1: In singles, 2: In doubles       */
   BOOL        GotComma  = FALSE;
   static BOOL CommaMode = FALSE;
   
   /* Check for special calls                                           */
   if(buffer == NULL && word == NULL)
   {
      CommaMode = !CommaMode;
      return(NULL);
   }

   /* Return a blank string if the input buffer is NULL                 */
   word[0] = '\0';
   if(buffer==NULL) return(NULL);

   /* Remove leading spaces                                             */
   KILLLEADSPACES(p, buffer);

   /* Copy up to next comma or white space                              */
   for(i=0; p[i]; i++)
   {
      GotComma = FALSE;  /* Assume this character is not an inv comma   */
      
      /* Set state machine if it is a ' or "                            */
      switch(p[i])
      {
      case '\'':                 /* Got a single inverted comma         */
         switch(Commas)
         {
         case 0:                 /* Not yet in commas, set state 1      */
            Commas   = 1;
            GotComma = TRUE;
            break;
         case 1:                 /* In single commas, set state 0       */
            Commas   = 0;
            GotComma = TRUE;
            break;
         case 2:                 /* In double commas, ignore single     */
            break;
         }
         break;
      case '"':                  /* Got a double inverted comma         */
         switch(Commas)
         {
         case 0:                 /* Not yet in commas, set state 2      */
            Commas = 2;
            GotComma = TRUE;
            break;
         case 1:                 /* In single commas, ignore double     */
            break;
         case 2:                 /* In double commas, set state 0       */
            Commas = 0;
            GotComma = TRUE;
            break;
         }
         break;
      default:                   /* Got some other character            */
         break;
      }
      
      /* If we're not in commas (state 0) and we get a , or white space
         then our word has ended.
      */
      if(!Commas && (p[i]==',' || p[i]==' ' || p[i]=='\t'))
         break;

      /* If this wasn't a state-changing ' or " then store it           */
      if(CommaMode || !GotComma)
         word[j++] = p[i];
   }

   /* Terminate output string                                           */
   word[j] = '\0';

   /* Move p on to the next word                                        */
   p += i;                  /* Move p onto the next character           */
   KILLLEADSPACES(p,p);     /* Strip any leading spaces                 */
   if(*p == ',') p++;       /* Kill a comma if found                    */

   if(*p == '\0') p = NULL;
   
   return(p);
}

/*************************************************************************

   Program:    
   File:       throne.c
   
   Version:    V1.1R
   Date:       11.03.94
   Function:   Convert between 1 and 3 letter aa codes
   
   Copyright:  (c) SciTech Software 1993
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0372) 275775
   EMail:      INTERNET: amartin@scitec.adsp.sub.org
                         martin@bsm.bioc.ucl.ac.uk
               UUCP:     ...{uunet|rutgers}!cbmehq!cbmuk!scitec!amartin
               JANET:    martin@uk.ac.ucl.bioc.bsm
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  29.09.92 Original    By: ACRM   
   V1.1  11.03.94 Added PCA, ASX and GLX to translation table.
                  PCA translates to E
                  Added routines to handle asx/glx

*************************************************************************/
/* Includes
*/
#include <string.h>

/************************************************************************/
/* Defines and macros
*/
#define NUMAAKNOWN 24

/************************************************************************/
/* Globals
*/

/* N.B. The order in sTab1[] and sTab3[] must be the same and they must
   end with X/UNK
*/
static char sTab1[]    = {'A','C','D','E','F',
                          'G','H','I','K','L',
                          'M','N','P','Q','R',
                          'S','T','V','W','Y',
                          'E','B','Z','X'
                         };
static char sTab3[][8] = {"ALA ","CYS ","ASP ","GLU ","PHE ",
                          "GLY ","HIS ","ILE ","LYS ","LEU ",
                          "MET ","ASN ","PRO ","GLN ","ARG ",
                          "SER ","THR ","VAL ","TRP ","TYR ",
                          "PCA ","ASX ","GLX ","UNK "
                         };

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>char throne(char *three)
   ------------------------
   Input:   char  *three    Three letter code
   Returns: char            One letter code

   Converts 3-letter code to 1-letter code.
   Handles ASX and GLX as X
   
   29.09.92 Original    By: ACRM
   11.03.94 Modified to handle ASX and GLX in the tables
*/
char throne(char *three)
{
   int j;

   if(three[2] == 'X')
      return('X');

   for(j=0;j<NUMAAKNOWN;j++)
      if(!strncmp(sTab3[j],three,3)) return(sTab1[j]);

   /* Only get here if the three letter code was not found              */
   return('X');
}

/************************************************************************/
/*>char thronex(char *three)
   -------------------------
   Input:   char  *three    Three letter code
   Returns: char            One letter code

   Converts 3-letter code to 1-letter code.
   Handles ASX and GLX as B and Z.
   
   29.09.92 Original    By: ACRM
*/
char thronex(char *three)
{
   int j;

   for(j=0;j<NUMAAKNOWN;j++)
      if(!strncmp(sTab3[j],three,3)) return(sTab1[j]);

   /* Only get here if the three letter code was not found              */
   return('X');
}

/************************************************************************/
/*>char *onethr(char one)
   ----------------------
   Input:   char  one     One letter code
   Returns: char  *       Three letter code (padded to 4 chars with a 
                          space)

   Converts 1-letter code to 3-letter code (actually as 4 chars).
   N.B. The last entry in the static conversion table *must* be UNK.
   
   07.06.93 Original    By: ACRM
*/
char *onethr(char one)
{
   int j;

   for(j=0;j<NUMAAKNOWN;j++)
      if(sTab1[j] == one) return(sTab3[j]);

   /* Only get here if the one letter code was not found                */
   return(sTab3[NUMAAKNOWN-1]);
}

/************************************************************************/
/*>char *QueryStrStr(char *string, char *substring)
   ------------------------------------------------
   This is like strstr() but allows a ? character in the substring
   which matches any character.

   15.12.95 Original    By: ACRM
*/
char *QueryStrStr(char *string, char *substring)
{
   int  i, j, 
        lenstr, 
        lensubstr;
   char *retval = NULL;

   /* Get lengths of the 2 strings                                      */
   lenstr    = strlen(string);
   lensubstr = strlen(substring);
   
   /* Walk through the main string                                      */
   for(i=0; i<=lenstr-lensubstr; i++)
   {
      /* Set return value to point to current position in main string   */
      retval = string + i;
      
      /* Walk through substring                                         */
      for(j=0; j<lensubstr; j++)
      {
         /* If there is a mismatch, set return to NULL and break out of
            the substring counting
         */
         if((substring[j] != '?') && 
            (substring[j] != '%') && 
            (string[i+j] != substring[j]))
         {
            retval = NULL;
            break;
         }
      }

      /* If retval hasn't been set to NULL, then we've got a match, so
         return pointer
      */
      if(retval)
         return(retval);
   }

   /* This should always be NULL...                                     */
   return(retval);
}
