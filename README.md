KabatMan
========

Version 2.26 - 7th October 2019
-------------------------------

*Copyright (c) 1994-2019*

*Prof. Andrew C.R. Martin / UCL / Reading*

*Department of Structural and Molecular Biology*

*University College London* 

Document created :   27th April, 1994


Introduction
------------

KabatMan is a simple database program for accessing the Kabat antibody
sequence database. Kabat, et al. have, for some years, published
collected antibody sequence data in the form of a book, but as the
number of sequences solved increases, the need for electronic
distribution becomes more and more necessary.

The Kabat data have now been made available as `dump' files by FTP
from the following site:

```
ncbi.nlm.nih.gov:repository/kabat/dump
```

As of V2.0, KabatMan now reads the new format introduced on 2nd April
and first made available during June, 1994. The old format files
may still be read by using the -o flag. In January 1995, the new format
changed somewhat and V2.4 reads this new format (and not the June 1994
format).

The new format files are available from:

```
        ftp://ncbi.nlm.nih.gov/repository/kabat/fixlen
```
or
```
        ftp://ftp.ebi.ac.uk/pub/databases/kabat/fixlen
```

The format of these data files, while easy to read by eye, is somewhat
difficult to scan by computer. Currently, no software is supplied to
derive data from the files, though it is the intention of the Kabat
database maintainers to provide a more computer-friendly environment
in the future. KabatMan (`Kabat Manager`) has been written to fill
the present gap. It assembles the antibody data from the Kabat dump
files and allows searches to be performed on the data. Standard PIR
sequence files may be created and other information (such as CDR
lengths) may be obtained while allowing many constraints to be placed
on the selection of sequences.

Building and installing
-----------------------

To build and install the program simply type
```
./install.sh
```

This will install the program in your `~/bin` directory. Simply edit
the first line of the install script if you wish to install elsewhere.

This will add an environment variable, `KABATDIR`, to your `~/.bashrc`
file:
```
   export KABATDIR=$HOME/bin/share/kabatman
```

If you are using *tcsh*, then you will need to add the
environment variable to your `~/.tcsh` file instead.
```
        setenv KABATDIR $HOME/bin/share/kabatman
```

Note that the program will first look in the current directory for the
data files and then in the directory specified by `KABATDIR`.

Command Syntax
--------------

KabatMan is normally invoked simply with the command `kabatman'.
Interaction is then interactive. Three options allow additional control
over the program:
```
        kabatman [-version] [-f] [-v[v...]] [-q] [-o]
```
(Square brackets indicate optional items; you don't type them!)

The `-version` (or `--version`) flag causes KabatMan to print its 
copyright/version message and exit.

The `-f` flag forces the reading of new raw Kabat data files even if the
kabat.dat file exists (see Section *The Data*)

The `-v` flag sets the information level to 2 or more (default is 1). Each
additional v increments the information level by one. A level of 2
causes the program to display the name of each sequence as it is
stored.

The `-q` flag sets the information level to 0 (default is 1). This
prevents the program from displaying the names of Kabat sequence files
as they are processed.

The `-o` flag causes the program to read old Kabat format files.

To leave the program, type QUIT or EXIT at the prompt.


The Data
--------

Since the data have been static since 2000, you probably won't want to
rebuild the data, but, if you do, then run the script
```
mirror/buildkabat.sh
```
to rebuild the Kabat database and then run
```
mirror/installkabat.sh
```
to install it.

Input data for KabatMan are the Kabat dump files. Either the old format
or the new format (available after June 1994) may be used. By default,
the new format is read.

The new format files unlike the old format dump files are not separated
into different species, protein type and chain types. It is necessary
to split the files up first in order for KabatMan to make associations
between appropriate light and heavy chains. The program SplitKabat
performs this function; it will read the new dump files and write them
out into separate files with names like human.ig.hc (as used for the
old format dump files). Any non-antibody sequences (T-cell receptors,
etc.) will be rejected during this procedure.

As of the 1995 release of the Kabat database, each entry in the 
database comes in a separate file. The database is supplied as a
set of compressed tar files which unpack into a set of directories:
000, 001, 002, 003, etc. containing files named by their KADBID entry 
code (000001, 000002, etc.) The directory names reflect the first
3 digits of the KADBID code.

All files must be merged into a single datafile for use by KabatMan.
This is achieved using a shell script of the form:

```
       touch kabat.jan95.dat
       for dir in 0*
       do
          cat $dir/* >>kabat.jan95.dat
       done
```

This must be run from the directory containing all the Kabat
sub-directories (000, 001, 002, etc.) and concatenates the data into a
file called `kabat.jan95.dat` (You should change the date as required).
You should ensure that there are no other files or directory starting
with the digit 0 and, when the database increases to include KADBID
entries starting with a 1, 2, 3, etc rather than just a 0 the script
should be suitably modified.

In addition, the Kabat database may be distributed with a `temp`
directory containing incremental additions and corrections to the
database. Again this is distributed as a compressed tar file
(`temp.tar.Z`) which unpacks to create a subdirectory called `temp`. These
files should be concatenated into a single file, `kabat.patch`, which is
applied to the main `kabat.jan95.dat` file created above using the
program `PatchKabat`. This program may also be used to delete entries
from the main `kabat.jan95.dat` file.

The `KADBID` numbers of entries to be deleted should be placed into a
file, one to a line. `PatchKabat` is then run using:
```
   patchkabat [-d delfile] [-p patchfile] kabatdatabasefile >patchedfile
              -d Specify file of deletions (entry numbers, one to a line)
              -p Specify patch file
```
Either (or both) `-d` or `-p` must be specified. For example, to apply the
additions and corrections in file `kabat.patch` to the database
`kabat.jan95.dat` and delete the entries in `kabat.del`, you would use the
command line:
```
   patchkabat -d kabat.del -p kabat.patch kabat.jan95.dat >kabat.jan95.p1.dat
```
Having prepared and merged the data, it must be split into separate
files for specifes and chain using `splitkabat`.

N.B. In common with KabatMan V2.4 itself, splitkabat V1.2 has been
modified to handle the 1995 Kabat database release.

The program is run by typing:
```
        splitkabat kabat1
```
This is repeated for kabat2, kabat3, etc. Data from each of the Kabat
files will be appended to the appropriate class file.

Also requried is a `file of files` which lists the Kabat dump files to
be read. The file of files has any number of lines of the format:
```
{ heavy | - } [light [light [light] ] ]
```
where heavy and light are names of matching heavy-chain and
light-chain Kabat dump files. (Square brackets indicate optional
items; you don't type them!)  If a light-chain dump file has no
matching heavy chain dump file, the heavy chain file must be specified
as `-`.

NOTE THIS FILE MUST BE CALLED `kabat.fof`

The following example file may be used for the Kabat dump files
available at the time of writing:
```
   cat.ig.hc
   chicken.ig.hc chicken.ig.lambda
   dog.ig.hc
   frog.ig.hc
   gopher.ig.hc
   - horse.ig.lambda
   human.ig.hc human.ig.kappa human.ig.lambda
   mouse.ig.hc mouse.ig.kappa mouse.ig.lambda
   rabbit.ig.hc rabbit.ig.kappa rabbit.ig.lambda
   rat.ig.hc rat.ig.kappa rat.ig.lambda
   shark.ig.hc
   - sheep.ig.lambda
   various.ig.hc various.ig.kappa various.ig.lambda
```

When the program is run by typing the command `kabatman`, the data 
contained in these files are processed, with correlations being made 
between heavy and light chains of the same name and the resulting data
are written to the file `kabat.dat` in the current directory.

If the file `kabat.dat` exists in the current directory, or in the
directory described by the environment variable `KABATDIR`, the file of
files and the Kabat dump files will be ignored. The `-f` flag on the
command line forces reading of the file of files and the generation of
`kabat.dat` from the Kabat dump files.

For the Chothia canonical information to be available, there must also
be a canonicals definition file in the current directory or in the
directory described by the environment variable `KABATDIR`. This file
must be called chothia.dat and has the format:
```
LOOP <loopname> <class> <length>
<residue> <types>
...
```
Where:
```
      <loopname> is L1...H3
      <class>    is the class label (1, 2, 1', etc)
      <length>   is the loop length (AbM definition)
      <residue>  is a residue specification (e.g. L25)
      <types>    is a list of allowed amino acid types (e.g. SYF)
```

Residue numbering is Kabat standard unless the keyword `CHOTHIANUMBERING`
appears at the beginning of the file (before any `LOOP` keywords).

By convention, class names not defined by Chothia (i.e. additional
classes defined on the basis of new crystal structures) are followed
by a * (e.g. Class 9* for CDR-L3).

Comments may be introduced with an exclamation mark (!) or a hash (#).
These characters must appear in the first column of the file.

The file may also contain one SOURCE record per loop record. This 
should follow the LOOP record and contains free-format text describing
the definition of the canonical. This information is not used by the
KabatMan program.


The Database Model
------------------

The database uses a simple functional data model (FDM). As well as
allowing data stored directly in the database tables to be accessed,
an FDM database allows derived data to be accessed. This reduces the
necessary storage size allowing information to be calculated `on the
fly'. For example, the CDRs of the antibodies may be extracted using
the KabatMan database. They are not stored explicitly as separate
columns in the database as would be necessary using a normal flat file
or relational database. Instead a function is stored which allows the
CDRs to be derived from the sequence data stored for the whole
antibody Fv. Similarly, the lengths of the CDRs are not stored
directly in the database tables; they are calculated using a function
in the database.

Clauses in the WHERE statement which filter the data are combined in
reverse polish notation. Selections resulting from each clause are
placed on a stack and the logical operators act on this stack to
combine the results. No attempt is made to optimise the searches (i.e.
successive components which are ANDed together only need to check the
items found in the previous hit list, but this is not done). Although
this incurs a time penalty, the current database takes no more than a
few seconds to search and the lack of optimisation simplifies the code.


The Query Language
------------------

Data in the database is accessed using a SQL-like query language. The
main deviation between the language used here and SQL is that clauses
within the WHERE statement are combined in reverse polish notation.
This is described in detail below. Also, since there is only one table
in the database, there is no FROM statement.

The query language uses three statements: SET, SELECT and WHERE. The
SELECT and WHERE commands take you into a mode where clauses are
specified with no introductionary SELECT or WHERE keyword; the prompt
changes from `KABATMAN>` to indicate the current mode. Conversely,
the SET command needs to be specified on each line where variables are
to be set and may only be given at the `KABATMAN>` prompt.

The SELECT and WHERE statements remain valid until a new statement is
given. This allows multiple SELECT statements (each followed by a
search execution) to obtain different data from a single WHERE
statement thus saving on typing. Similarly, a number of WHERE
statements may be tried to obtain the same data under different search
conditions. 

The QUIT or EXIT command is used to leave the program.

### 5.1 The SET Statement 

This allows variables to be set. These are used to control the
operation of the program and options set using this keyword remain
valid until the program is terminated or another SET command is given.
The SET keyword must be given at the start of each line containing
variables to be set and may only be given at the `KABATMAN>'
prompt. More than one variable/value pair may be specified on the
line. The syntax is thus:
```
        SET variable value [variable value ...]
```
(Square brackets indicate optional items; you don't type them!)

The following variables may be set:
```
   Variable     Values               Notes
   --------     ------               -----
   LOOPS        KABAT, ABM, CHOTHIA, Sets the definitions to use for 
                CONTACT              the CDRs
   INFO         integer              Sets the information level for 
                                     messages
   INSERTS      ON, OFF              Specifies whether CDRs should be 
                                     displayed with or without alignment
                                     dashes
   VARIABILITY  real (0...100)       Specifies a sequence identity
                                     threshold (see below).
   HTML         ON, OFF              Adds some HTML codes to the
                                     output to make the results
                                     look nicer over the web
   URL          string               Specify the URL
   CANONICAL    string               Specify the canonical class
                                     definition file
   DELIMiter    chararacter          Specify the field delimiter
                                     character
```

The VARIABILITY variable allows one to specify that only sequences
with less than a specified sequence homology be displayed. After
all WHERE sub-clauses have been executed, KabatMan will do a pairwise
comparison of all the remaining antibodies and reject one from each
pair which has a percentage sequence identity greater than this value.
For each sequence mismatch, a penalty of 100/MeanLength % is scored;
missing chains do not score a penalty; deletions in one sequence with
respect to the other are given a double penalty (i.e. 200/MeanLength).
By default, VARIABILITY is set to 0.0 which is used to indicate `off'
(this has the same effect as setting it to 100, but the program will
run faster since the variability code will not be called).

N.B. This option is SLOW if there are more than a couple of hundred
hits to be processed. For example, with 1124 hits (using the query
SELECT light WHERE complete = t;), takes 4 minutes on a 100MHz R4000
SG Indy. 


The option to specify the URL allows you to specfiy a different URL
string to be used when selecting the URLL or URLH fields. You may
also give the command SET URL DEFAULT to restore the default URL.
If you use this command, you must specify the format in the correct
fashion. It must start with href= and end with the HTML GET method
query. The actual accession code and anchor tags will be appended
to the string.

The default string is:
```
href=http://immuno.bme.nwu.edu/scripts/noninter.tcl?qt=
```

If this doesn't mean anything to you, then DO NOT use this option!


The SET CANONICAL command allows the canonical definitions (key residue)
file to be specified. The default file (which is read on startup and
may be re-read by giving the command SET CANONICAL DEFAULT) is called
`chothia.dat`. Any string which you specify with the SET CANONICAL
command (other than DEFAULT) will be appended with a `.` onto the name
`chothia.dat` and this file will be read. For example, specifying
SET CANONICAL STRICT will cause the program to read definitions from
a file called `chothia.dat.strict'


The SET DEMILITER command is fairly self-explanatory. The default is that
output fields are separated by a comma and a space. However, some of the 
names contain commas so one can change to any other character (e.g. |)
as a delimiter. A space is automatically appended to the delimiter.


### 5.2 The SELECT Statement

Specifies what data should be extracted from the database. Each field
may be separated by a space and/or a comma. Parameters are placed in
parentheses. In the case of the PIR command, the parameter is
optional; if not specified, output is to the standard output.

The select statement may extend over more than one line.

The following fields may be displayed:
```
   Property         Notes
   --------         -------------------------
   NAME             Name of antibody
   ANTIGEN          Antigen (or blank string)
   L1               Sequence of CDR L1
   L2               Sequence of CDR L2
   L3               Sequence of CDR L3
   H1               Sequence of CDR H1
   H2               Sequence of CDR H2
   H3               Sequence of CDR H3
   CLASS            Light chain class (derived from file or filename)
                    Note that there may also be heavy chain information
                    here (e.g. IGG-KAPPA)
   SOURCE           Animal source (derived from filename)
   REFerence        Publication reference
   LENgth(loop)     Length of specified CDR
   RESidue(pos)     Amino acid at specified residue
   PIRfile[(file)]  Create a PIR file (of given name if specified)
   LIGHT            The light chain sequence
   HEAVY            The heavy chain sequence
   CANonical(loop)  The Chothia canonical class for a loop
   IDLight          The accession code for the light chain
   IDHeavy          The accession code for the heavy chain
   URLLight         An HTML hypertext link to the light chain raw data
   URLHeavy         An HTML hypertext link to the heavy chain raw data
   SUBGroup(chain)  The subgroup (currently only for human sequences)
   DATe             The year of the earliest reference
   LFR1             Sequence of light framework 1
   LFR2             Sequence of light framework 2
   LFR3             Sequence of light framework 3
   LFR4             Sequence of light framework 4
   HFR1             Sequence of heavy framework 1
   HFR2             Sequence of heavy framework 2
   HFR3             Sequence of heavy framework 3
   HFR4             Sequence of heavy framework 4
   SEQuence[(file)] Create a numbered sequence file (of given name if 
                    specified)
```
(The required parts of field names are in capital letters)


### 5.3 The WHERE Statement

Specifies filters to be applied to the data in order to select a
subset which is to be displayed.

Each clause is specified in conventional order notation using the form:
```
        property  test  value
```
The following properties are supported:
```
   Property        Type     Notes
   --------        ----     -----
   NAME            string   The antibody name
   ANTIGEN         string   The antigen (if specified)
   L1              string   Sequence of CDR L1
   L2              string   Sequence of CDR L2
   L3              string   Sequence of CDR L3
   H1              string   Sequence of CDR H1
   H2              string   Sequence of CDR H2
   H3              string   Sequence of CDR H3
   CLASS           string   Light chain class (derived from file or name)
                            Note that there may also be heavy chain
                            information here (e.g. IGG-KAPPA)
   SOURCE          string   Animal source (derived from filename)
   REFerence       string   Publication reference
   LENgth(loop)    integer  Length of specified CDR
   RESidue(pos)    char     Amino acid at specified position
   COMPLETE        bool     True of both chains are present
   VARiability     integer  Not implemented
   LIGHT           string   The light chain sequence
   HEAVY           string   The heavy chain sequence
   CANonical(loop) string   Chothia class name or ? if unmatched
   IDLight         nnnnnn   The accession code for the light chain
   IDHeavy         nnnnnn   The accession code for the heavy chain
   SUBGroup(chain) string   The subgroup (currently only for human sequences)
   DATe            The year of the earliest reference
   LFR1            Sequence of light framework 1
   LFR2            Sequence of light framework 2
   LFR3            Sequence of light framework 3
   LFR4            Sequence of light framework 4
   HFR1            Sequence of heavy framework 1
   HFR2            Sequence of heavy framework 2
   HFR3            Sequence of heavy framework 3
   HFR4            Sequence of heavy framework 4
```
(The required parts of field names are in capital letters)
nnnnnn represents a 6-digit accession code.


The following tests may be made:
```
   Test                   Property type  Notes
   ----                   -------------  -----
   =, ==, eq              all            Exact match (NOTE 1)
   !=, <>, ne             all            Mismatch
   <, lt                  integer        Less than
   <=, le                 integer        Less than or equal
   >, gt                  integer        Greater than
   >=, ge                 integer        Greater than or equal
   cont,like,sim,inc,sub  string         Sub-string match (NOTES 2,3)
```

- NOTE 1: The exact match is not recommended for use with strings; it is
much safer to use the substring match.

- NOTE 2: These test words are all equivalent and are the minimum string
required. Thus `contains' may be used in place of `cont'; `includes' in
place of `inc', etc.

- NOTE 3: Sub-string matches on the CDRs, LIGHT and HEAVY properties ignore
dashes in the sequence allowing one to search for portions of sequence
while ignoring gaps left for insertions.


String values do not need to be placed in inverted commas unless they
include white space or inverted commas. Thus one can use the string
"HYHEL-5'CL" to include a single inverted comma or 'a "name"' to include
double inverted commas. Blank strings may be specified using '' or "".

You are recommended NOT to use exact matches for strings; use the
sub-string match operator.

Sub-string matches may contain ? wildcard characters to match single
characters. This is particularly useful for sequence matches, so you
can, for example, search for 10 residue L1s containing S?S??????Y.

Boolean values are entered as `true' or `false' (only the first letter, 
t or f, is actually considered). 

All comparisons ignore case. 


WHERE clauses are combined using Reverse Polish Notation. Each clause
generates a set which matches the specified test and places the set on
a stack. The logical operators AND, OR and NOT cause operations to be
performed on this stack of sets.

```
   Stack        Operator        Stack    Notes
   -----        --------        -----    -----
   1 2          AND             1        Performs a logical AND
   1 2          OR              1        Performs a logical OR
   1            NOT             1        Performs a logical NOT
```

The stack has a maximum depth of 10 items. This should be more than
enough for all likely queries providing that the logical operators are
not all placed at the end of the statement. (Should the stack depth
prove insufficient, the constant STACKDEPTH in kabatman.h may be
increased at compile time.)


### Running The Search

The search is run by typing a . or ; on a line by itself. These cause
output to be directed to the screen. Alternatively, output my be sent
to a file by typing a greater than sign followed by a filename (a
space may optionally be placed between the > and the filename). Output
from the search will then be sent to this file rather than to standard
output. Note that the PIR option in the SELECT statement also allows a
file to be specified for output.


### Examples


1. Find all complete antibodies where the antigen is known with loop
lengths:
```
SELECT name,antigen,length(l1),length(l2),length(l3),
       length(h1),length(h2),length(h3)
WHERE  antigen != '' complete = true and
;
```

2. Get the sequences of all complete mouse antibodies which bind to
lysozyme, placing the results in a file AbLysozyme.pir:
```
SELECT pir(AbLysozyme.pir)
WHERE  source   includes mouse
       antigen  includes lysozyme   AND
       complete =        true       AND
;
```

3. Find all antibodies with 11 residue CDR-L1s and a proline at the 
sixth position:
```
SELECT name, l1
WHERE  len(l1) = 11 res(L29) = P and
;
```

4. Find all complete antibodies with the sequence Ser-Ala-Ser-Ser-Ser 
in the light chain. (Note that there must be no spaces in the sequence:
```
SELECT name, light
WHERE  complete = t
       light includes SASSS AND
;
```

5. Find all Kappa-III light chains
```
SELECT name
WHERE  light <> ''
       class includes 'kappa' and
       subgroup(L) = 'III'
;
```


Loop Definitions
----------------

The `SET LOOPS` command allows one of the following four sets of loop
definitions to be selected:
```
   Loop      Kabat         AbM            Chothia       Contact
   ----      -----         ---            -------       -------
   L1        L24 -- L34    L24 -- L34     L24 -- L34    L30 -- L36
   L2        L50 -- L56    L50 -- L56     L50 -- L56    L46 -- L55
   L3        L89 -- L97    L89 -- L97     L89 -- L97    L89 -- L96
   H1        H31 -- H35B   H26 -- H35B    H26 -- H32    H30 -- H35B
   H2        H50 -- H65    H50 -- H58     H52 -- H56    H47 -- H58
   H3        H95 -- H102   H95 -- H102    H95 -- H102   H93 -- H101
```

Standard Kabat Numbering
------------------------

```
Light chain:    0     1     2     3     4     5     6     7     8     9   
               10    11    12    13    14    15    16    17    18    19 
               20    21    22    23    24    25    26    27  
               27A   27B   27C   27D   27E   27F               28    29 
               30    31    32    33    34    35    36    37    38    39 
               40    41    42    43    44    45    46    47    48    49 
               50    51    52    53    54    55    56    57    58    59 
               60    61    62    63    64    65    66    67    68    69 
               70    71    72    73    74    75    76    77    78    79
               80    81    82    83    84    85    86    87    88    89
               90    91    92    93    94    95   
               95A   95B   95C   95D   95E   95F   96    97    98    99 
              100   101   102   103   104   105   106
              106A                                      107   108   109

Heavy chain:    0     1     2     3     4     5     6     7     8     9 
               10    11    12    13    14    15    16    17    18    19 
               20    21    22    23    24    25    26    27    28    29 
               30    31    32    33    34    35 
               35A  35B                            36    37    38    39
               40    41    42    43    44    45    46    47    48    49
               50    51    52 
               52A   52B   52C   53    54    55    56    57    58    59 
               60    61    62    63    64    65    66    67    68    69 
               70    71    72    73    74    75    76    77    78    79 
               80    81    82 
               82A   82B   82C   83    84    85    86    87    88    89 
               90    91    92    93    94    95    96    97    98    99 
              100  
              100A  100B  100C  100D  100E  100F  100G  100H  100I  100J
              100K  101   102   103   104   105   106   107   108   109
              110   111   112   113
```

Deviations from the Standard Numbering
--------------------------------------

The Kabat standard numbering does not specify any naming convention for
insertions which occur in addition to the supplied insertion
lettering. This is particularly evident with very long CDR-H2
loops. The Kabat database may, for example, specify that additional
residues occur between H100G and H100H, but does not specify how these
should be named. KabatMan takes a simpler approach and just adds
letters sequentially throughout such an insertion region. This may
lead to differences in the naming of residues in such very long
insertion regions.


Acknowledgement
---------------

I should like to thank the maintainers of the Kabat sequence database
for their assistance in making the data available and in listening to
our suggestions about the new format. In particular, my thanks go to
George Johnson.


Revision History
----------------

V0.1  12.04.94 Development version
V1.0  27.04.94 Original release version
V1.1  11.05.94 Some internal changes for distribution purposes and
               fixed bug in command line parsing.
V1.2  11.05.94 Added new function for finding Chothia canonical classes.
               Blanks in the sequence are now read as ? rather than -
V2.0  30.06.94 Reads the new Kabat file format by default.
V2.1  11.07.94 Additional check made when matching light & heavy 
               chains; at least one author's name must match.
V2.2  21.07.94 New Kabat file reader takes strain information from
               ANNOTA STRN records. Now only takes the source from the
               filename if the filename is not `VARIOUS....'. Taking this
               information from the filename is not really required with
               the new format files.
               Also reads and compares all references rather than just
               the first.
V2.3  25.01.95 Added the SET VARIABILITY code.
V2.4  08.02.95 Modified for the January 1995 Kabat database format
V2.5  07.03.95 The Kabat reading code now places - instead of ? at the
               end of a sequence. This stops all the sequences which are
               truncated in CDR3 from giving maximum length CDR3s
V2.6  13.03.95 kabatman -version gives version information and exits.
               Fortran style comparisons (ne, eq, lt, etc.) now available.
               Fixed small bug in handling PIR option when a file is
               not specified. 
----  31.03.95 Documentation changes only.
V2.7  13.05.95 The Chothia data file may now contain SOURCE records which
               are ignored by KabatMan and comments introduced by a # as
               well as a !
---   16.06.95 Documentation changes only.
V2.8  22.06.95 Fixed bug when full chain sequences given in a WHERE clause
V2.9  23.06.95 Clean compiling under gcc -Wall
V2.10 27.06.95 Fixed bug in fuzzy string matching
---   27.11.95 Documentation changes only.
V2.11 15.12.95 Introduced single character wildcard matches
---   26.02.96 Documentation changes only.
V2.12 02.04.96 Added handling for IDs and URLs
V2.13 11.04.96 Added header lines to data file and datafile date
               Prints accession codes of skipped entries when reading
               Kabat data files. Allow --version as well as -version
V2.14 18.04.96 Added SET URL
V2.15 22.04.96 Fixed potential access violation in reading canonical
               data file.
               Fixed initialisation of URL string
V2.16 07.05.96 Chothia data file can now use Chothia or Kabat numbering
V2.17 29.05.96 Added SET CANONICAL
V2.18 10.09.97 Added SUBGROUP(chain) (Chain is either L or H)
V2.19 14.10.98 Added SET DELIMITER
V2.20 28.05.99 Fixed core-dump bug in EatInitials()
V2.21 13.07.00 Fixed long-standing bug preventing string equality 
               comparison when a string contained a ' or "
V2.22 31.07.00 Added Contact loop definitions
V2.23 03.04.02 Added DATE
V2.24 28.02.05 Added LFR1...HFR4
V2.25 24.08.06 Added SEQUENCE
V2.26 07.10.19 A maintenance release - all moved to GitHub and an install
               script added
