#!/usr/user/martin/bin/bash
#*************************************************************************
#
#   Program:    buildkabat
#   File:       buildkabat.sh
#   
#   Version:    V1.0
#   Date:       18.03.96
#   Function:   Build the KabatMan data (Similar to the script used
#               by the mirroring code, but doesn't tidy up completely).
#   
#   Copyright:  (c) Dr. Andrew C. R. Martin 1996
#   Author:     Dr. Andrew C. R. Martin
#   Address:    Biomolecular Structure & Modelling Unit,
#               Department of Biochemistry & Molecular Biology,
#               University College,
#               Gower Street,
#               London.
#               WC1E 6BT.
#   Phone:      (Home) +44 (0)1372 275775
#               (Work) +44 (0)171 419 3890
#   EMail:      INTERNET: martin@biochem.ucl.ac.uk
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#   This script is run from the mirror software (it pretends to be the
#   mail program).
#
#   The script performs the following functions:
#   1. Uncompresses the archives and builds into a single file
#   2. Moves any "update" files to required directories
#   3. Calls the splitkabat program to split the files by species/chain
#   4. Calls BuildFOF.perl to create a KabatMan file of files
#   5. Runs KabatMan to convert to the KabatMan data format
#
#   Note that this script will need modifying once we get to 100000
#   data files. We're still on < 30000, so we've got a way to go yet!
#
#*************************************************************************
#
#   Usage:
#   ======
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0  26.02.96 Original   By: ACRM
#   V1.1  28.02.96 kabat.stat is now kept in $KABATDIR rather than in
#                  $HTMLDIR
#   V1.2  18.03.96 Added handling of the updates.tar.Z file
#
#*************************************************************************
# The directory in which the mirrored database is kept
MIRRORDIR=/data/kabat
# The directory where we will keep the new processed data.
NEWKABAT=/data/kabat/newkabat
# The general bin directory where KabatMan is found
BINDIR=/usr/local/bin
# The bin directory where code for completing the mirroring is kept
MIRRORBIN=/home/andrew/mirror
# The directory where kabat.stat is kept
KABATDIR=/data/kabat/kabatman


#*************************************************************************
# Move to the appropriate directory and make temp copies of the files
cd $MIRRORDIR

# Create temp directory if needed
echo -n "Creating copies of tar files..."
if [ ! -e temp ]; then
  mkdir temp
fi
cp *.tar.Z temp
cd temp
echo "done"

# Uncompress all the tar archives and remove them
echo -n "Uncompressing and removing tar files..."
for file in *.tar.Z
do
zcat $file | tar -xf -
done
\rm -f *.tar.Z
echo "done"

# Now move in any updates files
echo -n "Installing any update files..."
if [ -e updates ]; then
  cd updates
  for file in 0*
  do
     $MIRRORBIN/InsertUpdate.perl $file
  done
  cd ..
  rm -rf updates
fi
echo "done"

# Remove all .bak backup files from revised entries
\rm -f */*.bak

# Combine all the data files into a single file
# Note the use of $dir/$dir* rather than $dir/* to ensure that only the
# proper files are used (other junk in the directories is skipped)
echo -n "Creating single raw data file..."
\rm -f kabat.raw.dat
touch kabat.raw.dat
for dir in 0*
do
   cat $dir/$dir* >>kabat.raw.dat
done
echo "done"

# Remove all the individual files
echo -n "Removing separate entry data files..."
\rm -rf 0*
echo "done"

# Run the splitkabat program to split these up by species
echo -n "Splitting raw data by species and chain..."
$BINDIR/splitkabat kabat.raw.dat
echo "done"

# Build the Kabat file of files
$MIRRORBIN/BuildFOF.perl >kabat.fof

# Run kabatman to build the new data file
echo -n "Running KabatMan to build Kabat data file..."
$BINDIR/kabatman -f >kabatman.log <<EOF
quit
EOF
echo "done"

