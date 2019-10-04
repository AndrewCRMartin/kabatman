#!/bin/sh
#*************************************************************************
#
#   Program:    installkabat
#   File:       installkabat.sh
#   
#   Version:    V1.2
#   Date:       22.04.96
#   Function:   Installs the updated Kabat data
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
#   V1.0  26.02.96 Original
#   V1.1  28.02.96 kabat.stat is now kept in $KABATDIR instead of $HTMLDIR
#                  Small changes to the way the new version of the data
#                  is copied to minimise time when data will be invalid
#   V1.2  22.04.96 Files in KABATDIR protected against writing
#
#*************************************************************************
NEWKABAT=/acrm/data/kabat/newkabat
KABATDIR=/acrm/data/kabat/kabatman
HTMLDIR=/acrm/www/html/abs

# Install the new kabat data; keep the old version 
cd $NEWKABAT
cp $KABATDIR/kabat.dat kabat.dat.old
cp kabat.dat $KABATDIR/kabat.dat.new
cp kabat.fof $KABATDIR
\mv -f $KABATDIR/kabat.dat.new $KABATDIR/kabat.dat

# Install the HTML
cp kabat.stat $KABATDIR
cp kabat_stats.tt $HTMLDIR
(cd $HTMLDIR; make)

# Make sure everything in the Kabat directory is protected from others
# writing to it (not sure where this problem came from...)
chmod og-w $KABATDIR/*

# Gzip the old version
gzip -f kabat.dat.old

# Making a copy for distribution to others
# /home/andrew/mirror/MakeFtpCopy.sh

