#!/bin/sh
#*************************************************************************
#
#   Program:    MakeFtpCopy
#   File:       MakeFtpCopy.sh
#   
#   Version:    V1.0
#   Date:       25.04.96
#   Function:   Installs a copy of the new Kabat data in the FTP area
#               and EMails interested parties
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
#   V1.0  25.04.96 Original
#
#*************************************************************************
#FTPUSERS="j-saldan@nimr.mrc.ac.uk"
FTPUSERS="andrew"
NEWKABAT=/data/kabat/newkabat
KABATDIR=/data/kabat/kabatman
HTMLDIR=/home/httpd/html/abs
FTPDIR=/home/httpd/html/abs/data

cp $KABATDIR/kabat.dat $FTPDIR
cd $FTPDIR
gzip -f kabat.dat


################ MAIL MESSAGE FOLLOWS ##################
mail -s "New Kabat data" $FTPUSERS <<EOF
The latest release of the KabatMan data file is now available by HTTP
from:

http://sapc34.rdg.ac.uk/abs/data/kabat.dat.gz

This mail was generated automatically, so if there are any problems,
please let me know.

Andrew

EOF
################ MAIL MESSAGE ENDS ##################
