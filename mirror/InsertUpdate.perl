#!/usr/bin/perl
#*************************************************************************
#
#   Program:    InsertUpdate
#   File:       InsertUpdate.perl
#   
#   Version:    V1.0
#   Date:       18.03.96
#   Function:   Move files from the Kabat updates directory into the
#               appropriate data directory
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
#
#*************************************************************************
if($#ARGV != 0 || $ARGV[0] eq "-h")
{
    print "Usage: InsertUpdate <filename>\n";
    exit 0;
}

$filename = $ARGV[0];

$outdir = "../".substr($filename,0,3);
system("mv $filename $outdir");


