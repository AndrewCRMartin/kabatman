#!/bin/perl
#*************************************************************************
#
#   Program:    PatachKabat
#   File:       patchkabat.perl
#   
#   Version:    V1.0
#   Date:       10.02.95
#   Function:   Patches the raw Kabat database file with a set of fixes
#               and additions
#   
#   Copyright:  (c) Dr. Andrew C. R. Martin 1995
#   Author:     Dr. Andrew C. R. Martin
#   Address:    Biomolecular Structure & Modelling Unit,
#               Department of Biochemistry & Molecular Biology,
#               University College,
#               Gower Street,
#               London.
#               WC1E 6BT.
#   Phone:      (Home) +44 (0)1372 275775
#               (Work) +44 (0)171 387 7050 X 3284
#   EMail:      INTERNET: martin@bsm.bioc.ucl.ac.uk
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
#   The 1995 format Kabat database comes as a set of separate files
#   stored in subdirectories. These must be merged into a single file
#   for processing by splitkabat and feeding into KabatMan.
#   This is done in the sh shell using:
# 
#      touch kabat.jan95.dat
#      for dir in 0*
#      do
#         cat $dir/* >>kabat.jan95.dat
#      done
#
#   A separate tar file temp.tar is also supplied with additions and
#   patches to the database. These must be merged into a single file
#   using:
#
#      cat temp/* >kabat.patch
# 
#   This program is used to patch these into the main file.
#
#   N.B. The program assumes that entries in the patch file are in 
#   numerically ascending order. The above cat command should assure
#   this is the case (but if you have a weird shell, you may have
#   to do some trickery to ensure the input to cat is in ascending
#   order)!
# 
#   In addition a file of deletions may be specified. This causes the
#   specified entries to be removed from the database.
# 
#*************************************************************************
#
#   Usage:
#   ======
#   patchkabat [-d delfile] [-p patchfile] kabatdatabasefile
#      -d Specify file of deletions (just entry numbers, one to a line)
#      -p Specify patch file
# 
#   Either (or both) -d or -p must be specified
#
#*************************************************************************
#
#   Revision History:
#   =================
#
#*************************************************************************
$dodel   = 0;
$dopatch = 0;

# Handle the command line
while(substr($ARGV[0],0,1) eq "-")
{
    if($ARGV[0] eq "-d")
    {
        $dodel   = 1;
        $delfile = $ARGV[1];
        shift;
    }
    elsif($ARGV[0] eq "-p")
    {
        $dopatch   = 1;
        $patchfile = $ARGV[1];
        shift;
    }
    else                        # Print usage message if -h given as parameter
    {
        print "\nPatchKabat V1.0 (c) 1995 Dr. Andrew C.R. Martin, UCL\n";
        print "\npatchkabat [-d delfile] [-p patchfile] kabatdatabasefile\n";
        print "   -d Specify file of deletions (entry numbers, one to a line)\n";
        print "   -p Specify patch file\n";
        print "Either (or both) -d or -p must be specified\n\n";

        exit;
    }

    shift;
}

# Check that either (or both) patch/delete files are specified
if($dodel == 0 && $dopatch == 0)
{
    print "\nYou must specify a patch file and/or a delete file!\n\n";
    exit;
}

# If a delete file has been specified, read it
if($dodel)
{
    &ReadDelFile($delfile);
}

# If doing patches, open the patch file
if($dopatch)
{
    open(PATCH,$patchfile) || die "Unable to open patch file: $patchfile\n";
}

# Process the main Kabat database file
$kadbid = 0;
$patchid = 0;
$patching = 0;
$maxpatchid = 0;

while(<>)
{
    if($dopatch)                # If doing patches
    {
        if($patchid == 0)       # No patch read, so get the ID
        {
            while($buffer = <PATCH>)
            {
                if($buffer =~ /KADBID/)
                {
                    ($junk, $patchid) = split(/[ \t\n]+/, $buffer);

                    if($patchid < $maxpatchid)
                    {
                        die "\nPatches must be in ascending order!\n\n";
                    }
                    if($patchid > $maxpatchid)
                    {
                        $maxpatchid = $patchid;
                    }
                    last;
                }
            }
        }
    }

    # Find out if this entry is to be patched or deleted
    if(/KADBID/)
    {
        $patching = 0;          # Clear flag for patching in progress
        ($junk, $kadbid) = split;

        if($dopatch)            # If doing a patch see if we've got an entry
        {
            if($kadbid == $patchid)
            {
                $patching = 1;
            }
        }

        if($dodel)
        {
            $skipping = &QueryDelete($kadbid);
        }
    }

    # If it's to be patched, read and write from patch file
    if($patching)
    {
        # Read through the record in the patch file and output it
        print $buffer;
        while($buffer = <PATCH>)
        {
            print $buffer;
            if($buffer =~ /RECEND/)
            {
                $patching = 0;
                $patchid  = 0;
                $skipping = 1;
                last;
            }
        }
    }

    # If we're skipping this record, read through it
    if($skipping)
    {
        while(<>)
        {
            if(/RECEND/)
            {
                last;
            }
        }
        $skipping = 0;
    }
    else                        # We're not skipping the record
    {
        print;
    }
}

# If there are any patches left, these simply need to be appended
if(!eof(PATCH))
{
    print $buffer;
    while(<PATCH>)
    {
        if(/KADBID/)
        {
            ($junk, $patchid) = split;

            if($patchid < $maxpatchid)
            {
                die "\nPatches must be in ascending order!\n\n";
            }
            if($patchid > $maxpatchid)
            {
                $maxpatchid = $patchid;
            }
        }
        print;
    }
}


#*************************************************************************
#   void ReadDelFile($delfile)
#   --------------------------
#   Read a file of entries to be deleted
#
#   10.02.95 Original    By: ACRM
#
sub ReadDelFile
{
    local ($delfile) = @_;
    $ndel = 0;
    open(DELFILE, $delfile) || die "Unable to open delete file: $delete\n";

    while(<DELFILE>)
    {
        chop;
        $DelIDs[$ndel] = $_;
        $ndel++;
    }
}

#*************************************************************************
#   BOOL QueryDelete($id)
#   ---------------------
#   Sees if a specified id is in the list to be deleted
#
#   10.02.95 Original    By: ACRM
#
sub QueryDelete
{
    local ($id) = @_;

    for($i=0; $i<$ndel; $i++)
    {
        if($id == $DelIDs[$i])
        {
            return(1);
        }
    }

    return(0);
}
    
