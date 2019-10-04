#!/usr/bin/perl
#*************************************************************************
#
#   Program:    
#   File:       
#   
#   Version:    
#   Date:       
#   Function:   
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

$pid         = $$;
$tmpfile     = "/tmp/fof.$pid";
$LastSpecies = "";

# Build a list of the *.ig.* files created by splitkabat
system("ls *.ig.* > $tmpfile");

# Read the temp file
open(TFILE, $tmpfile) || die "BuildFOF: Unable to read temp file";
$count = 0;
while(<TFILE>)
{
    chop;

    # Parse the input line
    ($InSpecies,$ig,$InChain) = split('\.');

    # If the species has changed
    if($InSpecies ne $LastSpecies)
    {
        # If we've got some (i.e. not first time around the loop
        if($count)
        {
            # If there's no H chain print a dash
            if($chain[0] ne "hc")
            {
                print "- ";
            }

            # Print the chains
            for($i=0; $i<$count; $i++)
            {
                print "$species[$i].ig.$chain[$i] ";
            }
            print "\n";
        }
        $count = 0;
        $LastSpecies = $InSpecies;
    }
    $species[$count] = $InSpecies;
    $chain[$count]   = $InChain;
    $count++;
}

# Do the last species group
if($count)
{
    # If there's no H chain print a dash
    if($chain[0] ne "hc")
    {
        print "- ";
    }
    
    # Print the chains
    for($i=0; $i<$count; $i++)
    {
        print "$species[$i].ig.$chain[$i] ";
    }
    print "\n";
}

unlink($tmpfile);


