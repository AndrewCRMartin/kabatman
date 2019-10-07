#!/usr/bin/perl -s
#*************************************************************************
#
#   Program:    KabatTest
#   File:       kabattest.pl
#   
#   Version:    V1.5
#   Date:       07.10.19
#   Function:   Test a sequence against the KabatMan database
#   
#   Copyright:  (c) UCL, Dr. Andrew C. R. Martin 1995-2019
#   Author:     Prof. Andrew C. R. Martin
#   Address:    Biomolecular Structure & Modelling Unit,
#               Department of Biochemistry & Molecular Biology,
#               University College,
#               Gower Street,
#               London.
#               WC1E 6BT.
#   EMail:      andrew@bioinf.org.uk
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
#   Tests a sequence read from stdin or from a file specified on the
#   command line and stored in the form of resnum resid, one to a
#   line (i.e. the format output by the kabatseq program) and tests
#   the sequence against the KabatMan database. Any unusual residues
#   are identified
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
#   V1.0  13.03.95 Original    By: ACRM
#   V1.1  14.03.95 All paths specified directly and files place in /tmp
#   V1.2  07.03.96 Tries to read NLight and NHeavy from the kabat.stat
#                  file
#   V1.3  29.09.99 Fixed for KabatMan which outputs a # in front of the
#                  number of hits
#   V1.4  21.10.99 Modified to use ACRMPerlVars module
#   V1.5  07.10.19 No longer uses ACRMPerlVars, but expects everything
#                  installed together with kabatman
# 
#*************************************************************************
use Cwd qw(abs_path);
use FindBin;
use lib $FindBin::Bin;

$kabatman  = abs_path("$FindBin::Bin") . "/kabatman";
$KabatStat = $ENV{'KABATDIR'} . "/kabat.stat"; # Stats on the Kabat d/b
$pid       = $$;
$Warning   = 0;
$NLight    = 3138;                # Default number of light chains in d/b
$NHeavy    = 3989;                # Default number of heavy chains in d/b
$Percent   = 1;                   # Percentage cutoff to show a results

# Print usage message if -h given as parameter
if($ARGV[0] eq "-h")
{
    print "\nKabatTest V1.1 (c) 1995 Andrew C.R. Martin, UCL\n";
    print "Usage: kabattest [sequence]\n\n";
    print "Takes a list of Kabat residue IDs / aa name pairs (one to a line)\n";
    print "and tests the sequence against the KabatMan database\n\n";

    print "Input from stdin or a file specified on the command line.\n";
    print "Typically input would be piped from the KabatSeq program.\n\n";

    exit;
}

open(RSEARCH, ">/tmp/runsearch.$pid") || die "Can't open temp file 1\n";
open(RID,     ">/tmp/resid.$pid")     || die "Can't open temp file 2\n";

# Attempt to open the $KabatStat file and read the number of light and
# heavy chains from there. We just read through the file and accept
# the last entry as being the true one
if(open(KABATSTAT, $KabatStat))
{
    while(<KABATSTAT>)
    {
        ($date,$temp1,$temp2,$NComplete)   = split;
        if($temp1 ne "" && $temp2 ne "")
        {
            $NLight = $temp1;
            $NHeavy = $temp2;
        }
    }
    close(KABATSTAT);
}

# Build a control file to run the search with KabatMan
# Also puts the resids in a file
while(<>)
{
    ($resnum, $code) = split;

    if($resnum eq "ERROR:")
    {
        print;

        # Set flag to say there was an error in this chain 
        # N.B. That $code is the second word in the ERROR: message and
        # must be either `Heavy' or `Light' since we use the first
        # letter of this word to indicate the chain of the error.
        # Any other letter indicates a fatal error.
        $chain = substr($code,0,1);
        if($chain ne "L" && $chain ne "H")
        {
            close(RSEARCH);
            close(RID);
            unlink("/tmp/runsearch.$pid");
            unlink("/tmp/resid.$pid");
            
            die "Fatal error in alignment";
        }
        $Error{$chain} = 1;
    }
    elsif($resnum eq "WARNING:")
    {
        $Warning = 1;
        print;
    }
    else
    {
        print RSEARCH "where res($resnum) = $code\n";
        print RSEARCH ".\n";
        print RID     "$resnum $code\n";
    }
}
close(RSEARCH);
close(RID);

# Run the KabatMan program placing the output in a file
system("$kabatman </tmp/runsearch.$pid >/tmp/results.$pid");

# Now run through the resid file and output file from the search
# identifying residues which are uncommon
open(RID,     "/tmp/resid.$pid")   || die "Can't open temp file 3\n";
open(RESULTS, "/tmp/results.$pid") || die "Can't open temp file 4\n";

$LastChain = ' ';
$AllNormal = 1;

while(<RESULTS>)
{
    ($j0,$key) = split;
    if($key eq "Number")
    {
        ($j0,$key,$j1,$j2,$j3,$num) = split;
        $buffer = <RID>;
        ($resid, $name) = split(' ',$buffer);
        
        $chain = substr($resid,0,1);

        # This just speeds things up by only recalculating the cutoff
        # when the chain actually changes
        if($chain ne $LastChain)
        {
            if($chain == 'L')
            {
                $cutoff = $Percent * $NLight / 100;
            }
            else
            {
                $cutoff = $Percent * $NHeavy / 100;
            }
            $LastChain = $chain;
        }

        # If we're below the cutoff display this residue with its
        # percentage and number of examples
        if($num < $cutoff)
        {
            $AllNormal = 0;

            if($chain eq "L" && !$Error{$chain})
            {
                printf "%5s = $name in %.3f%% of light chains ($num examples)\n", 
                       $resid, $num * 100 / $NLight;
            }
            elsif(!$Error{$chain})
            {
                printf "%5s = $name in %.3f%% of heavy chains ($num examples)\n", 
                       $resid, $num * 100 / $NHeavy;
            }
        }
    }
}

# Print a message if everything normal
if($AllNormal)
{
    printf "No unusual features: all residues seen in > %.1f%% of\n",$Percent;
    print "sequences in the database\n";
}

# Clean up, removing temp files
unlink("/tmp/runsearch.$pid");
unlink("/tmp/results.$pid");
unlink("/tmp/resid.$pid");




