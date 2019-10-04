#!/usr/bin/perl
#*************************************************************************
#
#   Program:    updatestats
#   File:       updatestats.perl
#   
#   Version:    V1.1
#   Date:       22.04.96
#   Function:   Update the HTML page containing stats on the KabatMan data
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
#   updatestats <statfile> <lighthits> <heavyhits> <completehits>
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0  26.02.96 Original   By; ACRM
#   V1.1  22.04.96 Modified to read first number out of hits files
#                  now that date stamps also appear there.
#                  
#
#*************************************************************************
if($#ARGV != 3)
{
    print "Usage: updatestats <statfile> <lighthits> <heavyhits> <completehits>\n";
    exit;
}

$statfile     = $ARGV[0];
$lighthits    = $ARGV[1];
$heavyhits    = $ARGV[2];
$completehits = $ARGV[3];

$newstatfile  = "$statfile.new";

# Generate the Main part of the html stats page
print <<__EOF;
[% INCLUDE "../header.tt" %]
[% INCLUDE "../main_menu.tt"
   abs = " id='mcurrent'"
%]
[% INCLUDE "abs_menu.tt"
   sequences = " id='current'"
%]


<h1>KabatMan - Kabat Sequence Database Statistics</h1>

<p>Statistics for the number of sequences contained in the Kabat
antibody sequence database.</p>

<p>Only immunoglobulin sequences are stored here; all T-cell receptor,
MHC sequences, etc. are rejected. In addition, sequences with fewer
than 75 residues are rejected so that the database contains only
essentially complete light or heavy chain sequences.</p>

<p><i>Complete antibodies</i> are those which have both light and
heavy chains. These are defined by ensuring the species, sequence name
and author details match between the chains. The raw data contains no
other pointers to link light and heavy chains to one another.</p>

<hr />
<h2>Previous Releases</h2>

__EOF
# Add details of previous releases

open(KABSTAT,    $statfile)       || die("Unable to read $statfile");
open(NEWKABSTAT, ">$newstatfile") || die("Unable to write $newstatfile");
while(<KABSTAT>)
{
    print NEWKABSTAT;           # Copy to the new file

    ($date, $nlightprev, $nheavyprev, $ncompleteprev) = split;
    print "<h3>$date</h3>\n<p>\n";
    print "Light chains: $nlightprev<br />\n";
    print "Heavy chains: $nheavyprev<br />\n";
    print "Complete antibodies formed: $ncompleteprev\n";
    print "</p>\n";
    print "\n";
    print "\n";
}


# Add details of current release
print "<hr />\n";
print "<h2>Current Release</h2>\n";

print "<p>Note that Kabat has been static since July 2000. Updates here reflect changes to the Kabatman code</p>\n";

# Get the current date into a string
$date = `date`;
($j,$month,$day,$j,$j,$year) = split(' ', $date);
$date =  "$day$month$year";

# Get the numbers of hits
$nlight    = &GetCount($lighthits);
$nheavy    = &GetCount($heavyhits);
$ncomplete = &GetCount($completehits);

# Write these data into the new stats file
print NEWKABSTAT "$date $nlight $nheavy $ncomplete\n";

# Write the entry for the Web page
print "<h3>$date</h3>\n<p>\n";
printf "Light chains: $nlight (%.1f%% change)<br />\n",
    100.0*($nlight-$nlightprev)/$nlightprev;
printf "Heavy chains: $nheavy (%.1f%% change)<br />\n",
    100.0*($nheavy-$nheavyprev)/$nheavyprev;
printf "Complete antibodies formed: $ncomplete (%.1f%% change)\n",
    100.0*($ncomplete-$ncompleteprev)/$ncompleteprev;
print "</p>\n<hr />\n";
print "\n";

print <<__EOF;
[% INCLUDE "../footer.tt" %]
__EOF


# Finally we close the files and rename the new stats file to the old one
close(KABSTAT);
close(NEWKABSTAT);
system("mv $newstatfile $statfile");

exit;


#*************************************************************************
# int GetCount($filename)
# -----------------------
# Extract the number of hits from a KabatMan output file
# 26.02.96 Original   By: ACRM
# 22.04.96 Modified split to cope with the date stamps now in the 
#          KabatMan output
# 28.05.99 Remodified it was picking up the = sign!
#
sub GetCount
{
    local ($file) = @_;

    open(FILE, $file) || die("Unable to read file $file\n");

    while(<FILE>)
    {
        chop;

        if(/Number of hits/)
        {
            ($j, $j, $j, $j, $j, $nhits, $j) = split;
            close(FILE);
            return $nhits;
        }
    }

    close(FILE);
    return 0;
}
