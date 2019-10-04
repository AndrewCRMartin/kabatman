#!/usr/bin/perl
# Simple script to take all the Kabat files from the current directory
# and extract the useful information from them (for antibodies) and
# write the data to a single file
#
# 10.10.13 Original   By: ACRM
#
use strict;

opendir(DIR, ".");
my @files = readdir(DIR);
closedir(DIR);

foreach my $file (@files)
{
    if($file =~ /^0/)
    {
        ProcessFile($file);
    }
}

sub ProcessFile
{
    my($file) = @_;
    my ($kadbid, $defini, $specie, $aaname);

    my $ok = 0;
    my $printed = 0;

    if(open(FILE, $file))
    {
        while(<FILE>)
        {
            if(/^KADBID/)
            {
                $kadbid = $_;
            }
            elsif(/^DEFINI/)
            {
                $defini = $_;
                if(($defini =~ /\s*IG\s+HEAVY/) ||
                   ($defini =~ /\s*IG\s+KAPPA/) ||
                   ($defini =~ /\s*IG\s+LAMBDA/))
                {
                    $ok = 1;
                }
            }
            elsif(/^SPECIE/)
            {
                $specie = $_;
            }
            elsif(/^AANAME/)
            {
                $aaname = $_;
            }
            elsif(/^SEQTPA/)
            {
                if($ok)
                {
                    if(!$printed)
                    {
                        print $kadbid;
                        print $defini;
                        print $specie;
                        print $aaname;

                        $printed = 1;
                    }
                    my $resnum = substr($_, 17, 6);
                    my $res    = substr($_, 28);

                    my $resCopy = $res;
                    $resCopy =~ s/\s//g;
                    if(length($resCopy))
                    {
                        print "$resnum $res";
                    }
                }
            }
        }
        close(FILE);
    }
}
