#!/usr/bin/perl -w
#Originally by AZ
#Updated by BAH to allow for command-line arguments, because why wouldn't you in the first place :)?
#this program takes a file with scaffold positions on chromosomes in amel4.5, as follows
#scaffold	scaffold:start	scaffold:stop	chr	start	stop	orientation
#Group1.1	1	1382403	chr1	1	1382403	+
#Group1.10	1	1405242	chr1	5928187	7333428	+
#Group1.11	1	2557	chr1	7383429	7385985	-
#Group1.12	1	4728	chr1	7435986	7440713	+
#Group1.13	1	5741	chr1	7490714	7496454	0
# it also takes a file with saffold positions to be translated, as follows
#LocusID Scaffold Pos1 Post2*optional
#it prints an ouput file with locusID and ChromosomeID and position 
use warnings;
# (1) quit unless we have the correct number of command-line args
$num_args = $#ARGV + 1;
if ($num_args != 2) {
    print "\nUsage: script.pl scaff_positions \n";
    exit;
}
$infile1=$ARGV[0];
chomp($infile1);
open(GETDB,$infile1) || die "can't open: $!";
$infile=$ARGV[1];
chomp($infile);
open(GET,$infile) || die "can't open: $!";
open(OUT,">".$infile."_onChr.txt") || die "can't open: $!";
$_=<GETDB>; #get the header line
chomp;
while(<GETDB>)
	{chomp;
         my @line=split;
         my $scaffold=shift(@line);
         $db{$scaffold}=$_;
    	 }
#for my $s (keys %db)
#    {print OUT "$s is $db{$s}\n";
#    }
close(GETDB);
#now do the analysis
#read the file to be converted
while(<GET>)
	{chomp;
        if ($_ =~/GroupUn/)
            {print OUT "$_\n";
            }
        else
        {my @line=split;
        $s=$line[1]; #reads in the scaffold of the position to be converted
        $sinfo=$db{$s}; #finds the coordinates for the scaffold in the scaffold_chromsome database
        my @coord = split(" ",$sinfo); # turns the string in the scaffold_chromsome database into an array
        #print OUT @coord."\n";
        my $c=$coord[3]; #c is the chromsome where the saffold is found
        my $dir=$coord[6]; #direction
        my $firstscaffold = $coord[1]; #first base of scaffold
        my $firstchr= $coord[4]; #first base of scaffold on chr
        my $lastchr= $coord [5]; #last base of scaffold on chr
        print OUT "$line[0]\t$c";
        if ($dir eq "+" || $dir eq "0")
            {for (my $i=2;$i<@line;$i++) 
                {$pos=$line[$i]; #this is the position to translate
                 #if orientation if + or 0, then translate as follows
                 # (position to translate - first base of scaffold) + first base of chr in table
                my $newpos= ($pos-$firstscaffold)+$firstchr;
                print OUT "\t$newpos";
                }
            print OUT "\n";
            }
        else
            {@results=();
                for (my $i=2;$i<@line;$i++) 
                {$pos=$line[$i]; #this is the position to translate
                 #if orientation if -, then translate as follows
                 #last base of chr in table - (position to translate - first base of scaffold)
                my $newpos= $lastchr - ($pos-$firstscaffold);
                unshift(@results,$newpos); # add this to the end of an array, so that positions are in ascending order in the report
                }
            for (my $i=0;$i<@results;$i++)
                {
                print OUT "\t".$results[$i];
                }
            print OUT "\n";
            }
            }
                
    	}
