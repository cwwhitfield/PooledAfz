#!/usr/bin/perl
use strict;
use warnings;

# Generates inputs for, runs, and graphs results from Ancestry_HMM for pooled sample(s).
# Pooled data forces us to do some weird things; some of these need to be examined
# more carefully. The functions get allele counts for ancestor refs (because we have them)
# but sequence depth for sample. The program doc asks for "pileup counts" and models
# sequence error, so using allele counts for ref may be a mistake.
#
# The functions are designed for our data where vcf files are scaffolds, but we may
# or may not want to examine whole chromosomes, so optionally converts to chromosome map.
# It should work if vcf is really chromosome, although args --scaff and --chr will be
# awkward.

my $help = q`
*******************************************************************************
Usage: perl pooledahmm.pl [function] --arg1 --arg2 ... -options

Extracts and formats data for use in Ancestry_HMM, runs it, and logs and
(optionally) displays results. Output and log files will be created in
directory "pooledahmm_data". The log contains input settings (eg, filters
used in vcftools extraction), sample(s) examined, scaffold or chromosome,
and summary Ancestry_HMM results (estimated admixture proportion, estimated
time in generations, and lnL). Use the --prefix input parameter to name a
project and keep it organized.

WARNING: Ancestry_HMM is not user friendly. If inputs are wrong (eg, missing
or misformated) it will generate output that looks somewhat like real data, or
hang while slowly consuming all server memory. (This programs intercepts some
possible problems.)

functions:
 [all]       Steps through functions below; skips the slow vcf-extrating
             steps (the first 3 below) if output already exists based on
             --prefix, --sample-name, and --scaff or --chr
 refcounts   Makes reference allele count files for --ref0-list and -ref1-list
             using --ref vcf
 sampledepth Makes sequence read depth files for --sample-name (or
             --sample-list) using --sample vcf
 fst         Makes FST file for --ref0-list and --ref1-list using --ref vcf
 mash        Combines ref count and sample depth data; with -f option picks
             best markers by FST in windows specified by --min-bp and --max-bp;
			 with -c option and --scaff-map converts scaffolds to chromosomes
 mash2ahmm   Calculates and appends recomb info from --recomb-map and formats
             input for Ancestry_HMM; can also convert to chromosome as above
 ahmm        Runs Ancestry_HMM and logs summary statistics
 ahmmpost    Extracts info from Ancestry_HMM .out and .stdout to generate
             multple data files (.maxpost, .lnL, .mean_ancestry, etc.)
 makeplots        Shows plot of maximum poster probability of ancestry by position
             (or use -v option when calling ahmm or "all")

input resource file names:
 --ref         vcf containing ref0 and ref1 data (needed by refcounts)
 --sample      vcf containing sample data (needed by sampledepth)
 --scaff-map   scaffold to chromosome conversion file (needed if -c option)
 --recomb-map  recombination map (needed by mash2ahmm)
 --recomb-by-chr recomb rates (cM/Mb) by chr; use in place of above
	   
other inputs:
 --ref0-list   file containg ref0 names (ancestor population) as rows
 --ref1-list   file containg ref1 names (ancestor population) as rows
 --sample-name name of single sample to be tested
 --sample-list file containing sample names as rows (in place of above)
 --scaff       scaffold restriction applied before conversion to chromosome;
               use commas only to specify >1 (e.g., 1.1,1.2,1.3)
 --chr         chromosome restriction applied after conversion to chromosome;
               use ,- to specify >1 (e.g., 1,3,5,9-16)
 
file naming control:
 --prefix    prefix name for all intermediate and output files (default "out")
 --suffix    suffix added only to output files; use to prevent overwriting if,
               for example, varying HMM parameters
 --work-dir  directory name for intermediate files (default "pooledahmm_files")
 --data-dir  directory name for Ancestry_HMM input/output and all script-
             processed logs and output (default "pooledahmm_data")
 
Note: The vcf extracting funcions (refcounts, sampledepth and fst) will check
for existing output and skip if found. You can save time by running these once
on whole genome (by not setting --scaff); all subsequent chromosome or scaffold
analyses using the same prefix will use existing extracted vcf data. In this
case, --ref, --ref0-list, --ref1-list and --sample can be omited. Later steps
(mash, mash2ahmm, ahmm) always overwrite existing files.			   
			   
parameters:
 --ref-filter  Quoted args to pass to vcftools when extracting from ref
 --samp-filter Quoted args to pass to vcftools when extracting from sample
 --min-bp    Defines window for marker pruning
 --max-bp    Togeter with min-bp defines window for marker pruning using best
             FST in each interval; will violate max if no marker in window
 --min-fst
 --ahmm-args Quoted arguments for Ancestry_HMM (eg, "-a 0.2"; don't use -i
             here, use --ploidy instead)
 --ploidy    Ploidy to use in Ancestry_HMM
 
options (you can group these after a single dash):
 -c          Convert scaffolds to chromosomes
 -u          Force mapping of unoriented scaffolds as if (+) orientation
 -v          View posterior probabilities in R plot (function "all" or ahmm) 
 
final output and analyses files:







See program file for example calls
   
*******************************************************************************

`;

#Example calls:



# INPUT HANDLER
die $help if ($#ARGV == -1 or $ARGV[0] eq 'help' or $ARGV[0] eq '-h');
my ($sub);

# Shared info used by one or more subroutines
my ($ref, $ref0_list, $ref1_list, $sample, $sample_list, $scaff_map, $recomb_map, $recomb_by_chr,
	$prefix, $suffix, $sample_name, $scaff_range, $chr_range, $samp_filter,
	$ref_filter, $ahmm_args, $min_bp, $max_bp, $min_fst, $work_dir, $data_dir, $ploidy,
	$chr_convert_option, $force_orientation_option, $view_option);
my %arguments = (
	'all', 			['sub', 	\&all			],
	'refcounts', 	['sub', 	\&refcounts		],
	'sampledepth', 	['sub', 	\&sampledepth	],
	'fst', 			['sub', 	\&fst			],
	'mash', 		['sub', 	\&mash			],
	'mash2ahmm', 	['sub', 	\&mash2ahmm		],
	'ahmm', 		['sub', 	\&ahmm			],
	'ahmmpost', 	['sub', 	\&ahmmpost		],
	'makeplots', 		['sub', 	\&makeplots			],
	
	'--ref', 		['file', 	\$ref			],
	'--ref0-list', 	['file', 	\$ref0_list		],
	'--ref1-list', 	['file', 	\$ref1_list		],
	'--sample', 	['file', 	\$sample		],
	'--sample-list',['file', 	\$sample_list	],
	'--scaff-map', 	['file', 	\$scaff_map		],
	'--recomb-map', ['file', 	\$recomb_map	],
	'--recomb-by-chr', ['file', \$recomb_by_chr	],
	
	
	
	'--prefix', 	['str', 	\$prefix,			'out'		],
	'--suffix', 	['str', 	\$suffix,			''			],
	'--sample-name',['str', 	\$sample_name,		''			],
	'--scaff',		['str', 	\$scaff_range,		''			],
	'--chr',		['str', 	\$chr_range,		''			],
	'--work-dir',	['str', 	\$work_dir,			'pooledahmm_files'],
	'--data-dir',	['str', 	\$data_dir,			'pooledahmm_data'],
	
	'--ref-filter',	['args', 	\$ref_filter,		''			],
	'--samp-filter',['args', 	\$samp_filter,		''			],
	'--ahmm-args', 	['args', 	\$ahmm_args,		''			],
	'--min-bp', 	['pos_int', \$min_bp,			0			],
	'--max-bp', 	['pos_int', \$max_bp,			0			],
	'--min-fst', 	['number', 	\$min_fst,			0			],
	'--ploidy', 	['pos_int', \$ploidy,			0			],
	
);

my %options = (			# entered singly as -v, or grouped -cdv
	# f => \$fst_option,
	c => \$chr_convert_option,
	u => \$force_orientation_option,
	v => \$view_option
);

# Init default values from table
foreach my $arginfo (values %arguments){					
	my $argtype = $arginfo->[0];
	if ($argtype eq 'file'){
		${$arginfo->[1]} = "";
	} elsif ($argtype ne 'sub'){
		${$arginfo->[1]} = $arginfo->[2];       # Behold my powers of dereferencing!
	}
}

# Process arguments
my @full_arg_list;
my $i = 0;
while($ARGV[$i]){
	my $arg = $ARGV[$i];
	if ($arg =~ /^\-(\w+)/) {	# process single dash single letter options, eg -vdcs
		foreach my $opt (split //, $1){
			my $optionref = $options{$opt};
			die "Did not recognize option \"$opt\" in \"$arg\"\n" if (not $optionref);
			$$optionref = 1;
			push @full_arg_list, "-$opt";
		}	
	} else {
		push @full_arg_list, $arg;
		my $arginfo = $arguments{$arg} or die "Unknown arg: \"$arg\"\n";
		my $argtype = $arginfo->[0];
		if ($argtype eq 'sub'){
			die "Expected function as first arg\n" if not ($i == 0);
			$sub = $arginfo->[1];
		} elsif ($argtype eq 'file'){
			$i++;
			my $filename = $ARGV[$i];
			push @full_arg_list, $filename;
			die "Expected filename after $arg\n" if (not $filename or $filename =~ /^-/);
			die "Could not find file \"$filename\" specified for $arg\n" if (not -f $filename);
			my $var_ref = $arginfo->[1];
			$$var_ref = $filename;
		} elsif ($argtype eq 'str'){
			$i++;
			my $value = $ARGV[$i];
			push @full_arg_list, $value;
			die "Expected string value after $arg\n" if (not $value or $value =~/^-/);
			my $var_ref = $arginfo->[1];
			$$var_ref = $value;
		} elsif ($argtype eq 'args'){
			$i++;
			my $value = $ARGV[$i];
			push @full_arg_list, "\"$value\"";  # the shell may have stripped quotes, so put them back
			# no die here; just take what user gives
			my $var_ref = $arginfo->[1];
			$$var_ref = $value;
		} elsif ($argtype eq 'pos_int'){
			$i++;
			my $value = $ARGV[$i];
			push @full_arg_list, $value;
			die "Expected positive integer after $arg\n" if (not $value or $value =~/\D/);
			my $var_ref = $arginfo->[1];
			$$var_ref = $value;
		} elsif ($argtype eq 'number'){
			$i++;
			my $value = $ARGV[$i];
			push @full_arg_list, $value;
			die "Expected number after $arg\n" if (not $value or $value !~/^[\d\.\-eE]+$/);
			my $var_ref = $arginfo->[1];
			$$var_ref = $value;
		}
	}
	$i++;
}
print "\nArguments as understood by pooledahmm.pl:\n\n";
print join (' ', @full_arg_list), "\n";

# General housekeeping
$sub = $sub || \&all;	# default subroutine
mkdir $work_dir if (not -e $work_dir);
mkdir $data_dir if (not -e $data_dir);

# File name constructors shared among functions
$work_dir = $work_dir . '/';
$data_dir = $data_dir . '/';
my $dotsuffix = $suffix && ".$suffix";	# either ".suffix" or ""
my $dotscaff = "";	# defined in chromosome_cycler or scaffold_cycler
my $dotloc = "";

# Capture vcftools stdout for this prefix to keep screen usefull
my $vcftools_stdout = $work_dir . "vcftools_stdout.$prefix.txt";

# Init log file
my $ahmm_log = $data_dir . "pooledahmm.log";
unless (-f $ahmm_log){
	print "Could not find log file $ahmm_log; opening new one\n";
	open LOG, "> $ahmm_log" or die "Could not open $ahmm_log for write";
	print LOG join (' ', "Time", "Prefix", "Suffix", "Sample","Chr","Scaff","Admix_Time", "Ancestry", "lnL", "Call"), "\n";
}
close LOG;

# Set up samples to test
die "Specify either --sample-name or --sample-list\n" if ($sample_name and $sample_list);
my @samples;
if ($sample_name){
	push @samples, $sample_name;
} else {
	open SAMPLELIST, $sample_list or die "Could not open $sample_list for read";
	while(<SAMPLELIST>){
		next if ($_ eq $/);
		s/\s//g;
		push @samples, $_;
	}
	close SAMPLELIST;
}
my $n_samples = $#samples + 1;
print "\n$n_samples samples to process: ";
print join (" ", @samples), "\n";

my $need_fsts = ($min_bp && $max_bp) || $min_fst;

my $restrict_chr = ""; 		# set in chromosome_cycler
my $restrict_scaff = ""; 	# set in scaffold_cycler
my $freindly_chr_scaff_name = "";

# These functions run the called function for each chromosome or scaffold
if ($chr_convert_option){
	chromosome_cycler();
} else {
	scaffold_cycler();
}

sub chromosome_cycler {
	print "\n************ running chromosome_cycler *************\n\n";

	# input testing
	die "Specify --chr when using -c option\n" unless ($chr_range);
	die "Don't specify --scaff when using -c option\n" if ($scaff_range);
	
	# set up array of chromosomes to test
	my @chromosomes;
	$chr_range =~ s/\s//g;
	my @prelim_chromosomes = split (/\,/, $chr_range);
	foreach (@prelim_chromosomes){
		if (/(\d+)\-(\d+)/){
			for (my $i = $1; $i <= $2; $i++){
				push @chromosomes, $i;
			}
		} else {
			die "Chromosomes must be digits only\n" if (/\D/);
			push @chromosomes, $_;
		}
	}
	my $n_chromosomes = $#chromosomes + 1;
	print "\n$n_chromosomes chromosomes to process: ";
	print join (" ", @chromosomes), "\n";
	
	# cycle through chromosomes doing requested function
	foreach my $chromosome (@chromosomes){
		$restrict_chr = $chromosome;
		$freindly_chr_scaff_name = "Chromosome $restrict_chr";
		print "Processing $freindly_chr_scaff_name\n";

		# File name constructors
		$dotscaff = "";
		$dotloc = ".$restrict_chr";

		
		
		
		# old
		# $dotscaff = $restrict_scaff && ".$restrict_scaff"; # either "" or ".scaffold" (can't do vcf extraction on chromosome)
		# $dotloc = $restrict_scaff || $restrict_chr;
		# $dotloc = $dotloc && ".$dotloc";		# now either ".chromosome", ".scaffold" or ""
		# $dotsuffix = $suffix && ".$suffix";	# either ".suffix" or ""



	
		# Log and call subroutine
		log_call();
		&$sub();	
	}
}

sub scaffold_cycler {
	print "\n************ running scaffold_cycler *************\n\n";
	# input testing
	die "Use -c option if specifying --chr\n" if ($chr_range);
	die "Specify --scaff or (with -c option) --chr\n" unless ($scaff_range);

	# set up array of scaffs to test
	die "" unless ($scaff_range);
	$scaff_range =~ s/\s//g;
	my @scaffolds = split (/\,/, $scaff_range);	# split on commas, assume no -'s
		
	my $n_scaffs = $#scaffolds + 1;
	print "\n$n_scaffs scaffolds to process: ";
	print join (" ", @scaffolds), "\n";
	foreach my $scaffold (@scaffolds){
		$restrict_scaff = $scaffold;
		$freindly_chr_scaff_name = "Scaffold $restrict_scaff";
		print "Processing $freindly_chr_scaff_name\n";

		# File name constructors
		$dotscaff = ".$restrict_scaff";
		$dotloc = ".$restrict_scaff";		

		# old
		# $dotscaff = $restrict_scaff && ".$restrict_scaff"; # either "" or ".scaffold" (can't do vcf extraction on chromosome)
		# $dotloc = $restrict_scaff || $restrict_chr;
		# $dotloc = $dotloc && ".$dotloc";		# now either ".chromosome", ".scaffold" or ""
		# $dotsuffix = $suffix && ".$suffix";	# either ".suffix" or ""

	
		# Log and call subroutine
		log_call();
		&$sub();
	}
}


sub log_call {
	# Log vcf extraction functions since a filter may have been used
	# skip for "all" or ahmm (will log at end of function) or others that don't use info that we need to record
	if ($sub eq \&refcounts or $sub eq \&sampledepth or $sub eq \&fst){
		open LOG, ">> $ahmm_log" or die "Could not open $ahmm_log for appending";
		my $time = localtime();
		my $full_arg_str = join (' ', @full_arg_list);
		print LOG join ("\t", $time, $prefix, $suffix, ($sample_name || $sample_list), $restrict_chr, $restrict_scaff, '---', '---', '---', $full_arg_str), "\n";
		close LOG;
	}
}




###############################    Subroutines    ######################################
sub all {
	# Do all (1st three will skip themselves if output already exists)
	print "\n************ running all functions for $freindly_chr_scaff_name *************\n\n";
	
	refcounts();
	sampledepth();
	if ($need_fsts){
		fst();
	}
	mash();
	mash2ahmm();
	ahmm();
	ahmmpost();
	if ($view_option){
		makeplots();
	}
}

sub refcounts {
	# Extracts allele counts for refs using vcftools
	print "\n************ running refcounts for $freindly_chr_scaff_name *************\n\n";
	
	# Skip?
	if (-f "$work_dir$prefix.ref01.frq.count.keyed"){
		print "Reference counts exist for $prefix (whole genome), skipping\n";
		return;
	} elsif ($restrict_scaff and -f "$work_dir$prefix.ref01$dotscaff.frq.count.keyed") {   
		print "Reference counts exist for $prefix$dotscaff, skipping\n";
		return;
	}
	
	die "Requires --ref\n" if (not $ref);
	die "Requires --ref0-list\n" if (not $ref0_list);
	die "Requires --ref1-list\n" if (not $ref1_list);
	
	if (not $restrict_scaff){
		print "Running refcounts for whole genome (--chr appplied after conversion)\n";
		print "This might take a while\n"
	}

	print "Appending vcftools stdout to $vcftools_stdout\n";
	
	# example call:
	# vcftools --vcf /media/data1/forty3/brock/align/AllBees_SNPs.raw.vcf --counts --max-alleles 2 --keep refA_names --chr 1.1 --out RefA_1.1	
	my @ref_call = ("vcftools --counts2 --max-alleles 2 $ref_filter");
	push @ref_call, ('--vcf', $ref);
	push @ref_call, ('--chr', $restrict_scaff) if ($restrict_scaff);
	my $ref0call = join (' ', @ref_call, '--keep', $ref0_list, '--out', "$work_dir$prefix.ref0$dotscaff", ">> $vcftools_stdout");
	my $ref1call = join (' ', @ref_call, '--keep', $ref1_list, '--out', "$work_dir$prefix.ref1$dotscaff", ">> $vcftools_stdout");
	print "\nSystem call:\n$ref0call\n";
	system ($ref0call);
	print "\nSystem call:\n$ref1call\n";
	system ($ref1call);
	
	# Outputs from vcftools will be, respectively:
	#	prefix.ref0[.scaff_or_chrom].frq.count
	#	prefix.ref1[.scaff_or_chrom].frq.count
	

	# Now combine and key these for mash
	my $ref0fc = "$work_dir$prefix.ref0$dotscaff.frq.count";
	my $ref1fc = "$work_dir$prefix.ref1$dotscaff.frq.count";	
	my $ref_counts_keyed = "$work_dir$prefix.ref01$dotscaff.frq.count.keyed";
	print "\nCreating combined keyed reference file\n";
	open REF0FC, $ref0fc or die "Could not open $ref0fc for read ";
	open REF1FC, $ref1fc or die "Could not open $ref1fc for read  ";
	open REFKEYED, '>', $ref_counts_keyed or die "Could not open $ref_counts_keyed for write ";
	
	my $printed = 0;
	my $missing = 0;
	$_ = <REF0FC>;				# skip headers
	my $r1_line = <REF1FC>;
	while (<REF0FC>){
		$r1_line = <REF1FC>;
		next if ($_ eq $/);
		# my ($scaff0, $pos0, $n_alleles0, $n_chr0, $count10, $count20) = split;  #faster to get what we want
		my ($scaff0, $pos0, $count10, $count20) = /(\S+)\t(\S+)\t\S+\t\S+\t(\S+)\t(\S+)/;
		my ($scaff1, $pos1, $count11, $count21) = $r1_line =~ /(\S+)\t(\S+)\t\S+\t\S+\t(\S+)\t(\S+)/;
		die "Mismatch in ref count files $scaff0 $pos0 $scaff1 $pos1 " if ($scaff0 ne $scaff1 or $pos0 ne $pos1);
		#Skip any missing data ("." appears in sample depth, not sure if ever in ref counts)
		if ($count10 !~ /\D/ and $count11 !~ /\D/){		# if count10 is good than count20 will be good
			my $key = "$scaff0-$pos0";
			print REFKEYED join ("\t", $key, $count10, $count20, $count11, $count21), "\n";
			$printed++;
		} else {
			$missing++;
		}
	}
	close REF0FC;
	close REF1FC;
	close REFKEYED;
	print "$printed lines printed; $missing missing data\n";
}

sub sampledepth {
	# Extracts sequence depth for pooled sample
	print "\n************ running sampledepth for $freindly_chr_scaff_name *************\n\n";

	foreach my $test_sample (@samples){
		# Skip?
		if (-f "$work_dir$prefix.$test_sample.RDAD.keyed"){
			print "Depth counts exist for $prefix.$test_sample (whole genome), skipping\n";
			next;
		} elsif ($restrict_scaff and -f "$work_dir$prefix.$test_sample$dotscaff.RDAD.keyed"){
			print "Depth counts exist for $prefix.$test_sample$dotscaff, skipping\n";
			next;
		}
		die "Requires --sample\n" if (not $sample);		
		
		print "\nRunning sampledepth for $test_sample\n";	
		if (not $restrict_scaff){
			print "Running sampledepth for whole genome (--chr appplied after conversion)\n";
			print "This might take a while\n"
		}
		print "Appending vcftools stdout to $vcftools_stdout\n";

		# example call:
		# vcftools --vcf prelim_AFZ.vcf --extract-FORMAT-info RD --max-alleles 2 --indv Sample1 --chr 1.1 --out Sample1_1_1
		my @sample_call = ("vcftools --max-alleles 2 $samp_filter");
		push @sample_call, ('--vcf', $sample);
		push @sample_call, ('--indv', $test_sample);
		push @sample_call, ('--chr', $restrict_scaff) if ($restrict_scaff);
		push @sample_call, ('--out', "$work_dir$prefix.$test_sample$dotscaff");
		my $sampleRDcall = join (' ', @sample_call, '--extract-FORMAT-info RD', ">> $vcftools_stdout");
		my $sampleADcall = join (' ', @sample_call, '--extract-FORMAT-info AD', ">> $vcftools_stdout");
		print "\nSystem call:\n$sampleRDcall\n";
		system ($sampleRDcall);
		print "\nSystem call:\n$sampleADcall\n";
		system ($sampleADcall);
		
		# Outputs from vcftools will be, respectively:
		#	prefix.sample_name[.scaff_or_chrom].RD.FORMAT
		#	prefix.sample_name[.scaff_or_chrom].AD.FORMAT
		
		# Now combine and key these for mash
		my $sampleRD = "$work_dir$prefix.$test_sample$dotscaff.RD.FORMAT";
		my $sampleAD = "$work_dir$prefix.$test_sample$dotscaff.AD.FORMAT";
		my $samp_depth_keyed = "$work_dir$prefix.$test_sample$dotscaff.RDAD.keyed";
		print "\nCreating combined keyed sample file for $test_sample\n";
		open SAMPRD, $sampleRD or die "Could not open $sampleRD for read ";
		open SAMPAD, $sampleAD or die "Could not open $sampleAD for read  ";
		open SAMPKEYED, '>', $samp_depth_keyed or die "Could not open $samp_depth_keyed for write ";
		
		my $printed = 0;
		my $missing = 0;
		$_ = <SAMPRD>;				# skip headers
		my $ad_line = <SAMPAD>;
		while (<SAMPRD>){
			$ad_line = <SAMPAD>;
			next if ($_ eq $/);
			chomp;
			my ($scaff_rd, $pos_rd, $count_rd) = split(/\t/);
			chomp($ad_line);
			my ($scaff_ad, $pos_ad, $count_ad) = split(/\t/, $ad_line);
			die "Mismatch in sample depth files $scaff_rd $pos_rd $scaff_ad $pos_ad " if ($scaff_rd ne $scaff_ad or $pos_rd ne $pos_ad);
			#Skip any missing data ("."s)
			if ($count_rd !~ /\D/ and $count_ad !~ /\D/){		# don't keep bad data ('.' used for missing data)
				my $key = "$scaff_rd-$pos_rd";
				print SAMPKEYED join ("\t", $key, $count_rd, $count_ad), "\n";
				$printed++;
			} else {
				$missing++;
			}
		}
		close SAMPRD;
		close SAMPAD;
		close SAMPKEYED;
		print "$printed lines printed; $missing missing data\n";	
	}
}

sub fst {
	print "\n************ running fst for $freindly_chr_scaff_name *************\n\n";
	
	# Skip?
	if (-f "$work_dir$prefix.ref01.weir.fst.sortedkeyed"){
		print "FST file exists for $prefix (whole genome), skipping\n";
		return;
	} elsif ($restrict_scaff and -f "$work_dir$prefix.ref01$dotscaff.weir.fst.sortedkeyed"){
		print "FST file exists for $prefix$dotscaff, skipping\n";
		return;
	}

	die "Requires --ref\n" if (not $ref);
	die "Requires --ref0-list\n" if (not $ref0_list);
	die "Requires --ref1-list\n" if (not $ref1_list);
	
	if (not $restrict_scaff){
		print "Running fst for whole genome (--chr appplied after conversion)\n";
		print "This might take a while\n"
	}

	print "Appending vcftools stdout to $vcftools_stdout\n";

	my @fst_call = ('vcftools --max-alleles 2');
	push @fst_call, ('--vcf', $ref);
	push @fst_call, ('--weir-fst-pop ', $ref0_list);
	push @fst_call, ('--weir-fst-pop ', $ref1_list);	# yes, double entry is the way it works
	if ($restrict_scaff){
		push @fst_call, ('--chr', $restrict_scaff);
	}
	push @fst_call, ('--out', "$work_dir$prefix.ref01$dotscaff");
	my $fst_call = join (' ', @fst_call, ">> $vcftools_stdout");
	print "\nSystem call:\n$fst_call\n";
	system ($fst_call);
	# Output of vcftools will be: prefix.ref01[.scaff_or_chrom].weir.fst
	
	# Create sorted and keyed version for mash	
	my $fst_file = "$work_dir$prefix.ref01$dotscaff.weir.fst";
	my $fst_sorted_keyed = "$fst_file.sortedkeyed";
	
	print "\nCreating sorted and keyed FST file...\n";
	open FST, $fst_file or die "Could not open $fst_file for read";
	my (%fsts, @ordered_keys);
	$_ = <FST>;
	while(<FST>){
		next if ($_ eq $/);
		chomp;
		my ($scaff, $pos, $fst) = split;
		next if ($fst eq '-nan');
		next if ($fst < 0.1);
		my $key = "$scaff-$pos";
		$fsts{$key} = $fst;
		push @ordered_keys, $key;
	}
	close FST;
		
	@ordered_keys = sort { $fsts{$b} <=> $fsts{$a} } @ordered_keys;
		# ties (specifically FST=1) should leave left-most marker first (we want this!)

	open FSTSORTKEY, '>', $fst_sorted_keyed or die "Could not open $fst_sorted_keyed for write";	
	foreach my $key (@ordered_keys){
		print FSTSORTKEY "$key\t$fsts{$key}\n";
	}
	close FSTSORTKEY;	
}

sub mash {
	# Combines allele counts (ref) and seq depth (sample) from combined keyed files above.
	# Skips markers not present in all samples.
	# Optionally picks most informative markers based on sorted FST_file using min_bp
	# to max_bp defined intervals.
	
	# This and later steps don't skip based on presence of previous output files because marker selection
	# will vary depending on input sample(s).
	
	print "\n*********** Running mash for $freindly_chr_scaff_name **********\n\n";
	
	my $ref_counts_keyed = "$work_dir$prefix.ref01$dotloc.frq.count.keyed";
	if (not -f $ref_counts_keyed) {
		if (-f "$work_dir$prefix.ref01.frq.count.keyed") {	# exists for whole genome?
			$ref_counts_keyed = "$work_dir$prefix.ref01.frq.count.keyed";
		} else {
			die "Did not find $prefix.ref01$dotloc.frq.count.keyed or $prefix.ref.frq.count.keyed in $work_dir";
		}
		# print "Using keyed ref allele count file $ref_counts_keyed\n";		
	}
	
	my $use_windows = ($min_bp and $max_bp);
	
	my $fst_sorted_keyed = "";
	if ($need_fsts){
		# die "Need --min_bp and --max-bp for FST marker pruning\n" if ($min_bp == 0 or $max_bp == 0);
		# die "Args must satisfy\: max_bp \> min_bp \> 0\n" if ($min_bp >= $max_bp or $min_bp < 1);
		if (-f "$work_dir$prefix.ref01$dotloc.weir.fst.sortedkeyed"){
			$fst_sorted_keyed = "$work_dir$prefix.ref01$dotloc.weir.fst.sortedkeyed";
		} elsif (-f "$work_dir$prefix.ref01.weir.fst.sortedkeyed"){
			$fst_sorted_keyed = "$work_dir$prefix.ref01.weir.fst.sortedkeyed";
		} else {
			die "Did not find sorted FST file $prefix.ref01$dotloc.weir.fst.sortedkeyed or $prefix.ref01.weir.fst.sortedkeyed in $work_dir\n";
		}	
		# print "Using sorted FST file $fst_sorted_keyed\n";
	}
	
	# Generate a list of common, non-missing, ordered positions and filter for --scaff or --chr
	print "Generating a common list of non-missing sample positions...\n";
	my (%filtered_keys, %alleles);
	# alleles indexted by $key; points to [ref0, ref0, ref1, ref1, s1, s1, s2, s2, ...] with n_samples * 2 + 4 columns
	my $this_sample_num = 0;
	foreach my $test_sample (@samples){
		$this_sample_num++;
		my $samp_depth_keyed = "$work_dir$prefix.$test_sample$dotloc.RDAD.keyed";
		if (not -f $samp_depth_keyed){
			if (-f "$work_dir$prefix.$test_sample.RDAD.keyed") {
				$samp_depth_keyed = "$work_dir$prefix.$test_sample.RDAD.keyed";
			} else {
				die "Did not find $prefix.$test_sample$dotloc.RDAD.keyed or $prefix.$test_sample.RDAD.keyed in $work_dir";
			}
		}
		open SAMPKEYED, $samp_depth_keyed or die "Could not open $samp_depth_keyed for read ";
		my $count = 0;
		
		# print "Debug: $restrict_scaff / $restrict_chr \n";
		
		while(<SAMPKEYED>){
			next if ($_ eq $/);
			if ($restrict_scaff){
				next unless (/^$restrict_scaff\-/);
			} else {	
				next unless (/^$restrict_chr\./);
			}			
			$count++;			
			my ($key, $rd, $ad) = /^(\S+)\t(\d+)\t(\d+)\n*$/;
			if (not $alleles{$key}){
				$alleles{$key} = [(-1) x ($n_samples * 2 + 4)]; #init ref to array with n*2+4 columns of -1
			}
			$alleles{$key}->[$this_sample_num * 2 + 4] = $rd;
			$alleles{$key}->[$this_sample_num * 2 + 5] = $ad;
			
			# $allele_counts{$key} = [$rd, $ad, -1, -1, -1, -1];  # will use -1 to skip if isn't found in ref
			
			# push @ordered_keys, $key;	#first sample will determine order
			
			$filtered_keys{$key} = $filtered_keys{$key} || 0;
			$filtered_keys{$key}++;		#  used to identify and examine only common snps
		}			
		close SAMPKEYED;	
		# print "$restrict_scaff\n"
		# print "$debug\n";
		print "$test_sample total positions: $count\n";		
	}
	
	# Figure out which snps are common to all samples
	my $n_common_nonmissing = 0;
	while (my ($key, $value) = each %filtered_keys){
		if ($value == $n_samples){
			$n_common_nonmissing++;
		} else {
			$filtered_keys{$key} = 0;	# evaluates false in tests below	
		}
	}
	print "Number of common non-missing positions: $n_common_nonmissing\n";
	
	# Filter for min_fst if there is one
	if ($min_fst){
		open FSTSORTKEY, "$fst_sorted_keyed" or die "Could not open $fst_sorted_keyed for reading";
		while(<FSTSORTKEY>){
			next if ($_ eq $/);
			my ($key, $fst) = /^(\S+)\t(\S+)\n*$/;
			last if ($fst < $min_fst);
			$filtered_keys{$key} = $filtered_keys{$key} && -999;	# mark it so all not found can be zeroed
		}
		close FSTSORTKEY;
		my $n_min_fst = 0;
		while (my ($key, $value) = each %filtered_keys){
			$value = $value || 0;
			if ($value == -999){
				$n_min_fst++;
			} else {
				$filtered_keys{$key} = 0;	# evaluates false in tests below	
			}
		}
		print "Number meeting min-fst: $n_min_fst\n";	
	}
	
	# Read ref data; prune now if simple min_bp; store alleles for later pruning if windowing
	print "\nReading reference ancestrial counts...\n";
	my @pruned_keys;
	my (@found_scaffolds, %scaff_ends);         # only used if windowing


	my $found = 0;
	my $kept = 0;
	my $pruned = 0;
	my $prev_scaff = 'Not this';
	my $prev_pos = -$min_bp;
	open REFKEYED, $ref_counts_keyed or die "Could not open $ref_counts_keyed for read ";
	while (<REFKEYED>){	# no header to skip
		my ($quickkey) = /^(\S+)\t/; 
		next if (not $filtered_keys{$quickkey});	# optimized given that the vast majority will fail here		
		next if ($_ eq $/);
		chomp;
		my ($key, $countRef0a1, $countRef0a2, $countRef1a1, $countRef1a2) = split;
		my ($scaff, $pos) = split(/\-/, $key);
		if ($use_windows or $min_bp <= ($pos - $prev_pos) or $scaff ne $prev_scaff){		# always if min_bp = 0
		# We might violate min_bp if scaffold assembly gaps < min_bp, oh well...
			my $alleles_pointer = $alleles{$key};
			$alleles_pointer->[0] = $countRef0a1;	# non -1 value here tells us that marker exists in ref & sample
			$alleles_pointer->[1] = $countRef0a2;
			$alleles_pointer->[2] = $countRef1a1;
			$alleles_pointer->[3] = $countRef1a2;
			if (not $use_windows){
				push @pruned_keys, $key;	# these are ready to use
			}
			if ($scaff ne $prev_scaff){					# only used if windowing, but overhead is small
				$scaff_ends{$prev_scaff} = $prev_pos;	# store the right-most position on each scaffold as we transition to next scaff
				push @found_scaffolds, $scaff;
				$prev_scaff = $scaff;
			}
			$kept++;
			$prev_pos = $pos;
		} else {
			$pruned++;
		}
		$found++;
	}
	$scaff_ends{$prev_scaff} = $prev_pos;
	close REFKEYED;	
		
	if (not $use_windows){
		# Everything in @pruned_keys is ready to use
		print "$found found in ref; $pruned pruned for min_bp; $kept kept\n";
		
	} else {
		print "$found found in ref; pruning for FST by scaffold...\n";
		
		# Need positions sorted by fst for each scaffold (input is already sorted high to low FST for possibly whole genome)
		my (%fst_ordered_pos_by_scaffold, %fsts);
		open FSTSORTKEY, "$fst_sorted_keyed" or die "Could not open $fst_sorted_keyed for reading";
		while(<FSTSORTKEY>){
			next if ($_ eq $/);
			# chomp;
			my ($key, $fst) = /^(\S+)\t(\S+)\n*$/;
			next if (not $alleles{$key});			# not present in all samples (the vast majority)
			next if ($alleles{$key}->[0] == -1);	# not present in ref
			my ($scaff, $pos) = split(/\-/, $key);
			$fst_ordered_pos_by_scaffold{$scaff} = $fst_ordered_pos_by_scaffold{$scaff} || []; # see below
			push @{$fst_ordered_pos_by_scaffold{$scaff}}, $pos;
			# %fst_ordered_pos_by_scaffold is indexed by scaffold
			# Each element is a pointer to an array
			# Each array is a set of positions ordered by FST (high to low)
			# The 1st position (ie, best) will be picked below consistent with min_bp to max_bp interval
			# Ties (in particular, FST = 1) leave left-most marker first, minimizing marker interval
		
			$fsts{$key} = $fst;   #informational only
		}
		close FSTSORTKEY;
		
		# We have excluded all markers that we can't use due to missing data; now pick soley from fst and position
		foreach my $scaff (@found_scaffolds){
			print "\nProcessing scaffold $scaff...\n";
			my $final_pos = $scaff_ends{$scaff};	# stop searching here	
			my $fst_ordered_pos = $fst_ordered_pos_by_scaffold{$scaff}; # ref to array holding fst-sorted positions on this scaffold
			$kept = 0;
			my $total_length = 0;
			my $total_fst = 0;
			my @used_fsts;
			# pick first pos in ordered_pos that is in allowed interval
			my $prev_pos = -$min_bp;		
			while ($prev_pos < $final_pos - $min_bp) {
				my $was_found = 0;
				foreach my $pos (@$fst_ordered_pos) {
					if ($pos >= $prev_pos + $min_bp and $pos <= $prev_pos + $max_bp){ # find 1st (ie, best) position in interval
						my $key = "$scaff-$pos";
						push @pruned_keys, $key;	# this is what we need
						$kept++;	
						$total_length += $pos - $prev_pos;
						$total_fst += $fsts{$key};
						push @used_fsts, $fsts{$key}; # feedback information only
						$prev_pos = $pos;
						$was_found = 1;
						last;
					}
				}
				if (not $was_found){	# no snps in allowed window; push to right and violate max_bp
					print " ! No SNPs found in interval ", ($prev_pos + $min_bp), " to ", ($prev_pos + $max_bp), "\n";
					$prev_pos += $min_bp;
				}
			}
			$pruned = $found - $kept;
			print "$pruned pruned for min_bp max_bp intervals; $kept kept\n";
			if ($kept > 0){
				my $ave_interval = $total_length / $kept;
				my $ave_fst = $total_fst / $kept;
				@used_fsts = sort {$b <=> $a} @used_fsts;
				my $q1 = $used_fsts[int ($kept / 4)];
				my $med = $used_fsts[int ($kept / 2)];
				my $q3 = $used_fsts[int ($kept * 0.75)];
				print "Average marker interval: $ave_interval\n";
				print "Average marker FST: $ave_fst\n";		
				print "Upper quartile:     $q1\n";		
				print "Median:             $med\n";		
				print "Lower quartile:     $q3\n";	
			} else {
				print "Kept no markers for this scaffold!\n";
			}		
		}
	} # end of fst pruning
	
	# Print mash file for each sample
	print "\nPrinting mash file(s) for:";
	$this_sample_num = 0;
	foreach my $test_sample (@samples){
		$this_sample_num++;
		print " $test_sample";		
		my $mash_file = "$work_dir$prefix.$test_sample$dotloc.mash";     # output name
		open MASH, '>', $mash_file or die "Could not open $mash_file for writing";
		foreach my $key (@pruned_keys){
			my ($scaff, $pos) = split(/\-/, $key);
			my $alleles_pointer = $alleles{$key};
			# output is scaff, pos, ref0, ref0, ref1, ref1, samp, samp
			print MASH join ("\t", $scaff, $pos, $alleles_pointer->[0], $alleles_pointer->[1], $alleles_pointer->[2], $alleles_pointer->[3], $alleles_pointer->[$this_sample_num * 2 + 4], $alleles_pointer->[$this_sample_num * 2 + 5]), "\n";

		}
		close MASH;
	}
	print "\n";
	# Convert to chromosome
	if ($chr_convert_option){
		undef %alleles;			# we are done but conversion will need more memory before these go out of scope,
		undef %filtered_keys;	#  so garbage collect now
		undef @pruned_keys;
		
		print "\nConverting scaffolds to chromosomes...\n";
		foreach my $test_sample (@samples){	
			my $mash_file = "$work_dir$prefix.$test_sample$dotloc.mash";		
			scaff2chrom($mash_file, $scaff_map, "$mash_file.onChr", 0, $force_orientation_option); # this will generate new file with .onChr extention
		}
	}	
}	

sub mash2ahmm {
	# Generates input file for Ancestry_HMM from this program's *.mash or *.mash.onChr.
	# Appends recomb rate (in Morgans) as right hand column and strips 1st chrom column.
	# Recomb map must follow this format (last column is Morgans):
	# 1	50001	60000	+	0
	# 1	60001	70000	+	0.023255814
	# 1	70001	80000	+	0
	# 1	80001	90000	+	0
	# 1	90001	100000	+	0
	# 1	100001	110000	+	0.023255814
	# 1	110001	120000	+	0
	#
	# (strand polarity used in case this is a scaffold rather than chrom map)
	
	# Assumes that input is already filtered for --scaff or --chr; therefore, there will
	# be consistancy of direction whithin this function call
	
	print "\n*********** Running mash2ahmm for $freindly_chr_scaff_name **********\n\n";
	die "Provide --recomb-map or --recomb-by-chr\n" if ($recomb_map and $recomb_by_chr);
	die "-c option must be used with --recomb-by-chr" if ($recomb_by_chr and not $chr_convert_option);

	# Assumes that input is already filtered for --scaff or --chr	
	
	# Do scaffold to chromosome conversion if not done already
	# Also set mash file to 1st sample (to be used for recomb rate determination for all samples)
	my $mash_file = "$work_dir$prefix.$samples[0]$dotloc.mash";	
	if ($chr_convert_option){
		if (not -f "$mash_file.onChr"){		# May have it already from end of last step
			die "Did not find file $mash_file or $mash_file.onChr in $work_dir" if (not -f $mash_file);
			# Do it for all samples even though only tested first sample
			foreach my $conv_sample (@samples){
				my $conv_mash_file = "$work_dir$prefix.$conv_sample$dotloc.mash";
				scaff2chrom($conv_mash_file, $scaff_map, "$conv_mash_file.onChr", 0, $force_orientation_option);
			}
		}
		$mash_file = "$mash_file.onChr";
	} else {
		die "Did not find file $work_dir$mash_file" if (not -f $mash_file);
	}
	
	my @interval_recomb_rates;
	my $n_intervals = 0;
	if ($recomb_by_chr){
		# Get recomb by chr and convert from cM/Mb to Morgans/bp:  1/(100*1000000) = 1/100000000
		my %morgans_per_bp_by_chr;
		open RECOMBBYCHR, $recomb_by_chr or die "Could not open $recomb_by_chr for reading";
		while (<RECOMBBYCHR>){
			next if ($_ eq $/);
			chomp;
			my ($chr, $cMperMb) = split;
			$morgans_per_bp_by_chr{$chr} = $cMperMb / 100000000;
		}
		close RECOMBBYCHR;
		
		# Assume that all input sample files have same markers (but test line number for safety)
		# So we can calculate interval rates once for $sample[0] and then construct all sample outputs
		print "\nProcessing $samples[0] for marker interval recomb rates\n";
		print "These will be appended to all samples to create Ancestry_HMM input files...\n";
		my $total_intervals_length = 0;
		my $total_morgans = 0;
		my $left_snp = 0;
		my $this_is_first_one = 1;
		
		open MASH, $mash_file or die "Could not open $mash_file for reading";		
		while (<MASH>){
			# We need to generate recomb rate (Morgans) between previous and current position
			next if ($_ eq $/);
			chomp;
			my ($chr, $right_snp) = /(\S+)\t(\d+)/;
			if ($this_is_first_one){
				push @interval_recomb_rates, 0.5;	# 1st value is ignored		
				$left_snp = $right_snp;		
				$this_is_first_one = 0;
				next;
			}
			$n_intervals++;	# will be number lines - 1
			my $length = $right_snp - $left_snp;
			my $morgans = $morgans_per_bp_by_chr{$chr} * $length;
			push @interval_recomb_rates, $morgans;
			$total_intervals_length += $length;
			$total_morgans += $morgans;
			$left_snp = $right_snp;	# set up for next interval
		}
		close MASH;
		print "Map length in bp (w/ recomb data): $total_intervals_length\n";
		print "Morgans: $total_morgans\n";
		my $ave_rate_infile = 'div/0';
		if($total_intervals_length > 0){
			$ave_rate_infile = 100000000 * $total_morgans / $total_intervals_length;
		}
		print "Ave rate (cM/Mb): $ave_rate_infile\n";
		
	} else {
		# Process recomb map once into memory
		my $use_restriction = $restrict_scaff;
		if ($restrict_chr){
			$use_restriction = $restrict_chr;
		}
		my (@left_window, @right_window, @rate);
		my $window_dir = 'init';
		my $total_map_length = 0;
		my $total_map_morgans = 0;
		my $n_windows = 0;
		open RECOMBMAP, $recomb_map or die "Could not open $recomb_map for reading";
		while (<RECOMBMAP>){	# Assumes no header!
			next if ($_ eq $/);
			chomp;
			my ($chr, $left_window, $right_window, $this_window_dir, $rate) = split;	
			# $chr could be either scaffold or chromosome depending on inputs)
			next if ($use_restriction and $chr ne $use_restriction);
			push @left_window, $left_window;
			push @right_window, $right_window;
			push @rate, $rate;
			if ($window_dir eq 'init'){
				$window_dir = $this_window_dir;
			} 
			die "What!? Scaffold switched directions!" if ($window_dir ne $this_window_dir);
			$n_windows++;
			$total_map_length += $right_window - $left_window;
			$total_map_morgans += $rate;
		}
		close RECOMBMAP;
		if ($window_dir eq 'init') {
			die "\nScaffold $use_restriction appears to be missing from recomb map; can't determine rate\n\n";
		} elsif ($window_dir ne '+' and $window_dir ne '-') {
			# if chromosome analysis, all should be (+) so we never get here
			# if scaffold analysis, then user has suppled an unoriented scaffold
			die "\nScaffold $use_restriction appears to be unoriented; can't determine recomb rates\n\n";		
		}	
		
		print "Processed data from recomb map $freindly_chr_scaff_name ";	
		print "($window_dir):\n";

		
		print "Map length in bp: $total_map_length\n";
		print "Morgans: $total_map_morgans\n";
		my $ave_rate = 100000000 * $total_map_morgans / $total_map_length;
		print "Ave rate (cM/Mb): $ave_rate\n";
		print "Expect small difference below due to interval boundaries\n";
		
		# Assume that all input sample files have same markers (but test line number for safety)
		# So we can calculate interval rates once for $sample[0] and then construct all sample outputs

		print "\nProcessing $samples[0] for marker interval recomb rates\n";
		print "These will be appended to all samples to create Ancestry_HMM input files...\n";

		my $total_intervals_length = 0;
		my $total_morgans = 0;
		my $gaps = 0;
		my $left_snp = 0; 
		# my $printed_lines = 0;
		my $this_is_first_one = 1;
		my $ran_off_the_end = 0;
		
		# Slide map windows to right, accounting for orientation reversal in (-)
		my $i = 0;		# starting window index for (+) orientation
		my $incr = 1;	# increment forward through map if (+)
		if ($window_dir eq '-'){
			$i = $n_windows - 1;
			$incr = -1;	# increment backward through map if (-)
		}
	
			
		open MASH, $mash_file or die "Could not open $mash_file for reading";		
		MASHLOOP: while (<MASH>){
			# We need to generate recomb rate (Morgans) between previous and current position
			# Use 0.5 if any portion of interval falls into gaps in recomb map (are these assembly gaps???)
			# If interval includes 10% of a window then we calculate 10% of that window's rate.
			# We could have markers at the end that fall outside map windows, in which case bail out without printing it
			next if ($_ eq $/);
			chomp;
			my ($chr, $right_snp) = /(\S+)\t(\d+)/;
			if ($this_is_first_one){
				push @interval_recomb_rates, 0.5;
				# $printed_lines++;			
				$left_snp = $right_snp;		
				$this_is_first_one = 0;
				next;
			}
			$n_intervals++;	# will always be #lines - 1
			if ($ran_off_the_end){	# this is a permanent condition; all remaining intervals get 0.5
				push @interval_recomb_rates, 0.5;
				$gaps++;
				next MASHLOOP;	
			}
			
			# We now have $left_snp and $right_snp defining the interval we want to measure
			my $morgans = 0;		
			while ($left_snp > $right_window[$i]){	# slide map windows right until right window reaches left_snp
				$i += $incr;						# advance sliding windows to right
				if ($i >= $n_windows or $i < 0){	#  off the end of map
					$ran_off_the_end = 1;	# write 0.5 for remaining lines and don't try again
					print " ! Ran off the end of map at snp interval $left_snp $right_snp (using 0.5)\n";
					next MASHLOOP;
				}
			}
			if ($left_snp < $left_window[$i]){	# Left snp is in a map gap; use 0.5 since we have no idea
				print " ! Gap in recomb map at left snp $left_snp (using 0.5)\n";
				push @interval_recomb_rates, 0.5;
				$gaps++;
				$left_snp = $right_snp;		# setup for next interval
				next MASHLOOP;
			}
			if ($right_snp <= $right_window[$i]){	# marker interval wholey contained within one window
				# portion of window covered by entire interval
				$morgans = $rate[$i] * (($right_snp - $left_snp) / ($right_window[$i] - $left_window[$i]));		
			} else {	# found left side of interval; we'll need to step through windows to find right side
				# portion of right side of window window covered by left side of interval
				$morgans = $rate[$i] * (($right_window[$i] - $left_snp) / ($right_window[$i] - $left_window[$i]));
				my $gap_test_pos = $right_window[$i];	# test this against next window to identify gap
				$i += $incr;
				if ($i >= $n_windows or $i < 0){   #  off the end of map
					$ran_off_the_end = 1;
					print " ! Ran off the end of map at snp interval $left_snp $right_snp (using 0.5)\n";
					next MASHLOOP;
				}
				while ($right_snp > $right_window[$i]){	# move windows right to find right side of interval; watch for map gaps
					if ($gap_test_pos < $left_window[$i] - 1){	# Map gap; use 0.5 since we have no idea
						print " ! Gap in recomb map $gap_test_pos to $left_window[$i] (using 0.5)\n";
						push @interval_recomb_rates, 0.5;
						$gaps++;	
						$left_snp = $right_snp;		# setup for next interval
						next MASHLOOP;
					}
					$morgans += $rate[$i]; # add entire window (this is a big interval)
					$gap_test_pos = $right_window[$i];
					$i += $incr;
					if ($i >= $n_windows or $i < 0){   #  off the end of map
						$ran_off_the_end = 1;
						print " ! Ran off the end of map at snp interval $left_snp to $right_snp (using 0.5)\n";
						next MASHLOOP;
					}
				}
				if ($right_snp < $left_window[$i]){
					# Right snp is in a map gap; use 0.5 since we have no idea
					print " ! Gap in recomb map at right snp $right_snp (using 0.5)\n";
					push @interval_recomb_rates, 0.5;
					$gaps++;	
					$left_snp = $right_snp;		# setup for next interval
					next MASHLOOP;					
				}
				# portion of left side of window covered by right side of interval
				$morgans += $rate[$i] * (($right_snp - $left_window[$i]) / ($right_window[$i] - $left_window[$i]));	
			}
			# If we didn't bail for some reason with 0.5 value, here is real calculated rate
			push @interval_recomb_rates, $morgans;		
			$total_intervals_length += $right_snp - $left_snp;
			$total_morgans += $morgans;		
			$left_snp = $right_snp;		# setup for next interval
		}
		close MASH;	
		if ($ran_off_the_end){	# need one more
			push @interval_recomb_rates, 0.5;
			$gaps++;	
		}	
		# close OUTFILE;
		
		#debug
		if ($n_intervals != $#interval_recomb_rates){
			die "$n_intervals  $#interval_recomb_rates";
		}

		
		print "\nIntervals without recomb data: $gaps (used 0.5 in file)\n";
		print "Number of intervals with recomb data: ", $n_intervals - $gaps, "\n";
		print "Map length in bp (w/ recomb data): $total_intervals_length\n";
		print "Morgans: $total_morgans\n";
		my $ave_rate_infile = 'div/0';
		if($total_intervals_length > 0){
			$ave_rate_infile = 100000000 * $total_morgans / $total_intervals_length;
		}
		print "Ave rate (cM/Mb): $ave_rate_infile\n";
	}


	
	# Open sample mash files and re-write as Ancestry_HMM input files with recomb rate (and no 1st chr column)
	
	print "\nWriting Ancestry_HMM input files for:";
	foreach my $test_sample (@samples){
		print " $test_sample";
		my $mash_file = "$work_dir$prefix.$test_sample$dotloc.mash";	
		$mash_file = "$mash_file.onChr" if ($chr_convert_option);
		my $ahmm_in = "$work_dir$prefix.$test_sample$dotloc.ahmm.in";	# output name

		open MASH, $mash_file or die "Could not open $mash_file for reading";		
		open OUTFILE, '>', $ahmm_in or die "Could not open $ahmm_in for writing";
		
		# Number lines in all files had better match $n_intervals
		my $line_number = 0;
		my $print_line_break = 0;
		while(<MASH>){
			next if ($_ eq $/);
			
			## If this line AND THE NEXT are 0.5, then don't print this line
			my $next_line_recomb_rate = $interval_recomb_rates[$line_number + 1] || 0.5;
			if ($interval_recomb_rates[$line_number] < 0.5 or $next_line_recomb_rate < 0.5){
			
				chomp;
				s/\S+\t//;	#remove 1st column (chromosome) but keep all the rest
				if ($print_line_break) {
					print OUTFILE "\n" if ($line_number > 0);
				} else {
					$print_line_break = 1;
				}
				print OUTFILE "$_\t$interval_recomb_rates[$line_number]";	# print with appended recomb rate
				# Ancestry_HMM is buggy if there is newline at eof; hence the awkward \n construction				
			}
			$line_number++;

		}
		if ($line_number != $n_intervals + 1){
			die "$line_number  $n_intervals ";
		}
		close MASH;
		close OUTFILE;
	}
	print "\n";	
}

sub ahmm {
	# WARNING! Ancestry_HMM provides real(ish) looking results if inputs missing, or goes
	# into infinite memory leak if input line is misformated
	# BUG! Ancestry_HMM doubles up the last snp in the file if input has a final line return
	print "\n************ running ahmm for $freindly_chr_scaff_name *************\n\n";

	# Test input file integrity and build the -i argument for Ancestry_HMM
	my $i_input = "";
	foreach my $test_sample (@samples){
		my $ahmm_in = "$work_dir$prefix.$test_sample$dotloc.ahmm.in";				# input name
		my $ahmm_out = "$data_dir$prefix.$test_sample$dotloc$dotsuffix.ahmm.out";	# output name
		$i_input .= " -i $ploidy $ahmm_in $ahmm_out";											

		# Check input since Ancestry_HMM generates real-looking results from a blank file
		open INFILE, $ahmm_in or die "Could not open $ahmm_in for read";
		# Looks like:
		# 121496	0	2	4	0	24	3	0
		# 140811	21	1	0	18	5	5	0.022776696319832
		# 151623	0	22	17	1	12	7	0.023734931680168
		# 171534	8	14	18	0	20	2	0
		my $j = 0;
		while (<INFILE>){
			chomp;
			my @test = split(/\t/);
			foreach my $test (@test[0..6]){
				if (not defined($test) or $test =~ /\D/){
					die "$ahmm_in appears to be empty, misformated or corrupted";
				}
			}
			if ($test[7] !~ /^[\d\.\-eE]+$/){
				die "$ahmm_in appears to be empty, misformated or corrupted";
			}
			$j++;
			last if ($j > 10);
		}
		die "$ahmm_in appears to be empty, misformated or corrupted" if ($j < 5);
		close INFILE;
		
		# We could get confused if there is already and outfile and ancestry_hmm fails during run
		if (unlink $ahmm_out){	# deletes existing file if there was one
			print "WARNING! overwriting $ahmm_out\n";
		}
	}
	
	# Name the stdout file (w/ warning and safety delete if overwriting)
	my $ahmm_stdout = "$data_dir$prefix$dotloc$dotsuffix.ahmm.stdout";				# stdout name
	if (unlink $ahmm_stdout){	# deletes existing file if there was one
		print "WARNING! overwriting $ahmm_stdout\n";
	}	
		
	# Call Ancestry_HMM
	# Example: ancestry_hmm -a 0.2 -i 20 Sample1_1.1_10-20kbAHMM.in Sample1_1.1_10-20kbAHMM.out [repeat -i for each sample]
	my $ahmm_call = "ancestry_hmm $ahmm_args $i_input > $ahmm_stdout 2>&1";	# redirect stdout and capture stderr
	print("\nHere is the system call...\n");
	print ("$ahmm_call\n");
	system ($ahmm_call);
	
	# Open stdout and get summary info
	my $ad_time = 'not found';
	my $ancestry = 'not found';
	my $lnL = 'not found';
	open AHMMSTDOUT, $ahmm_stdout or die "Could not open $ahmm_stdout for read";
	while(<AHMMSTDOUT>){
		if(/estimated_admixture time\:\s*([\d\.]+)/){
			$ad_time = $1;
		}
		if(/estimated ancestry\:\s*([\d\.]+)/){
			$ancestry = $1;
		}
		if(/\tlnL\:\s*([\-\d\.]+)/){
			$lnL = $1;
		}
	}
	close AHMMSTDOUT;
		
	# Log it
	open LOG, '>>', $ahmm_log or die "Could not open $ahmm_log for appending";
	my $time = localtime();
	my $full_arg_str = join (' ', @full_arg_list);
	print LOG join ("\t", $time, $prefix, $suffix, ($sample_name || $sample_list), $restrict_chr, $restrict_scaff, $ad_time, $ancestry, $lnL, $full_arg_str), "\n";
	close LOG;
	
	print "\nInfo from stdout writen to log file...\n";
	print "Estimated admixture time : $ad_time\n";
	print "Estimated ancestry       : $ancestry\n";
	print "lnL                      : $lnL\n";
}

sub ahmmpost {
	print "\n************ running ahmmpost for $freindly_chr_scaff_name *************\n\n";
	
	my $ahmm_stdout = "$data_dir$prefix$dotloc$dotsuffix.ahmm.stdout";				# stdout from Ancestry_HMM
	my $ahmm_lnL = "$data_dir$prefix$dotloc.ahmm.lnL";								# appends time-lnL pairs from stdout for all runs with prefix
	my $ahmm_maxpost = "$data_dir$prefix$dotloc$dotsuffix.ahmm.maxpost";			# files of sample max postior proportions
	my $ahmm_maxpost95ci = "$data_dir$prefix$dotloc$dotsuffix.ahmm.maxpost95ci";	# above with 95% CIs
	my $ahmm_mean_ancestry = "$data_dir$prefix$dotloc$dotsuffix.ahmm.mean_ancestry"; # contains mean maxpost +- SD

	# Open stdout and get run info including golden section search
	my $ad_time = 'not found';
	my $ancestry = 'not found';
	my $lnL = 'not found';
	my %all_lnLs;
	open AHMMSTDOUT, $ahmm_stdout or die "Could not open $ahmm_stdout for read";
	while(<AHMMSTDOUT>){
		if(/estimated_admixture time\:\s*([\d\.]+)/){
			$ad_time = $1;
		}
		if(/estimated ancestry\:\s*([\d\.]+)/){
			$ancestry = $1;
		}
		if(/\tlnL\:\s*([\-\d\.]+)/){
			$lnL = $1;
		}
		if(/(\d+)\t([\d\.]+)\t([\d\.]+)\t([\d\.]+)\t([\d\.]+)\t([\-\d\.]+)\t([\-\d\.]+)/){ # golden section search
			my $low_test = $3;
			my $high_test = $4;
			my $lnL_low = $6;
			my $lnL_high = $7;
			$all_lnLs{$low_test} = $lnL_low;
			$all_lnLs{$high_test} = $lnL_high;
		}
	}
	close AHMMSTDOUT;
	if ($ad_time ne 'not found'){
		# Ancestry_HMM calculates one more at the midpoint of the last section (the summary value)
		$all_lnLs{$ad_time} = $lnL;	
	}
	
	# Append all Time-lnL pairs
	if (-f $ahmm_lnL){
		open LNL, '>>', $ahmm_lnL or die "Could not open $ahmm_lnL for appending";
	} else {
		open LNL, '>', $ahmm_lnL or die "Could not open $ahmm_lnL for writing";
		print LNL join ("\t", "Prefix", "Ancestry", "Suffix", "Time", "lnL"), "\n";
	}
	my @lnL_keys = sort {$a <=> $b} (keys %all_lnLs);	
	foreach my $key (@lnL_keys){
		print LNL join ("\t", $prefix, $ancestry, $suffix, $key, $all_lnLs{$key}), "\n";
	}
	close LNL;

	# Calculate max posterior probability proportions (w/ and w/out CI) and store value for sample mean & SD calculation
	my @positions;
	my (%maxpost, %low95, %high95);	# indexed by test_sample, holds array references (each array pushed in order of positions)
	my %maxposts_by_pos; # hash of array refs, structured for mean & SD calculation
	# my (@lnL_mean, @lnL_SD_plus, @lnL_SD_minus);
	foreach my $test_sample (@samples){
		my $ahmm_out = "$data_dir$prefix.$test_sample$dotloc$dotsuffix.ahmm.out";			# output name
		open AHMMOUT, $ahmm_out or die "Could not open $ahmm_out for read";
		while(<AHMMOUT>){
			chomp;
			my @line = split;
			my $haplotypes = $#line - 1;	# position plus n+1 columns
			my $max_post_value = 0;
			my $sum_posteriors = 0;
			my $max_post_index = 0;
			my $low95_index = 0;
			my $high95_index = 0;
			for (my $i = 1; $i < $haplotypes; $i++){
				my $posterior = $line[$i];
				if ($posterior > $max_post_value or ($posterior == $max_post_value and rand(2) < 1)){
					$max_post_value = $posterior;
					$max_post_index = $i;
				}
				$sum_posteriors += $posterior;
				if (not $low95_index and $sum_posteriors > 0.025){	#0 is false
					$low95_index = $i;	#ci includes this proportion
				}
				if (not $high95_index and $sum_posteriors > 0.975){
					$high95_index = $i;	#ci includes this proportion
				}
			}
			my $max_post_proportion = ($max_post_index - 1) / $haplotypes;  #eg, 2 haplotypes: index 1,2,3 -> 0%, 50%, 100%	
			my $low95ci = ($low95_index - 1) / $haplotypes;
			my $high95ci = ($high95_index - 1) / $haplotypes;
			
			if ($test_sample eq $samples[0]){
				push @positions, $line[0];	# same for all samples so only do once
				# push @maxpost_array_refs, [];	# init a reference to an array
			}
			# arrays pushed in the same order as @positions, indexed by sample
			push @{$maxpost{$test_sample}}, $max_post_proportion;
			push @{$low95{$test_sample}}, $low95ci;
			push @{$high95{$test_sample}}, $high95ci;
			# these arrays will hold all values for a given position; used to calculate position mean and SD
			push @{$maxposts_by_pos{$line[0]}}, $max_post_proportion;
		}
		close AHMMOUT;
	}
	
	# Print mean ancestry +- SD
	open MEANMAXPOST, '>', $ahmm_mean_ancestry or die "Could not open $ahmm_mean_ancestry for write\n";
	print MEANMAXPOST "\tMean\tSD\n";
	foreach my $pos (@positions){
		my ($mean, $sd) = mean_and_sd ($maxposts_by_pos{$pos});	# expects array ref
		print MEANMAXPOST join ("\t", $pos, $mean, $sd), "\n";
	}
	close MEANMAXPOST;
	
	# Print max posterior proportiona and CIs column-wise, and mean
	open MAXPOST, '>', $ahmm_maxpost or die "Could not open $ahmm_maxpost for write\n";
	open MAXPOST95, '>', $ahmm_maxpost95ci or die "Could not open $ahmm_maxpost95ci for write\n";
	print MAXPOST join ("\t", "", @samples), "\n";					# print header with all sample names
	print MAXPOST95 "\t";
	foreach my $test_sample (@samples){			
		print MAXPOST95 "$test_sample\t\t";							# 3 columns per sample
	}
	print MAXPOST95 "\n";
	for (my $i = 0; $i <= $#positions; $i++){
		print MAXPOST $positions[$i];								# print position
		print MAXPOST95 $positions[$i];		
		foreach my $test_sample (@samples){
			my $max_post_proportion = $maxpost{$test_sample}->[$i];
			my $low95ci = $low95{$test_sample}->[$i];
			my $high95ci = $high95{$test_sample}->[$i];
			print MAXPOST "\t$max_post_proportion";					# print values column-wise
			print MAXPOST95 "\t$max_post_proportion\t$low95ci\t$high95ci";				
		}
		print MAXPOST "\n";								
		print MAXPOST95 "\n";		
	}	
	close MAXPOST;
	close MAXPOST95;	
}

sub makeplots {
	print "\n************ running makeplots for $freindly_chr_scaff_name *************\n\n";
	
	# TO DO: make work for scaffold
	die "makeplots only works on chromosomes; fix is on my TO DO list" if (not $restrict_chr);

	my $ahmm_maxpost = "$data_dir$prefix$dotloc$dotsuffix.ahmm.maxpost";	# files of sample max postior proportions
	my $outsuffix = "$data_dir$prefix$dotloc$dotsuffix.mean_ancestry";		# output plot (R adds .png)
	
	# R
	my $r_script = q`#!/usr/bin/env Rscript

###
# Plot Ancestry Lineages
###

#takes in an output from CW's ancestry plotting 
#as of right now, only takes a single chromosome for 1-30 samples

#to do:
	#Generalize a bit more....
	#add in errors for missing packages
	#colour by scaffold?
	
#CW usage: Rscript ancestryPlots.R <FILE1 .ahmm.maxpost> <chromosome> <scaffold map> <outfile>
#e.g: Rscript ancestryPlots.R nofilter.16.Ne5k.ahmm.maxpost  16

#Load packages --------------------
require(ggplot2)
require(ggthemes)
require(plyr)
require(wesanderson)
require(reshape)
require(scales)
require(grid)

#load data.frames --------------
args = commandArgs(trailingOnly=TRUE)
anc =  read.table(file= args[1] ,header=F, skip=1)
chr = args[2]
samps = ncol(anc) -1
names(anc) = c("POS", paste("SAMPLE",c(1:samps),sep=""))
scaffs = read.table(file=args[3],header=T)	# CW, no longer hard-coded

#anc$relPOS = seq(1,nrow(anc),1)
#anc$POS = NULL
#anc.lineage = 1-anc
#sudo Rscript ancestryPlots.R /home/charlie/PooledAfz/pooledahmm_data/nofilter.16.Ne5kHiRes.ahmm.maxpost 16

#Calculate MN and SD for lineage -------------------------
mean.anc = apply(anc[-1],1 , mean,na.rm=T)
sd.anc = apply(anc[-1],1 , sd,na.rm=T)

#create dataframes -------------------------------------
mean.anc.df = data.frame(cbind(mn = mean.anc,sd = sd.anc,pos = anc$POS))
scaffs = scaffs[scaffs$chromosome==chr,]



# plot --------------------------------


mean.anc.plot = ggplot(mean.anc.df,aes(x=pos,y=mn)) +
				geom_rect(data = scaffs, 
							aes(xmin = start, xmax = stop, ymin = -0.05, ymax = 0.65) , #add in rectangles for scaffolds
							fill="grey", alpha=0.5, inherit.aes = FALSE) +
				geom_errorbar(aes(ymin = mn - sd, ymax = mn + sd), color = wes_palette("GrandBudapest")[1]) +				
				geom_point(aes(x=pos,y=mn), color = wes_palette("GrandBudapest")[3],  size = 2.3) +
				
				theme_bw() +
				scale_x_continuous(labels = comma) +
				scale_y_continuous(limits = c(-0.05, .65)) +
					theme(strip.background = element_blank(),
							axis.line.x = element_line(size = 1, colour = "black"),
							axis.line.y = element_line(size = 1, colour = "black"),
							axis.line = element_line(size = 1, colour = "black"),
							text = element_text(size=18),
							axis.ticks = element_line(size = 1), 
							panel.grid.major = element_blank(),
							panel.grid.minor = element_blank(),
							panel.background = element_blank(), 
							panel.border = element_blank()
							) +						
				labs(x = "Position (bp)", y = "Proportion M Ancestry") 

				
				
				
#mean.anc.plot



ggsave(filename = paste(args[4], ".png",sep=""), mean.anc.plot, width = 12, height = 6, units = "in" )	


	`;

	print "Creating temp.R file and calling it...\n";	
	open RSCRIPT, '>temp.R' or die "Could not open temp.R for write";
	print RSCRIPT $r_script;
	close RSCRIPT;
	
	#CW usage: Rscript ancestryPlots.R <FILE1 .ahmm.maxpost> <chromosome> <scaffold map> <outfile>
	my $sys_call = "Rscript temp.R $ahmm_maxpost $restrict_chr $scaff_map $outsuffix";
	print $sys_call, "\n\n";
	system ($sys_call);
}


############  UTILITY FUNCTIONS  #############

sub mean_and_sd {
	# expects array ref rather than array (faster this way) so just pass \@array rather than @array
	my $array_ref = $_[0];
	my $n = 0;
	my $sum = 0;
	foreach (@$array_ref){
		$sum += $_;
		$n++;
	}
	if ($n == 1){
		return ($sum, 0);
	}
	my $mean = $sum / $n;
	my $sqtotal = 0;
	foreach (@$array_ref) {
			$sqtotal += ($mean - $_) ** 2;
	}
	my $std = ($sqtotal / ($n - 1)) ** 0.5;
	return ($mean, $std);
}


sub scaff2chrom {
	# Born as AZ script, now mangled beyond recognition
	# Takes any input file with scaffold, pos at first 2 columns and outputs
	# same file contents with chrom, pos at first 2 columns (retaining all
	# additional columns). Sorts output by sensible chromosome/position order.
	
	# Pass 4th true argument if infile has a headerline
	# Pass 5th true argument to force unoriented scaffolds (0) as if they were (+);
	#
	# Skips all data on unoriented scaffolds! (unless arg 5 true)

	# Takes a file with scaffold positions on chromosomes in amel4.5, as follows
	#	scaffold	scaffold:start	scaffold:stop	chr	start	stop	orientation
	#	1.1	1	1382403	1	1	1382403	+
	#	1.10	1	1405242	1	5928187	7333428	+
	#	1.11	1	2557	1	7383429	7385985	-
	#	1.12	1	4728	1	7435986	7440713	+
	#	1.13	1	5741	1	7490714	7496454	0

	# Be sure to strip "Group", "chr", etc. from map file; expects only digits and '.'
	
	# Expects input file names: file_to_conv, scaff_map, outfile
	print "\nConverting $_[0]\n";
	# print "Using map file: $_[1]\n";
	# print "Outfile: $_[2]\n";
		
	open INFILE, $_[0] or die "Could not open $_[0] for read";
	open MAP, $_[1] or die "Could not open $_[1] for read";	
	open OUTFILE, '>', $_[2] or die "Could not open $_[2] for write"; 	
	
	my $has_headerline = $_[3];
	my $force_orientation = $_[4];
	
	my %map;	# will hold elements [scaff_start, scaff_stop, chr, start, stop, orientation]	
	$_ = <MAP>; # skip the header line
	while(<MAP>){
		chomp;
		my @line = split;
		$map{$line[0]} = [@line[1..6]];
	}
	close MAP;

	# Read the file to be converted and store everything in @outdata
	# with elements [chr, pos, remaining_line_to_print]
	my @outdata;
	my $n_unfound_scaffolds = 0;
	my $n_discarded_unoriented = 0;
	
	if ($has_headerline){
		$_ = <INFILE>;
	}
	while(<INFILE>){
		next if ($_ eq $/);

		chomp;
		my ($test_scaff, $pos, $extra) = /(\S+)\t(\d+)\t(.*)/;
		my $mapref = $map{$test_scaff};
		if (not $mapref){
			$n_unfound_scaffolds++;
			print "Didn't find: $test_scaff  $pos\n";
			next;
		}
		# these indexes shifted down 1 from file columns (now 0-5)
		my $chr = $mapref->[2]; #chr is the chromsome where the saffold is found
		my $dir = $mapref->[5]; #direction
		my $scaff_start = $mapref->[0]; #first base of scaffold
		my $chr_start = $mapref->[3]; #first base of scaffold on chr
		my $chr_stop = $mapref->[4]; #last base of scaffold on chr 
		my $newpos = 0;
		if ($dir eq '-'){
			$newpos = $chr_stop - ($pos - $scaff_start);
		} elsif ($dir eq "+" or $force_orientation) {
			$newpos = ($pos - $scaff_start) + $chr_start;
		} else {
			$n_discarded_unoriented++;
			next;
		}
		push @outdata, [$chr, $newpos, $extra];
	}
	close INFILE;
	print "$n_unfound_scaffolds not found; $n_discarded_unoriented discarded in unoriented scaffolds; $#outdata converted\n";
	
	# Reorder in a sensible way
	@outdata = sort {		# $a and $b are elements of the sorting array (and are array references)
							# return -1 if $a should be before $b, otherwise 1
		if ($a->[0] eq $b->[0]){			# equal chr
			if ($a->[1] < $b->[1]) { -1 }	# sort on pos
			else { 1 }
		}
		elsif ($a->[0] < $b->[0]) { -1 }	# sort on chr
		else { 1 }	
	} @outdata;
	
	# Print
	foreach my $outref (@outdata){
		my $chr = $outref->[0];
		my $pos = $outref->[1];
		my $extra = $outref->[2];
		print OUTFILE "$chr\t$pos\t$extra\n";
	}
	close OUTFILE;
}

# These are no longer used but here in case I need them elsewhere

sub join_intersect {
	my $has_header = $_[3];
	open FILE1, $_[0] or die "Could not open $_[0] for read";
	open FILE2, $_[1] or die "Could not open $_[1] for read";
	open OUT, '>', $_[2] or die "Could not open $_[2] for write";
	my $f1_line = <FILE1>;
	my $f2_line = <FILE2>;
	chomp($f1_line, $f2_line);
	my @f1_elements = split(/\t/, $f1_line);
	my @f2_elements = split(/\t/, $f2_line);
	my $f1_n = $#f1_elements;	# not counting joining [0] column
	my $f2_n = $#f2_elements;
	if ($has_header){
		print OUT join ("\t", @f1_elements, @f2_elements[1 .. $f2_n]), "\n";	#1st row headers
		$f1_line = <FILE1>;
		chomp($f1_line) if ($f1_line);
		$f2_line = <FILE2>;
		chomp($f2_line) if ($f2_line);
	}
	while($f1_line and $f2_line){
		@f1_elements = split(/\t/, $f1_line);
		@f2_elements = split(/\t/, $f2_line);
		if ($f1_elements[0] == $f1_elements[0]){
			print OUT join ("\t", @f1_elements, @f2_elements[1 .. $f2_n]), "\n";
			$f1_line = <FILE1>;
			chomp($f1_line) if ($f1_line);
			$f2_line = <FILE2>;			
			chomp($f2_line) if ($f2_line);
		} elsif ($f1_elements[0] < $f1_elements[0]){
			$f1_line = <FILE1>;			
			chomp($f1_line) if ($f1_line);
		} else {
			$f2_line = <FILE2>;
			chomp($f2_line) if ($f2_line);			
		}
		
	}
}

sub join_union {
	my $has_header = $_[3];
	# Replacement for linux (doesn't work with system):  
	# 	join -t $'\t' -a 1 -a 2 <(sort file1) <(sort file2) > file3
	# 	sort -g file3 
	# But... this assumes tab diliniation and numerical order on join (1st) column
	#  except for the first row which are joined without test
	
	open FILE1, $_[0] or die "Could not open $_[0] for read";
	open FILE2, $_[1] or die "Could not open $_[1] for read";
	open OUT, '>', $_[2] or die "Could not open $_[2] for write";
	my $f1_line = <FILE1>;
	my $f2_line = <FILE2>;
	chomp($f1_line, $f2_line);
	my @f1_elements = split(/\t/, $f1_line);
	my @f2_elements = split(/\t/, $f2_line);
	my $f1_n = $#f1_elements;	# not counting joining [0] column
	my $f2_n = $#f2_elements;
	my @f1_dummy = ('---') x $f1_n;
	my @f2_dummy = ('---') x $f2_n;
	if ($has_header){
		print OUT join ("\t", @f1_elements, @f2_elements[1 .. $f2_n]), "\n";	#1st row headers
		$f1_line = <FILE1>;
		chomp($f1_line) if ($f1_line);
		$f2_line = <FILE2>;
		chomp($f2_line) if ($f2_line);
	}
	while($f1_line or $f2_line){
		if ($f1_line and $f2_line){
			@f1_elements = split(/\t/, $f1_line);
			@f2_elements = split(/\t/, $f2_line);
			if ($f1_elements[0] == $f1_elements[0]){
				print OUT join ("\t", @f1_elements, @f2_elements[1 .. $f2_n]), "\n";
				$f1_line = <FILE1>;
				chomp($f1_line) if ($f1_line);
				$f2_line = <FILE2>;			
				chomp($f2_line) if ($f2_line);
			} elsif ($f1_elements[0] < $f1_elements[0]){
				print OUT join ("\t", @f1_elements, @f2_dummy), "\n";
				$f1_line = <FILE1>;			
				chomp($f1_line) if ($f1_line);
			} else {
				print OUT join ("\t", $f2_elements[0], @f1_dummy, @f2_elements[1 .. $f2_n]), "\n";
				$f2_line = <FILE2>;
				chomp($f2_line) if ($f2_line);			
			}
		} elsif ($f1_line) {	#eof file2
			@f1_elements = split(/\t/, $f1_line);
			print OUT join ("\t", @f1_elements, @f2_dummy), "\n";
			$f1_line = <FILE1>;
			chomp($f1_line) if ($f1_line);
		} else {				#eof file1
			@f2_elements = split(/\t/, $f2_line);
			print OUT join ("\t", $f2_elements[0], @f1_dummy, @f2_elements[1 .. $f2_n]), "\n";
			$f2_line = <FILE1>;
			chomp($f2_line) if ($f2_line);		
		}
	}
}


