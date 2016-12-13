#!/usr/bin/perl
use strict;
use warnings;

# Ad hoc program used to merge recomb map windows avoiding any 0 values.
# Any window(s) with 0 value after initial merging is (are) then merged
# with neighboring windows (both sides) recursively until all have non-0 values.

# Checks for gaps in map and leaves as is (ie, doesn't merge across missing region).

my $help = q`
usage: perl merge_windows.pl <min_bp> <infile> <outfile> <min_recomb_mean_percent>

min_recomb_percent is percent of mean recomb value to be imputed as minimum for 0-value
windows that cannot be merged with neighbors to acheive non-0 (default = 5)

Expeted infile format (no header):
1	1	10000	+	0
1	10001	20000	+	0
1	20001	30000	+	0
1	30001	40000	+	0
1	40001	50000	+	0
1	50001	60000	+	0
1	60001	70000	+	0.015503876
1	70001	80000	+	0
1	80001	90000	+	0

`;

die $help if ($#ARGV != 2 and $#ARGV != 3);
my $min_bp = $ARGV[0];
my $infile = $ARGV[1];
my $outfile =  $ARGV[2];
my $min_recomb_mean_percent = $ARGV[3];
$min_recomb_mean_percent = $min_recomb_mean_percent || 5;

open INFILE, $infile or die "Could not open $infile for read";
open OUTFILE, '>', $outfile or die "Could not open $outfile for write";

# Read lines summing rates to achieve $min_bp
my @data;	#array of array refs [0] chr, [1] start, [2] stop, [3] orientation, [4] rate

my $current_chr = 'not this';
my $prev_end = 0;
my $i = -1;
my $running_length = 0;
my $rate_sum = 0;
my $length_sum = 0;
my $n_original_windows = 0;
while (<INFILE>){
	next if ($_ eq $/);
	chomp;
	my ($chr, $start, $end, $orient, $rate) = split;
	$rate_sum = $rate_sum += $rate;
	$length_sum = $length_sum += $end - $start;
	$n_original_windows++;
	if ($chr ne $current_chr or $start > ($prev_end + 1) or $running_length > $min_bp){ # new window
		$i++;
		my $new_data = [$chr, $start, $end, $orient, $rate];
		push @data, $new_data;
		$running_length = $end - $start;
	} else { # add to current window
		my $current_data = $data[$i];
		$current_data->[2] = $end;
		$current_data->[4] += $rate;
		die "Inconsistant orientation" if ($current_data->[3] ne $orient);
		$running_length += $end - $start;
	}
	$prev_end = $end;
	$current_chr = $chr;
}	# right-most windows on chromosomes will be < min_bp, we'll fix that below
my $n_windows = $i + 1;
my $min_impute_rate = ($rate_sum / $length_sum) * $min_recomb_mean_percent / 100;
print "Input map length: $length_sum; morgans: $rate_sum; windows: $n_original_windows\n";
print "Windows after initial min_bp trimming: $n_windows\n";


open DEBUG, '>', 'debug.txt' or die "Could not open debug.txt for write";
foreach my $data_line (@data){
	print DEBUG join ("\t", $data_line->[0], $data_line->[1], $data_line->[2], $data_line->[3], $data_line->[4]), "\n";
}

my $n_merges = 0;
my $n_zero_merges = 0;
my $n_left_merges = 0;
my $n_right_merges = 0;
my $n_min_imputed = 0;
$i = 0;
while ($i < $n_windows){
	my $current_data = $data[$i];
	if ($current_data->[4] > 0 and ($current_data->[2] - $current_data->[1]) >= $min_bp){ # non-0 rate and window >= min_bp
		$i++;
	} else {
		$n_merges++;
		if ($current_data->[4] == 0){
			$n_zero_merges++;
		}
		# see if we can merge left and/or right
		my ($merge_left, $merge_right);
		if ($i > 0){
			my $left_data = $data[$i - 1];
			if ($current_data->[0] eq $left_data->[0] and $current_data->[1] == $left_data->[2] + 1){
				$merge_left = 1;
			}
		}
		if ($i < $n_windows - 1){
			my $right_data = $data[$i + 1];
			if ($current_data->[0] eq $right_data->[0] and $current_data->[2] + 1 == $right_data->[1]){
				$merge_right = 1;
			}
		}
		if ($merge_left and $merge_right) { # pick left or right randomly
			if (rand(2) < 1) {
				$merge_left = 0;
			} else {
				$merge_right = 0;
			}
		}
		
		if ($merge_left) {
			my $left_data = $data[$i - 1];
			$left_data->[2] = $current_data->[2];	# end
			$left_data->[4] += $current_data->[4];	# rate
			splice @data, $i, 1;
			$n_windows--;
			$n_left_merges++;
		} elsif ($merge_right) {
			my $right_data = $data[$i + 1];
			$right_data->[1] = $current_data->[1];	# start
			$right_data->[4] += $current_data->[4];	# rate
			splice @data, $i, 1;
			$n_windows--;
			$n_right_merges++;			
		} else {
			$i++;
			my $length = $current_data->[2] - $current_data->[1];
			$current_data->[4] = $length * $min_impute_rate;
			$n_min_imputed++;
		}		
	}
}

# Still need minimum value
my $n_additional_merges = $n_merges - $n_zero_merges;
print "Merges for zero rate: $n_zero_merges\n";
print "Additional merges for short ends: $n_additional_merges\n";
print "Left merges: $n_left_merges; right merges: $n_right_merges\n";
print "Windows using min_impute_rate: $n_min_imputed\n";

$rate_sum = 0;
$length_sum = 0;
foreach my $data_line (@data){
	$rate_sum += $data_line->[4];
	$length_sum += $data_line->[2] - $data_line->[1];
	print OUTFILE join ("\t", $data_line->[0], $data_line->[1], $data_line->[2], $data_line->[3], $data_line->[4]), "\n";
}
print "Final map length: $length_sum; morgans: $rate_sum; windows: $n_windows\n";

