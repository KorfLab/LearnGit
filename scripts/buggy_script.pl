#!/usr/bin/perl
# 
# buggy_script.pl
# The 'fix-the-typo-and-add-a-new-typo-challenge'!
#
# Remember...confirm you have added a bug by running 'perl -c buggy_script.pl' 

use warnings;
use strict;

# Fill array with 5–20 random integers (with values ranging from 1–10)
# Calculates mean, population variance and std deviation of numbers in array
my @numbers;

my $random = int(rand(15))+5;
print "Generating $random random numbers: ";   

for(my $i = 1; $i <= $random; $i++){
	my $rand_int = int(rand(10)) + 1;
	print "$rand_int ";
    push @numbers, $rand_int;  #push row onto array
}
print "\n\n";


# Calculate Mean and population variance and population std deviation in One Pass
my $sum = 0;
my $square_sum = 0;
my $n = scalar @numbers;

# loop over each number
for(my $i = 0; $i < scalar @numbers; $i++){
    $sum += $numbers[$i];
    $square_sum += ($numbers[$i] * $numbers[$i]);
}

#  calculate basic stats and round down to 2 d.p.
my $mean     = sprintf("%.2f", $sum / $n);
my $variance = sprintf("%.2f", $square_sum / $n - $mean * $mean);
my $sd       = sprintf("%.2f", sqrt($variance));


# print final output
print "n = $n\n";
print "Sum = $sum\n";
print "Mean = $mean\n";
print "Variance = $variance\n";
print "Standard deviation = $sd\n\n\n";

exit;

__END__
Some people asked if variance could actually be calculated in a single pass.  Yes, see the links below for one pass algorithms:
http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
http://www.strchr.com/standard_deviation_in_one_pass