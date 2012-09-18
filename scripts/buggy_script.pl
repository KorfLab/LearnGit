#!/usr/bin/perl
use warnings;
use strict;

#Fill table with numbers 1-16
#Calculates mean, population variance and std deviation of array

my @row;
my @array;

print "Filling array with 1..16\n";   

for(my $i=0;$i<16;$i+=4){
    push @row , ($i+1,$i+2,$i+3,$i+4);  #fill row with value
    push @array, @row;  #push row onto array
}

#Calculate Mean and Population variance and std deviation in One Pass
my $sum = 0;
my $square_sum = 0;
my $n = scalar @array;

for(my $i=0; $i< scalar @array; $i++){
    $sum += $array[$i];
    $square_sum+=$array[$i] * $array[$i];
}

my $mean = $sum/$n;
my $variance = $square_sum/$n - $mean * $mean;
my $sd = sqrt($variance);

(abs($mean-8.5) < 0.0001) ? printf "Mean:\t%f\n",$mean :
    printf "Incorrect Mean:\t%f\n", $mean;

(abs($variance - 21.25) < 0.0001) ? printf "Population Variance:\t%f\n", $variance :
    printf "Incorrect Variance:\t%f\n",$variance;

(abs($sd - 4.609772)<0.0001) ? printf "Population Std. Deviation:\t%f\n", $sd :
    printf "Incorrect Std Dev:\t%f\n",$sd;
