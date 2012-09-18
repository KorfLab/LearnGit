#!/usr/bin/perl
use warnings;
use strict;

#Fill table with 4 rows, numbers 1-16
#Prints the table while adding and after adding rows

my @row;
my @array;

print "Filling in array\n";   

for(my $i=0;$i<16;$i+=4){
    @row=($i+1,$i+2,$i+3,$i+4);  #fill row with value
    push @array, \@row;  #push row onto array
    my $line=join "\t" , @row; #join the row - just for printing on next line
    print $line . "\n";  #print the row
}

#Print the rows from $array_ref
print "After filled in:\n";

for(my $row_number=0;$row_number<4;$row_number++){
    @row= @{$array[$row_number]};
    my $line=join "\t", @row;
    print $line . "\n";
}