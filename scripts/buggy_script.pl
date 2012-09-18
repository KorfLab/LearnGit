#!/usr/bin/perl
use warnings;
use strict;

#This script has one bug.  Fix the bug and make sure that you get the correct
#answer.   Commit the fix.   Then create a new bug.   Commit that bug and Push
#it to Github

print "Usage: $0 \n This script is pretty much useless" unless @ARGV==0;


my $array_ref=fillArray();  #get reference to 4x4 array with numbers 1-16
                            #Also prints the rows pushed onto the array

#Print the rows from $array_ref
print "After filled in:\n";

for(my $row=0;$row<4;$row++){
    my @row= @{$array_ref->[$row]};
    my $line=join "\t", @row;
    print $line . "\n";
}


#creates an 4x4 array with numbers from 1 to 16.
#returns an array-reference to the array
sub fillArray{
    my @array;
    my @row;
    
    print "Filling in array\n";
    
    for(my $i=0;$i<16;$i+=4){
        @row=($i+1,$i+2,$i+3,$i+4);  #fill row with value
        
        my $line=join "\t" , @row; #join the row - just for printing on next line
        print $line . "\n";  #print the row
        
        push @array, \@row;  #push row onto array
    }
    
    print "\n";
    
    return \@array;  #return reference to @array
}