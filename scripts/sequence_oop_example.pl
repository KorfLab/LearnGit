#!/usr/bin/perl
# 
# Name:        <script_name.pl>
# Description: <one line description>
# Author:      <your name>
# License:     <BSD, MIT, GPL, CC, PD, choose one>

use strict;
use warnings FATAL => 'all';
use Getopt::Std;

################
# command line #
################

die "
usage: <program_name> [options] <arguments...>
options:
  -h
  -p <parameter>
" unless @ARGV; # adjust as necessary

my %opt;
getopts('hp:', \%opt);

my ($var1, $var2) = @ARGV;

# global variables and other initializations

#############
# main loop #
#############

############
# clean up #
############

###############
# subroutines #
###############

sub foo {}

BEGIN {} # remove if you don't need
END {} # remove if you don't need

=head1 Documentation

Include extra documentation here if you like

=cut


__END__

free form comments and thoughts
