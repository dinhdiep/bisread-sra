#!/usr/bin/perl

use strict;

my $prefix = $ARGV[0];
open(OUT1, ">$prefix.R1.fq") || die("Error writing.\n");
open(OUT2, ">$prefix.R2.fq") || die("Error writing.\n");

my $last_line = "NA";

while(my $line = <STDIN>){
    chomp($line);
    my @tmp1 = split "\t", $last_line;
    my @tmp2 = split "\t", $line;
    my ($substr1, $substr2) = split " ", $tmp2[0];
    $tmp2[0] = $substr1;
    ($substr1, $substr2) = split " ", $tmp1[0];
    $tmp1[0] = $substr1;
    if($tmp1[0] eq $tmp2[0]){
        print OUT1 "$tmp1[0]\n$tmp1[1]\n$tmp1[2]\n$tmp1[3]\n";
        print OUT2 "$tmp2[0]\n$tmp2[1]\n$tmp2[2]\n$tmp2[3]\n";
    }
    $last_line = $line;
}
close(OUT1);
close(OUT2);
