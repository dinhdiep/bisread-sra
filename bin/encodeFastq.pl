#!/usr/bin/perl -w

use strict;


sub encodeFastq {
    my $fileName = shift;
    my $outFileName = shift;
    my $trimMode = shift;
    my $fivep = 0;
    my $threep = 0;

    $fivep = 27 if($trimMode eq 'BSPP');
    $fivep = 5 if($trimMode eq 'WGBS' or $trimMode eq 'RRBS');

    if($fileName =~ /\.gz$/){
        open(FQ, "gunzip -c $fileName |") || die("Can't open pipe to $fileName!");
    }else{
        open(FQ, "$fileName") || die ("Error in opening file $fileName!");
    }
    open(FQ_OUT, ">$outFileName") || die ("Error in opening file $outFileName!");
    while(my $line1 = <FQ>){
        chomp($line1);
        my @tmp = split /[\s+\t]/, $line1;
        # first 6 fields separated by : are cluster ID and by _ are UMI
        $line1 = $tmp[0];
        my $line2 = <FQ>;
        $line2 =~ tr/\./N/;
        chomp($line2);
        my $line3 = <FQ>;
        chomp($line3);
        my $line4 = <FQ>;
        chomp($line4);
        next if(!$line1 || !$line2 || !$line3 || !$line4);
        if($trimMode eq 'RRBS-ND' and length($line2) > 2){
            if($line2 =~ /^CAA/ or $line2 =~ /^CGA/){
                $line2 = substr($line2, 2, length($line2) - 2);
                $line4 = substr($line4, 2, length($line4) - 2);
            }
        }
        if($threep || $fivep){
            my $start = $fivep;
            my $total = length($line2) - $threep - $fivep;
            $line2 = substr($line2, $start, $total);
            $line4 = substr($line4, $start, $total);
        }
        my $guess = guess_strand($line2);
        my $prefix = $line1;
                
        my $s = $line2;
        if($guess eq "R"){
            $line2 =~ s/G/a/g;
        }else{
            $line2 =~ s/C/t/g;
        }
        $line1 = $prefix . "\\" . $s . "|" . $guess;
        print FQ_OUT "$line1\n$line2\n$line3\n$line4\n";
    }
    close(FQ);
    close(FQ_OUT);
}

sub guess_strand{
    my $seq = shift;
    my %baseCounts;
    $baseCounts{'A'}=0.001;
    $baseCounts{'T'}=0.001;
    $baseCounts{'G'}=0.001;
    $baseCounts{'C'}=0.001;
    while(my $base = chop($seq)){
        $baseCounts{$base}++;
    }
    if($baseCounts{'T'}/$baseCounts{'C'} > $baseCounts{'A'}/$baseCounts{'G'}) {
        return "F";
    }else{
        return "R";
    }
}

encodeFastq($ARGV[0], $ARGV[1], $ARGV[2]);

