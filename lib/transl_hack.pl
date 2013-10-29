#!/usr/bin/perl -w

@files = split "\n", `ls Best_MSA/*afa`;

$ctr = 0;

foreach $file (@files){


    print $file, "\n";
    $new_file = $file;
    $new_file =~ s/afa/msf/;
    $new_file =~ s/MSA/MSF/;
    print $new_file, "\n";

    $cmd = "/home/ivanam/perlscr/translation/afa2msf.pl ".
	" $file > $new_file ";
    print $cmd, "\n";
    `$cmd`;

    $ctr++;

    #last if ($ctr==3);

}
