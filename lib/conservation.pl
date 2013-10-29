#!/usr/bin/perl -w 



($afa_file) = @ARGV;


#########
# afa 2msf

########
# restrict
$cmd = "$restrict  $input_fpath $query_name > $restricted_msf";


open(FH, "<$cmd_template") || diehard("specs", "cno:$cmd_template,$!");
undef $/;
$prms_string = <FH>;
$/="\n";
close(FH);

$prms_string =~ s/align/$align_msf/;
$prms_string =~ s/query/$qry/;
$prms_string =~ s/outn/$specs_out/;
$prms_string =~ s/method/$method_string/;

if($structure){
    $prms_string =~ s/!pdbf/$pdbf/;
}
$cmd_specs = "$jobdir/specs_cmd";

open(FH, ">$cmd_specs") || diehard("specs","Cno:$cmd_specs, $!\n");
undef $/;
print FH $prms_string;
$/="\n";
	
close(FH);

my $stdout = "$jobdir/specs.stdout";
my $stderr = "$jobdir/specs.stderr";

$cmd = "$specs $cmd_specs 1>$stdout 2>$stderr";
system($cmd);
########################################
#generate picture
    
`awk \'\$1 \!= \"\%\" \{print \$5 \"  \" \$3 \"   \" \$1\}\' $score_f > $score_eachmethod`; 
	

$cmd = "java -jar $seqReport $score_eachmethod $jobdir/$input_fn";

my @resi_cnt = split(/\s/,`wc -l $score_eachmethod`);
my $num_resi = $resi_cnt[0];
my $range = 400;
($errmsg, my $png_ref) = getPngFile($num_resi,$png_f,$score_eachmethod, $range);

($errmsg eq "") || diehard("specs", $errmsg);
    
############################################################################
#generate excel spreadsheet
$cmd = "$specs2excel $jobdir/specs_out.score $jobdir/$input_fn";

system($cmd) && diehard("specs", "Error running $cmd:$!");

$excel_out .= ".xls";

###########################################################################
#zip the how directory
my $dirzipfile = "$jobdir.zip";
$cmd = "$zip -r $dirzipfile $jobdir > /dev/null";
system($cmd) && diehard("SPECS", "Error running $cmd");

($errmsg, my $html) = getHtmlBody($score_f, $excel_out, $zipfile, $dirzipfile,$png_ref);

print $html;

sub getHtmlBody(@){
    my $html="";
    my @downloading = @_;
    
    $html .= "<tr><td>";
    $html .= "<table width = \"720\" border=\"0\" align=\"center\" cellpadding=\"0\" cellspacing=\"0\">\n";
    $html .= "<tr><td> Downloads:\n";
    $html .= "<ul>\n";
   
    foreach my $elem(@downloading){
	
	$html .= "<tr>\n";
	
	$html .= "<td>\n";
	if(ref($elem) eq "ARRAY"){
	    foreach my $png(@{$elem}){
		my @tmp = split(/\./, $png);
		(my $from, my $to) = split(/_/, $tmp[$#tmp-1]);
		$html .= "<br>\nFrom $from to $to<br>\n";
		$html .= "<img src=\"./showImg.cgi?params=$png\" align=\"center\">\n";
		
	    }
	    $html .= "</td></tr>\n";
	}
	else{
	    if($elem =~ /(xls)$|(score)$|(zip)$/){
		
		if($elem =~ /(xls)$/){
            
		    $html .= "<li><a href=\"http://epsf.bmad.bii.a-star.edu.sg/cgi-bin/".
			"struct_server/download.cgi?ID=$elem\">output in excel format</a></li>\n";
		}
		if($elem =~ /(score)$/){
		    $html .= "<li><a href=\"http://epsf.bmad.bii.a-star.edu.sg/cgi-bin/struct_server/".
			"download.cgi?ID=$elem\">download score file</a></li>\n";
		}
		if($elem =~ /(zip)$/ && -e $elem){ #when user upload a sequence file $elem file not exist
		   
		    
		    if($elem =~ /(\d+\.zip)$/){
		    
			$html .= "<li><a href=\"http://epsf.bmad.bii.a-star.edu.sg/cgi-bin/".
			    "struct_server/download.cgi?ID=$elem\">directory</a></li>\n";
		    }
		    elsif($elem =~ /\.afa/){
			$html .= "<li><a href=\"http://epsf.bmad.bii.a-star.edu.sg/cgi-bin/".
			    "struct_server/download.cgi?ID=$elem\">the alignements(sorted by species)</a></li>\n";
		    }
		    elsif($elem =~ /\.pse/){
			$html .= "<li><a href=\"http://epsf.bmad.bii.a-star.edu.sg/cgi-bin/".
			    "struct_server/download.cgi?ID=$elem\">pymol session file</a></li>\n";
		    }
		    elsif($elem =~ /\_specs.py/){
			$html .= "<li><a href=\"http://epsf.bmad.bii.a-star.edu.sg/cgi-bin/".
			    "struct_server/download.cgi?ID=$elem\">chimera session file for specialization</a></li>\n";
		    }
		    elsif($elem =~ /\_cons.py/){
			$html .= "<li><a href=\"http://epsf.bmad.bii.a-star.edu.sg/cgi-bin/".
			    "struct_server/download.cgi?ID=$elem\">chimera session file for conservation</a></li>\n";
		    }

		}
	    }
	    else{
		$html .= "<tr><td>\n$elem";
	    }
	    
	}
	

    }
    $html .= "</td>\n</tr>\n";
    $html .= "</table>";
    
    return("",$html);
}

