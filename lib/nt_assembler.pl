#! /usr/bin/perl -w
use DBI;
#use id2species;

sub formatted_sequence ( @);
sub set_id2species();


# da li da hardcodiramo path do databaze?

( @ARGV >= 2 )  || die "Usage: $0  <db_name> <ensembl gene id> <outfile>.\n"; 

($dbname, $gene_id,  $outname) =   @ARGV; 



($dbh = DBI->connect("DBI:mysql:$dbname", "marioot", "tooiram")) ||
    die "Cannot connect to $dbname.\n";

%id2species = ();
set_id2species();

# get all exons

$species = $id2species{substr $gene_id,0,6};
$query  = "SELECT dna_seq, left_flank, right_flank, exon_key, start_in_gene, is_canonical, ensembl_exon_id, source ";
$query .= "FROM exon_$species WHERE ensembl_gene_id = '$gene_id'";

%exon         = ();
%is_canonical = ();
%start        = ();
@header       = ();
$sth= prepare_and_exec( $query);

$all_strand = 0;
while ( ($dna_seq, $left_flank, $right_flank,  $key,  $start, $is_canon, $ensembl_id, $source) = $sth->fetchrow_array) {

    $dna_seq || next;

    $len      = length ($left_flank);
    $dna_str  = '';
    $dna_str .= lc substr $left_flank, $len-5;
    #$dna_str .= 'Z';
    $dna_str .= uc $dna_seq;
    #$dna_str .= 'Z';
    $dna_str .= lc substr $right_flank,0, 5;
    #$dna_str .= 'z';
    $exon{$key} =   $dna_str;
    $is_canonical{$key} = $is_canon;
    $start{$key}        = $start;
    $end{$key}          = $end;
    $source{$key}       = $source;

    if ($ensembl_id ne "anon") {
	$header = ">$ensembl_id\n "

    } else {

#	@aux = split "_", $key;
#	$header = ">".$aux[1];
#	if ($aux[2] == 2) {
#	    $header .= "_sw#";
#	} elsif ($aux[2] == 1) {
#	    $header .= "_known";
#	} else {
#	    $header .= "_pred";
#	}
#	$header .= "\n";

        $header = ">".$source."_".$start."\n"
    }

    $header{$key}      =  $header;    
}

@keys_sorted    = sort { $start{$a} <=>  $start{$b} } keys %start;


@exon_sorted         = ();
@is_canonical_sorted = ();

foreach $key (@keys_sorted) {
    push @exon_sorted,  $exon{$key};
    push @is_canonical_sorted, $is_canonical{$key};
    push @header_sorted, $header{$key};
}

#if ($all_strand < 0) {
#    @exon         = reverse(@exon_sorted);
#    @is_canonical = reverse(@is_canonical_sorted);
#
#} else {
#
@exon         =  @exon_sorted;
@is_canonical =  @is_canonical_sorted;
@header       =  @header_sorted;
#}



@gaps = ();
for $exon (@exon) {

    $gaps = "-"x length($exon);
    push @gaps, $gaps;
}

@exon_seqs = ();
foreach $i (0 .. $#exon) {

    @seq     = @gaps;
    $seq[$i] = $exon[$i];
    
    push @exon_seq, join "", @seq;
}


@canon = @gaps;

foreach $i (0 .. $#exon) {
    if ( $is_canonical[$i]) {
	$canon[$i] = $exon[$i];
    }
}

$canonical = join "", @canon;

#####################################################################
# finally, output the whole thing;



open (OUT, ">$outname") || die "Cno $outname: $!.\n";

print  OUT ">canonical\n";
print  OUT formatted_sequence($canonical);
print  OUT  "\n";


foreach $i (0 .. $#exon) {

    print  OUT  $header[$i];
    print  OUT  formatted_sequence($exon_seq[$i]);
    print  OUT "\n";
}

close OUT;






#untie the databse
$rc=$dbh->disconnect;

############################################################################3


sub prepare_and_exec
{
        my $query=shift;
        my $sth;
        my $rv;
        my $rc;
        $sth=$dbh->prepare($query);
        if (!$sth)
        {
               print "Can't prepare $query\n";
               $rc=$dbh->disconnect;
               die;
        }
        #$sth->trace(2);
        $rv=$sth->execute;
        if (!$rv)
        {
               print "Can't execute $query\n";
               $rc=$dbh->disconnect;
               die;;
        }
        return $sth;
}



sub one_line_query {


    my $sth;

    print "\n*******************************************\n";
    print "QUERY:         $query\n";
    print "RESULTS:             \n";
    $sth= prepare_and_exec( $query);
    while(@row = $sth->fetchrow_array)
    {
        foreach $row (@row)
        {
            print "\t\t $row\n";
        }
    }
}
#############################################################
#############################################################
#############################################################
sub formatted_sequence ( @) {

    my $ctr, 
    my $sequence = $_[0];
    ( defined $sequence) || die "Error: Undefined sequence in formatted_sequence()";
    $ctr = 50;
    while ( $ctr < length $sequence ) { 
	substr ($sequence, $ctr, 0) = "\n";
	$ctr += 51; 
    } 
    
    return $sequence; 
} 

#############################################################
#############################################################
#############################################################
sub set_id2species() {

    $id2species{'ENSAME'} = 'ailuropoda_melanoleuca';
    $id2species{'ENSAPL'} = 'anas_platyrhynchos';
    $id2species{'ENSACA'} = 'anolis_carolinensis';
    $id2species{'ENSBTA'} = 'bos_taurus';
    $id2species{'ENSCJA'} = 'callithrix_jacchus';
    $id2species{'ENSCAF'} = 'canis_familiaris';
    $id2species{'ENSCPO'} = 'cavia_porcellus';
    $id2species{'ENSCHO'} = 'choloepus_hoffmanni';
    $id2species{'ENSDAR'} = 'danio_rerio';
    $id2species{'ENSDNO'} = 'dasypus_novemcinctus';
    $id2species{'ENSDOR'} = 'dipodomys_ordii';
    $id2species{'ENSETE'} = 'echinops_telfairi';
    $id2species{'ENSECA'} = 'equus_caballus';
    $id2species{'ENSEEU'} = 'erinaceus_europaeus';
    $id2species{'ENSFCA'} = 'felis_catus';
    $id2species{'ENSFAL'} = 'ficedula_albicollis';
    $id2species{'ENSGMO'} = 'gadus_morhua';
    $id2species{'ENSGAL'} = 'gallus_gallus';
    $id2species{'ENSGAC'} = 'gasterosteus_aculeatus';
    $id2species{'ENSGGO'} = 'gorilla_gorilla';
    $id2species{'ENSG00'} = 'homo_sapiens';
    $id2species{'ENSSTO'} = 'ictidomys_tridecemlineatus';
    $id2species{'ENSLAC'} = 'latimeria_chalumnae';
    $id2species{'ENSLAF'} = 'loxodonta_africana';
    $id2species{'ENSMMU'} = 'macaca_mulatta';
    $id2species{'ENSMEU'} = 'macropus_eugenii';
    $id2species{'ENSMGA'} = 'meleagris_gallopavo';
    $id2species{'ENSMIC'} = 'microcebus_murinus';
    $id2species{'ENSMOD'} = 'monodelphis_domestica';
    $id2species{'ENSMUS'} = 'mus_musculus';
    $id2species{'ENSMPU'} = 'mustela_putorius_furo';
    $id2species{'ENSMLU'} = 'myotis_lucifugus';
    $id2species{'ENSNLE'} = 'nomascus_leucogenys';
    $id2species{'ENSOPR'} = 'ochotona_princeps';
    $id2species{'ENSONI'} = 'oreochromis_niloticus';
    $id2species{'ENSOAN'} = 'ornithorhynchus_anatinus';
    $id2species{'ENSOCU'} = 'oryctolagus_cuniculus';
    $id2species{'ENSORL'} = 'oryzias_latipes';
    $id2species{'ENSOGA'} = 'otolemur_garnettii';
    $id2species{'ENSPTR'} = 'pan_troglodytes';
    $id2species{'ENSPSI'} = 'pelodiscus_sinensis';
    $id2species{'ENSPMA'} = 'petromyzon_marinus';
    $id2species{'ENSPPY'} = 'pongo_abelii';
    $id2species{'ENSPCA'} = 'procavia_capensis';
    $id2species{'ENSPVA'} = 'pteropus_vampyrus';
    $id2species{'ENSRNO'} = 'rattus_norvegicus';
    $id2species{'ENSSHA'} = 'sarcophilus_harrisii';
    $id2species{'ENSSAR'} = 'sorex_araneus';
    $id2species{'ENSSSC'} = 'sus_scrofa';
    $id2species{'ENSTGU'} = 'taeniopygia_guttata';
    $id2species{'ENSTRU'} = 'takifugu_rubripes';
    $id2species{'ENSTSY'} = 'tarsius_syrichta';
    $id2species{'ENSTNI'} = 'tetraodon_nigroviridis';
    $id2species{'ENSTBE'} = 'tupaia_belangeri';
    $id2species{'ENSTTR'} = 'tursiops_truncatus';
    $id2species{'ENSVPA'} = 'vicugna_pacos';
    $id2species{'ENSXET'} = 'xenopus_tropicalis';
    $id2species{'ENSXMA'} = 'xiphophorus_maculatus';
}


