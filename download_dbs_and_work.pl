#!/usr/bin/perl
use strict;
use warnings;
use Bio::Index::Fasta;
use Bio::SearchIO;
use lib "/home/remo/src/RS9";
use Biosub::Utils;
use Bio::Index::Fasta;
use DBI;
use DBD::mysql;
use Data::Dumper;
use Bio::SeqIO;
use Bio::Tools::CodonTable;
use Bio::Tools::SeqStats;

BEGIN {

  $ENV{BLASTPLUSDIR} = "/home/remo/src/ncbi-blast-2.2.24+/bin";
  $ENV{BLASTDBDIR} = "/media/LOCAL_DATA/BLASTDB";
  $ENV{UNIPROTDIR} = '/media/LOCAL_DATA/UNIPROT';
  $ENV{EMBLDIR} = '/media/LOCAL_DATA/EMBL';
  $ENV{PAL2NALDIR} = "/home/remo/src/pal2nal.v13";
  $ENV{ANNSRCDIR} = "/home/remo/Dropbox/SEQ/ANNOTATION_PIPELINE";

  chdir("$ENV{BLASTPLUSDIR}") or die "\nError in checking dir $ENV{BLASTPLUSDIR}: $!\n\n";
  chdir("$ENV{BLASTDBDIR}") or die "\nError in checking dir $ENV{BLASTDB}: $!\n\n";
  chdir("$ENV{UNIPROTDIR}") or die "\nError in checking dir $ENV{UNIPROTDIR}: $!\n\n";
  chdir("$ENV{EMBLDIR}") or die "\nError in checking dir $ENV{EMBLDIR}: $!\n\n";
  chdir("$ENV{PAL2NALDIR}") or die "\nError in checking dir $ENV{PAL2NALDIR}: $!\n\n";
  chdir("$ENV{ANNSRCDIR}") or die "\nError in checking dir $ENV{ANNSRCDIR}: $!\n\n";

#  $ENV{BLASTPLUSDIR} = "/Users/remo/src/ncbi-blast-2.2.24+/bin";
#  $ENV{BLASTDBDIR} = "/Volumes/EXTERNAL_DATA/BLASTDB";
#  $ENV{UNIPROTDIR} = '/Volumes/EXTERNAL_DATA/UNIPROT';
#  $ENV{EMBLDIR} = '/Volumes/EXTERNAL_DATA/EMBL';
#  $ENV{PAL2NALDIR} = "/Users/remo/src/pal2nal.v13";
#  $ENV{ANNSRCDIR} = "/Users/remo/Dropbox/SEQ/ANNOTATION_PIPELINE";

#  chdir("$ENV{BLASTPLUSDIR}") or die "\nError in checking dir $ENV{BLASTPLUSDIR}: $!\n\n";
#  chdir("$ENV{BLASTDBDIR}") or die "\nError in checking dir $ENV{BLASTDB}: $!\n\n";
#  chdir("$ENV{UNIPROTDIR}") or die "\nError in checking dir $ENV{UNIPROTDIR}: $!\n\n";
#  chdir("$ENV{EMBLDIR}") or die "\nError in checking dir $ENV{EMBLDIR}: $!\n\n";
#  chdir("$ENV{PAL2NALDIR}") or die "\nError in checking dir $ENV{PAL2NALDIR}: $!\n\n";
#  chdir("$ENV{ANNSRCDIR}") or die "\nError in checking dir $ENV{ANNSRCDIR}: $!\n\n";

}

# ACE AND ESTFILE MUST TO BE IN THE SAME DIRECTORY
my $ESTFILE = 'EYE.fa';
my $ESTDIR = "/home/remo/Dropbox/SEQ";
#my $ACEFILE = "clone.fa.cap.ace";
#my $ESTDIR = "/Users/remo/Dropbox/SEQ";

# DATABASE SPECIFIC SETTINGS
my $DBuniprot = 'uniprot_idmapping';
my $DBest = 'OCTOPUS_EST';
my $USR = 'mysql_dev';
my $PWD = 'riiGbs';
#my $PWD = 'dEvEl0pEr';
my $HOST;

my $DOWNLOAD_UNIPROT = 0;
my $CONCATENATE = 0;
my $CREATE_BLASTDB = 0;
my $CREATE_MYSQL_ID_MAPPING = 0;
my $CREATE_UNIPROT_FASTA_INDEX = 0;
my $CREATE_UNIPROT_SPECIES_FASTA = 0;
my $CREATE_UNIPROT_SPECIES_BLASTDB = 0;
my $RUN_UNIPROT_SPECIES_BLASTX = 0;
my $DOWNLOAD_EMBL = 0;
my $CREATE_EMBL_SPECIES_FASTA = 0;
my $CREATE_EMBL_SPECIES_FASTA_INDEX = 0;
my $CALCULATE_KA_KS = 1;
my $POPULATE_RATIO = 1;
my $POPULATE_EST = 0;

my %SPECIES = (
               'Homo sapiens'            =>   'homo_sapiens',
               'Drosophila melanogaster' =>   'drosophila_melanogaster',
               'Danio rerio'             =>   'danio_rerio',
               'Aplysia californica'     =>   'aplysia_californica',
               'Lottia gigantea'         =>   'lottia_gigantea',
               'Ciona intestinalis'      =>   'ciona_intestinalis',
               'Caenorhabditis elegans'  =>   'caenorhabditis_elegans'
              );

# FROM HERE WE WORK IN THE UNIPROT FOLDER
chdir("$ENV{UNIPROTDIR}") or die "Error in changing dir: $!\n";
system("cp $ESTDIR\/$ESTFILE $ENV{UNIPROTDIR}\/");
system("cp $ESTDIR\/$ESTFILE $ENV{EMBLDIR}\/");

if($DOWNLOAD_UNIPROT) {
  # Download swissprot
  system('wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz');
  system('gunzip -c uniprot_sprot.fasta.gz > uniprot_sprot.fasta');
  # Download trembl
  system('wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz');
  system('gunzip -c uniprot_trembl.fasta.gz > uniprot_trembl.fasta');
  # Download id_mapping
  system('wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz');
  system('gunzip -c idmapping_selected.tab.gz > idmapping_selected.tab');
}

if($CONCATENATE) {
  # Concatenate the files
  system('cat uniprot_trembl.fasta >> uniprot_sprot.fasta');
  system('mv uniprot_sprot.fasta uniprot_sprot_trembl.fasta');
  system('rm uniprot_trembl.fasta');
}

if($CREATE_BLASTDB) {
  # Create blast db
  system("$ENV{BLASTPLUSDIR}/makeblastdb -in uniprot_sprot_trembl.fasta -dbtype prot -out uniprot_sprot_trembl -title uniprot_sprot_trembl -logfile uniprot_sprot_trembl.log");
  system("mv uniprot_sprot_trembl.0* $ENV{BLASTDBDIR}/");
  system("mv uniprot_sprot_trembl.p* $ENV{BLASTDBDIR}/");
  system("mv uniprot_sprot_trembl.log $ENV{BLASTDBDIR}/");
}

if($CREATE_MYSQL_ID_MAPPING) {
  # Need to create a db named 'uniprot_idmapping' before
  # And need to have the program in the proper folder
  system("perl $ENV{ANNSRC}/parse_idmapping_and_populate_db.pl");
}

if($CREATE_UNIPROT_FASTA_INDEX) {
  my $inx = Bio::Index::Fasta->new(-filename => 'uniprot_sprot_trembl.index',
                                   -write_flag => 1);
  $inx->make_index('uniprot_sprot_trembl.fasta');
}

if($CREATE_UNIPROT_SPECIES_FASTA) {
  my $in = Bio::SeqIO->new(-file => 'uniprot_sprot_trembl.fasta',
                           -format => 'fasta');
  my $file = {};
  foreach my $sn(keys %SPECIES) {
    my $cn = $SPECIES{$sn};
    my $out = Bio::SeqIO->new(-file => ">uniprot_$cn\.fa",
                              -format => 'fasta');
    $file->{$cn} = $out;
  }
  while(my $seq = $in->next_seq) {
    if($seq->description =~ /OS=(\w+ \w+) /) {
      my $s = $1;
      if(exists $SPECIES{$s}) {
        $file->{$SPECIES{$s}}->write_seq($seq);
      }
    }
  }
}

if($CREATE_UNIPROT_SPECIES_BLASTDB) {
  foreach my $sn(keys %SPECIES) {
    my $cn = $SPECIES{$sn};
    # Create blast db
    system("$ENV{BLASTPLUSDIR}/makeblastdb -in uniprot_$cn\.fa -dbtype prot -out uniprot_$cn -title uniprot_$cn -logfile uniprot_$cn\.log");
    #system("mv uniprot_$cn\.0* $ENV{BLASTDBDIR}/"); only neede when you have big files which makeblastdb splits
    system("mv uniprot_$cn\.p* $ENV{BLASTDBDIR}/");
    system("mv uniprot_$cn\.log $ENV{BLASTDBDIR}/");
  }
}

if($RUN_UNIPROT_SPECIES_BLASTX) {
  foreach my $sn(keys %SPECIES) {
    my $cn = $SPECIES{$sn};
    my $dbname = "uniprot_$cn";
    my $out = run_blastplus('blastx',$dbname);
    my $tab = parse_blast_results($out);
  }
}

# FROM HERE WE START TO WORK IN THE EMBL DIRECTORY
chdir("$ENV{EMBLDIR}") or die "Error in changing dir: $!\n";

if($DOWNLOAD_EMBL) {
  # Download EMBL_CDS
  system('wget ftp://ftp.ebi.ac.uk/pub/databases/embl/cds/cds.fasta.gz');
  system('gunzip cds.fasta.gz');
}

if($CREATE_EMBL_SPECIES_FASTA) {
  my $in = Bio::SeqIO->new(-file => 'cds.fasta',
                           -format => 'fasta');
  my $file = {};
  foreach my $sn(keys %SPECIES) {
    my $cn = $SPECIES{$sn};
    my $out = Bio::SeqIO->new(-file => ">cds_$cn\.fa",
                              -format => 'fasta');
    $file->{$cn} = $out;
  }
  while(my $seq = $in->next_seq) {
    # EMBL_CDS:CAA00063 CAA00063.1 Homo sapiens
    if($seq->description =~ /\S+ (\w+ \w+) /) {
      my $s = $1;
      if(exists $SPECIES{$s}) {
        my $id = $seq->id;
        $id =~ s/EMBL_CDS://;
        $seq->id($id);
        $file->{$SPECIES{$s}}->write_seq($seq);
      }
    }
  }
}

if($CREATE_EMBL_SPECIES_FASTA_INDEX) {
  my @fasta;
  mkdir('species_fasta');
  foreach my $sn(keys %SPECIES) {
    my $cn = $SPECIES{$sn};
    my $fasta = "cds_$cn\.fa";
    system("mv $fasta species_fasta\/");
    push(@fasta,"species_fasta\/$fasta");
  }
  my $index = Bio::Index::Fasta->new(-filename => 'species_index.idx',
                                     -write_flag => 1);
  $index->make_index(@fasta);
}

if($CALCULATE_KA_KS) {
  my $DBHuniprot = connect_to_db($DBuniprot,$USR,$PWD,$HOST);
  my $index = 'species_index.idx';
  my $inx = Bio::Index::Fasta->new(-filename => $index);
  my $pre_cnt_string = get_cnt();

  my $ratiofile = "$ENV{UNIPROTDIR}\/$ESTFILE\.codeml";
  open(RATIOS,">$ratiofile");
  print RATIOS join("\t",qw(organism clone hit cds dn_ds dn ds))."\n";

  foreach my $sn(keys %SPECIES) {
    my $cn = $SPECIES{$sn};
    # Can be dangerous if you change names of blast output and tables
    # EST.blastx_uniprot_aplysia_californica.table
    my $parsed_blastx = "$ENV{UNIPROTDIR}\/$ESTFILE\.blastx\_uniprot\_$cn\.table";
    my $sigla = lc("$cn"); $sigla =~ s/^(\w)\w+\_(\w)\w+$/$1$2/;
    $parsed_blastx =~ s/\.fasta//;
    $parsed_blastx =~ s/\.fa//;
    my $table = get_table_in_href($parsed_blastx);

    my $seqio = Bio::SeqIO->new(-file => $ESTFILE,
                                -format => 'fasta');
    my $dir = "$ESTFILE\_$sigla\_PAL2NAL";
    $dir =~ s/\.fa//;
    mkdir($dir);

    while(my $pre_est_seq = $seqio->next_seq) {
      if(exists $table->{$pre_est_seq->id}) {

        my $res = $table->{$pre_est_seq->id};
        my $embl_cds_id = get_embl_cds($DBHuniprot,$res->{hit});
        next unless $embl_cds_id;

        my $embl_cds = $inx->fetch($embl_cds_id);
        # We need to check why we need to do this...
        next unless $embl_cds;
        # A check to define the codon table http://doc.bioperl.org/releases/bioperl-1.6.1/Bio/Tools/CodonTable.html
        my($embl_protein,$codontable_id) = get_embl_translation($embl_cds);
        my($est_seq,$est_protein) = get_translation($pre_est_seq,$res->{strand},$res->{frame},$codontable_id);

        # At this point we have:
        # CDS seq objects         $est_seq         $embl_cds
        # Protein seq objects     $est_protein     $embl_protein

        my $prot = "$dir/".$est_seq->id.'_prot.fa';
        my $nucl = "$dir/".$est_seq->id.'_nucl.fa';
        my $aln = "$dir/".$est_seq->id.'_prot.aln';
        my $dnd = "$dir/".$est_seq->id.'_prot.dnd';
        my $pal2nal = "$dir/".$est_seq->id.'.pal2nal';
        my $codeml = "$dir/".$est_seq->id.'.codeml';

        my $tree = "$dir/".$est_seq->id.'.tree';
        open(TREE,">$tree");
        print TREE '('.$est_seq->id.', '.$embl_cds->id.');';
        close(TREE);
 
        my $cnt = "$dir/".$est_seq->id.'.cnt';
        my $cnt_string = "$pre_cnt_string";
        $cnt_string =~ s/seqfile \= test.codon/seqfile \= $pal2nal/;
        $cnt_string =~ s/treefile \= test.tree/treefile \= $tree/;
        $cnt_string =~ s/outfile \= test.codeml/outfile \= $codeml/;
        #$cnt_string =~ s/icode \= 0/icode \= 1/ if $codontable_id == 2;
        open(CNT,">$cnt");
        print CNT $cnt_string;
        close(CNT);

        my $out_nt = Bio::SeqIO->new(-file => ">$nucl",
                                     -format => 'fasta');    
        $out_nt->write_seq($est_seq) if $res->{strand} == 1;
        $out_nt->write_seq($est_seq->revcom) if $res->{strand} == -1;
        $out_nt->write_seq($embl_cds);    

        my $out_pr = Bio::SeqIO->new(-file => ">$prot",
                                     -format => 'fasta');
        $out_pr->write_seq($est_protein);
        $out_pr->write_seq($embl_protein);
        system("muscle -quiet -clwstrict -in $prot -out $aln");
#        system("perl $ENV{PAL2NALDIR}\/pal2nal.pl $aln $nucl -output paml -nomismatch -nogap > $pal2nal");
        system("perl $ENV{PAL2NALDIR}\/pal2nal.pl $aln $nucl -output paml -nogap > $pal2nal");
        system("codeml $cnt") unless -z $pal2nal;
        my $ratios = get_codeml_res($codeml);
        print RATIOS join("\t",$sigla,$est_seq->id,$res->{hit},$embl_cds->id,@$ratios)."\n";
      }
    }
  }
}

if($POPULATE_RATIO) {
  my $ratiofile = "$ENV{UNIPROTDIR}\/$ESTFILE\.codeml";
  my $ratiotable = "$ESTFILE\.codeml";
  $ratiotable =~ s/.fa//;
  $ratiotable =~ s/\./\_/g;
  my $DBHest = connect_to_db($DBest,$USR,$PWD,$HOST);
  my $create = "CREATE TABLE $ratiotable (
                  organism VARCHAR(3),
                  clone VARCHAR(256),
                  hit VARCHAR(256),
                  cds VARCHAR(256),
                  dnds DOUBLE,
                  dn DOUBLE,
                  ds DOUBLE,
                  KEY organism(organism),
                  KEY clone(clone),
                  KEY hit(hit),
                  KEY cds(cds),
                  KEY dnds(dnds),
                  KEY dn(dn),
                  KEY ds(ds)
                )";
  $DBHest->do($create);
  open(IN,$ratiofile);
  my $header = <IN>;
  while(my $row = <IN>) {
    chomp($row);
    $row =~ s/-nan//;
    my @field = split(/\t/,$row);
    my $dnds = $field[4];
    my $dn = $field[5];
    my $ds = $field[6];
    my $string = "INSERT INTO $ratiotable SET organism = ".$DBHest->quote($field[0]).
                 ',clone = '.$DBHest->quote($field[1]).',hit = '.$DBHest->quote($field[2]).
                 ',cds = '.$DBHest->quote($field[3]);
    $string .= ",dnds = $dnds" if $dnds;
    $string .= ",dn = $dn" if $dn;
    $string .= ",ds = $ds" if $ds;
    my $sth = $DBHest->prepare($string);
    $sth->execute;
  }
}

if($POPULATE_EST) {
  my $DBHest = connect_to_db($DBest,$USR,$PWD,$HOST);
  my $seqtable = "$ESTFILE";
  $seqtable =~ s/\.fasta//;
  $seqtable =~ s/\.fa//;
  $seqtable .= '_seq';
  my $create = "CREATE TABLE $seqtable (
                  clone VARCHAR(256) PRIMARY KEY,
                  length INT(6),
                  A INT(5), T INT(5), C INT(5), G INT(5), N INT(5),
                  seq TEXT,
                  KEY length(length),
                  KEY A(A), KEY T(T), KEY C(C), KEY G(G), KEY N(N))";
  $DBHest->do($create);
  my $seqio = Bio::SeqIO->new(-file => "$ESTDIR\/$ESTFILE",
                              -format => 'fasta');
  while(my $seq = $seqio->next_seq) {
    my $seqstring = uc($seq->seq);
    $seq->seq($seqstring);
    my $href = Bio::Tools::SeqStats->count_monomers($seq);
    my $string = "INSERT INTO $seqtable SET ".
                 "clone = ".$DBHest->quote($seq->id).
                 ", length = ".$seq->length.
                 ", seq = ".$DBHest->quote($seq->seq).
                 ", A = ".($href->{A}||'0').", T = ".($href->{T}||'0').
                 ", C = ".($href->{C}||'0').", G = ".($href->{G}||'0');
    $string .= ", N = ".($href->{N}||'0') if exists $href->{N};
    my $sth = $DBHest->prepare($string);
    $sth->execute;
  }
}

sub get_codeml_res {
  my $file = shift;
  my $string = `tail -1 "$file"`;
  # dN/dS= 0.0053  dN= 0.3757  dS=70.7432
  $string =~ /.+dN\/dS\=\s*(\S+)\s*dN\=\s*(\S+)\s*dS\=\s*(\S+)$/;
  my @ratio = ($1,$2,$3);
  print "$string\n@ratio\n";
  return \@ratio;
}

sub get_embl_translation {
  my $seq = shift;
  # Test the vertebrate default codon table
  my $tvert = $seq->translate(-orf => 1, -codontable_id => 1);
  # Test the vertebrate mitocondrial codon table
  my $tvmit = $seq->translate(-orf => 1, -codontable_id => 2);
  return($tvert,1) if $tvert->length >= $tvmit->length;
  return($tvmit,2) if $tvert->length < $tvmit->length;
}

sub get_cnt {
  my $string;
  while(<DATA>) {
    $string .= $_;
  }
  return $string;
}

# Primary key should be in the first column
sub get_table_in_href {
  my $file = shift;
  my $res = {};
  open(IN,$file);
  my $head = <IN>;
  chomp($head);
  my @rowname = split(/\t/,$head); # qw(clone hit description evalue coverage identity strand frame)
  while(my $row = <IN>) {
    chomp($row);
    my @field = split(/\t/,$row);
    for(my $c = 1; $c <= $#rowname; $c++) {
      $res->{$field[0]}->{$rowname[$c]} = $field[$c];
    }
  }
  return $res;
}

sub get_embl_cds {
  my $DBH = shift;
  my $id = shift;
  my $sth = $DBH->prepare('SELECT EMBL_CDS from idmapping WHERE UniprotKB_AC = '.$DBH->quote($id));
  $sth->execute;
  my $res = $sth->fetchrow_hashref;
  my $field = $res->{EMBL_CDS};
  my @ids = split(/\; /,$field);
  foreach my $el(@ids) {
    $el =~ s/\.\d$//;
    return $el if $el =~ /[\w\d]/;
  }
}

sub get_translation {
  my $seq = shift;
  my $strand = shift;
  my $frame = shift;
  my $codontable_id = shift;
  my $pre_nucl;
  my $nucl;
  my $prot;
  if($strand == -1) {
    $pre_nucl = $seq->revcom;
    $prot = $pre_nucl->translate(-frame => $frame, -terminator => 'X', -codontable_id => $codontable_id);
  }
  elsif($strand == 1) {
    $prot = $seq->translate(-frame => $frame, -terminator => 'X', -codontable_id => $codontable_id);
    $pre_nucl = $seq;
  }
  $nucl = $pre_nucl->subseq(1,$pre_nucl->length) if $frame == 0;
  $nucl = $pre_nucl->subseq(2,$pre_nucl->length) if $frame == 1;
  $nucl = $pre_nucl->subseq(3,$pre_nucl->length) if $frame == 2;

  my $shorter = $prot->subseq(1,$prot->length-1);
  $prot->seq($shorter) if $prot->length x 3 > $seq->length;
  $seq->seq($nucl);
  return($seq,$prot);
}

sub run_blastplus {
  my $program = shift;
  my $dbname = shift;

  my $analysis = "$program\_$dbname";
  my $query = $ESTFILE;
  my $name = "$ESTFILE";
  $name =~ s/\.fasta$//;
  $name =~ s/\.fa$//;

  my $exe = $ENV{BLASTPLUSDIR}.'/'.$program;
  my $db = $ENV{BLASTDBDIR}.'/'.$dbname;
  my $out = $name.'.'.$analysis;

  my $command = "$exe -query $query -db $db -out $out -num_threads 10 -best_hit_overhang 0.1 -evalue 0.01";

  print "$command\n";
  print "Executing $analysis\...\n";
  print "$analysis... ";

  system($command);
  sleep(3);

  if(-e $out) { print "OK\n\n" }
  else { print "\nERROR: problems running $analysis\n\n" }
  return $out;
}

sub parse_blast_results {
  my $out = shift;
  my $out_table = "$out\.table";
  open(OUT,">$out_table");
  print OUT join("\t",qw(clone hit description evalue coverage identity strand frame))."\n";

  my $in = new Bio::SearchIO(-format => 'blast', 
                             -file => $out);

  while(my $result = $in->next_result) {
    my $program = lc($result->algorithm);
    my $candidate = {};
    while(my $hit = $result->next_hit) {
      while(my $hsp = $hit->next_hsp) {
        my $strand = $hsp->strand;
        my $frame = $hsp->frame('query');
        my $bad_description = 1;
        my $identity = $hsp->percent_identity;
        my $coverage = $hsp->length('total') / $result->query_length * 100;
        if($program eq 'blastx') {
          $coverage = $hsp->length('total') / ($result->query_length / 3) * 100;
          $bad_description = check_description($hit);
        }
        $coverage =~ s/^(\d+)\.\d+$/$1/;
        $identity =~ s/^(\d+)\.\d+$/$1/;

        # HARD CODED CUTOFFS!!!
        next unless ($identity >= 30 && $coverage >= 20);

        if(
        (! exists $candidate->{hsp}) || 
        ($hsp->evalue < $candidate->{hsp}->evalue) || 
        ($candidate->{bad_description} > $bad_description)) {
          $candidate = candidate_this($hit,$hsp,$coverage,$identity,$strand,$frame,$bad_description);
        }
      }
    }

    if(exists $candidate->{hsp}) {
    print OUT $result->query_name ."\t".
              strip_id($candidate->{hit}->accession) ."\t".
              $candidate->{hit}->description ."\t".
              $candidate->{hsp}->evalue ."\t".
              $candidate->{coverage} ."\t".
              $candidate->{identity}, "\t".
              $candidate->{strand}."\t".
              $candidate->{frame}."\n";
    }
  }
  close(OUT);
  return $out_table;
}

sub check_description {
  # SPECIAL EVALUATION FOR UNIPROT RESULTS IN ORDER
  # TO AVOID TO GET HIT WIT UNINFORMATIVE DESCRIPTIONS
  # perl -e '$string="quel giorno cadevo";@a=qw(quel uno tre);$c=map{$string=~/$_/}@a;print"$c\n"'
  my $hit = shift;
  my @bad_descriptors = (
    'Predicted protein OS=',
    'Predicted protein (Fragment) OS=',
    'Putative uncharacterized protein',
    'Uncharacterized protein',
    'whole genome shotgun sequence'
  );
  my $bad_description = map { $hit->description =~ /$_/ } @bad_descriptors;
  $bad_description ++ if $hit->description =~ /^\w+ OS\=Drosophila/;
  $bad_description ++ if $hit->description =~ /\d\d\d\d OS\=/;
  $bad_description ++ if $hit->description =~ /\d\d\d\d protein OS\=/;
  $bad_description ++ if $hit->description =~ /\d\d\d\d protein (Fragment) OS\=/;
  return $bad_description;
}

sub candidate_this {
  my $candidate = {};
  $candidate->{hit} = shift;
  $candidate->{hsp} = shift;
  $candidate->{coverage} = shift;
  $candidate->{identity} = shift;
  $candidate->{strand} = shift;
  $candidate->{frame} = shift;
  $candidate->{bad_description} = shift;
  return $candidate;
}

sub strip_id {
  my $id = shift;
  if($id =~ /\|/) {
    my @field = split(/\|/,$id);
    return $field[1];
  }
  else {
    return $id;
  }
}

sub connect_to_db {
  my $db = shift;
  my $usr = shift;
  my $pwd = shift;
  my $host = shift;
  my $dsn = 'dbi:mysql:'.$db;
  $dsn .= ':'.$host if $host; # IN THE CURRENT DBI POD VERSION THERE IS THE '@' IN THE PLACE OF ':'
  my $dbh = DBI->connect($dsn,$usr,$pwd,{PrintWarn=>1,PrintError=>1,RaiseError=>1}) or die $DBI::errstr;
  return $dbh;
}

__DATA__
      seqfile = test.codon
     treefile = test.tree
      outfile = test.codeml

        noisy = 0   * 0,1,2,3,9: how much rubbish on the screen
      verbose = 0   * 1: detailed output, 0: concise output
      runmode = -2  * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

    cleandata = 1   * "I added on 07/07/2004" Mikita Suyama

      seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        model = 2
                    * models for codons:
                        * 0:one, 1:b, 2:2 or more dN/dS ratios for branches

      NSsites = 0   * dN/dS among sites. 0:no variation, 1:neutral, 2:positive
        icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
        Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all

    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2   * initial or fixed kappa
    fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate
        omega = 1   * initial or fixed omega, for codons or codon-transltd AAs

    fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = .0  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * different alphas for genes
        ncatG = 4   * # of categories in the dG or AdG models of rates

        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0   * (1/0): rates (alpha>0) or ancestral states (alpha=0)
       method = 0   * 0: simultaneous; 1: one branch at a time

