#!/usr/bin/perl
use strict;
use warnings;
use lib "$ENV{HOME}/src/bioperl-live";
use lib "$ENV{HOME}/src/bioperl-run/lib";
use lib "$ENV{HOME}/src/RS9";
use Biosub::Utils;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::Index::Fasta;
use Data::Dumper;

# Load the config file
my $CONF = Biosub::Utils->get_conf('CONF_BLAST_PIPELINE');
# Make directory to use for outputs and create the log
# The '0' is when you not need to create the directory
# when all the output have already been created
my($RES_DIR,$LOG) = Biosub::Utils->start_log($CONF,1);
# Check for the existence of all files and dirs
# in the config file
Biosub::Utils->check_conf_files_folders($CONF);

open(LOG,">>$LOG") or die "\nCannot open log file $LOG\: $!\n\n";

my $blastx_uniprot = run_blastplus('blastx','uniprot_sprot_trembl');
#my $blastn_rfam = run_blastplus('blastn','Rfam');
#my $blastn_nt = run_blastplus('blastn','nt');

my $blastx_uniprot_tab = parse_blast_results($blastx_uniprot);
#my $blastn_rfam_tab = parse_blast_results($blastn_rfam);
#my $blastn_nt_tab = parse_blast_results($blastn_nt);

sub run_blastplus {
  my $program = shift;
  my $dbname = shift;

  my $analysis = "$program\_$dbname";
  my $query = $CONF->{INPUT_DIR}.'/'.$CONF->{INPUT_FILE};
  my $name = $CONF->{INPUT_FILE};
  $name =~ s/\.fa$//;
  $name =~ s/\.fasta$//;

  my $exe = $CONF->{BLASTPLUS_DIR}.'/'.$program;
  my $db = $CONF->{BLASTDB_DIR}.'/'.$dbname;
  my $out = $RES_DIR.'/'.$name.'.'.$analysis;

  my $command = "$exe -query $query -db $db -out $out -num_threads 12 -best_hit_overhang 0.1 -evalue 0.01";

  print LOG "$command\n";
  print LOG "Executing $analysis\...\n";
  print LOG "$analysis... ";

  system($command);
  sleep(3);

  if(-e $out) { print LOG "OK\n\n" }
  else { print LOG "\nERROR: problems running $analysis\n\n" }
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
              Biosub::Utils->strip_id($candidate->{hit}->accession) ."\t".
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
