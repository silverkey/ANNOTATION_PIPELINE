#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use DBI;
use DBD::mysql;

# The idmapping_selected.tab table is downloaded from UNIPROT
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
# The definition of the fields is in the file
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README

# DATABASE SPECIFIC SETTINGS
my $DB = 'uniprot_idmapping_GI';
my $USR = 'mysql_dev';
my $PWD = 'riiGbs';
my $HOST;

my $DBH = connect_to_db($DB,$USR,$PWD,$HOST);

my $FILE = 'idmapping_selected.tab';

my @ALL_KEYS = qw(
  UniProtKB-AC
  UniProtKB-ID
  GeneID
  RefSeq
  GI
  PDB
  GO
  IPI
  UniRef100
  UniRef90
  UniRef50
  UniParc
  PIR
  NCBI-taxon
  MIM
  UniGene
  PubMed
  EMBL
  EMBL-CDS
  Ensembl
  Ensembl_TRS
  Ensembl_PRO
  Additional PubMed
);

my %TO_USE_KEYS = qw(
  UniProtKB-AC    VARCHAR(256)
  UniProtKB-ID    VARCHAR(256)
  RefSeq          VARCHAR(256)
  GO              TEXT
  MIM             VARCHAR(256)
  UniGene         VARCHAR(256)
  EMBL            VARCHAR(256)
  EMBL-CDS        VARCHAR(256)
  Ensembl         VARCHAR(256)
  Ensembl_TRS     VARCHAR(256)
  Ensembl_PRO     VARCHAR(256)
);

open(IN,$FILE);
open(OUT,">$FILE\.parsed");
print OUT join("\t",keys %TO_USE_KEYS)."\n";
open(OUTGO,">$FILE\.parsed.go");
print OUTGO join("\t",'uniprot','go')."\n";

my $CREATE = 'CREATE TABLE idmapping (UniProtKB-AC VARCHAR(256) PRIMARY KEY, ';
for(my $c = 0; $c <= $#ALL_KEYS; $c ++) {
  if(exists $TO_USE_KEYS{$ALL_KEYS[$c]}) {
    $CREATE .= "$ALL_KEYS[$c] $TO_USE_KEYS{$ALL_KEYS[$c]}\, " unless $ALL_KEYS[$c] eq 'UniProtKB-AC';
  }
}
$CREATE =~ s/\-/\_/g;
$CREATE =~ s/\, $/\)/;

$DBH->do($CREATE);

while(my $row = <IN>) {
  chomp($row);
  my $string;
  my @field = split(/\t/,$row);
  for(my $c = 0; $c <= $#ALL_KEYS; $c ++) {
    if(exists $TO_USE_KEYS{$ALL_KEYS[$c]}) {
      $string .= ($field[$c] || '\N')."\t";
    }
    if($ALL_KEYS[$c] eq 'GO') {
      print_go_field($field[0],$field[$c]);
    }
  }
  $string =~ s/\t$/\n/;
  print OUT $string;
}

close(IN);
close(OUT);

my $LOADER = "LOAD DATA LOCAL INFILE \'$FILE\.parsed\' INTO TABLE idmapping IGNORE 1 LINES";

$DBH->do($LOADER);

sub print_go_field {
  my $id = shift;
  my $field = shift || '\N';
  my @go = split(/\; /,$field);
  print OUTGO join("\t",$id,$_)."\n" foreach @go;
  print join("\t",$id,$_)."\n" foreach @go;
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

__END__

Cut and paste from the current README

1. UniProtKB-AC
2. UniProtKB-ID
3. GeneID (EntrezGene)
4. RefSeq
5. GI
6. PDB
7. GO
8. IPI
9. UniRef100
10. UniRef90
11. UniRef50
12. UniParc
13. PIR
14. NCBI-taxon
15. MIM
16. UniGene
17. PubMed
18. EMBL
19. EMBL-CDS
20. Ensembl
21. Ensembl_TRS
22. Ensembl_PRO
23. Additional PubMed

