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

# The file with the GO definitions coming from
# http://www.geneontology.org/doc/GO.terms_alt_ids

# DATABASE SPECIFIC SETTINGS
my $DB = 'uniprot_idmapping';
my $USR = 'mysql_dev';
my $PWD = 'riiGbs';
my $HOST;

# FILE SETTINGS
my $UNIPROT = 'idmapping_selected.tab';
my $GO = 'GO.terms_alt_ids';

my $DBH = connect_to_db($DB,$USR,$PWD,$HOST);

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
  GI              TEXT
);

open(IN,$UNIPROT);
open(OUT,">$UNIPROT\.parsed");
print OUT join("\t",keys %TO_USE_KEYS)."\n";
open(OUTGO,">$UNIPROT\.parsed.go");
print OUTGO join("\t",'uniprot','go')."\n";
open(OUTGI,">$UNIPROT\.parsed.gi");
print OUTGI join("\t",'uniprot','gi')."\n";

my $CREATE = 'CREATE TABLE idmapping (UniProtKB-AC VARCHAR(256) PRIMARY KEY, ';
for(my $c = 0; $c <= $#ALL_KEYS; $c ++) {
  if(exists $TO_USE_KEYS{$ALL_KEYS[$c]}) {
    $CREATE .= "$ALL_KEYS[$c] $TO_USE_KEYS{$ALL_KEYS[$c]}\, " unless $ALL_KEYS[$c] eq 'UniProtKB-AC';
  }
}
$CREATE =~ s/\-/\_/g;
$CREATE =~ s/\, $/\)/;
$DBH->do($CREATE);

my $CREATEgo = 'CREATE TABLE idgo (id VARCHAR(256), go VARCHAR(256), CONSTRAINT idgo UNIQUE(id,go))';
$DBH->do($CREATEgo);

my $CREATEgi = 'CREATE TABLE idgi (id VARCHAR(256), gi VARCHAR(256), CONSTRAINT idgi UNIQUE(id,gi))';
$DBH->do($CREATEgi);

while(my $row = <IN>) {
  chomp($row);
  my $string;
  my @field = split(/\t/,$row);
  for(my $c = 0; $c <= $#ALL_KEYS; $c ++) {
    if(exists $TO_USE_KEYS{$ALL_KEYS[$c]}) {
      my $toadd = $field[$c] || '\N';
      $string .= '"'.$toadd.'"'."\t";
    }
    if($ALL_KEYS[$c] eq 'GO') {
      print_go_field($field[0],$field[$c]);
    }
    if($ALL_KEYS[$c] eq 'GI') {
      print_gi_field($field[0],$field[$c]);
    }
  }
  $string =~ s/\t$/\n/;
  print OUT $string;
}

close(IN);
close(OUT);
close(OUTGO);
close(OUTGI);

my $LOADER = qq( LOAD DATA LOCAL INFILE '$UNIPROT\.parsed'
                 INTO TABLE idmapping 
                 FIELDS TERMINATED BY '\\t' ENCLOSED BY '"'
                 LINES TERMINATED BY '\\n'
                 IGNORE 1 LINES );
$DBH->do($LOADER);

my $LOADERgo = qq( LOAD DATA LOCAL INFILE '$UNIPROT\.parsed.go'
                   INTO TABLE idgo 
                   FIELDS TERMINATED BY '\\t' ENCLOSED BY '"'
                   LINES TERMINATED BY '\\n'
                   IGNORE 1 LINES );
$DBH->do($LOADERgo);

my $LOADERgi = qq( LOAD DATA LOCAL INFILE '$UNIPROT\.parsed.gi'
                   INTO TABLE idgi 
                   FIELDS TERMINATED BY '\\t' ENCLOSED BY '"'
                   LINES TERMINATED BY '\\n'
                   IGNORE 1 LINES );
$DBH->do($LOADERgi);

populate_go_term_definition($GO);

sub print_go_field {
  my $id = shift;
  my $field = shift || '\N';
  my @go = split(/\; /,$field);
  print OUTGO join("\t", '"'.$id.'"', '"'.$_.'"')."\n" foreach @go;
}

sub print_gi_field {
  my $id = shift;
  my $field = shift || '\N';
  my @go = split(/\; /,$field);
  print OUTGI join("\t", '"'.$id.'"', '"'.$_.'"')."\n" foreach @go;
}

sub populate_go_term_definition {
  my $file = shift;
  system("grep -v \'\\\!\' $file \> $file\.mod");
  my $create = "CREATE TABLE GO (
    go VARCHAR(256) PRIMARY KEY,
    secondary VARCHAR(256),
    definition VARCHAR(256),
    class ENUM('F','P','C')
  )";
  $DBH->do($create);
  my $load = "LOAD DATA LOCAL INFILE \'$file\.mod\' INTO TABLE GO";
  $DBH->do($load);
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

