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
# select e.*,definition,class from EST_blastx_uniprot e, EST_blastx_uniprot_GO eg, GO g where e.hit=eg.hit and eg.go=g.go and definition like '%neuro%';
my $load_go = 0;
my $go_term_definition = '/media/LOCAL_DATA/GO/GO.terms_alt_ids';

# DATABASE SPECIFIC SETTINGS
my $DBuniprot = 'uniprot_idmapping';
my $DBest = 'OCTOPUS_EST';
my $USR = 'mysql_dev';
my $PWD = 'riiGbs';
#my $PWD = 'dEvEl0pEr';
my $HOST;
my $DBHuniprot = connect_to_db($DBuniprot,$USR,$PWD,$HOST);
my $DBHest = connect_to_db($DBest,$USR,$PWD,$HOST);

my $file = '/media/EXTERNAL_DATA/BLAST_PIPELINE_EYE/EYE.blastx_uniprot_sprot_trembl.table';
my ($filename,$dir) = get_filename_dir($file);

my $tablename = "$filename";
$tablename =~ s/\.table$//;
$tablename =~ s/\_trembl$//;
$tablename =~ s/\_sprot$//;
$tablename =~ s/\./\_/g;
my $gotablename = "$tablename\_GO";
my $gofile = "$dir\/$gotablename";

my $create_bl = "CREATE TABLE $tablename (
  clone VARCHAR(256) PRIMARY KEY,
  hit VARCHAR(256),
  description TEXT,
  evalue DOUBLE,
  coverage INT(3),
  identity INT(3),
  strand ENUM('1','-1'),
  frame ENUM('0','1','2'),
  KEY hit(hit),
  KEY evalue(evalue),
  KEY coverage(coverage),
  KEY identity(identity)
)";
$DBHest->do($create_bl);
my $loader_bl = "LOAD DATA LOCAL INFILE \'$file\' INTO TABLE $tablename IGNORE 1 LINES";
$DBHest->do($loader_bl);

my $create_go = "CREATE TABLE $gotablename (
  hit VARCHAR(256),
  go VARCHAR(256),
  KEY hit(hit),
  KEY go(go)
)";
$DBHest->do($create_go);
populate_go_table($file,$gofile,$gotablename);

populate_go_term_definition($go_term_definition) if $load_go;

sub populate_go_term_definition {
  my $file = shift;
  system("grep -v \'\\\!\' $file \> $file\.mod");
  my $create = "CREATE TABLE GO (
    go VARCHAR(256) PRIMARY KEY,
    secondary VARCHAR(256),
    definition VARCHAR(256),
    class ENUM('F','P','C')
  )";
  $DBHest->do($create);
  my $load = "LOAD DATA LOCAL INFILE \'$file\.mod\' INTO TABLE GO";
  $DBHest->do($load);
}

sub populate_go_table {
  my $blastxtable = shift;
  my $file = shift;
  my $gotablename = shift;
  my $table = get_table_in_href($blastxtable);
  foreach my $clone(keys %$table) {
    my $hit = $table->{$clone}->{hit};
    my $go = get_go($hit);
    foreach my $goid(@$go) {
      my $sth = $DBHest->prepare("INSERT INTO $gotablename SET".
                                 " hit = ".$DBHest->quote($hit).
                                 ", go = ".$DBHest->quote($goid));
      $sth->execute;
    }
  }
}

sub get_go {
  my $hit = shift;
  my $sth = $DBHuniprot->prepare('SELECT GO FROM idmapping WHERE UniProtKB_AC = '.$DBHuniprot->quote($hit));
  $sth->execute;
  my $res = $sth->fetchrow_hashref;
  return undef unless $res->{GO};
  my @go = split(/\; /,$res->{GO});
  return \@go;
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

sub get_filename_dir {
  my $file = shift;
  my @path = split(/\//,$file);
  return (pop(@path),join('/',@path));
}
