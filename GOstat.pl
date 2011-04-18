#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use DBI;
use DBD::mysql;

my $CREATE_UNIPROT_COUNTS = 0;

# DATABASE SPECIFIC SETTINGS
my $DBuniprot = 'uniprot_idmapping_2';
my $DBest = 'OCTOPUS_EST';
my $USR = 'mysql_dev';
my $PWD = 'riiGbs';
#my $PWD = 'dEvEl0pEr';
my $HOST;
my $DBHu = connect_to_db($DBuniprot,$USR,$PWD,$HOST);
my $DBHe = connect_to_db($DBest,$USR,$PWD,$HOST);

my $gotable = 'EYE_blastx_uniprot_GO';
my $esttable = 'EYE_blastx_uniprot';

# select distinct clone, go from EST_blastx_uniprot e, EST_blastx_uniprot_GO eg where e.hit=eg.hit;

if($CREATE_UNIPROT_COUNTS) {
  my $create = 'CREATE TABLE go_counts(go VARCHAR(256) PRIMARY KEY, n INT(5), KEY n(n))';
  $DBHu->do($create);
  my $insert = 'INSERT INTO go_counts SELECT DISTINCT go, count(*) n FROM GO WHERE GO IS NOT NULL GROUP BY go';
  $DBHu->do($insert);
}

my $ipu = get_uniprot_universe();
my $estu = get_population_universe($esttable);

print join("\t",qw(id in_est in_ip tot_est tot_ip))."\n";
work_on_population($gotable);

sub get_uniprot_universe {
  my $sth = $DBHu->prepare('SELECT COUNT(DISTINCT UniProtKB_AC) n FROM idmapping');
  $sth->execute;
  my $res = $sth->fetchrow_hashref;
  return $res->{n};
}

sub get_population_universe {
  my $sth = $DBHe->prepare("SELECT COUNT(DISTINCT hit) n FROM $esttable");
  $sth->execute;
  my $res = $sth->fetchrow_hashref;
  return $res->{n};
}

sub work_on_population {
  my $gotable = shift;
  my $sth = $DBHe->prepare("SELECT DISTINCT go FROM $gotable");
  $sth->execute;
  while(my $res = $sth->fetchrow_hashref) {
    my $goid = $res->{go};
    my $goid_ipr = count_this_uniprot($goid);
    my $goid_est = count_this_est($goid,$gotable);
    print join("\t",$goid,$goid_est,$goid_ipr,$estu,$ipu)."\n";
  }
}

sub count_this_uniprot {
  my $id = shift;
  my $sth = $DBHu->prepare("SELECT COUNT(DISTINCT uniprot) n FROM GO WHERE go = ".$DBHu->quote($id));
  $sth->execute;
  my $res = $sth->fetchrow_hashref;
  return $res->{n};
}

sub count_this_est {
  my $id = shift;
  my $gotable = shift;
  my $sth = $DBHe->prepare("SELECT COUNT(DISTINCT hit) n FROM $gotable WHERE go = ".$DBHe->quote($id));
  $sth->execute;
  my $res = $sth->fetchrow_hashref;
  return $res->{n};
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

