#!/usr/local/bin/perl

# POD docs at end of file

use strict;
use Getopt::Long;
use FileHandle;
use GO::AppHandle;

$|=1;


if (!@ARGV) {
    system("perldoc $0");
    exit 1;
}

use Getopt::Long;

my $apph = GO::AppHandle->connect(\@ARGV);
my $opt = {
           field=>"acc",
           cutoff=>0.1,
           min_gps=>2,
          };
GetOptions($opt,
	   "help|h",
           "input|i=s%",
           "field=s",
           "filter=s%",
           "min_gps=s",
           "cutoff=s",
           "speciesdb=s@",
           "evcode|e=s@",
	  );

$apph->filters($opt->{filter} || {});
if ($opt->{evcode}) {
    $apph->filters->{evcodes} = $opt->{evcode};
}
if ($opt->{speciesdb}) {
    $apph->filters->{speciesdbs} = $opt->{speciesdb};
}


if ($opt->{help}) {
    system("perldoc $0");
    exit;
}

my @ids = @ARGV;
my $input = $opt->{input};
if ($input) {
    open(F,$input);
    while(<F>) {
        chomp;
        push(@ids,$_);
    }
    close(F);
}
my $field = $opt->{field};

my $eh = $apph->get_enriched_term_hash([ map { {$field=>$_} } @ids]);
my @erows =
  sort {
      $a->{p_value} <=> $b->{p_value}
  } values %$eh;
foreach (@erows) {
    next unless $_->{p_value} <= $opt->{cutoff};
    next if $_->{n_gps_in_sample_annotated} < $opt->{min_gps};
    
    printf("%s sample:%d/%d database:%d/%d P-value:%s Corrected:%s \"%s\" Genes: %s\n",
           $_->{term}->acc,
           $_->{n_gps_in_sample_annotated},
           $_->{n_gps_in_sample},
           $_->{n_gps_in_database_annotated},
           $_->{n_gps_in_database},
           $_->{p_value},
           $_->{corrected_p_value},
           $_->{term}->name,
           join('; ',map {sprintf("%s[%s:%s]", $_->symbol, $_->xref->xref_dbname, $_->acc)} @{$_->{gps_in_sample_annotated}}))
}
exit 0;


__END__

=head1 NAME

go-find-enriched-terms.pl

=head1 SYNOPSIS

  go-find-enriched-terms.pl -d go -h localhost -field synonym YNL116W YNL030W YNL126W

  go-find-enriched-terms.pl -d go -h localhost -field acc -i gene-ids.txt

=head1 DESCRIPTION

Performs a term enrichment analysis. Uses hypergeometric distribution,
takes entire DAG into account.

First the database will be queried for matching gene products. Any
filters in place will be applied (or you can pass in a list of gene
products previously fetched, eg using $apph->get_products).

The matching products count as the *sample*. This is compared against
the gene products in the database that match any pre-set filters
(statistics may be more meaningful when a filter is set to a
particular taxon or speciesdb-source).

We then examine terms that have been used to annotate these gene
products. Filters are taken into account (ie if !IEA is set, then no
IEA associations will count). The DAG is also taken into account - so
anything annotated to a process will count as being annotated to
biological_process. This means the fake root "all" will always have
p-val=1. Currently the entire DAG is traversed, relationship types are
ignored (in future it may be possible to specify deduction rules -
this will be useful when the number of relations in GO progresses
beyond 2, or when this code is used with other ontologies)

Results are returned as a hash-of-hashes, outer hash keyed by term
acc, inner hash specifying the fields:


=head1 ARGUMENTS

Arguments are either connection arguments, generic to all go-db-perl
scripts, and arguments specific to this script

=head2 CONNECTION ARGUMENTS

Specify db connect arguments first; these are:

=over

=item -dbname [or -d]

name of database; usually "go" but you may want point at a test/dvlp database

=item -dbuser 

name of user to connect to database as; optional

=item -dbauth

password to connect to database with; optional

=item -dbhost [or -h]

name of server where the database server lives; see
http://www.godatabase.org/dev/database for details of which servers
are active. or you can just specify "localhost" if you have go-mysql
installed locally

=back

=head2 SCRIPT ARGUMENTS

=over

=item -field FIELDNAME

May be: acc, name, synonym

=item -input [or -i] FILE

a file of ids or symbols to search with; newline separated

=item -filter FILTER=VAL

see L<GO::AppHandle> for explanation of filters

multiple args can be passed:

  -filter taxid=7227 -filter 'evcode=!IEA'

=item -speciesdb SPECIESDB

filter by source database

multiple args can be passed

  -e SGD -e FB

=item -evcode [or -e] EVCODE

filter by evidence code

negative arguments can be passed

  -e '!IEA'

this opt can be passed multiple times:

  -e ISS -s IDA -e IMP

=item -cutoff P-VAL

p-value report threshold

=back

=head1 OUTPUT

The default output produces tab-delimited rows with the following data:

=head2 DOCUMENTATION

L<http://www.godatabase.org/dev>

=head2 SEE ALSO

L<GO::AppHandle>

=cut

