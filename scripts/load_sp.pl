#!/usr/local/bin/perl -w

# requires bioperl

BEGIN {
    if (defined($ENV{GO_ROOT})) {
	use lib "$ENV{GO_ROOT}/perl-api";
    }
}
use strict;
use FileHandle;
use GO::AppHandle;
use GO::Model::Xref;
use Bio::DB::SwissProt;
use Bio::SeqIO;
use Bio::Index::Swissprot;

if (!@ARGV) {
    print usage();
    exit;
}

use Getopt::Long;

my $apph = GO::AppHandle->connect(\@ARGV);
my $h = {};
GetOptions($h,
           "help=s",
           "speciesdb=s",
           "store",
           "detailed",
	   "swissdir=s",
           "skipnoterms",
           "location=s",
           "evcode|e=s@",
           "out=s",
          );

my $outfh;
if ($h->{out}) {
    $outfh = FileHandle->new(">".$h->{out});
}

if ($h->{help}) {
    print usage();
    exit;
}

if ($h->{evcode}) {
    $apph->filters->{evcodes} = $h->{evcode};
}

my $cachedir = $h->{cachedir} || ".";
`mkdir $cachedir` unless -d $cachedir;

my $swissdir = $h->{swissdir} || "proteomes";

my $taxa_id_lookup =
  $apph->get_taxa_id_lookup;
my %spindexh = ();

my @bad = ();
while (my $fn = shift @ARGV) {

    my $speciesdb;
    if ($h->{speciesdb}) {
        $speciesdb = $h->{speciesdb};
    }
    elsif ($fn =~ /gp2protein\.(\w+)\.gz/)  {
        $speciesdb = $1;
    }
    elsif ($fn =~ /gp2protein\.(\w+)/)  {
        $speciesdb = $1;
    }
    else {
        die("$fn not right format; must be gp2protein.SPECIESDB");
    }
    my $fh;
    if ($fn =~ /\.gz/) {
        $fh = FileHandle->new("gzip -dc $fn |") || die "could not gzip -dc $fn";
    }
    else {
        $fh = FileHandle->new($fn) || die "cannot open $fn";
    }
    while(<$fh>) {
        chomp;
        next if /^\!/;
        next unless $_;
        my ($modacc, $spxref) = split(/\t/, $_);
        if (!($spxref =~ /(\S+):(\S+)/)) {
            warn("seq dbxref must be DB:ACC format");
            next;
        }
        my $spid = $2 || die($_);
        my $curr_speciesdb = $speciesdb;
        if ($modacc =~ /(\w+)\:(\S+)/i) {
            $curr_speciesdb = $1;
            $modacc = $2;
	    if ($curr_speciesdb eq 'MGI' && int($modacc)) {
		$modacc = "MGI:$modacc";
	    }
        }
	my $taxa_id = $taxa_id_lookup->{lc("$curr_speciesdb:$modacc")};
	if (!$taxa_id) {
	  print STDERR "No taxa id for $curr_speciesdb:$modacc!\n";
	  next;
	}
	my $f = "$swissdir/$taxa_id.SPC";
	if (! -f $f) {
            if (-f "$f.gz") {
                if (system("gzip -d $f.gz")) {
                    print STDERR "cannot uncompress $f.gz\n";
                    next;
                }
            }
            else {
                # hmm, should we fetch over the network?
                print STDERR "No proteome file for taxa $taxa_id\n";
                next;
            }
	}

	my $spindex =
	  $spindexh{$taxa_id};
	
	if (!$spindex) {
	    $spindex =
	      Bio::Index::Swissprot->new(-file=>"$f.idx",
					 -write_flag=>"WRITE");
	    $spindex->make_index("$f");
	    $spindexh{$taxa_id} = $spindex;
	}
        print STDERR "Fetching: $spid in taxa $taxa_id\n";
        my $bseq;
        eval {
            # Fetch seq;
            # try from cached swissprot file first:
            if ($spindex) {
                $bseq =
                  $spindex->fetch($spid);
            }
#            if (!$bseq) {
#                print STDERR "Trying web query....\n";
		#                $bseq =
#                  $db->get_Seq_by_id($spid);
#            }
        };
        if ($@) {
            bad("err fetching $spid = $@");
            next;
        }
        if (!$bseq) {
            bad("Can't find $spid!!!");
            next;
        }
#        $outfa->write_seq($bseq);
#        $outswiss->write_seq($bseq);

	my $seq = GO::Model::Seq->new;
	$seq->pseq($bseq);

        my $xref =
          GO::Model::Xref->new;
        $xref->xref_key($seq->accession);
        $xref->xref_dbname("UniProt");
        $seq->add_xref($xref);
        $seq->display_id($seq->accession);

        my $gps =
          $apph->get_products({xref=>{xref_key=>$modacc, xref_dbname=>$curr_speciesdb}});
        if (!@$gps) {
            warn("$modacc $curr_speciesdb not in db");
            next;
        }
        elsif (@$gps > 1) {
            warn(@$gps." >1 match $modacc $curr_speciesdb");
        }
	foreach my $gp (@$gps) {
            if ($h->{store}) {
                $apph->set_product_seq($gp, $seq);
            }
            if ($outfh) {
                my $terms = $apph->get_terms({product=>$gp});
                if (!@$terms && $h->{skipnoterms}) {
                    next;
                }
                my @h_elts = ();
                foreach my $term (@$terms) {
                    my $al = $term->selected_association_list;
                    my %codes = ();
                    map { $codes{$_->code} = 1 } map { @{$_->evidence_list} } @$al;
                    push(@h_elts,
                         sprintf("[%s%s evidence=%s]",
                                 $term->public_acc,
                                 $h->{detailed} ? " \"".$term->name."\"" : "",
                                 join(";", keys %codes),
                                )
                        );
                }
                my $hdr =
                  sprintf("%s|%s SPTR:%s symbol:%s%s taxa:$taxa_id %s %s",
                          $curr_speciesdb,
                          $modacc,
                          $spid,
                          $gp->symbol,
                          $gp->full_name ? " \"".$gp->full_name."\"" : "",
                          join(" ", @h_elts),
                          join(" ", map {$_->xref_dbname.":".$_->xref_key } @{$seq->xref_list}),
                         );
                $seq->description($hdr);
                print $outfh $seq->to_fasta;
            }
        }
    }
    $fh->close;
}
$outfh->close if $outfh;
exit 0;

sub bad {
    push(@bad, join("", @_));
    warn($bad[-1]);
}

sub usage {
    print "fetch_sp.pl  [-d dbname] [-h dbhost]  [-dbms dbms] [-detailed] [-location country-code] [-cachedir dir] [-evcode code] <gp2protein file>\n";
    print <<EOM;

This script will produce a fasta file of sequences, from a gp2protein
file, by downloading the sequences from an expasy server.

As a side effect it will produce swiss and fasta files, mimicing the
source data (eg with no GO information)

If you specify -store, it will also store the sequences in the go database specified.

If you specify

it relies on you having downloaded a swissprot format file of
sequences. this should preferably be filtered for the organism of
interest (the script cycles through all entries in the swissprot file
outputtinga fasta entry every time it sees a mapping in the gp2protein
file)

use -detailed to output the GO name as well as the GO accession on
every line

 REQUIREMENTS: You must have bioperl installed; see http://www.bioperl.org

EOM
}
