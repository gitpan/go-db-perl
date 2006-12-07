#!/usr/local/bin/perl -w

# All tests must be run from the software directory;
# make sure we are getting the modules from here:
use lib '.';
use strict;
use GO::TestHarness;
use Set::Scalar;

n_tests(2);

my $apph = get_readonly_apph();
stmt_ok;

# lets check we got stuff

$apph->filters({speciesdb=>"SGD"});
my $h = $apph->get_enriched_term_hash([ map { {synonym=>$_} } qw(YNL116W YNL030W) ]);
#my $h2 = $apph->get_enriched_terms([{symbol=>"Acadm"}]);

$apph->disconnect;
stmt_ok;

