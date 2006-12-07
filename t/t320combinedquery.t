#!/usr/local/bin/perl -w

# All tests must be run from the software directory;
# make sure we are getting the modules from here:

use lib '.';
use strict;
use Test;
BEGIN { plan tests => 1 }
use GO::TestHarness;
use GO::AppHandle;

# ----- REQUIREMENTS -----

# we want to be able to do combined queries;
# eg fetch me every product that is
# transmembrane receptor and membrane

# ------------------------

my $apph = get_readonly_apph(@ARGV);
my $pl = $apph->get_deep_products({terms=>[4888,5615], operator=>"and"});
foreach my $p (@$pl) {
    stmt_note($p->symbol);
}
my $pl1 = $apph->get_deep_products({terms=>[4888]});
my %ph1 = map {$_->id => $_ } @$pl1;
my $p12 = $apph->get_deep_products({terms=>[5615]});
my %ph2 = map {$_->id => $_ } @$pl1;
stmt_check(!
           grep {
               !($ph1{$_->id} && $ph2{$_->id})
           } @$pl
          );
