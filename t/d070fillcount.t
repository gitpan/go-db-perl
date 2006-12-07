#!/usr/local/bin/perl -w

# All tests must be run from the software directory;
# make sure we are getting the modules from here:
use lib '.';
use strict;
use GO::TestHarness;

n_tests(5);


use GO::Parser;

 # Get args

create_test_database("go_path");
my $apph = getapph();

my $parser = new GO::Parser ({handler=>'db'});
$parser->xslt('oboxml_to_godb_prestore');
$parser->handler->apph($apph);
stmt_ok;

$parser->parse ("./t/data/baby-function.dat");
$apph->add_root;
$parser->cache_errors(1);
$parser->parse_assocs("./t/data/mini-fb-assocs.fb");

stmt_ok;

# lets check we got stuff

$apph->fill_path_table;

#lets check count after filtering out certain relationship type
$apph->fill_count_table(undef, ['is_a']); #want only is_a
my $t1 = $apph->get_term({acc=>"GO:0030427"});
stmt_check($t1->n_deep_products == 2);

$apph->fill_count_table;

$t1 = $apph->get_term({acc=>"GO:0030427"});
stmt_check($t1->n_deep_products == 3);

my $g = $apph->get_graph(3673, -1);
#$g->to_text_output;
my $t= $g->get_term("GO:0003673");
my $pc = $t->n_deep_products;
stmt_note($pc);
my $pl = $apph->get_products({deep=>1, term=>$t});
stmt_note(scalar(@$pl));
stmt_check($pc == 563);
