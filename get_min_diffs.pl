#!/usr/bin/perl
##Eric Smith

use strict;
use warnings;

my $base_dir = "/Research/postdoc/p_apium/CRISPR";
my $spacers = "${base_dir}/CRISPR_array_spacers.fasta";
my $spacers_out = "${base_dir}/all_spacers_closest_match.txt";
my %spacers;
my %species_min_diffs;

open SPACE, "$spacers";
while (my $line = <SPACE>) {
  chomp $line;
  my $space_name = substr($line, 1);
  my $spacer_seq = <SPACE>;
  chomp $spacer_seq;
  $spacers{$space_name}{seq} = $spacer_seq;
  $spacers{$space_name}{min_diff} = 99;
  $spacers{$space_name}{match} = "x";
  $spacers{$space_name}{match_rev} = "N";
  $spacers{$space_name}{per_match} = 0;
}
close SPACE;
foreach my $spacer_a (keys %spacers) {
  $spacer_a =~ /(\S+)_spacer/;
  my $species_a = $1;
  my $seq_a = $spacers{$spacer_a}{seq};
  unless ($species_min_diffs{$species_a}) {
    $species_min_diffs{$species_a}{min_diff} = 99;
    $species_min_diffs{$species_a}{percent_match} = 0; 
  }
  foreach my $spacer_b (keys %spacers) {
    $spacer_b =~ /(\S+)_spacer/;
    my $species_b = $1;
    next if $species_a eq $species_b;
    my $seq_b = $spacers{$spacer_b}{seq};
    my $min_diff = get_diffs($seq_a, $seq_b);
    $seq_b = reverse_comp($seq_b);
    my $rev_min_diff = get_diffs($seq_a, $seq_b);
    if ($min_diff < $spacers{$spacer_a}{min_diff}) {
      my $percent_match = ((length($seq_a) - $min_diff) / length($seq_a)) * 100;
      $spacers{$spacer_a}{min_diff} = $min_diff;
      $spacers{$spacer_a}{match} = $spacer_b;
      $spacers{$spacer_a}{match_rev} = "N";
      $spacers{$spacer_a}{per_match} = $percent_match;
      if ($min_diff < $species_min_diffs{$species_a}{min_diff}) {
        $species_min_diffs{$species_a}{min_diff} = $min_diff;
        $species_min_diffs{$species_a}{percent_match} = $percent_match;
      }
    } elsif ($rev_min_diff < $spacers{$spacer_a}{min_diff}) {
      my $percent_match = ((length($seq_a) - $rev_min_diff) / length($seq_a)) * 100;
      $spacers{$spacer_a}{min_diff} = $rev_min_diff;
      $spacers{$spacer_a}{match} = $spacer_b;
      $spacers{$spacer_a}{match_rev} = "Y";
      $spacers{$spacer_a}{per_match} = $percent_match;
      if ($rev_min_diff < $species_min_diffs{$species_a}{min_diff}) {
        $species_min_diffs{$species_a}{min_diff} = $rev_min_diff;
        $species_min_diffs{$species_a}{percent_match} = $percent_match;
      }
    }
  }
}

open OUT, ">$spacers_out";
foreach my $print_spacer (sort keys %spacers) {
  print OUT "$print_spacer\n\tclosest_match: $spacers{$print_spacer}{match}\n\tdifferences: $spacers{$print_spacer}{min_diff}\n\tpercent_match: $spacers{$print_spacer}{per_match}\n\treverse_comp: $spacers{$print_spacer}{match_rev}\n\n";
}
close OUT;

print STDOUT "Absolute minimum differences (percent match) across all spacers within each genome:\n";
foreach my $print_species (sort keys %species_min_diffs) {
  print STDOUT "${print_species}: $species_min_diffs{$print_species}{min_diff} ($species_min_diffs{$print_species}{percent_match})\n";
}

sub reverse_comp {
  my ($seq) = @_;
  $seq =~ tr/ACTGactg/TGACtgac/;
  $seq = reverse($seq);
  return $seq;
}

sub get_diffs {
  my ($seq_a, $seq_b) = @_;
  my $short_seq = $seq_a;
  my $long_seq = $seq_b;
  my $min_diff = 99;
  my $len = length($seq_a);
  if ($len > length($seq_b)) {
    $len = length($seq_b);
    $short_seq = $seq_b;
    $long_seq = $seq_a;
  }
  
  for (my $i = 0; $i < length($long_seq) - $len; $i++) {
    my $long_seq_sub = substr($long_seq, $i, $len);
    my $x = $long_seq_sub ^ $short_seq;
    my $diff = $x =~ tr/\0//c;
    $min_diff = $diff if $diff < $min_diff;
  }
  return $min_diff;
}