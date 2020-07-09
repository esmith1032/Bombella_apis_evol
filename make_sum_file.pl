#!/usr/bin/perl
##Eric Smith

use strict;
use warnings;

my $base_dir = "/Research/postdoc/p_apium/orthoMCL/glox_outgroup/all_genomes_new_paras_added_bombella";
my $int_prots = "${base_dir}/interesting_gene_types.txt";
my $prot_pos = "${base_dir}/para_sacc_prot_locs.txt";
my $chrom_info = "${base_dir}/chrom_len_offset.txt";
my $outfile = "${base_dir}/int_gene_locs.txt";
my %chroms;
my %prot_info;
my $needed_prots = '';

open CHROM, "$chrom_info";
while (my $line = <CHROM>) {
  chomp $line;
  my @values = split /\t/, $line;
  $chroms{$values[0]}{strand} = $values[1];
  $chroms{$values[0]}{length} = $values[2];
  $chroms{$values[0]}{offset} = $values[3];
}
close CHROM;

open INTS1, "$int_prots";
while (my $line = <INTS1>) {
  chomp $line;
  next if $line =~ /^species/;
  my @values = split /\t/, $line;
  $needed_prots .= "$values[1],";
}
close INTS1;

open PROTS, "$prot_pos";
while (my $line = <PROTS>) {
  chomp $line;
  my @values = split /\t/, $line;
  next unless $needed_prots =~ /$values[1],/;
  $prot_info{$values[1]}{chrom} = $values[2];
  $prot_info{$values[1]}{start} = $values[3];
  $prot_info{$values[1]}{end} = $values[4];
}
close PROTS;

open INTS2, "$int_prots";
open OUT, ">$outfile";
print OUT "species\tgene\ttype\tlocation\n";
while (my $line = <INTS2>) {
  chomp $line;
  next if $line =~ /^species/;
#  my ($species, $gene, $type) = split /\t/, $line;
  my @values = split /\t/, $line;
  my $chrom = $prot_info{$values[1]}{chrom};
  my $pos;
  if ($chroms{$chrom}{strand} eq "forward") {
    $pos = $prot_info{$values[1]}{start} + $chroms{$chrom}{offset};
  } else {
    my $new_start = $chroms{$chrom}{length} - $prot_info{$values[1]}{end};
    $pos = $new_start + $chroms{$chrom}{offset};
  }
  print OUT "$values[0]\t$values[1]\t$values[2]\t$pos\n";
}
close OUT;
close INTS2;