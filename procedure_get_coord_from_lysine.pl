#!/usr/bin/perl
#Find coordinates for PLP-bound lysine in a PDB file
#Read Orf_Pred from identificazione_genomi_hmmsearch.pl (input file from pdb2fasta)
use strict;

my $usage="$0 Orf_Pred_file pdb_file_directory_path\n    extract coordinates of catalytic lysine\n";

if (@ARGV<2){
    print $usage;
    exit();
}


my($tab_file,$pdb_dir)=@ARGV;

open(fp,"<$tab_file") or die();
my %coord;
while(my $line=<fp>){
 my @f=split /\s+/, $line;
 my $model=$f[6];
 my $K_resn=$f[-2]; #penultimate field in tab rows
 $model="$model.pdb";
 
 open (pdb,"<$pdb_dir/$model") or die();
 while(my $pdb_line=<pdb>){
    $coord{$model}="$1" if ($pdb_line=~/^ATOM.*NZ\s.*[A-Z]{3}\s[\sA]\s*$K_resn\s+(\S+\s+\S+\s+\S+\s+)/);
 }
  printf "$model $K_resn $coord{$model}\n";
 close(pdb);
 
}
close(fp);
