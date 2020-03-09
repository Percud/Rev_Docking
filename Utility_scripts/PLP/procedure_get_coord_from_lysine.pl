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
 printf "M: $model\n";
 open (pdb,"<$pdb_dir/$model") or die();
 my $found=0;
 while(my $pdb_line=<pdb>){
    if ($pdb_line=~/^ATOM.*NZ\s.*[A-Z]{3}\s[\s\S]\s*$K_resn\s+(\S+\s+\S+\s+\S+\s+)/){
        $coord{$model}="$1";
        printf "$model $K_resn $coord{$model}\n";
        $found=1;
    }
 }
 printf "#Coordinates not found for $model\n" if (!$found); 
 close(pdb);
 
}
close(fp);
