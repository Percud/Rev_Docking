#!/usr/bin/perl
=head1 NAME

procedure_get_cood_from_resn.pl find coordinates for a specific pdb_id, chain, residue, and atom

=head1 description

procedure_get_cood_from_rens.pl resn_file pdb_files_directory_path

input_file format (rnsn_file)
#pdbid	chain	resn
Q8N5Z0_1_425_5eun_pdb   A       263
...
=head1
=cut

use strict;
use Getopt::Long;

my $usage="\n$0 -a [CA|CB|NZ..] resn_file pdb_files_directory_path\n 
extract coordinates for a specific pdbid chain, residue, and atom\n\n";

my $Atom="CA"; #default atom: c alpha 
GetOptions('a|atom=s' => \$Atom);

die ("$usage") if (@ARGV <=> 2);

my($resn_file,$pdb_dir)=@ARGV;

open(fp,"<$resn_file") or die();

while(my $line=<fp>){
  next if (/^\s*#/); #skip comments
  my ($model,$chain,$K_resn)=split /\s+/, $line;

  $model="$model.pdb";
  open (pdb,"<$pdb_dir/$model") or die();
  
  my $found=0;
  while(my $pdb_line=<pdb>){
    if ($pdb_line=~/^ATOM.*$Atom\s.*[A-Z]{3}\s$chain\s*$K_resn\s+(-*\d+\.\d+)\s*(-*\d+\.\d+)\s*(-*\d+\.\d+)/){
        printf "$model $K_resn $chain $1 $2 $3\n";
        $found=1;
    }
 }
 printf "#Coordinates not found for $model\n" if (!$found); 
 close(pdb);
}
close(fp);

exit 0;
