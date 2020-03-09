#!/usr/bin/perl
#from the internet pdb2fasta.pl - Chengxin Zhang pdb2fasta in perl
#modified to add missing residues
use strict;
use File::Basename;

if (@ARGV<1){
    print "pdb2fasta.pl pdb.pdb > seq.fasta\n";
    print "    convert PDB file pdb.pdb to FASTA sequence file seq.fasta\n";
    exit();
}

my %aa3to1=(
   'ALA'=>'A', 'VAL'=>'V', 'PHE'=>'F', 'PRO'=>'P', 'MET'=>'M',
   'ILE'=>'I', 'LEU'=>'L', 'ASP'=>'D', 'GLU'=>'E', 'LYS'=>'K',
   'ARG'=>'R', 'SER'=>'S', 'THR'=>'T', 'TYR'=>'Y', 'HIS'=>'H',
   'CYS'=>'C', 'ASN'=>'N', 'GLN'=>'Q', 'TRP'=>'W', 'GLY'=>'G',
   'MSE'=>'M', 'LLP' => 'K',
);

foreach my $pdb_file(@ARGV){
    my $filename=basename("$pdb_file");
    if ($filename=~/(\S+)\.\S*/){$filename="$1";}
    my %chain_hash;

    open(fp,"<$pdb_file");
    my %resn; #residue number;
    my %presn; #previous resnum;
    while(my $line=<fp>){
        chomp($line);
	last if ($line=~/^ENDMDL/); # just the first model
        if ($line=~/^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])\s*(\d{1,4})/
         or $line=~/^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE|LLP)\s([\s\w])\s*(\d{1,4})/){
            my $res="$1";
            my $chain="$2";
            $resn{$chain}="$3";
            
            #add missing residues as x
            my $diff=$resn{$chain}-$presn{$chain}-1;
            if($diff){my $missing='x' x $diff; $chain_hash{$chain}="$chain_hash{$chain}$missing"};
            
            $chain_hash{$chain}="$chain_hash{$chain}$aa3to1{$res}";
            
            $presn{$chain}=$resn{$chain};
        }
    }
    close(fp);

    foreach my $chain(sort keys %chain_hash){
        printf ">$filename:$chain\n$chain_hash{$chain}\n";
    }
}
