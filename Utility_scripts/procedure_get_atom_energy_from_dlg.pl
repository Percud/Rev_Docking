#!/usr/bin/perl
#
# get vdW and electrostatic energy from dlg 
#

use strict;

my $usage="$0 ATOM_RES directory_dlg_file\n";

my($atom_id,$d)=(shift,shift);
die ($usage) unless ($d);

opendir(DIR,$d) or die;
my @F = grep /\.dlg/, readdir (DIR);
closedir(DIR);


printf "Name\tBest\tLARGE\n";


foreach my $f (@F) {
	my ($clu,@A);
	my ($large_rank,$large_num,@large_runs);
    my %E;    
    open T, "<$d/$f" or die $!;
    foreach (<T>){
    	#CLUSTERING HISTOGRAM
        if (/^\s+(\d+).*\|.*\|\s+(\d+)\s+.*\s+(\S+)\s+\|\#/){
            my ($rank,$run,$num)=($1,$2,$3);
            ($large_rank,$large_num)=($rank,$num) if (!defined $large_num or $num > $large_num);
        }

	#LOWEST ENERGY DOCKED CONFORMATION from EACH CLUSTER
        $clu=$1 if (/Rank = (\d+)/);
        my @A=split /\s+/ if (/^ATOM.*$atom_id/);
        $E{$clu}=0 if (!defined $E{$clu});
	#Sum vdW and Elec energies for ATOM of ATOM_RES
        $E{$clu}+=$A[-4]+$A[-3] if (defined $A[0]);    
    }     
    printf "$f\t%s\t%s\t(%s)\n", $E{1},$E{$large_rank},$large_num;
}


exit;
