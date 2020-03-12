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


printf "Name\tBest (run)\tLarge (run)\tBest_Large (run)\t(num)\n";


foreach my $f (@F) {
	my ($run);
	my ($large_rank,$large_num);
	my ($best_run, $best_large_run, $best_large_run_atom);
    my %E; #energy for each run
    my %C; #cluster for each run
    open T, "<$d/$f" or die $!;
    foreach (<T>){
    	#FINAL GENETIC ALGORITHM DOCKED STATE
        $run=$1 if (/Run = (\d+)/);
        if (/^DOCKED: ATOM.*$atom_id/){
		my @A=split /\s+/ ;
        	$E{$run}=0 if (!defined $E{$run});
		#Sum vdW and Elec energies for ATOM of ATOM_RES
        	$E{$run}+=$A[-4]+$A[-3] if (defined $A[0]);
		}
	
    	#CLUSTERING HISTOGRAM
        if (/^\s+(\d+).*\|.*\|\s+(\d+)\s+.*\s+(\S+)\s+\|\#/){
            my ($rank,$run,$num)=($1,$2,$3);
	    $best_run=$run if (!defined $best_run); #the first listed is best cluster 
            ($large_rank,$best_large_run,$large_num)=($rank,$run,$num) if (!defined $large_num or $num > $large_num);
        }
	#RMSD TABLE - cluster for each run
	if (/^\s+(\d+)\s+\d+\s+(\d+).*RANKING/){
	$best_large_run_atom=$best_large_run if (!defined $best_large_run_atom);
		my ($c,$r) = ($1,$2);
		$best_large_run_atom = $r if ($c==$large_rank and $E{$r}<$E{$best_large_run_atom});
	}
	last if (/INFORMATION ENTROPY ANALYSIS FOR THIS CLUSTERING/);

    }
        	
    printf "%-50s\t%s ($best_run)\t%s ($best_large_run)\t%s ($best_large_run_atom)\t#%s\n", $f,$E{$best_run},$E{$best_large_run},$E{$best_large_run_atom},$large_num;
}


exit;
