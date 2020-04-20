#!/usr/bin/perl
# Docks a ligand into several pdb structures
# with autodock - 

use strict;
use File::Basename;
use File::Spec;

my $usage="$0 configuration_file\n";

if (@ARGV<1){
    print $usage;
    exit();
}

my($conf)=@ARGV;

my %par= (
	 MGL_ROOT=>'/usr/lib/mgltools_x86_64Linux2_1.5.6',
	 receptor_dir=>'./pdb',        
         ligand      =>'ligand.pdb',  
         coord_file  =>'coord_tab',   
         script_path =>'./script',
         gpf_par     =>'./reference.gpf',
         dpf_par     =>'./reference.dpf',
         outdir      =>'./outdir'
);

### read configuration_file
my $conf_log;
open(fp,"<$conf") or die();
while (<fp>){
    next if (/^\s*\#/);
    $conf_log.=$_;
    $par{$1}=$2 if (/^\s*(\S+)\s+(\S+)/);
    }
close fp;

#Set MGL_ROOT environmental variable
$ENV{'MGL_ROOT'}=$par{'MGL_ROOT'};


mkdir("$par{outdir}") or die "Can't mkdir $par{outdir} $!";

for my $k (keys %par){
    $par{$k}=File::Spec->rel2abs( $par{$k} ) ; #absolute_path
    die("$par{$k} not found") unless (-e $par{$k}); 
    }
   
my $opath=`pwd`;
chdir("$par{outdir}") or die "Can't chdir $par{outdir} $!";

printf "Results written in: $par{outdir}\nLog written in: $par{outdir}/out_log\n";
open (outlog,">out_log") or die(); #log file
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
printf outlog "command: %s %s\n",$0, join(' ',@ARGV);  
printf outlog "job started at %02d:%02d:%02d on %02d %02d %04d ", $hour,$min,$sec, $mday, $mon, $year+1900;  
printf outlog "\n\n#### configuration parameters ####\n%s\n",$conf_log;


my $script=$par{script_path};
############### calls to scripts #####
my $prep_ligand="$script/pythonsh $script/prepare_ligand4.py -A 'bonds_hydrogens' -U 'nphs' -l $par{ligand}";
my $prep_receptor="$script/pythonsh $script/prepare_receptor4.py -A 'checkhydrogens' -e ";
my $prep_gpf4="$script/pythonsh $script/prepare_gpf4.py -i reference.gpf";
my $prep_dpf4="$script/pythonsh $script/prepare_dpf42.py -i reference.dpf";
my $prep_flexreceptor4="$script/pythonsh $script/prep_flexreceptor4 -P CA_CB_CG_CD";
#######################################


############### read coordinates #####
open(fp,"<$par{coord_file}") or die();
my %coord;
while(my $line=<fp>){
    next if ($line=~/^#/); #skip comments
 
    my($name,$K_resn,$chain,$x,$y,$z)=split /\s+/, $line;
    push @{$coord{$name}}, "$x $y $z" if($name);
}

my $cmd;
####### prepare ligand if needed ######
printf outlog "\n\n#### ligand ####\n",$conf_log;
my ($ligand_name,$ligand_dir,$ligand_ext)=fileparse($par{ligand}, qr/\.[^.]*/);
    if ($ligand_ext eq ".pdbqt"){
	$cmd="cp $par{ligand} $par{outdir}"; 
	print outlog "\n$cmd\n";
        my $error=system ($cmd); if ($error) {  die "Copy failed: $!"};
        }
    else {
	$cmd="$prep_ligand -o $ligand_name.pdbqt";
	print outlog "\n$cmd\n";
         my $error=system ($cmd); if ($error) { die "failed: $!" }; 
        }    
    my $ligand="$ligand_name.pdbqt";

    
foreach my $pdb (sort keys %coord){
    for my $i (0 .. $#{ $coord{$pdb} }){
		
	########## prepare receptor ##########
    	my ($receptor_name,$receptor_ext)=($1,$2) if ($pdb=~/(\S+)\.(\S+)/);
    	my $receptor_name_chain=sprintf "%s.%s",$receptor_name,$i;
	my $receptor_chain="$receptor_name_chain.pdbqt";
	my $chain;
	my $K_resn;
	#my $receptor_chain_flex="$receptor_name_chain_flex.pdbqt";
	#my $receptor_chain_rigid="$receptor_name_chain_rigid.pdbqt";
    	printf outlog "\n\n#### processing: %s ####\n",$receptor_chain;
    	printf stderr "\n\n#### processing: %s ####\n",$receptor_chain;
    	next if (!-e  "$par{receptor_dir}/$pdb"); 
    	$cmd="$prep_receptor -r $par{receptor_dir}/$pdb -o $receptor_chain -e";
    	print outlog "\n$cmd\n";
    	my $error=system ($cmd); if ($error) { print outlog "**failed: $!"; warn "**failed: $!" };
	
	########## prepare flex_receptor ##########
    	$cmd="$prep_receptor -r $receptor_chain -s :$chain:LYS$K_resn -g $receptor_chain_rigid -x $receptor_chain_flex";

	##########   prepare gpf   ##########

	    $cmd="cp $par{gpf_par} ./reference.gpf";
	    if (-e $par{gpf_par}){
	      print outlog "\n$cmd\n";
	       my $error=system ($cmd); if ($error) { die "Copy failed: $!" };
	      }
	    open(fp,">>reference.gpf") or die();
	    printf fp "\ngridcenter %s\n", $coord{$pdb}[$i];
	    close fp;
	    $cmd="$prep_gpf4 -r $receptor_chain_rigid -l $ligand -o $receptor_name_chain$ligand_name.gpf -x $receptor_chain_flex";
	    print outlog "\n$cmd\n";
	    my $error=system ($cmd); if ($error) { print outlog "**failed: $!"; warn "**failed: $!" };        
	##########       autogrid4   ##########
	    $cmd="autogrid4 -p $receptor_name_chain$ligand_name.gpf -l $receptor_name_chain$ligand_name.glg";
	    print outlog "\n$cmd\n";
	    my $error=system ($cmd); if ($error) { print outlog "**failed: $!"; warn "**failed: $!" };  
	##########     prepare dpf   ##########
	    $cmd="cp $par{dpf_par} ./reference.dpf";
	    if (-e $par{dpf_par}){
	      print outlog "\n$cmd\n";
	       my $error=system ($cmd); if ($error) { die "Copy failed: $!" };
	      }
	    open(fp,">>reference.dpf") or die();
	    close fp;
	    $cmd="$prep_dpf4 -r $receptor_chain_rigid -x $receptor_chain_flex -l $ligand -o $receptor_name_chain$ligand_name.dpf";
	    print outlog "\n$cmd\n";
	    my $error=system ($cmd); if ($error) { print outlog "**failed: $!"; warn "**failed: $!"  };
	    $cmd="sed -i 's/\#/\ \#/g' $receptor_name_chain$ligand_name.dpf"; #to compensate for a prepare_dpf4 bug
	    my $error=system ($cmd); if ($error) { print outlog "**failed: $!"; warn "**failed: $!"  };
	##########       autodock4 [fork]  ##########
	    $cmd="autodock4 -p $receptor_name_chain$ligand_name.dpf -l $receptor_name_chain$ligand_name.dlg";
	    my $pid=fork();
	    if ($pid == -1) {
	       print outlog "**failed: $!"; warn "**failed: $!";
	   } elsif ($pid == 0) {
	      exec $cmd;
	   }
	    print outlog "\n$cmd\n";
     } #end for my $i (0 .. $#{ $coord{$pdb} })
     
#######################################
}
while (wait() != -1) {} # wait for children termination;
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
printf outlog "\n\ncommand: %s %s\n",$0, join(' ',@ARGV);  
printf outlog "job ended at %02d:%02d:%02d on %02d %02d %04d\n\n", $hour,$min,$sec, $mday, $mon, $year+1900;  
close outlog;
