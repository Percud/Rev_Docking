#!/usr/bin/perl

use strict;

my $usage="$0 csv_file\n";

my($file)=(shift);
die ($usage) unless ($file);
open T, "<$file" or die $!;


foreach my $line(<T>){
  my @a = split /\t/, $line;
  my ($acc,$uniprot,$gene,$desc)=@a[0..3];
  my @h = grep /^http/, @a; #http addresses

  foreach (@h){
    my ($name,$smodel);
    my($from,$to,$template,$provider)=($1,$2,$3,$4) if (/from=(\d+).to=(\d+).template=(\S+).provider=(\S+)/);
    $name=sprintf "%s_%s_%s_%s_%s.pdb", $uniprot,$from,$to,$template,$provider;
    system ("wget --content-disposition \"$_\" -O $name\n");
    if ($provider eq "swissmodel"){
      open(S,"<$name") or die($!);
      my @F=<S>;
      my $f=join("\n",@F);
      my $GMQE=$1 if ($f=~/GMQE|QSPRD\s+(\S+)/m);
      my $QMN4=$1 if ($f=~/QMN4\s+(\S+)/m);
      my $OSTAT=$1 if ($f=~/OSTAT\s+(\S+)/m);
      $smodel="$GMQE\t$QMN4\t$OSTAT";
      }      
    $smodel="" if (!$smodel);
    printf "$uniprot\t$acc\t$gene\t$desc\t$name\t$smodel\n";
    }
 }

