#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my $file = "ExecutionFileDiversity";
my $fout;
open($fout, '>' ,$file);
my $Instance=0;
my $Sed=0;
my $PATH=  `cd ..; pwd`;
chomp $PATH;
#my @configuration2 = split ' ', $configuration;#~ s/ /_/g; 

#foreach my $time(("1440", "2880")){
foreach my $time(("10")){
  foreach my $days(("20", "40", "60")){
   foreach my $px(("0.2", "0.6")){
     foreach my $ITV(("0.2", "0.4")){
         for(my $seed=1; $seed<=35; $seed++){
           print $fout "~/$PATH/main ".$PATH."/entradas/Platillos.csv $PATH/entradas/Restricciones.csv $days $time $PATH/salidas/dias_".$days."_minutos_".$time."_seed_".$seed.".csv ".$seed." $PATH/History/populationAndFitness_days_".$days."_min_".$time."_seed_".$seed."_Px_".$px."_".$ITV." $px $ITV";
           print $fout "\n";
         }
       }
     }
  }
}
