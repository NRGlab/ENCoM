#!/usr/bin/perl

use strict;

my $foldx = "NA";
my @covf = ();
my $wtf = "NA";

my $factor = 3;

for my $i (0..@ARGV-1) {
	if ($ARGV[$i] eq "-FoldX") {$foldx = $ARGV[$i+1];}
	if ($ARGV[$i] eq "-mutl") {
		my $j = 0;
		while(1) {
			++$j;
			if ($ARGV[$j+$i] =~ m/-/) {last;}
			unless(-e $ARGV[$j+$i]) {last;}
			push(@covf,$ARGV[$i+$j]);
		}
	}
	if ($ARGV[$i] eq "-wt") {$wtf = $ARGV[$i+1];}
	if ($ARGV[$i] eq "-factor") {$factor = $ARGV[$i+1];}
}

print "FoldX Results File:$foldx\n";
print "Wild Type Covariance:$wtf\n";
print "Mutants Covariance:@covf\n";

# Load WT Covariance and data

my $wt_ener;
my %wtcov;
my %wtres;
open IN, "$wtf" or die "Cannot find WT covariance -> $wtf\n";

while(my $l = <IN>) {
	if ($l =~ m/Energy:(\S+)/) {$wt_ener = $1;}
	if ($l =~ m/COV\:\s+(\d+)\s+(\d+)\s+(\S+)\s+\S+\s+(\S+)/) {
		my $id = $1;
		my $id2 = $2;
		my $res = $3;
		my $val = $4;
		
		if ($id != $id2) {next;}
		$wtcov{$id} = $val;
		$wtres{$id} = $res;
		
	}
}
close IN;
if ($wt_ener eq undef) {warn "Cannot find WT energy\n";}

my %mut;
my %resn;
my %all_diff;
my @ener = ();
my @files = ();

my @abs_diff = ();

foreach my $f (@covf) {
	open IN, "$f" or die "Cannot find mutant covariance -> $f\n";
	my %temp_diff;
	while(my $l = <IN>) {
		if ($l =~ m/Energy:(\S+)/) {push(@ener,$1);}
		if ($l =~ m/File:(\S+)/) {push(@files,$1);}
		if ($l =~ m/COV\:\s+(\d+)\s+(\d+)\s+(\S+)\s+\S+\s+(\S+)/) {
			my $id = $1;
			my $id2 = $2;
			my $res = $3;
			my $val = $4;
		
			if ($id != $id2) {next;}
			$mut{$f}{$id} = $val;
			$resn{$f}{$id} = $res;
			push(@abs_diff,abs($wtcov{$id}-$val));
			#if ($wtcov{$id} eq undef) {next;}
			$temp_diff{$res} = $val-$wtcov{$id};
			$all_diff{$f}{$id} = $val-$wtcov{$id};
		}
	}
	close IN;
	
	# Look most pertub res
	
	my $c = 0;
	print "$f\n";
	print "	Most affected residue (Delta S to WT)\n";
	my $mean = 0;
	my $mean_abs = 0;
	foreach my $k (sort {abs($temp_diff{$b}) <=> abs($temp_diff{$a})} keys %temp_diff) {
		if ($c < 5) {printf "		%10s : %.4f\n",$k,$temp_diff{$k};}
		++$c;
		$mean += $temp_diff{$k};
		$mean_abs += abs($temp_diff{$k});
	}
	$mean /= $c;
	$mean_abs /= $c;
	my $sd = 0;
	my $sd_abs = 0;
	foreach my $k (sort {abs($temp_diff{$b}) <=> abs($temp_diff{$a})} keys %temp_diff) {
		$sd += ($mean - $temp_diff{$k})**2;
		$sd_abs += ($mean_abs - abs($temp_diff{$k}))**2;
	}
	printf "	Mean difference         :%.4f +/- %.4f\n",$mean,sqrt($sd/$c);
	printf "	Mean absolute difference:%.4f +/- %.4f\n",$mean_abs,sqrt($sd_abs/$c);
	
}

# Calculate maximum difference, mean, SD

@abs_diff = sort {$a <=> $b} @abs_diff;
my $mean = 0;
foreach my $val (@abs_diff) {$mean += $val;}
$mean /= scalar(@abs_diff);
my $sd = 0;
foreach my $val (@abs_diff) {$sd += ($mean - $val)**2;}
$sd = sqrt($sd/scalar(@abs_diff));

printf "\nAll absolute mean difference:%.4f +/- %.4f\n",$mean,$sd;
printf("Max absolute difference:%.4f\n",$abs_diff[scalar(@abs_diff)-1]);

my $max = $abs_diff[scalar(@abs_diff)-1];

print "The color values are scaled in fonction of $factor times the standard deviation or the maximum absolute difference, whichever is smaller\n";

if ($max > $sd * $factor) {
	printf "%.4f is bigger than $factor times the SD (%.4f), will use the sd for scaling\n",$max,$sd * $factor;
	$max = $sd * $factor;
}

open OUT, ">Diff.pml" or die;

for my $i (0..50) {
	my $id = sprintf("%.2f",$i/50);
	print OUT "set_color colordef$i = [$id,$id,1.00]\n";
}
for my $i (51..100) {
	my $id = sprintf("%.2f",1-($i-50)/50);
	print OUT "set_color colordef$i = [1.00,$id,$id]\n";
}

print "@files\n";
my @allname = ();
for my $i (0..@files-1) {
  #Get name???
  my $name = "";
  if ($files[$i] =~ m/([A-Za-z0-9-_]+)\./) {$name = "_$1";}
  push(@allname,"mutant_$i$name");
  print OUT "load $files[$i], mutant_$i$name\n";
}

print OUT "hide everything\n";
print OUT "show cartoon\n";

for my $i (0..@files-1) {
	
	foreach my $id (sort {$a <=> $b} keys %{$mut{$covf[$i]}}) {
		my $val = $all_diff{$covf[$i]}{$id};
		my $res = $resn{$covf[$i]}{$id};
		my $diff = 0;
		if ($val < 0) {
		#	print "$diff = int(-$val/$max*50)\n";
			$diff = int($val/$max*50)+50;
		} else {
			$diff = int($val/$max*50)+50;
		}
		
		if ($diff > 100) {$diff = 100;}
		if ($diff < 0) {$diff = 0;}
		if ($res =~ m/^(\D{3})(\d+)(\D)/) {
			my $buf = "\/\/\/$3\/$2\/ & $allname[$i]";
			printf OUT "color colordef%d,$buf\n",$diff;
			if ($res ne $wtres{$id}) {
				printf OUT "show stick, $buf\n";
			}
			#print "$res $val $diff\n";
		} else {
			print "Failed to parse this residue:$res\n";
		}
	}
}



open IN, "$foldx";
my $record = -1;
my @foldx;
while(my $l = <IN>) {
	
	if ($l =~ m/Pdb	total energy	Backbone Hbond	Sidechain Hbond	Van der Waals/) {
		++$record;
		next;
	}
	if ($record == -1) {next;}
	#print "$l";
	my @sp = split(/\s+/,$l);
	push(@foldx,$sp[1]);
	#printf "%.4f %.4f %.4f\n",$ener[$record]-$wt_ener,$sp[1],($ener[$record]-$wt_ener)*-0.8 + $sp[1]*0.35;
	++$record;
}
close IN;
print "Prediction of DDG\n";
print "ENCoM_DS FoldX Combined\n";
if (scalar(@foldx) != scalar(@ener)) {print "******************************\nWarning no Value for Foldx -> All set to 0\n******************************\n";}
for my $i (0..@files-1) {
  
  printf "$files[$i] %.4f %.4f %.4f\n",$ener[$i]-$wt_ener,$foldx[$i],($ener[$record]-$wt_ener)*-0.8 + $foldx[$i]*0.35;

}
















