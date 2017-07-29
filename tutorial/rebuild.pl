#!/usr/bin/perl


use strict;
my $ca_flag = 0;

#Search for modeller version
my $pdb = "NA";
my $script = "NA";
for my $i (0..@ARGV-1) {
	if ($ARGV[$i] eq "-i") {$pdb = $ARGV[$i+1];}
	if ($ARGV[$i] eq "-script") {$script = $ARGV[$i+1];}
}



my $modeller = "9.14";

for my $i (8..16) {
	my $out = `whereis mod9.$i`;
	#print "$out\n";
	if ($out =~ m/(mod9\.\d+)/) {
		#print "I detech this version of Modeller -> $1\n";
		$modeller = $1;
	}

}

open IN, "$pdb" or die;
my $buf = "";
my $count = 0;
my @in = <IN>;

write_align();

mkdir "models";
foreach my $l (@in) {
	if ($l =~m/^Model\s*\d+/) {$buf = "";print "$l";next;}
#	print "$l\n";
	if ($ca_flag == 1) {
		if ($l =~ m/^ATOM.........CA/ ) {$buf .= "$l";}
	} else {
		$buf .= "$l";
	}
	if ($l =~ m/^ENDMDL/) {
		
		++$count;
		
		open OUT, ">./temp.pdb" or die;
		print OUT "$buf";
		close OUT;
		
		$buf = "";
		
		my @seq = @{get_seq("./temp.pdb")};
		my $convert = convert_seq(\@seq);
		
		my $out = `$modeller ./align_seq.py temp.pdb temp`;
		
		open SEQ,"all_aligne.pir" or die;
		my @pir = <SEQ>;
		close SEQ;
		open SEQ,">all_aligne.pir" or die;
		foreach my $iseq (0..@pir-1) {
			my $lseq = $pir[$iseq];
			if ($lseq =~ m/^>/) {
				print SEQ "$lseq";
				++$iseq;
				$lseq = $pir[$iseq];
				print SEQ "$lseq";
				print SEQ "$convert*\n";
			} else {
				next;
			}
		}
		close SEQ;
		system "$modeller $script temp.pdb temp";
		if (-e "tobuild.B99990001.pdb") {} else {
			print LOG "Failed to model\n";
			next;
		}
		rename_seq("tobuild.B99990001.pdb",\@seq);
		system "rm tobuild.B99990001.pdb";
		system "mv fix_seq.pdb ./models/model_$count.pdb";
	}
}
system "rm align_seq.py";

sub convert_seq() {
	my @seq = @{$_[0]};
	my $convert = "";
	my @last = ();
	my $c = -1;
	foreach my $lseq (@seq) {
		++$c;
		my @sp = split(/\s+/,$lseq);
		printf "%d -> @sp et @last\n",abs($sp[0] - $last[0]);
		if ($c != 0 && $sp[1] eq $last[1] && abs($sp[0] - $last[0]) > 1) {
			$convert .= "/";
		}
		if ($c != 0 && $sp[1] ne $last[1]) {$convert .= "/";}
		$convert .= mod_aa($sp[2]);		
		@last = @sp;

	}
	
	return($convert);
}

sub get_seq() {
  my @seq = ();

  open IN, "$_[0]";
  my $ener;
  while(<IN>){
    my $line = $_;
    if ($line =~ m/REMARK ENCoM Generated energie=(.*)/) {
      my $temp = $1;
      if ($ener eq undef) {$ener = $temp;}
    }

    #	ATOM    880  CA  THR A 131
    if ($line =~ m/^ATOM.........CA..(...).(.)\s+(\d+)/) {
    #  printf "%d $line",scalar(@seq);
      my $ami = $1;
      my $num = $3;
      my $cha = $2;

      if ($cha eq " ") {$cha = "Z";}
      $num =~ s/\s+//;
      #$ami = mod_aa($ami);
      push(@seq,"$num $cha $ami");

      } else {
      #	print "$line";
    }
  }
  close IN;




  return(\@seq);
}

sub mod_aa_reverse() {
  my $d = $_[0];

  if ($d eq "K") {return("LYS");}
  if ($d eq "R") {return("ARG");}
  if ($d eq "F") {return("PHE");}
  if ($d eq "L") {return("LEU");}
  if ($d eq "I") {return("ILE");}
  if ($d eq "Y") {return("TYR");}
  if ($d eq "E") {return("GLU");}
  if ($d eq "D") {return("ASP");}
  if ($d eq "G") {return("GLY");}
  if ($d eq "W") {return("TRP");}
  if ($d eq "A") {return("ALA");}
  if ($d eq "V") {return("VAL");}
  if ($d eq "T") {return("THR");}
  if ($d eq "H") {return("HIS");}
  if ($d eq "C") {return("CYS");}
  if ($d eq "S") {return("SER");}
  if ($d eq "N") {return("ASN");}
  if ($d eq "Q") {return("GLN");}
  if ($d eq "P") {return("PRO");}
  if ($d eq "M") {return("MET");}

  print "Cannot find:-$d-\n";
  return("X");
}

sub mod_aa() {
  my $d = $_[0];

  if ($d eq "LYS") {return("K");}
  if ($d eq "ARG") {return("R");}
  if ($d eq "PHE") {return("F");}
  if ($d eq "LEU") {return("L");}
  if ($d eq "ILE") {return("I");}
  if ($d eq "TYR") {return("Y");}
  if ($d eq "GLU") {return("E");}
  if ($d eq "ASP") {return("D");}
  if ($d eq "GLY") {return("G");}
  if ($d eq "TRP") {return("W");}
  if ($d eq "ALA") {return("A");}
  if ($d eq "VAL") {return("V");}
  if ($d eq "THR") {return("T");}
  if ($d eq "HIS") {return("H");}
  if ($d eq "CYS") {return("C");}
  if ($d eq "SER") {return("S");}
  if ($d eq "ASN") {return("N");}
  if ($d eq "GLN") {return("Q");}
  if ($d eq "PRO") {return("P");}
  if ($d eq "MET") {return("M");}

  print "Cannot convert amino acid to one letter Code, check it -> $d\n";
  return("X");
}

sub rename_seq() {
	my $file = $_[0];
	my @seq = @{$_[1]};
	
	open PDB, "$file" or die;
	open OUT, ">fix_seq.pdb" or die;
	my $last = "UNK";
	my $count = -1;
	while(my $l = <PDB>) {
	#	print "$l";
		if ($l =~ m/^ATOM.............(...).(.)(....)/) {
			my $aa = $1;
			my $ch = $2;
			my $num = $3;
			$num =~ s/\s+//gi;
			my $now = "$aa $ch $num";
			
			if ($now ne $last) {++$count;}
			my @sp = split(/\s+/,$seq[$count]);
			$sp[0] = sprintf("%4d",$sp[0]);
			$l =~ s/^(ATOM.............)(...)(.)(.)(....)/$1$sp[2]$3$sp[1]$sp[0]/;
			#print "\t$count:$now -> $seq[$count]\n";
			print OUT "$l";
			if ($seq[$count] eq undef) {last;}
			$last = $now;
		}
	}
	close PDB;
	close OUT;
	

}


sub write_align() {


	open ALN, ">align_seq.py";
	print ALN "from modeller import *
from modeller.automodel import *

import sys

prot = sys.argv[1]
prot_file = sys.argv[2]
env = environ()
aln = alignment(env)
env.io.hetatm = True
mdl = model(env, file=prot)
aln.append_model(mdl, align_codes=prot, atom_files=prot_file)
aln.append_model(mdl, align_codes='tobuild', atom_files=prot_file)
aln.align()
aln.write(file='all_aligne.pir', alignment_format='PIR')\n";
close ALN;


}
