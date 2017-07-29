#!/usr/bin/perl

use strict;

my $init = "UNDEF";
my $target = "UNDEF";
my $chain = "all";
my $remark = 0;
my $order = 0;
my $renum = 0;
my $break = 1;
for my $i (0..@ARGV-1) {
	if ($ARGV[$i] eq "-i") {$init = $ARGV[$i+1];}
	if ($ARGV[$i] eq "-t") {$target = $ARGV[$i+1];}
	if ($ARGV[$i] eq "-c") {$chain = $ARGV[$i+1];}
	if ($ARGV[$i] eq "-r") {$remark = 1;}
	if ($ARGV[$i] eq "-or") {$order = 1;}
	if ($ARGV[$i] eq "-renum") {$renum = 1;}
	

}

unless (-e $init) {
	print "Please specify a PDB file to clean\n";
	exit;
}

if ($target eq "UNDEF") {
	print "Please specify a destination name for the cleaned PDB\n";
	exit;
}

my @rstr = ();
if ($remark == 1) {@rstr = get_remark($init);}
my $random = int(rand(10000));
system "cp $init ./$random.pdb";
if ($renum == 1) {
  renum("$random.pdb","$target");
  system "mv $target $random\.pdb";
}

clean("$random\.pdb",$target,$chain);

system "cp $target $random.pdb";
fix_insertion("$random.pdb","$target");
system "rm $random.pdb";
if ($order == 1) {
  system "cp $target $random.pdb";
  reorder_chain("$random.pdb","$target");
  system "rm $random.pdb";
}

if ($break == 1) {
  system "cp $target $random.pdb";
  fix_chain_break("$random.pdb","$target");
  system "rm $random.pdb";
}


if ($remark == 1) {
  system "cp $target $random.pdb";
  system "rm $target";
  open OUT, ">$target" or die;
  print OUT "$rstr[0]";
  close OUT;
  system "cat $random.pdb >> $target";
  open OUT, ">>$target" or die;
  print OUT "$rstr[1]";
  close OUT;
  system "rm $random.pdb";

}

sub fix_chain_break() {
    my $init = $_[0];
    my $targ = $_[1];
    my %atm = load_pdb_atm($init,1,"NA");
    my %seq = load_pdb_seq($init,1,"NA");
    my @keys = sort {$seq{$a} <=> $seq{$b}} keys %seq;
    open OUT, ">$targ" or die;
    my $add = 0;
    my @spkl;
    foreach my $i (0..@keys-1) {
        my $k = $keys[$i];
        my @spk = split(/\s+/,$k);
        if ($i > 0) {
            # If chain change, reset the add
            if (($spk[2] ne $spkl[2])) {
                $add = 0;
            } else {
                # if they are residue that are next to each other
                if (abs($spk[1] - $spkl[1]) == 1) {
                    my $kl = $keys[$i-1];
                    my $dist = min_dist($atm{$k},$atm{$kl}," CA ");
                    #$dist = sqrt($dist);
                    if ($dist > 16) {
                        print "I found this chain Break -> $k $kl $dist\n";
                        $add += 1;
                    }
                }
            }
        }
        @spkl = @spk;
        
        #Apply the add
        foreach my $l (split(/\n/,$atm{$k})) {
            if ($l =~ m/^.................(....)(.)(....)(.)/) {
              my $num = $3;
			  my $strc = sprintf("%4d",$num+$add);
			  $l =~ s/^(......................)(.....)/$1$strc /;
			  print OUT "$l\n";
		  }
        }
    }
}

sub reorder_chain() {
  open IN, "$_[0]";
	open OUT, ">$_[1]";
	my @toprint = <IN>;
  close IN;
  my %chain = get_chain($_[0]);
  my @sortchain = ("A".."Z");
#  print "@sortchain\n";
  my %nchain; 
   my $counter = 0;
    foreach my $ch (sort {$chain{$b} <=> $chain{$a}} keys  %chain) {
    	print "	I have $chain{$ch } atoms in chain $ch  -> and now it will be chain $sortchain[$counter]\n";
    	$nchain{$ch} = $sortchain[$counter];
    	++$counter;
    }
    my $num = 1;
  foreach my $ch (sort {$chain{$b} <=> $chain{$a}} keys  %chain) {
    foreach my $i (0..@toprint-1) {
      my $l = $toprint[$i];
      if ($l =~ m/^.....................($ch)/ ) {} else {next;}
    	 
    	
    	if ($l =~ m/^.....................(.)/ ) {
    	  my $tch = $1;
    	  $l =~ s/^(.....................)(.)/$1$nchain{$tch}/;
    	}
    	my $snum = sprintf("%5d",$num);
    	 $l =~ s/^(......)(.....)/$1$snum/;
    	 #print "$l";
    	print OUT "$l";
    	++$num;
    	#die;
    }
  }
}

sub renum() {
    open IN, "$_[0]";
	open OUT, ">$_[1]";
	my $c = 0;
	my $last = -1;
	my $lastch = "NA";
	while(my $l = <IN>) {
		if ($l =~ m/^ATOM/ || $l =~ m/^HETATM/) {
			if ($l =~ m/^.................(....)(.)(....)(.)/) {
			  #print "$1 $2 $3 $4  -> $c\n";
			  if ($2 ne $lastch) {$c= 0;}
			  $lastch = $2;
			  my $now = "$1 $2 $3 $4";
			  if ($now ne $last) {++$c;}
			  my $strc = sprintf("%4d",$c);
			  $l =~ s/^(......................)(.....)/$1$strc /;
			  $last = $now;
			  print OUT "$l";
			}
	  }
	}
	close IN;
  close OUT;
    print "I'm done re-numbering\n";
}

sub get_chain() {
  my $het = $_[1];
 	my $mask = $_[2];
  open FILE, "$_[0]" or die;

  my %count;
  my $c = 0;
  # Get ResnumC

  while(my $l = <FILE>) {
          unless($l =~ m/$mask/) {next;}
          if ($het == 0 && $l =~ m/^HET/) {next;}
          if ($l =~ m/^.................(...) (.)(....)/) {
                  if ($2 eq " ") {$2 = "-";}

                  

                  $count{$2} += 1;
                  ++$c;
								#	print "$l";


          }
  }
  close FILE;
  return(%count);
}

sub fix_insertion() {
	open IN, "$_[0]";
	open OUT, ">$_[1]";
	my @in = <IN>;
  close IN;
  
  # Get Resnumc
 	my %uni;
 	my %replace;
 	my %toadd;
 	my $last = "NA";
 	my $lastnum = "NA";
  foreach my $l (@in) {
 #   print "INI:$l";
    
    # If last resnum is different and was a insert, add to chain
    
    # En fait, si resnumcALT est different mais num est le meme, ++num
    if ($l =~ m/^ATOM.............(....)(.)(....)(.)/) {
       my $all = "$1$2$3$4";
       my $ch = $2;
       my $num = $3;
       if ($last ne "NA" && $all ne $last && $num eq $lastnum) {
        
        $toadd{$ch} += 1;
        print "\t\tLAST:$last All:$all Ins:$lastnum Ch:$ch ToAdd:$toadd{$ch}\n";
      }
    }
    # Get Latest
     if ($l =~ m/^ATOM.............(....)(.)(....)(.)/) {
      my $all = "$1$2$3$4";
       my $ins = $4;
      $last = $all;
      $lastnum = $3;
     }
    # Add everything that need to be 
    
    if ($l =~ m/^ATOM.................(.)(\s*\d+)(.)/) {
     my $num = $2;
     my $ch = $1;
    
     my $new_num = sprintf("%4d",$num + $toadd{$ch});
  #   print "old:$num New:$new_num Toadd:\$toadd{$ch} = $toadd{$ch}\n";
     $l =~ s/(^ATOM..................)(\s*\d+)/$1$new_num/;
    }
    # Find insert et trouve un remplacement
    if ($l =~ m/^ATOM.............(.........)(.)/) {
      my $resnumc = $1;
      my $ins = $2;
      my $all = "$1$2";
      if ($ins ne " " && $replace{$all} eq undef) {
        $replace{$all} = "$resnumc ";
      }
      if ($replace{$all} ne undef) {
        $l =~ s/$all/$replace{$all}/;
      }
    }
    
    
    
  #  print "FIN:$l"; 
    print OUT "$l";
	}
	

}







sub clean() {
 # Script perl qui enleve les HETATM et les acides aminée exotiques 1X

  open IN, "$_[0]"or die "Cannot open $_[0]";

  open OUT, ">$_[1]" or die "Cannot open $_[1]";
  my @in = <IN>;
  close IN;
	my $kchain = $_[2];
	print "I keep $kchain chain(s)\n";
	my %ctk; # Chain to keep
	unless($kchain eq "all") {
		foreach my $tc (split(/,/,$kchain)) {
			if ($tc eq undef) {next;}
			if (length($tc) != 1) {next;}
			$ctk{$tc} += 1;
		}
	}
	#die;
  my %chain;

  my $done = 0;
  my $long = 0;
  my $log = "";
  my @rem = (0,0,0,0);
  my %dupli;
  my @toprint = ();
  my %cresnumc;
  my %cresnumc_ca;
  
  #Negative number
  
  my %toadd;
  
  my %trans;
  
  foreach my $l (@in) {
		if ($l =~ m/^ATOM/) {++$long;}
		
		# Transforme SCH to cysteine
		if ($l =~ m/^...........(....)..(SCH......)/) {
			my $atom = $1;
			my $trans_resnumc = $2;
			if ($trans{$trans_resnumc} eq undef) {
				print "Converting SCH to CYS -> $trans_resnumc\n";
				++$trans{$trans_resnumc};
			}
			if ($atom eq "  N " || $atom eq "  CA" || $atom eq "  CB" || $atom eq "  SG" || $atom eq "  C " || $atom eq "  O ") {
			
				
				$l =~ s/HETATM/ATOM  /;
				$l =~ s/SCH/CYS/;
				
				
			} else {next;}
			
			#die;
		}  
		# Transforme MSE to MET
		if ($l =~ m/^...........(....)..(MSE......)/ ) {
			my $atom = $1;
			my $trans_resnumc = $2;
			if ($trans{$trans_resnumc} eq undef) {
				print "Converting MSE to MET -> $trans_resnumc\n";
				++$trans{$trans_resnumc};
			}
			$l =~ s/HETATM/ATOM  /;
			$l =~ s/MSE/MET/;
			$l =~ s/SE / SD/;
		}	
		
    if ($l =~ m/^HETATM/) {
    	++$long;
   #   print "$l";
      ++$rem[0];
      next;
    }
  
    
    if ($l =~ m/^ATOM.........H......../) {
    	++$rem[1];
      next;
    }
    if ($l =~ m/H\s*$/) {
      #print "$l";
      ++$rem[1];
      next;
    }
    
    # Add a chain
    if ($l =~ s/^(ATOM.................)( )/$1Z/) {
    }
  	my $nchain = " ";
   if ($l =~ m/^ATOM.......(.....).(.........)(.)/) {
    	my $tresnumc = $2;
    	my $insert = $3;
    	if ($tresnumc =~ m/(...).(.)(.*)/) {
       	my $tchain = $2;
       	$nchain = $2;
       	my $num = $3;
       	$num =~ s/\s+//gi;
       	if ($num < 0) {
       		if ($toadd{$tchain} eq undef) {
       		  printf "Found negative residu numbering in chain $tchain, will add %d to all residue number\n",-$num;
       		  $toadd{$tchain} = -$num;
       		 }
       		
       	}
       	
       	my $new_num = sprintf("%4d",$num + $toadd{$tchain});
      # 	print "Old:$num New:$new_num \$toadd{$tchain} = $toadd{$tchain} ->   $tresnumc \n";
       	$l =~ s/$tchain\s*$num/$tchain$new_num/;
       }
   
   } 
    my $duplik = "NA";
    my $resnumc = "NA";
   
    if ($l =~ m/^ATOM/) {
        if ($l =~ m/^ATOM.......(.....).(..........)/) {
            $duplik = "$1 $2";
            $resnumc = $2;
          
           ++$cresnumc{$resnumc};
         } else {print "$l\n";next;}
    }
    
    if ($dupli{$duplik} ne undef && $duplik ne "NA") {print "Alternate conformations $duplik\n";++$rem[2];next;} else {$dupli{$duplik} = 1;}
    $l =~ s/^(ATOM............)\S/$1 /;
   
    
    
   
    if ($l =~ m/^MODEL/) {
      #print "$l";
    }
     if ($done > 100 && $l =~ m/^MODEL\s+\d+/) {
    	$log .= "	I have multiple model, I kept only the first one !\n";
    	last;
    }
    unless($l =~ m/^ATOM/) {next;}
    
    # Regarder s'il s'agit d'un des 20 aa
    if ($l =~ m/^ATOM.............(...)/) {
    	my $taa = $1;
    	my $caa = mod_aa($taa);
    	if ($caa eq "X") {next;}
    }
    
    
    if ($l =~ m/^.....................(.)/ ) {
    	my $tc = $1;
    	$tc =~ s/\s+//;
      $chain{$tc} += 1;
     
      if ($kchain ne "all") {
      	if ($ctk{$tc} eq undef) {$rem[3] += 1;next;}
      }
    }
    ++$done;
   	push(@toprint,$l);
   

  }
  foreach my $l (@toprint) {
  	my $resnumc;
  	if ($l =~ m/^ATOM.......(.....).(..........)/) {
  		 $resnumc = $2;
  		 if ($cresnumc{$resnumc} < 4) {
  		 	print "$l";
  		 	print "Removed this residue because there is not enough atoms -> $resnumc = $cresnumc{$resnumc} atom(s)\n";
  		# 	next;
  		 }
  	}
  	print OUT "$l";
  }
  close OUT;
  my $tot = 0;
  my @all_chain = keys  %chain;
  #$log .= "I have all those chain:@all_chain\n";
  foreach(keys  %chain) {
    ++$tot;
  }
  if ($tot != 1) {
   # print "$done/$long\n";
    foreach(keys  %chain) {
    	$log .= "	I have $chain{$_} atoms in chain $_\n";
    }

  }
  $log .= "	I had initially $long atoms, kept only $done, because I removed $rem[0] HETATM, $rem[1] Hydrogen Atoms, $rem[2] alternative conformations and $rem[3] chain atom were not kept\n";
  print "$log\n";
  #die;
	return($log)
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

sub get_remark() {
  open IN, "$_[0]" or die;
  my @str = ();
  my $id = 0;
  while (my $l = <IN>) {
    if ($l =~ m/^ATOM/ || $l =~ m/^HETA/) {$id = 1;next;}
    
      $str[$id] .= "$l";
  }
  close IN;
  return(@str);
}




sub load_pdb_seq {
	my $het = $_[1];
	my $mask = $_[2];
	open FILE, "$_[0]" or die "Can't open $_[0]";
	if ($mask eq undef) {$mask = ".";}
	my %count;
	my $c = 0;
	# Get ResnumC
	
	while(my $l = <FILE>) {
		if ($het == 0 && $l =~ m/^HET/) {next;}
		if ($mask ne "NA") {unless($l =~ m/$mask/) {next;}}
		if ($l =~ m/^.................(...) (.)(....)/) {
		  my $two = $2;
		  my $one = $1;
		  my $tre = $3;
			if ($two eq " ") {$two = "-";}
			
			my $resnumc = "$one $tre $two";
			$resnumc =~ s/\s+/ /gi;
			if ($count{$resnumc} ne undef) {next;}
			$count{$resnumc} = $c;
			++$c;
			
			
			
		}
	}
	close FILE;
	return(%count);
}

sub min_dist {
        my $str1 = $_[0];
        my $str2 = $_[1];
        my $mask = $_[2];
        my $tomask = 1;
        if ($mask eq undef) {$tomask = 0;}
        my $min = 999999999999;
        foreach my $l1 (split(/\n/,$str1)) {
            if ($tomask == 1) {if ($l1 =~ m/$mask/) {}else {next;}}
                foreach my $l2 (split(/\n/,$str2)) {
                    if ($tomask == 1) {if ($l2 =~ m/$mask/) {}else {next;}}
                        my @cor1;
                        my @cor2;
                        if ($l1 =~ m/^............................\s+(\S+)\s+(\S+)\s+(\S+)/) {@cor1 = ($1,$2,$3);}
                        if ($l2 =~ m/^............................\s+(\S+)\s+(\S+)\s+(\S+)/) {@cor2 = ($1,$2,$3);}
                        my $dist = 0;
                        
                        for my $i (0..2) {
                                $dist += ($cor1[$i]-$cor2[$i])**2
                        }
                        if($dist < $min) {$min = $dist;}
                }
        }
	return($min);
}

sub load_pdb_atm {
  my $het = $_[1];
  my $mask = $_[2];
  open FILE, "$_[0]" or die "Can't open $_[0]\n" ;
  if ($mask eq undef) {$mask = ".";}
  my %count;
  my $c = 0;
  # Get ResnumC

  while(my $l = <FILE>) {
    if ($mask ne "NA") {unless($l =~ m/$mask/) {next;}}
          if ($het == 0 && $l =~ m/^HET/) {next;}
          if ($l =~ m/^.................(...) (.)(....)/) {
            my $two = $2;
            my $one = $1;
            my $tre = $3;
                  if ($two eq " ") {$two = "-";}

                  my $resnumc = "$one $tre $two";
                  $resnumc =~ s/\s+/ /gi;

                  $count{$resnumc} .= $l;
                  ++$c;
								#	print "$l";


          }
  }
  close FILE;
  return(%count);
}

