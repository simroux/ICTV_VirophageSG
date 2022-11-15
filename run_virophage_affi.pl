#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h='';
my $dir='';
my $input_file="";
my $out_dir="";
my $n_cpu=4;
my $full_out=0;
my $db_dir="Db/";
GetOptions ('help' => \$h, 'h' => \$h, 'i=s'=>\$input_file, 'o=s'=>\$out_dir, 'd=s'=>\$db_dir, 't=s'=>\$n_cpu, 'f'=>\$full_out);
if ($h==1 || $input_file eq "" || $out_dir eq ""){ # If asked for help or did not set up any argument
	print "# Script to affi a fasta file
# Arguments :
# -i: input file
# -o: output directory
# Optional:
# -d: path to the database folder (if different from Db/)
# -t: number of threads (default = 4)
# -f: include all contigs in the output file, even if not identified as either virophage or PLV
";
	die "\n";
}
## Cutoff and taxa names stored in centralized files
my $cutoff_file=$db_dir."/Cutoffs.tsv";
my $db_name=$db_dir."/names.tsv";
##
my $db_hmm=$db_dir."/All_markers.hmm";
my $blast_db=$db_dir."/MCP_blast_db";
my $db_hmm_plv=$db_dir."/PLV_PC_054.hmm";

$input_file=~/(.*)\.[^\.]+/;
my $tmp=$1;
if ($tmp=~/\/([^\/]+)/){
      $tmp=$1;
}
my $basename=$tmp;
print "$input_file -- $basename\n";

my $wdir=$out_dir;
if (!(-d $wdir)){
	print "Please create the directory $out_dir before running the affiliation script\n";
	die("\n");
}
if (!($wdir=~/\//)){$wdir.="/";}
## First step, prodigal
my $faa_file=$wdir.$basename.".faa";
if (!(-e $faa_file)){
      my $gff_file=$wdir.$basename.".gff";
      &run_cmd("prodigal -p meta -i $input_file -f gff -o $gff_file -a $faa_file");
}
## Then, match to markers to know if these are virophages
my $hmm_step_1=$wdir.$basename."_vs_all_markers-step_1.tsv";
if (!(-e $hmm_step_1)){
      &run_cmd("hmmsearch -o /dev/null --noali -E 1.0e-05 --tblout $hmm_step_1 --cpu $n_cpu $db_hmm $faa_file");
}
## Parse it into a tsv file
my $best_marker_result=$wdir.$basename."_best_markers-step_1.tsv";
if (!(-e $best_marker_result)){
      &parse_step_1($faa_file,$hmm_step_1,$best_marker_result);
}
## Get the list of virophages from there
my $virophage_list=$wdir.$basename."_virophages_list.tsv";
my $virophage_faa=$wdir.$basename."_virophages_cds.faa";
if (!(-e $virophage_list)){
      &select_virophages($best_marker_result,$virophage_list,$faa_file,$virophage_faa,$cutoff_file);
}
## Step 2: Blast virophages to MCPs
my $blast_result=$wdir.$basename."_virophage_cds-vs-marker.tsv";
if (!(-e $blast_result)){
      &run_cmd("blastp -query $virophage_faa -db $blast_db -out $blast_result -outfmt 6 -num_threads $n_cpu");
}
## Parse and affiliate to families (or not)
my $family_affi=$wdir.$basename."_virophage_cds_to_family.tsv";
if (!(-e $family_affi)){
      &parse_step_2($virophage_faa,$blast_result,$family_affi,$cutoff_file);
}
## Step 3: bonus hmm search of the PgVV marker
my $hmm_step_3=$wdir.$basename."_vs_PLV_Marker.tsv";
if (!(-e $hmm_step_3)){
      &run_cmd("hmmsearch -o /dev/null --noali -E 1.0e-05 --tblout $hmm_step_3 --cpu $n_cpu $db_hmm_plv $faa_file");
}
## Parse it into a tsv file
my $PLV_result=$wdir.$basename."_best_markers-PLV.tsv";
if (!(-e $PLV_result)){
      &parse_step_1($faa_file,$hmm_step_3,$PLV_result);
}
## Make a final file with affiliation of everything
my $summary_file=$wdir.$basename."_final_affiliation.tsv";
&summarize($input_file,$virophage_list,$best_marker_result,$family_affi,$PLV_result,$db_name,$summary_file,$cutoff_file,$full_out);

sub summarize{
      my $in_fna=$_[0];
      my $list_viro=$_[1];
      my $marker_result=$_[2];
      my $family_affi=$_[3];
      my $plv_result=$_[4];
      my $file_name=$_[5];
      my $out_file=$_[6];
	my $cutoff_file=$_[7];
	my $tag_full=$_[8];
      my %translate_name;
      open my $tsv,"<",$file_name;
      while(<$tsv>){
            chomp($_);
            my @tab=split("\t",$_);
            $translate_name{$tab[0]}=$tab[1];
      }
      close $tsv;
      my %cutoffs;
      open my $tsv,"<",$cutoff_file;
      while(<$tsv>){
            chomp($_);
            my @tab=split("\t",$_);
            $cutoffs{$tab[0]}{$tab[1]}=$tab[2];
      }
      close $tsv;
      my %info;
      open my $fa,"<",$in_fna;
      while(<$fa>){
            chomp($_);
            if ($_=~/^>(\S+)/){
                  my $id=$1;
                  $info{$id}{"selected"}=1;
            }
      }
      close $fa;
      open my $tsv,"<",$list_viro;
      while(<$tsv>){
            chomp($_);
            my @tab=split("\t",$_);
            if (defined($info{$tab[0]})){
                  $info{$tab[0]}{"Maveriviricetes"}=1;
            }
            elsif($tab[0] ne "Genome"){
                  die("unexpected id $tab[0]\n");
            }
      }
      close $tsv;
      open my $tsv,"<",$marker_result;
      while(<$tsv>){
            chomp($_);
            my @tab=split("\t",$_);
            if (defined($info{$tab[0]})){
                  if ($tab[1] eq "MCP" && $tab[2]>=$cutoffs{"MCP"}{"HMM"}){
                        $info{$tab[0]}{"markers"}{$tab[1]}=$tab[2];
                  }
                  elsif($tab[2]>=$cutoffs{$tab[1]}{"HMM"}){
                        $info{$tab[0]}{"markers"}{$tab[1]}=$tab[2];
                  }
            }
      }
      close $tsv;
      open my $tsv,"<",$family_affi;
      while(<$tsv>){
            chomp($_);
            my @tab=split("\t",$_);
            if ($tab[0]=~/(.*)_\d+$/){$tab[0]=$1;}
            if (defined($info{$tab[0]})){
                  if ($tab[1] ne "NA"){
                        $tab[1]=$translate_name{$tab[1]};
                        $info{$tab[0]}{"family"}=$tab[1];
                  }
            }
            elsif($tab[0] ne "Prot"){
                  die("unexpected id $tab[0]\n");
            }
      }
      close $tsv;
      open my $tsv,"<",$plv_result;
      while(<$tsv>){
            chomp($_);
            my @tab=split("\t",$_);
            if (defined($info{$tab[0]})){
                  if ($tab[2]>=$cutoffs{"PLV"}{"HMM"}){
                        $info{$tab[0]}{"PLV_hit"}=1;
                  }
            }
            elsif($tab[0] ne "Genome"){
                  die("unexpected id $tab[0]\n");
            }
      }
      close $tsv;
	print "##################### FINAL RESULTS ###########################\n";
      open my $s1,">",$out_file;
      print $s1 "Genome\tAssignation (Class;Family)\tMaveriviricetes marker list\n";
	my $tag=0;
      foreach my $genome (sort keys %info){
            my $line=$genome;
		$tag=0;
            if ($info{$genome}{"Maveriviricetes"}==1){
                  $line.="\tMaveriviricetes";
                  if (defined($info{$genome}{"family"})){$line.=";".$info{$genome}{"family"};}
                  else{$line.=";Unclassified";}
                  my @list=sort keys %{$info{$genome}{"markers"}};
                  $line.="\t".join(" ",@list);
			$tag=1;
                  # if ($info{$genome}{"PLV_hit"}==1){$line.="\tAdditional hit to PLV";}
            }
            else{
                  if ($info{$genome}{"PLV_hit"}==1){$line.="\tPossible_PLV"; $tag=1;}
                  else{$line.="\tNA";}
                  $line.="\tNA";
            }
		if ($tag==1 || $tag_full==1){
            	print $s1 $line."\n";
			print $line."\n";
		}
      }
      close $s1;

	print "#################################################################\n";
	print "Results stored in $out_file\n";
}



sub parse_step_2{
      my $in_faa=$_[0];
      my $blast_result=$_[1];
      my $out_file=$_[2];
      my $cutoff_file=$_[3];
      my %check;
      open my $fa,"<",$in_faa;
      while(<$fa>){
            chomp($_);
            if ($_=~/^>(\S+)/){
                  my $id=$1;
                  $check{$id}=1;
            }
      }
      close $fa;
      my %cutoffs;
      open my $tsv,"<",$cutoff_file;
      while(<$tsv>){
            chomp($_);
            my @tab=split("\t",$_);
            $cutoffs{$tab[0]}{$tab[1]}=$tab[2];
      }
      close $tsv;
      my %store;
      open my $tsv,"<",$blast_result;
      while(<$tsv>){
            chomp($_);
            my @tab=split("\t",$_);
            my $query=$tab[0];
            my $hit=$tab[1];
            my @t=split(/\|/,$hit);
            my $clade=$t[0];
            my $hit_g=$t[1];
            my $score=$tab[11];
            if ($score>=$cutoffs{$clade}{"Blast"}){
                  if (!defined($store{$query}{"hits"}{$hit}) || ($store{$query}{"hits"}{$hit}{"score"} < $score)){
                        $store{$query}{"hits"}{$hit}{"score"}=$score;
                        $store{$query}{"hits"}{$hit}{"clade"}=$clade;
                  }
            }
      }
      close $tsv;
      open my $s1,">",$out_file;
      print $s1 "Prot\tClade\tScore\n";
      foreach my $prot (sort keys %check){
            if (!defined($store{$prot})){
                  print $s1 $prot."\tNA\tNA\n";
            }
            else{
                  my @t=sort {$store{$prot}{"hits"}{$b}{"score"} <=> $store{$prot}{"hits"}{$a}{"score"} or $a cmp $b} keys %{$store{$prot}{"hits"}};
                  my $best=$t[0];
                  print $prot."\t".$store{$prot}{"hits"}{$best}{"clade"}."\t".$store{$prot}{"hits"}{$best}{"score"}."\n";
                  print $s1 $prot."\t".$store{$prot}{"hits"}{$best}{"clade"}."\t".$store{$prot}{"hits"}{$best}{"score"}."\n";
            }
      }
      close $s1;
}

sub select_virophages{
      my $tsv_in=$_[0];
      my $out_file=$_[1];
      my $in_faa=$_[2];
      my $out_faa=$_[3];
      my %check;
      my %check_prot;
      my %cutoffs;
      open my $tsv,"<",$cutoff_file;
      while(<$tsv>){
            chomp($_);
            my @tab=split("\t",$_);
            $cutoffs{$tab[0]}{$tab[1]}=$tab[2];
      }
      close $tsv;
      open my $tsv,"<",$tsv_in;
      while(<$tsv>){
            chomp($_);
            my @tab=split("\t",$_);
            if ($tab[1] eq "MCP" && $tab[2]>=$cutoffs{"MCP"}{"HMM"}){
                  $check{$tab[0]}=$tab[2];
                  $check_prot{$tab[3]}=1;
            }
      }
      close $tsv;
      open my $s1,">",$out_file;
      print $s1 "Genome\tMCP score\n";
      foreach my $genome (sort keys %check){
            print $s1 $genome."\t".$check{$genome}."\n";
      }
      close $s1;
      my $tag=0;
      open my $s1,">",$out_faa;
      open my $fa,"<",$in_faa;
      while(<$fa>){
            chomp($_);
            if ($_=~/^>(\S+)/){
                  my $id=$1;
                  $tag=0;
                  if ($check_prot{$id}==1){
                        $tag=1;
                        print $s1 $_."\n";
                  }
            }
            elsif($tag==1){
                  print $s1 $_."\n";
            }
      }
      close $fa;
      close $s1;
}


sub parse_step_1{
      my $in_faa=$_[0];
      my $hmm_result=$_[1];
      my $out_file=$_[2];
      my %check;
      open my $fa,"<",$in_faa;
      while(<$fa>){
            chomp($_);
            if ($_=~/^>(\S+)/){
                  my $id=$1;
                  if ($id=~/(.*)_\d+$/){
                        my $genome=$1;
                        $check{$genome}{$id}=1;
                  }
                  else{
                        die("Pblm with protein id $id\n");
                  }
            }
      }
      close $fa;
      my %store;
      open my $tsv,"<",$hmm_result;
      while(<$tsv>){
            chomp($_);
            if ($_=~/^#/){next;}
            my @tab=split(" ",$_);
            if ($tab[2]=~/(.*)_\d+$/){
                  my $marker=$1;
                  if ($tab[0]=~/(.*)_\d+$/){
                        my $genome=$1;
                        my $code=$tab[0].";".$tab[2];
                        if (!defined($store{$genome}{$marker}{$code}{"score"}) || ($tab[5]>$store{$genome}{$marker}{$code}{"score"})){
                              $store{$genome}{$marker}{$code}{"score"}=$tab[5];
                              $store{$genome}{$marker}{$code}{"prot"}=$tab[0];
                        }
                  }
                  else{
                        print "Pblm with genome $tab[0]\n";
                  }
            }
            else{
                  print "Pblm with marker $tab[1]\n";
            }
      }
      close $tsv;
      open my $s1,">",$out_file;
      print $s1 "Genome\tMarker\tScore\tProtein\n";
      foreach my $genome (sort keys %check){
            if (defined($store{$genome})){
                  foreach my $marker (sort keys %{$store{$genome}}){
                        my @list_code=sort { $store{$genome}{$marker}{$b}{"score"} <=> $store{$genome}{$marker}{$a}{"score"}  } keys %{$store{$genome}{$marker}};
                        print $s1 $genome."\t".$marker."\t".$store{$genome}{$marker}{$list_code[0]}{"score"}."\t".$store{$genome}{$marker}{$list_code[0]}{"prot"}."\n";
                  }
            }
            else{
                  print $s1 $genome."\tNA\tNA\tNA\n";
            }
      }
      close $s1;
}


sub run_cmd{
	my $cmd=$_[0];
	if ($_[1] ne "veryquiet"){print "$cmd\n";}
	my $out=`$cmd`;
	if ($_[1] ne "quiet" && $_[1] ne "veryquiet"){
		if ($_[1] eq "stderr"){print STDERR "$out\n";}
		else{print "$out\n";}
	}
	return($out);
}
