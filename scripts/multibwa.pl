#!/usr/bin/perl
use strict; use warnings;
use threads;


# 1.3 Reference
our $TOP_DIR="/mnt/backup/DetoxS/UMITest";
our $REF_DIR="$TOP_DIR/References/Broad_UMI";
our $SPECIES_DIR="$REF_DIR/Human_RefSeq";
our $REF_SEQ_FILE="$SPECIES_DIR/refMrna_ERCC_polyAstrip.hg19.fa";

my($numThreads)=@ARGV;
system("mkdir -p /tmp/locks.$$");
my @threads = initThreads($numThreads);
print "@threads\n";
our @cmds,
our @done;

my @dirs=split(' ',`ls output`);
foreach my $dir (@dirs){
 system("rm output/$dir/*.sam");
	push(@cmds,"output/$dir");	
}
foreach(@threads){
	$_ = threads->create(\&doOperation);
}

foreach(@threads){
 $_->join();
}
system("rm -rf /tmp/locks*");
				

sub doOperation{
	my $id = threads->tid();
	foreach my $i (0..$#cmds){
		my $lockfile="/tmp/locks.$$/lock.$i";
 	my $donefile="/tmp/locks.$$/done.$i";
		if (mkdir($lockfile)){
	  printf stderr "%d %s\n",$id,$cmds[$i];
	  my $dir=$cmds[$i];
	  system("cat $dir/*.fq >  $dir/all.fastq");
	  my $cmd="(bwa aln -l 24 -t 1 $REF_SEQ_FILE $dir/all.fastq | bwa samse -n 20 $REF_SEQ_FILE - $dir/all.fastq | grep -v '^\@' > $dir/all.sam) >& /dev/null";
	  system("$cmd");
	  system("rm $dir/all.fastq\n");
	  mkdir($donefile);
		}
		else{
			#printf stderr "%d %s $i exists\n",$id,$lockfile;
		}	
	}
	threads->exit();
}
sub initThreads{
	my($numThreads)=@_;
	my @initThreads;
	for(my $i = 1;$i<=$numThreads;$i++){
		push(@initThreads,$i);
	}
	return @initThreads;
}
