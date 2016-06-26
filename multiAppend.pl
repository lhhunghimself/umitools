#!/usr/bin/perl
use strict; use warnings;
use threads;

my($file,$numThreads)=@ARGV;
system("mkdir -p /tmp/locks.$$");
my @threads = initThreads($numThreads);
print "@threads\n";
our @cmds,
our @done;
my @lines=split(/\n/,`cat $file`);
foreach my $cmd (@lines){
	if($cmd =~/\S/){
	 #check if output has been defined
	 if ($cmd !~ /\-o /){
	 	my (@parts)=split(' ',$cmd);
	 	foreach my $part(@parts){
				if($part =~/\_R1\./){
					$part =~ s/\_R1//;
					if($part =~/\.gz/){
						$part=substr($part,0,-3);
					}	
					$cmd.= " -o $part";
				 last;
				}
			}
			push(@cmds,"umiappend $cmd");	 
		}
	}
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
	  system("$cmds[$i]");
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
