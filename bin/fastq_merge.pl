#!/usr/bin/perl
use Parallel::ForkManager;

#make directory with projectname and subdirectories for various analysis
system("mkdir ../fastq_merged") unless (-d fastq_merged);

open FILE, '../data/03feb2017fastqlist.csv' or die
<FILE>;

my @keys = map{chomp;$_} split(',',<FILE>);

#create an array samples
my @samples;

#splitting each line based on ',' and storing in an array @r
#pushing the reference of this array in another array @samples
while(<FILE>){
	chomp;
	#my @r = split(',');
	my %h;
	@h{@keys}=split(',');
	$h{'seqname'}=$h{'SAMP_name'};
	push(@samples,\%h);
}

#run parallel jobs
my $pm=new Parallel::ForkManager(15);
my $ext1 ="_mate1.fastq.gz";
my $ext2="_mate2.fastq.gz";
foreach (@samples)
{
	$pm->start and next;
		my $cmd1 = "cat $_->{'mate1'} > ../fastq_merged/$_->{'SAMP_name'}$ext1";
		my $cmd2 = "cat $_->{'mate2'} > ../fastq_merged/$_->{'SAMP_name'}$ext2";
		print $cmd1,"\n";
		print $cmd2,"\n";
        system($cmd1);
        system($cmd2);
        $pm->finish;
}
$pm->wait_all_children;

