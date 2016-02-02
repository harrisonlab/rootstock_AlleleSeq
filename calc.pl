#! /usr/bin/perl -w -s

open(IN,$ARGV[0]);
my @inp=<IN>;

foreach(@inp){
	my @split = split(/\t/,$_);
	if ($split[1] >= $ARGV[1] and $split[1] <= $ARGV[2]) {
		print "$_";
	}
}
