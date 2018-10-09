#!/usr/bin/env perl

while(<>) {
	chomp;
	if (/ROOT\s(\w+)/ || /JOINT\s(\w+)/) {
		$joint = $1;
	}
	if (/CHANNELS\s\d+\s(.+)/) {
		@chans = split ' ', $1;
		print "$joint $_\n" for @chans;
	}
}
