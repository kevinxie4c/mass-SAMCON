#!/usr/bin/env perl

while(<>) {
	chomp;
	if (/ROOT\s(\w+)/ || /JOINT\s(\w+)/) {
		if (!$printed{$1}) {
			print "$1\n";
			$printed{$1} = 1;
		}
	}
}
