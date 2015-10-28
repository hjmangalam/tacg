#!/usr/bin/env perl

# this script runs tacg thru several tests to see if it at least passes a light
# stress test without bonking

# still to be tested:
#     --HTML (this needs some more work to add HTML support for those
#             stanzas that have been added since I stopped work on the Web
#             stuff
#     -R (alt pattern file)
#     different input formats
#     -x
#
#
require 5.005;
my $md5found;
BEGIN {eval "use Digest::MD5"; $md5found = $@ ? 0 : 1}
# $md5found = 0;  # uncomment to simulate no MD5
if ($md5found == 0) {
    print "\nOoops! Digest::MD5 module not found - continuing with simpler error checking\n\n" ;
} else {
    print "\nGood! Digest::MD5 module found ... continuing with test\n\n";
}

use strict;
use vars qw( $cmdline $cmdline1 $MD5_hash  %result_hash $fqsd
);
use Getopt::Long; # GNU-style getopt
#use Digest::MD5;
$ENV{TACGLIB} = "../Data";
# or read in the data to verify that the analysis ran as expected using
# an eval

my $tacg =       "../tacg";
#my $tacg =       "/usr/local/bin/tacg41";

my $seq =        "../Seqs/hlef.seq";
my $fastaseq =   "../Seqs/hsp.cluster.1000.fasta";
my $REBASE =     "../Data/rebase.data";
my $regex =      "../Data/regex.data";
my $rulefile =   "../Data/rules.data";
my $matrixfile = "../Data/matrix.data";
my $GCC = `gcc --version |grep -i gcc`;
my $out = "";
my $diff_app = "/usr/bin/diff -b "; # /usr/bin/kompare or /usr/bin/tkdiff
my $save = 0;
my $diff_test = "";
my $test_name = "";

&GetOptions(
	'save'         => \$save,       #  save all the output to 'golden' dir.
	'diff=s'       => \$diff_test,  # string is the test to diff the output against
	'test=s'       => \$test_name,  # explicit run to test
);


#my $hostname = `hostname`;
#print "local host name = $hostname\n";
# generate the hashes below by grepping thru test.log:
# grep '#wc' test.log
my %wcTable = (
"raw" => " 98 317 7348 out", #wc
"raw_infile_fasta" => " 896 3444 65013 out", #wc
"rev" => " 489 1777 39019 out", #wc
"comp" => " 497 1765 39667 out", #wc
"revcomp" => " 482 1775 38452 out", #wc
"notics_strands_numstart" => " 3437 61073 481015 out", #wc
"notics_strands_numstart_infile_fasta" => " 3010 23745 385618 out", #wc
"overhang_5" => " 95 416 4015 out", #wc
"overhang_3" => " 63 272 2553 out", #wc
"overhang_0" => " 89 379 3822 out", #wc
"Fragments_0" => " 43796 148482 3469376 out", #wc
"Fragments_1" => " 47708 181453 3634675 out", #wc
"Fragments_2" => " 47708 181453 3634673 out", #wc
"Fragments_3" => " 51620 214424 3799972 out", #wc
"FASTA_Fragments_0" => " 7213 26154 389268 out", #wc
"FASTA_Fragments_1" => " 10599 42000 485415 out", #wc
"FASTA_Fragments_2" => " 10599 42000 485389 out", #wc
"FASTA_Fragments_3" => " 13985 57846 581536 out", #wc
"begin_end_0" => " 43584 146266 3452822 out", #wc
"begin_end_1" => " 19891 79809 1842073 out", #wc
"begin_end_2" => " 16993 76296 1813014 out", #wc
"begin_end_3" => " 1993 9081 202928 out", #wc
"begin_end_4" => " 7086 38174 933768 out", #wc
"wide_ladder_map" => " 188 1008 23502 out", #wc
"extract_hindiii_bamhi_0" => " 205 2267 16101 out", #wc
"extract_regex_0" => " 179 1968 14785 out", #wc
"extract_hindiii_bamhi_1" => " 149 1604 11325 out", #wc
"extract_regex_1" => " 97 1006 7454 out", #wc
"extract_hindiii_bamhi_2" => " 161 1742 12205 out", #wc
"extract_regex_2" => " 87 891 6601 out", #wc
"regex_alt1" => " 335 3762 28521 out", #wc
"regex_FILE" => " 138 987 7883 out", #wc
"patterns_errors_0" => " 17 98 612 out", #wc
"patterns_errors_1" => " 22 128 826 out", #wc
"patterns_errors_2" => " 61 480 3300 out", #wc
"matrix_test_76" => " 834 7241 50088 out", #wc
"matrix_test_79" => " 571 4900 33896 out", #wc
"matrix_test_81" => " 483 4108 28436 out", #wc
"matrix_test_85" => " 332 2719 18819 out", #wc
"matrix_test_86" => " 302 2455 16984 out", #wc
"matrix_test_87" => " 290 2347 16250 out", #wc
"matrix_test_88" => " 235 1859 12877 out", #wc
"matrix_test_89" => " 205 1635 11258 out", #wc
"matrix_test_90" => " 190 1511 10403 out", #wc
"matrix_test_91" => " 187 1476 10161 out", #wc
"matrix_test_92" => " 168 1342 9186 out", #wc
"matrix_test_93" => " 66 451 3041 out", #wc
"matrix_test_95" => " 34 179 1172 out", #wc
"matrix_test_98" => " 25 130 817 out", #wc
"matrix_test_99" => " 25 128 803 out", #wc
"dam_dcm_test" => " 1714 13637 98117 out", #wc
"dcm_test" => " 1916 15850 116942 out", #wc
"clone_1" => " 258 1147 14191 out", #wc
"clone_2" => " 19 113 702 out", #wc
"clone_3" => " 4 23 141 out", #wc
"cost_10" => " 67 319 1993 out", #wc
"cost_25" => " 28 222 1891 out", #wc
"cost_75" => " 42 213 1289 out", #wc
"graphics_test_X" => " 21 388 3258 out", #wc
"graphics_test_Y" => " 184 467 4171 out", #wc
"graphics_test_L" => " 172 389 3864 out", #wc
"Codon_tables_0" => " 197 2342 15405 out", #wc
"Codon_tables_1" => " 150 1765 11602 out", #wc
"Codon_tables_2" => " 490 6262 42605 out", #wc
"Codon_tables_3" => " 338 4185 28481 out", #wc
"Codon_tables_4" => " 490 6262 42604 out", #wc
"Codon_tables_5" => " 418 4929 34348 out", #wc
"Codon_tables_6" => " 338 4185 28481 out", #wc
"Codon_tables_7" => " 333 4106 27973 out", #wc
"Codon_tables_8" => " 197 2342 15403 out", #wc
"Codon_tables_9" => " 197 2342 15403 out", #wc
"Codon_tables_10" => " 485 6183 42094 out", #wc
"Codon_tables_11" => " 587 6727 48861 out", #wc
"Codon_tables_12" => " 266 3225 21445 out", #wc
"Codon_tables_13" => " 16 86 527 out", #wc
"Codon_tables_14" => " 16 86 527 out", #wc
"Codon_tables_15" => " 16 86 527 out", #wc
"Codon_tables_16" => " 16 86 527 out", #wc
"Proximity_1" => " 34 142 1260 out", #wc
"Proximity_2" => " 37 156 1564 out", #wc
"Proximity_3" => " 34 149 1246 out", #wc
"silent_short" => " 399 1553 31724 out", #wc
"silent_long" => " 1311 5437 105599 out", #wc
"ORF1_short_orfmap" => " 557 7331 48370 out", #wc
"ORF2_medium" => " 72 638 4206 out", #wc
"ORF3_long_orfmap" => " 93 629 5279 out", #wc
"ORF6xS_44_orfmap" => " 218 2516 17595 out", #wc
"ORF6xL_44_orfmap" => " 1796 27374 182043 out", #wc
"ORF6xL_150_orfmap" => " 84 672 6086 out", #wc
);

# generate the hashes below by grepping thru test.log:
# grep '#MD5' test.log
my %MD5_table = (
"raw" => "ccf790625524591e0e77cac5335db7c6", #MD5
"raw_infile_fasta" => "75a08fae811fc2fc265007070e44ec6e", #MD5
"rev" => "79adc8c64ad485190d0c9a09ae5b1377", #MD5
"comp" => "0ec08b476e1c2a9364690c429b360dd8", #MD5
"revcomp" => "398f7da07aea851f5aa3d4840aa06878", #MD5
"notics_strands_numstart" => "cb636e1d06231b73d3ea1d7d8adcfb80", #MD5
"notics_strands_numstart_infile_fasta" => "135ce7cd9c70af84ccc376b7f4f2d89b", #MD5
"overhang_5" => "5a10dc5985d72039cbc19dc6d5c87192", #MD5
"overhang_3" => "944cfdac7bf357f59b385b0754a5a2b5", #MD5
"overhang_0" => "92b8d0bd29a7b89ef5680a2d476c4b41", #MD5
"Fragments_0" => "595974048c9a3821321ea6e335df8790", #MD5
"Fragments_1" => "22f07db2ba32f58e89e718c097f7b97c", #MD5
"Fragments_2" => "f350364ce43f5bf721fd8bb67c9ecd42", #MD5
"Fragments_3" => "cb51ed6e809c58f83bd5c57fa41972c0", #MD5
"FASTA_Fragments_0" => "a3cfb4148c62f548c6b750738c5cd8bb", #MD5
"FASTA_Fragments_1" => "fbae0fa47caeae6233fbd921344116bb", #MD5
"FASTA_Fragments_2" => "a7ce5120c9edf0598432374c18854754", #MD5
"FASTA_Fragments_3" => "e4aaf46fefe76872a0132abc1bbb6a8f", #MD5
"begin_end_0" => "51a84f1dc26594b713c2e9966a7f973a", #MD5
"begin_end_1" => "bf071906113bbcb2d126850eb64b9091", #MD5
"begin_end_2" => "44a7ec8ac8bc20d1601d69e5625536f5", #MD5
"begin_end_3" => "28b836d89483c9df28e06632f6e0770e", #MD5
"begin_end_4" => "6c39e6a48475818b6f6508a42a949d2e", #MD5
"wide_ladder_map" => "05664efc3d8447a0375b897b0190cdd8", #MD5
"extract_hindiii_bamhi_0" => "ad0f69cfaabe0d71023426a2e3435f73", #MD5
"extract_regex_0" => "4c612d528c67bc64e421d95197eaf355", #MD5
"extract_hindiii_bamhi_1" => "536a5ae4ac9e5c3a9b9326a254b818a8", #MD5
"extract_regex_1" => "aaa754f165fb8805ad8308b2eaa32f91", #MD5
"extract_hindiii_bamhi_2" => "dfb45bdd8d0fd0e3305984c1819e1245", #MD5
"extract_regex_2" => "6fe8515225a74cd9dda6d82508d603c4", #MD5
"regex_alt1" => "7f6e3d3439a3ac9dbe04b06b475711f2", #MD5
"regex_FILE" => "5fafba7b4017704d0e92628c3822270c", #MD5
"patterns_errors_0" => "0728ff10829cfb625ff4e99a030b1618", #MD5
"patterns_errors_1" => "57d1e43272be6d938b6df4caef68cd16", #MD5
"patterns_errors_2" => "766e98c424d602f7cab690ec66b27c35", #MD5
"matrix_test_76" => "01bc20c346601a4f97889b0339dd495b", #MD5
"matrix_test_79" => "d737c75442a4761aae393bcf72cb3d18", #MD5
"matrix_test_81" => "c0eeb2cf7ccf171a430a1c1fa57c0f7d", #MD5
"matrix_test_85" => "2b32e03892ee1cf85a161f00c1dfdc93", #MD5
"matrix_test_86" => "e36201b7115d076af86cf61471bf40b4", #MD5
"matrix_test_87" => "d0c23dd33d1c7bb2d3cc6325c3464de2", #MD5
"matrix_test_88" => "f61c3417bd90285b061b05018a96adc1", #MD5
"matrix_test_89" => "5b4b8eb09d13cd2c835f7f2d4327236b", #MD5
"matrix_test_90" => "93ac892552a61420fd3a0f4c0b3cb6ed", #MD5
"matrix_test_91" => "49b7323b487186c69c45158749453bf1", #MD5
"matrix_test_92" => "403d47c0f8294a374ba59ff8561885a0", #MD5
"matrix_test_93" => "b77eb8e1de241d71005d88f398a3b60b", #MD5
"matrix_test_95" => "ea476f4a1ff145e9f68c66a8b7fc372e", #MD5
"matrix_test_98" => "e128eaf68ae065d0ab7a375e956bad03", #MD5
"matrix_test_99" => "4b6e0586c6ee4b4df6cb47c4d107834d", #MD5
"dam_dcm_test" => "25a181b83cbb6868272cd6e971e9b581", #MD5
"dcm_test" => "62c971ccc085eb5f843363f0f4443ab3", #MD5
"clone_1" => "357918e7ae58fc7b3b9db05a0956a438", #MD5
"clone_2" => "836d044b01a84daad58aeb08796a0baa", #MD5
"clone_3" => "19f58dc34575e5c6d299fa5f9cfa52e1", #MD5
"cost_10" => "cee8c8cf390de12136fcf94bb94a20a9", #MD5
"cost_25" => "27e72eeb074499013a2174ae43ceebca", #MD5
"cost_75" => "ef12f48964b49f24c71b4be78c3ab429", #MD5
"graphics_test_X" => "94efb51ae1715ae7afdd9e512014c709", #MD5
"graphics_test_Y" => "f62a6a7de2b07641afcdb81105e4188b", #MD5
"graphics_test_L" => "d02c075f57fb4b3b493a062f8ba4300d", #MD5
"Codon_tables_0" => "77252a1d6ef467a066a848a9a542fcb7", #MD5
"Codon_tables_1" => "a8b3bf86438092d6424f3355c195a1bf", #MD5
"Codon_tables_2" => "d966d145001eefe0fea8e648f1fe8d56", #MD5
"Codon_tables_3" => "95562e78dbe9d4fbc9b8fa231f1a48d0", #MD5
"Codon_tables_4" => "6ca0ad6aa8837af5ab5e51b7d766292e", #MD5
"Codon_tables_5" => "69850a9b7abed8666b303376afc93eb6", #MD5
"Codon_tables_6" => "a0b141ab235a658827f8594f775bde49", #MD5
"Codon_tables_7" => "51825edef1256603b469d0feaba71915", #MD5
"Codon_tables_8" => "a39d42f06224a3593c9cb7ac594d06d5", #MD5
"Codon_tables_9" => "8076fa49d57e360f4c3ff76100784287", #MD5
"Codon_tables_10" => "ed3c631d058b5755c3933c558b75ed5e", #MD5
"Codon_tables_11" => "0a0c5d93a1977dcb23f2d023c86f5cfa", #MD5
"Codon_tables_12" => "974477a288c193e86565f85dd59c8ae8", #MD5
"Codon_tables_13" => "c339784829281e0a1e64ab40b4cdc92f", #MD5
"Codon_tables_14" => "cf66f9804be78a493b52c16b4d2f63e3", #MD5
"Codon_tables_15" => "024bb411784a0aa984c32ff8386bad77", #MD5
"Codon_tables_16" => "19400305f6665bb0d4d0fc9b4005f937", #MD5
"Proximity_1" => "5fe719d0fb492dac500b9c45226dac08", #MD5
"Proximity_2" => "5d86a7fc9281b1f1806986752f69d4fd", #MD5
"Proximity_3" => "d3a1f217d51397c912826288a8a2f34f", #MD5
"silent_short" => "95e71bd281f38445cebf538e9d1ed3a8", #MD5
"silent_long" => "dee7c8cda4eff567c54583646c7bd80b", #MD5
"ORF1_short_orfmap" => "5c0e8c212cfe30088cf9c159adae614f", #MD5
"ORF2_medium" => "6c5db0d97f8882aafb9035bc50ce1109", #MD5
"ORF3_long_orfmap" => "85910b6b7128d8cf23eee26f55e9653a", #MD5
"ORF6xS_44_orfmap" => "288c9adb14d6b78446e5294b42e484d9", #MD5
"ORF6xL_44_orfmap" => "040c57bb35fc407b33dfaf023e5ff477", #MD5
"ORF6xL_150_orfmap" => "0e4dfab27310370e66738c914b6bbd23", #MD5
);

print << "POI";

Now running thru some automated tests of the tacg executable that
you built.. This will take as long as a few minutes to run, depending on the
speed of your machine.  If it fails, please send the logfile back to me:
hjm\@tacgi.com

The signatures for these tests were generated with:
gcc (GCC) 4.0.3 (Ubuntu 4.0.3-1ubuntu5)
on:
Linux bodi 2.6.16.18 #3 PREEMPT Tue Jun 6 19:10:29 PDT 2006 i686 GNU/Linux

The gcc you used to build tacg was probably:
$GCC
If they're different, certain tests may fail due to compiler differences
and there are also output differences between the same compiler/libs on
different architectures.

Please check not only the MD5 signatures, but also the wc sigs if
you're compiling on different compilers and different platforms.

Nevertheless, I'd like to know if things appear to bugger up.

FYI, this test is using [$tacg] as the target executable.

<Hit a key to continue>

POI

my $tmp = <STDIN>;

open(LOG, ">test.log") or die "Can't open logfile 'test.log': $!";


my @begin = (1,   765,   22989,  77673,  101798);
my @end   = (0, 77673,  101798,  86735,       0);
my @matrix_cutoffs = (76, 79, 81, 85, 86, 87, 88, 89, 90, 91, 92, 93, 95, 98, 99);

print LOG "from:  ", `uname -a`, "\n";
print LOG "using:  ", `../tacg -v`, "\n";

if ($save){
	my $fqsd = `pwd`; chomp $fqsd;
   $fqsd .= "/golden";
   print "dirname = [$fqsd]\n";
	if (-e $fqsd && -d $fqsd) {
		print "$fqsd exists.  Continuing.\n";

	} else {
		print "Making dir [$fqsd] for 'golden data'. OK? [Ny]\n";
		my $tmp = <STDIN>;
		if ($tmp !~ "[Yy]") {die "OK - run this again to say where you want the dir.\n";}
		system("mkdir $fqsd"); # or die "Can't make [$fqsd] - check permissions\n";
	}
	$out = $fqsd;
} else {$out = "out";}
	my $test = "";


# # Test --raw
# 	$test = "raw";
# 	if ($out ne "out") {$out = "$out/$test";}
#    $cmdline = "$tacg --raw  -L -e 333 < $seq > ";
#    $MD5_hash = process_tacg_results($out, $cmdline, $test);


# Test --raw
	$test = "raw";
   $cmdline = "$tacg --raw  -L -e 333 < $seq > ";
   $MD5_hash = process_tacg_results($out, $cmdline, $test);

#$tmp = <STDIN>;
#print "...\n";

# Test -- input w/ infile with a multi-line fasta file
	$test = "raw_infile_fasta";
   $cmdline = "$tacg  -L -e 333 --infile $fastaseq > ";
   $MD5_hash = process_tacg_results($out, $cmdline, $test);
#$tmp = <STDIN>;

# Test --rev
	$test = "rev";
   $cmdline = "$tacg --rev -L -e3333   < $seq > ";
   $MD5_hash = process_tacg_results($out, $cmdline, $test);

# Test --comp
	$test = "comp";
   $cmdline = "$tacg --comp -L -e3333   < $seq > ";
   $MD5_hash = process_tacg_results($out, $cmdline, $test);

# Test --revcomp
	$test = "revcomp";
   $cmdline = "$tacg --revcomp -L -e3333   < $seq > ";
   $MD5_hash = process_tacg_results($out, $cmdline, $test);

# Test --notics & --strands 1 & numstart at a negative (also tests -T)
	$test = "notics_strands_numstart";
   $cmdline = "$tacg -M22 -lc -T3,1 -L --strands 1 --notics --numstart -3446  -b44444 -e99999 -w 120  < $seq > ";
   $MD5_hash = process_tacg_results($out, $cmdline, $test);

 # Test --notics & --strands 1 & numstart at a negative (also tests -T) --infile & FASTA
	$test = "notics_strands_numstart_infile_fasta";
   $cmdline = "$tacg -M22 -lc -T3,1 -L --strands 1 --notics --numstart -3446  -b44 -e888 -w 120  --infile $fastaseq > ";
   $MD5_hash = process_tacg_results($out, $cmdline, $test);


# Test -o (overhang)
	$test = "overhang_5";
   $cmdline = "$tacg  -M4 -S -lc -o5 -b44444 -e99999  < $seq > ";
   $MD5_hash = process_tacg_results($out, $cmdline, $test);

	$test = "overhang_3";
   $cmdline = "$tacg  -M4 -S -lc -o3 -b44444 -e99999  < $seq > ";
   $MD5_hash = process_tacg_results($out, $cmdline, $test);

	$test = "overhang_0";
   $cmdline = "$tacg  -M4 -S -lc -o0 -b44444 -e99999  < $seq > ";
   $MD5_hash = process_tacg_results($out, $cmdline, $test);

#if (0) { # loop over for short tests

# Test Fragments
for (my $f=0; $f<4; $f++) {
	$test = "Fragments_$f";
   $cmdline = "$tacg -n 4  -Lc -S  -T 6,3 -F $f -O'123456x,100' < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);
}

# Test Fragments & FASTA
for (my $f=0; $f<4; $f++) {
	$test = "FASTA_Fragments_$f";
   $cmdline = "$tacg -n 4  -Lc -S  -T 6,3 -F $f -O'123456x,100' < $fastaseq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);
}

# Test begin/end with several options
for (my $i=0; $i<5; $i++) {
	my $width = 60 + ($i * 15);
	$test = "begin_end_$i";
   $cmdline = "$tacg -b $begin[$i] -e $end[$i] -n 4 -Lc -S -T 6,3  -w $width < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);

#}


}
#Test ladder Map (-l) - some occassional hiccups here on some platforms
$test = "wide_ladder_map";
$cmdline = "$tacg -w 120 -n6 -l < $seq > ";
$MD5_hash = process_tacg_results($out, $cmdline, $test);

# Test extract (-X) option
for (my $i=0; $i<3; $i++) {
	$test = "extract_hindiii_bamhi_$i";
   $cmdline = "$tacg -b $begin[$i] -e $end[$i] -x'hindiii,bamhi' -X 11,17,0 < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);

	$test = "extract_regex_$i";
   $cmdline = "$tacg -b $begin[$i] -e $end[$i] --regex 'Rxtest:gc(cc|tac)nr{2,4}yt' -X 15,23,1 < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);
}

# Test other regex ops
$test = "regex_alt1";
$cmdline = "$tacg  --regex 'Rxtest2:gyt{3,6}n{2,5}ymk(ty|ghy)' -X 15,23,1 < $seq > ";
$MD5_hash = process_tacg_results($out, $cmdline, $test);

$test = "regex_FILE";
$cmdline = "$tacg  --regex 'FILE:$regex' -sl -S2 < $seq > ";
$MD5_hash = process_tacg_results($out, $cmdline, $test);


# Test cmdline pattern flag (-p), with errors
for (my $i=0; $i<3; $i++) {
	$test = "patterns_errors_$i";
   $cmdline = "$tacg -b $begin[$i] -e $end[$i] -p 'testpat,gcwwgtgatr,$i' -S < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);
}

if (0){

# Test --rule
	$test = "rule";
   $cmdline = "$tacg  --rule 'Rule1, \
   ((MwoI:1:8 & NaeI:1:8 & NarI:1:8) ^ ( NciI:1:8 & NcoI:1:8) & \
   (NdeI:1:8 ^ ((NgoMIV:1:8 & NheI:1:8)))),866' < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);

# Test --rulefile
	$test = "rulefile";
   $cmdline = "$tacg  --rulefile '$rulefile'  < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);
}


# Test matrix
   for (my $i=0; $i<15; $i++){
      my $co = $matrix_cutoffs[$i];
		$test = "matrix_test_$co";
      $cmdline = "$tacg -# $co -R $matrixfile -S  -b 43747 -e 94321 < $seq > ";
    $MD5_hash = process_tacg_results($out, $cmdline, $test);
   }
# Test dam & dcm
	$test = "dam_dcm_test";
   $cmdline = "$tacg --dam --dcm -n5 -S  -c < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);

# Test -g (gel)
	$test = "dcm_test";
   $cmdline = "$tacg  --dcm -n5 -S2 -c -O123x,88 < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);

# Test --clone
	$test = "clone_1";
   $cmdline = "$tacg  --clone '222_4444,5555_6666,44567_55555,55999_65888,4445x5554,6667x9999,11111x22222,88888x111111' < $seq > ";
   $MD5_hash = process_tacg_results($out, $cmdline, $test);

	$test = "clone_2";
   $cmdline = "$tacg  --clone '222_4444,5555_6666,44567_55555,55999_65888,4445x5554,6667x9999,11111x22222,88888x111111'  -e 22221 < $seq > ";
   $MD5_hash = process_tacg_results($out, $cmdline, $test);

	$test = "clone_3";
   $cmdline = "$tacg  --clone '222_4444,5555_6666,5557x6566,55999_65888,4445x5554,6667x9999,11111x22222,88888x111111' < $seq >& ";
   $MD5_hash = process_tacg_results($out, $cmdline, $test);

# Test --cost
	$test = "cost_10";
   $cmdline = "$tacg  -S -M4 --cost 10 -b44444 -e99999  < $seq > ";
   $MD5_hash = process_tacg_results($out, $cmdline, $test);

	$test = "cost_25";
   $cmdline = "$tacg  -s -m2 --cost 25 -b4344 -e9999  -w 120 < $seq > ";
   $MD5_hash = process_tacg_results($out, $cmdline, $test);

	$test = "cost_75";
   $cmdline = "$tacg  -S 2 --cost 75 -b4344 -e9999  < $seq > ";
   $MD5_hash = process_tacg_results($out, $cmdline, $test);

# Test -G (graphics)
my @opts = ("X", "Y", "L");
for (my $i=0; $i<3; $i++) {
	$test = "graphics_test_$opts[$i]";
   $cmdline = "$tacg -b 60000 -e 90000 -pjjj,mmmmmm,0 -G 200,$opts[$i]   < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);
}

# Test -C (Codons)
for (my $i=0; $i<17; $i++) {
	my $start = $i +1;
	$test = "Codon_tables_$i";
   $cmdline = "$tacg  -b $start -e 10000 -C $i -O 123456x,44  < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);
}

# Test topology - this one buggers up - diffs in 'diff'?
print "skipping topology test for now..\n";
# 	$test = "topology";
#    $cmdline  = "$tacg -f0 -n5 -e 33333  -S1 < $seq  > out1; ";
#    $cmdline .= "$tacg -f1 -n5 -e 33333  -S1 < $seq  > out2; ";
#    $cmdline .= "diff out1 out2 > $out/$test ";
# #   system('$cmdline; $cmdline1; diff out out1 > diff.f0.f1; mv diff.f0.f1 out');
#   $MD5_hash = process_tacg_results($out, $cmdline, $test);

# Test -P (Proximity)
	$test = "Proximity_1";
   $cmdline = "$tacg -b 5000 -e 25000 -P'hindiii,-l500,bamhi'  < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);

	$test = "Proximity_2";
   $cmdline = "$tacg -b 26000 -e 56000 -P'hindiii,+200-500,bamhi'  < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);

	$test = "Proximity_3";
   $cmdline = "$tacg -b 40000 -e 55000 -P'hindiii,-200-700,bamhi'  < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);

# Test --silent
	$test = "silent_short";
   $cmdline = "$tacg -e 550 --silent -LT1,1  < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);

	$test = "silent_long";
   $cmdline = "$tacg -b 550 -e 2535 --silent -LT1,1  < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);

# Test --ps & --pdf
# if (-x "/usr/bin/gs" || -x $ENV{TACG_GS_BIN}) {

#    $cmdline = "$tacg  -n6 -M2 --ps -e 10000 < $seq  > $out; mv tacg_Map.ps out";
# 	$MD5_hash = process_tacg_results($cmdline, "postscript_short");
#
#    $cmdline = "$tacg  -n6 -M2 --ps -e 50000 < $seq  > $out; mv tacg_Map.ps out";
# 	$MD5_hash = process_tacg_results($cmdline, "postscript_long");
#
#    $cmdline = "$tacg  -n6 -M2 --pdf -e 10000 < $seq  > $out; mv tacg_Map.pdf out";
# 	$MD5_hash = process_tacg_results($cmdline, "PDF_short");
#
#    $cmdline = "$tacg  -n5 -M2 --pdf -e 50000 -O123456,55 < $seq  > $out; mv tacg_Map.pdf out";
# 	$MD5_hash = process_tacg_results($cmdline, "PDF_long");
# }	else {
# 	print "You don't have an executable ghostscript installed at '/usr/bin/gs.'\n",
#    "Skipping the postscript & PDF test.  If it's installed elsewhere,\n",
#    "symlink it to '/usr/bin/gs'\n";
# }

# Test -O

	$test = "ORF1_short_orfmap";
   $cmdline = "$tacg -b 5000 -e 25000 -O123456x,30  --orfmap < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);

	$test = "ORF2_medium";
   $cmdline = "$tacg  -e 5000 -O123456x,60  < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);

	$test = "ORF3_long_orfmap";
   $cmdline = "$tacg  -O123x,150  --orfmap < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);

	$test = "ORF6xS_44_orfmap";
   $cmdline = "$tacg   -b 14 -e 10000 -w 90 -O 123456x,44  --orfmap  < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);

	$test = "ORF6xL_44_orfmap";
   $cmdline = "$tacg   -b 14  -w 134 -O 123456x,44  --orfmap  < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);

	$test = "ORF6xL_150_orfmap";
   $cmdline = "$tacg  -O123x,150 -w 120 --orfmap < $seq > ";
 $MD5_hash = process_tacg_results($out, $cmdline, $test);

# print out Result Summary at end, of all results
	print "\n\n   ======================= Summary of Results =======================\n";
	print "                             Test             Status\n";
	print "                -------------------------------------------\n";


	foreach my $key (sort(keys %result_hash)) {
		printf "%40s    ", $key;
		if ($result_hash{$key} eq "EXEC"){print "Exec err\n";}
		elsif ($result_hash{$key} eq "MD5") {print "md5 ok\n";}
		elsif ($result_hash{$key} eq "WC")  {print "wc ok\n";}
		elsif ($result_hash{$key} eq "FAIL"){print "! Different Output !\n";}
	}


   print << "CLEANUP";

   The output for --ps and --pdf are not tested here as they both
   reference the date which will change and the PDF conversion seems to be
   quite variable.

      If you have problems with the postscript output, remember that the
   new output is APPENDED to any current file named 'tacg_Map.ps'.
      Also you may have font incompatibilities or missing fonts which may
   result in the ghostscript interpreter crashing as well.
      Also, the postscript routines are the least well-debugged so you
   might just be hitting a true bug..
      To verify the PDF and postscript output, take a look at the resulting
   pages with 'gv' or the acrobat reader.


   Don't forget to clean out the test directory after running the tests
   as they may generate a huge amount of data that you probably will never
   need again.  The only thing you might need again is the 'testtacg.pl'
   script, and the '.err' files - the output of the run if it doesn't
   match with what I've determinied to the REAL output on my dev machine.
   If you're concerned about the mismatch, please mail me (hjm\@tacgi.com)
   to ask if I need the output and we can arrange the exchange.

   Incidentally, this test took ~32s to run my dev machine (IBM thinkpad a22p
   - 900MHz PentiumIII/Coppermine); ~14 on a logy T60 (1GHz).  You can time
   yours by invoking it as:
   'touch tacg; time make test'

CLEANUP


sub process_tacg_results {
	 my ($out) = shift;
	 my ($cmdline) = shift;
    my ($testtype) = shift;
    if ($diff_test ne "" and $testtype ne $diff_test) {
    	#print " $testtype skipped...\n";
    	return "NULL";
    }
    if ($out ne "out") {$out = "$out/$test";}
 	 $cmdline .= $out;
#    print "about to exec:  $cmdline \n";
	 my $return_value = system("$cmdline");; # execute the command
	 $return_value /= 256;
#	 print "return value = $return_value\n";
    my $file = $out;
    my $hash = "";
    my $wc = "";
    my $errfile = "";


   if ($return_value > 1) {
   	$result_hash{$testtype} = "EXEC FAILURE"; # execution failure
      print LOG "Failed commandline = $cmdline\n";
      print "EXEC FAILURE:     $testtype\n";
      print "cmd: $cmdline \n";
#      $hash = "UNDEFINED";
#      return($hash);
   } else {
      open(FILE, $out) or die "Can't open $out: $!";
      binmode(FILE);
      $wc = `wc out`;
#      print "wc before: [$wc]\n";
      $wc =~ s/\s+/\ /g;

      chomp $wc;chop $wc;
#      print " wc after: [$wc]\n\n";

   #      print LOG "wc of out = $wc\n";
      if ($md5found == 1) { $hash = Digest::MD5->new->addfile(*FILE)->hexdigest ;}
      else                { $hash = "UNDEFINED"; }
      print LOG "\"$testtype\" => ", "\"$hash\", #MD5\n\"$testtype\" => ", "\"$wc\", #wc\n$cmdline\n\n";
   }
   if ($md5found == 1) {
      if ( $MD5_table{$testtype} eq $hash ) {
         print "passed test (MD5): $testtype\n";
         $result_hash{$testtype} = "MD5"; # passed MD5 check
      } else {
         print "\n!!FAILED!! test: $testtype\n";
         print "the MD5 sig: $hash \n  should be: $MD5_table{$testtype}\n";
         if ($wc eq $wcTable{$testtype}) {
             print "however the wordcount matches:\n";
             $result_hash{$testtype} = "WC"; # passed WC check
         } else {
             print "and the wordcount is off as well:\n";
             $result_hash{$testtype} = "FAIL"; # failed both MD5 & WC check
         }
         print "the wc sig: $wc \n should be: $wcTable{$testtype}\n";
         print "cmd: $cmdline \n\n";
         $errfile = "$testtype" . ".MD5.err";
         if (!$save) { system("mv out $errfile");}
      }
   } else {
      if ($wcTable{$testtype} eq $wc) {
         print "passed test (wc): $testtype\n";
         $result_hash{$testtype} = "WC"; # passed WC check
      } else {
         print "\n!!FAILED!! test: $testtype\n";
         $result_hash{$testtype} = "FAIL"; # failed both MD5 & WC check
         print "the wc sig: [$wc] \n should be: [$wcTable{$testtype}]\n";
         print "cmd: $cmdline \n\n";
         $errfile = "$testtype" . ".wc.err";
         if (!$save) { system("mv out $errfile"); }
      }
   }
	# then match the test and diff the current output (should be 'out') and
	# the previous 'golden' output
	# $diff_test has to be something and the only result to test that makes sense is if it's
	# equal to 'WC' or 'FAIL' (not a perfect match and not a seqfault)
	if ($diff_test ne "" && ($result_hash{$diff_test} eq "WC" || $result_hash{$diff_test} eq "FAIL")) {
		# then just diff|kompare|tkdiff it
#		print "[$diff_app $errfile golden/$diff_test]\n";
		system("$diff_app $errfile golden/$diff_test ");
		print "Want to start a kompare of the files? [yN]";
		$tmp = <STDIN>;
		if ($tmp =~ "[Yy]") {
			exec("/usr/bin/kompare $errfile golden/$diff_test ") or die "Can't kompare...\n";
		}
	}

   return $hash;
} # end process_tacg_results()

