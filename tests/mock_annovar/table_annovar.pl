#!/usr/bin/env perl
# Mock of ANNOVAR's table_annovar.pl for tests and CI (no licence needed).
# $Date: 2099-01-01 mock $
#
# Reproduces the output layout of `table_annovar.pl --vcfinput --otherinfo`:
#   <outfile>.<buildver>_multianno.txt  Chr Start End Ref Alt <annotations> Otherinfo1..N
#   <outfile>.<buildver>_multianno.vcf  input VCF with ANNOVAR fields added to INFO
#   <outfile>.avinput
# Annotation values are deterministic placeholders derived from the position.

use strict;
use warnings;
use Getopt::Long;

my ($buildver, $outfile, $protocol, $operation, $nastring) = ('hg19', 'out', '', '', '.');
my ($vcfinput, $otherinfo, $remove, $onetranscript, $thread);
GetOptions(
    'buildver=s'  => \$buildver,  'outfile=s'   => \$outfile,
    'protocol=s'  => \$protocol,  'operation=s' => \$operation,
    'nastring=s'  => \$nastring,  'vcfinput'    => \$vcfinput,
    'otherinfo'   => \$otherinfo, 'remove'      => \$remove,
    'onetranscript' => \$onetranscript, 'thread=i' => \$thread,
) or die "mock table_annovar.pl: bad arguments\n";
my ($query, $dbdir) = @ARGV;
die "mock table_annovar.pl: usage: table_annovar.pl <query> <humandb> [options]\n" unless $query && $dbdir;
die "mock table_annovar.pl: --vcfinput is required\n" unless $vcfinput;

my @dbs = split /,/, $protocol;
my @ops = split /,/, $operation;
die "mock table_annovar.pl: protocol/operation length mismatch\n" unless @dbs == @ops;
for my $db (@dbs) {
    -e "$dbdir/${buildver}_$db.txt" or die "Error: the database file $dbdir/${buildver}_$db.txt does not exist\n";
}

# Annotation column names per database
my @anno_cols;
for my $i (0 .. $#dbs) {
    my ($db, $op) = ($dbs[$i], $ops[$i]);
    if ($op =~ /^g/) {
        push @anno_cols, map { "$_.$db" } qw(Func Gene GeneDetail ExonicFunc AAChange);
    } elsif ($db =~ /^clinvar/) {
        push @anno_cols, qw(CLNALLELEID CLNDN CLNDISDB CLNREVSTAT CLNSIG);
    } elsif ($db eq 'revel') {
        push @anno_cols, 'REVEL';
    } elsif ($db =~ /^pubmind/) {
        push @anno_cols, qw(PVID pathogenicity_sum paper_level_pathogenicity_score confidence);
    } else {
        push @anno_cols, $db;
    }
}

my @funcs  = qw(exonic intronic intergenic UTR3 splicing exonic);
my @exfunc = ('nonsynonymous SNV', '.', '.', '.', '.', 'synonymous SNV');
my @clnsig = qw(Pathogenic Benign Likely_pathogenic Uncertain_significance .);
my @pubmind = ('Pathogenic', 'Uncertain significance', 'Conflicting', 'Likely benign');

open my $in,  '<', $query or die "cannot read $query: $!\n";
open my $txt, '>', "$outfile.${buildver}_multianno.txt" or die $!;
open my $vcf, '>', "$outfile.${buildver}_multianno.vcf" or die $!;
open my $av,  '>', "$outfile.avinput" or die $!;

my $header_done = 0;
while (my $line = <$in>) {
    chomp $line;
    if ($line =~ /^##/) { print $vcf "$line\n"; next; }
    if ($line =~ /^#CHROM/) {
        print $vcf qq(##INFO=<ID=ANNOVAR_DATE,Number=1,Type=String,Description="Flag the start of ANNOVAR annotation for one alternative allele">\n);
        print $vcf "##INFO=<ID=$_,Number=.,Type=String,Description=\"$_ annotation provided by ANNOVAR\">\n" for @anno_cols;
        print $vcf qq(##INFO=<ID=ALLELE_END,Number=0,Type=Flag,Description="Flag the end of ANNOVAR annotation for one alternative allele">\n);
        print $vcf "$line\n";
        next;
    }
    my @f = split /\t/, $line;
    my ($chr, $pos, $ref, $alt) = @f[0, 1, 3, 4];

    # VCF -> ANNOVAR coordinates (bi-allelic input)
    my ($start, $end, $aref, $aalt) = ($pos, $pos + length($ref) - 1, $ref, $alt);
    if (length($ref) > 1 && length($alt) == 1 && substr($ref, 0, 1) eq $alt) {        # deletion
        ($start, $end, $aref, $aalt) = ($pos + 1, $pos + length($ref) - 1, substr($ref, 1), '-');
    } elsif (length($alt) > 1 && length($ref) == 1 && substr($alt, 0, 1) eq $ref) {  # insertion
        ($start, $end, $aref, $aalt) = ($pos, $pos, '-', substr($alt, 1));
    }

    my @vals;
    my $k = $pos % 6;
    for my $col (@anno_cols) {
        if    ($col =~ /^Func\./)       { push @vals, $funcs[$k] }
        elsif ($col =~ /^Gene\./)       { push @vals, "GENE$k" }
        elsif ($col =~ /^ExonicFunc\./) { push @vals, $exfunc[$k] }
        elsif ($col eq 'CLNSIG')        { push @vals, $clnsig[$pos % 5] }
        elsif ($col eq 'REVEL')         { push @vals, sprintf('%.3f', ($pos % 100) / 100) }
        elsif ($pos % 3)                 { push @vals, $nastring }            # PubMind: every third position
        elsif ($col eq 'PVID')          { push @vals, "PVID$pos" }
        elsif ($col eq 'pathogenicity_sum') { push @vals, $pubmind[$pos % 4] }
        elsif ($col eq 'paper_level_pathogenicity_score') { push @vals, $pos % 4 == 1 ? '' : sprintf('%.2f', ($pos % 7) / 7) }
        elsif ($col eq 'confidence')    { push @vals, $pos % 4 }
        else                            { push @vals, $nastring }
    }

    my @other = ('0.5', $f[5], '10', @f);   # zygosity/AF, QUAL, DP, full VCF line
    if (!$header_done) {
        print $txt join("\t", qw(Chr Start End Ref Alt), @anno_cols, map { "Otherinfo$_" } 1 .. scalar @other), "\n";
        $header_done = 1;
    }
    print $txt join("\t", $chr, $start, $end, $aref, $aalt, @vals, @other), "\n";
    print $av  join("\t", $chr, $start, $end, $aref, $aalt, @other), "\n";

    my $info = join(';', 'ANNOVAR_DATE=2099-01-01', (map { "$anno_cols[$_]=$vals[$_]" } 0 .. $#anno_cols), 'ALLELE_END');
    $f[7] = ($f[7] eq '.' ? '' : "$f[7];") . $info;
    $f[7] =~ s/ /_/g;
    print $vcf join("\t", @f), "\n";
}
print $txt join("\t", qw(Chr Start End Ref Alt), @anno_cols, 'Otherinfo1'), "\n" unless $header_done;
close $_ for ($in, $txt, $vcf, $av);
print STDERR "NOTICE: mock table_annovar.pl finished ($outfile.${buildver}_multianno.txt)\n";
