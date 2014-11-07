#!/usr/bin/perl
use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta; #makes a searchable db from my fasta file
###############################################################################
##########       Make ChLG FASTA from scaffolds to ChLG AGP           #########
###############################################################################
# Make DB of scaffolds
my $fasta_in_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/Tcas5.2_scaffolds.fasta";
#my $fasta_in_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/test.fa";
my $db = Bio::DB::Fasta->new("$fasta_in_file");
# Make ChLG FASTA
my $fasta_out_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/Tcas5.2_chlg_pre_header.fasta";
my $seq_out = Bio::SeqIO->new('-file' => ">$fasta_out_file",'-format' => 'fasta');		#Create new fasta outfile object.
# open final contigs to scaffolds AGP
my $chlg_from_scaffolds_in_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/Tcas5.2_chlg_from_scaffolds.agp";
open (my $chlg_from_scaffolds_in, "<", $chlg_from_scaffolds_in_file) or die "can't open $chlg_from_scaffolds_in_file: $!";
my $first=1;
my $old_mol="X";
my $scaffold_id;
my $new_seq  = ''; ### initialize first superscaffold

my $seq_in = Bio::SeqIO->new(-file => "<$fasta_in_file", -format => 'fasta');

while(<$chlg_from_scaffolds_in>)
{
    unless (/^#/)
    {
        my @columns = (split(/\t/));
        my $new_mol=$columns[0];
        if ($columns[4] eq "W")
        {
            ###################################################################
            #################    starting/changing ChLGs   ####################
            ###################################################################
            if ($new_mol ne $old_mol) ## if we are not on the same molecule
            {
                unless ($first==1) # unless this is the first molecule map write the old one
                {
                    my $scaffold_obj = Bio::Seq->new( -display_id =>  $scaffold_id, -seq => $new_seq, -alphabet => 'dna');
                    $seq_out->write_seq($scaffold_obj); ## write the finished superscaffold
                }
                $scaffold_id = "$columns[0]"; ## initialize new superscaffold
                $new_seq = '';
                $first=0;
            }
            ###################################################################
            ####    Continue building superscaffolds: append gaps    ##########
            ###################################################################
            if ($old_mol eq $new_mol)
            {
                $new_seq = "$new_seq"."n" x 100; ## add n's if the last contig is on the same molecule
            }
            ###################################################################
            ####              append scaffolds to ChLGs              ##########
            ###################################################################
            my ($start,$stop) = ($columns[6], $columns[7]);
            my $id = "$columns[5]";
            $new_seq = "$new_seq".$db->seq("$id:$start,$stop"); ## add the new sequence to the growing superscaffold
            $old_mol=$new_mol; ## now the current molecule will be listed as the last molecule we have seen
            if (eof) ## if this is the last row in the stitchmap table
            {
                my $scaffold_obj = Bio::Seq->new( -display_id =>  $scaffold_id, -seq => $new_seq, -alphabet => 'dna');
                $seq_out->write_seq($scaffold_obj); ## Write the final sequence object
            }
        }
    }
}

#my $fasta_out_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/Tcas5.2_chlg_pre_header.fasta";
open (my $fasta_out, "<", $fasta_out_file) or die "can't open $fasta_out_file: $!";

my $new_fasta_out_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/Tcas5.2_chlg.fasta";

open (my $new_fasta_out, ">", $new_fasta_out_file) or die "can't open $new_fasta_out_file: $!";

while(<$fasta_out>)
{
    if (/^>/)
    {
        chomp;
        print $new_fasta_out "$_ [organism=Tribolium castaneum] [strain=Georgia GA2] [country=USA: Kansas] [collection-date=Apr-2003]\n";
    }
    else
    {
        print $new_fasta_out "$_";
    }
}
unlink $fasta_out_file;

