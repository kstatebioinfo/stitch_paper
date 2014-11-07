#!/usr/bin/perl
use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta; #makes a searchable db from my fasta file
use Text::Wrap;

## 3)Find placed by header
### 3b) Reverse : Super_scaffolds and Flip AGPs where needed (e.g. stitch built backwards)
## 4) Find rest if not placed
## 5) Cat both
## 6) Make new raw AGP
## 7) Get type of gap for those in new AGP that are from BNG
## 8) Write new ChLG with 100 bp spacers and genetic map as the type

## ISSUE with collapse AGP?
#Super_scaffold_63	1	847447	1	W	scaffold_99	1	847447	+
#Super_scaffold_63	847448	851048	2	N	3601	scaffold	yes	map
#Super_scaffold_63	851049	2184743	3	W	scaffold_100	1	1333695	-
#Super_scaffold_63	2184744	2265261	4	N	80518	scaffold	yes	map
#Super_scaffold_63	2265262	2911381	5	W	scaffold_101	1	646120	+
#Super_scaffold_63	2911382	2944542	6	N	33161	scaffold	yes	map
#Super_scaffold_63	2944543	3336913	7	W	scaffold_102	1	392371	-
#Super_scaffold_63	3336914	3782898	8	N	445985	scaffold	yes	map
#Super_scaffold_63	3782899	4456720	9	W	scaffold_103	1	673822	-
#
#Super_scaffold_55	1	324174	1	W	scaffold_1203	1	324174	+
#Super_scaffold_55	324175	572800	2	N	248626	scaffold	yes	map
#Super_scaffold_55	572801	883697	3	W	scaffold_3	1	310897	-
#Super_scaffold_55	883698	917448	4	N	33751	scaffold	yes	map
#Super_scaffold_55	917449	1131550	5	W	scaffold_4	1	214102	-
#Super_scaffold_55	1131551	1171393	6	N	39843	scaffold	yes	map
#Super_scaffold_55	1171394	1634415	7	W	scaffold_5	1	463022	+
#Super_scaffold_55	1634416	1675674	8	N	41259	scaffold	yes	map
#Super_scaffold_55	1675675	1897757	9	W	scaffold_6	1	222083	-
#Super_scaffold_55	1897758	2230651	10	N	332894	scaffold	yes	map
#Super_scaffold_55	2230652	3389387	11	W	scaffold_7	1	1158736	+
#Super_scaffold_55	3389388	3497099	12	N	107712	scaffold	yes	map
#Super_scaffold_55	3497100	3749323	13	W	scaffold_977	1	252224	+
#Super_scaffold_55	3749324	3786090	14	N	36767	scaffold	yes	map
#Super_scaffold_55	3786091	3887837	15	W	scaffold_8	1	101747	+
#Super_scaffold_55	3887838	3913286	16	N	25449	scaffold	yes	map
#Super_scaffold_55	3913287	4286656	17	W	scaffold_9	1	373370	-
#Super_scaffold_55	4286657	4303960	18	N	17304	scaffold	yes	map
#Super_scaffold_55	4303961	4461328	19	W	scaffold_10	1	157368	-
#Super_scaffold_55	4461329	4506698	20	N	45370	scaffold	yes	map
#Super_scaffold_55	4506699	5870923	21	W	scaffold_11	1	1364225	-
#Super_scaffold_55	5870924	5945417	22	N	74494	scaffold	yes	map
#Super_scaffold_55	5945418	7260238	23	W	scaffold_12	1	1314821	-
#Super_scaffold_55	7260239	7260338	24	U	100	scaffold	yes	map
#Super_scaffold_55	7260339	8421036	25	W	scaffold_13	1	1160698	-
my $app_in_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/pre_Tcas5.2_1.agp";
open (my $agp_in , "<", $app_in_file) or die "can't open $app_in_file: $!";
my $agp_out_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/pre_Tcas5.2_1.5.agp";
open (my $agp_out , ">", $agp_out_file) or die "can't open $agp_out_file: $!";
print $agp_out "##agp-version   2.0\n";

my %agp;
while (<$agp_in>)
{
#    chomp;
    my $object = (split(/\t/))[0];
    if ($agp{$object})
    {
         $agp{$object} .= $_;
    }
    else
    {
        $agp{$object} = $_;
    }
}

my $fasta_in_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/pre_Tcas5.2_1.fasta";
my $db = Bio::DB::Fasta->new("$fasta_in_file");
my $fasta_out_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/pre_Tcas5.2_2.fasta";
my $seq_out = Bio::SeqIO->new('-file' => ">$fasta_out_file",'-format' => 'fasta');		#Create new fasta outfile object.

#$reversed_obj = $seq_obj->revcom;

my %do_reverse = ("Super_scaffold_61"=>1, "Super_scaffold_26"=>1, "Super_scaffold_21"=>1, "Super_scaffold_34"=>1,"Super_scaffold_37"=>1, "Super_scaffold_36"=>1);

my $chlg_in_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/ChLG_order.txt";
open (my $chlg_in, "<", $chlg_in_file) or die "can't open $chlg_in_file: $!";
my %chlg;
###############################################################################
######     Make AGP and FASTA for placed or superscaffoled unplaced   #########
###############################################################################

while (<$chlg_in>)
{
    unless(/^#/)
    {
        chomp;
        $chlg{$_} =1; #make list of placed
        if ($do_reverse{$_})
        {
            my @agps = split(/\n/,$agp{$_});
            for my $line (@agps) #flip orientations
            {
                if ($line =~ /\+/)
                {
                    $line =~ s/\+/-/;
                }
                elsif ($line =~ /-/)
                {
                    $line =~ s/-/\+/;
                }
                
            }
            $agp{$_} = join("\n",@agps);
            my $agp_lines = (join("\n",(reverse(split(/\n/,$agp{$_})))));#reverse agp
            print $agp_out "$agp_lines\n"; # print to agp
            
            my $seq_obj = $db->get_Seq_by_id($_); #get fasta record
            my $reversed_obj = $seq_obj->revcom; #reverse fasta record
            $seq_out->write_seq($reversed_obj); # print fasta record
        }
        else
        {
            if ($agp{$_})
            {
                print $agp_out "$agp{$_}"; # print to agp
            }
            else
            {
                my $seq_obj = $db->get_Seq_by_id($_); #get fasta record
                my $length = $seq_obj->length;
                print $agp_out "$_\t1\t$length\t1\tW\t$_\t1\t$length\t+\n";
            }
            my $seq_obj = $db->get_Seq_by_id($_); #get fasta record
            $seq_out->write_seq($seq_obj); # print fasta record
        }
    }
}
close $agp_out;
close $chlg_in;
###############################################################################
##############            Repair AGP coordinates             ##################
###############################################################################
open (my $new_input , "<", $agp_out_file) or die "can't open $agp_out_file: $!";
#/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/pre_Tcas5.2_1.5.agp
my $output_agp = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/pre_Tcas5.2_1.6.agp";
open (my $new_output, ">", $output_agp) or die "Can't open $output_agp: $!";
print $new_output "##agp-version   2.0\n";
my $previous_object = "first";
my $current_object_beg = 1;
my $current_part_number = 1;
while (<$new_input>)
{
    unless ((/^#/)||(/^\s*$/))
    {
        chomp;
        my ($object,$object_beg,$object_end,$part_number,$component_type,$gap_length  ,$gap_type     ,$linkage      ,$Linkage_evidence,$component_id,$component_beg,$component_end,$orientation);
        $object = (split(/\t/))[0];
        $component_type = (split(/\t/))[4];
        unless ($object eq "start")
        {
            if ($previous_object ne $object)
            {
                $current_object_beg = 1;
                $current_part_number = 1;
            }
        }
        if ($component_type eq "W")
        {
            ($object,$object_beg,$object_end,$part_number,$component_type,$component_id,$component_beg,$component_end,$orientation     ) = split(/\t/);
            $object_end = $current_object_beg + $component_end - 1;
            print $new_output "$object\t$current_object_beg\t$object_end\t$current_part_number\t$component_type\t$component_id\t$component_beg\t$component_end\t$orientation\n";
        }
        else
        {
            ($object,$object_beg,$object_end,$part_number,$component_type,$gap_length  ,$gap_type     ,$linkage      ,$Linkage_evidence) = split(/\t/);
            $object_end = $current_object_beg + $gap_length - 1;
            print $new_output "$object\t$current_object_beg\t$object_end\t$current_part_number\t$component_type\t$gap_length\t$gap_type\t$linkage\t$Linkage_evidence\n";
        }
        $previous_object = $object; #reset previous
        $current_object_beg = $object_end + 1;
        $current_part_number = $current_part_number + 1;
    }
    
}

###############################################################################
##############       Add unplaced objects back to FASTA      ##################
###############################################################################
my $seq_in = Bio::SeqIO->new(-file => "<$fasta_in_file", -format => 'fasta');

while (my $seq = $seq_in->next_seq)
{
    my $id = $seq->display_id();
    #    print "$id\n";
    
    unless ($chlg{$id})
    {
        $seq_out->write_seq($seq);
        my $seq_obj = $db->get_Seq_by_id($id); #get fasta record
        my $length = $seq_obj->length;
        print $new_output "$id\t1\t$length\t1\tW\t$id\t1\t$length\t+\n";
    }
    
    
}
close $new_output;
#my $output_agp = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/pre_Tcas5.2_1.6.agp";
open (my $new_input_agp, "<", $output_agp) or die "Can't open $output_agp: $!";
###############################################################################
##############           Make raw contig AGP and fasta       ##################
###############################################################################
my $make_contigs = `perl /homes/bioinfo/bioinfo_software/bionano/Irys-scaffolding/KSU_bioinfo_lab/stitch/make_contigs_from_fasta.pl /homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/pre_Tcas5.2_2.fasta`;
print "$make_contigs\n";
###############################################################################
##############      add BNG map info to raw contig AGP       ##################
###############################################################################
open (my $bng_agp, "<", $output_agp) or die "can't open $output_agp: $!";
my %mapped_gaps;
while (<$bng_agp>)
{
    unless(/^#/)
    {
        my @mapped_gap = (split(/\t/))[0,1,2];
        my $component_type = (split(/\t/))[4];
        my $start = join("\t",@mapped_gap);
        $start .= "\t";
        $mapped_gaps{$start} = $component_type;
    }
}
close $bng_agp;
my $raw_contig_agp_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/pre_Tcas5.2_2.fasta_contig.agp";
open (my $raw_contig_agp, "<", $raw_contig_agp_file) or die "can't open $raw_contig_agp_file:$!";
my $bng_contig_agp_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/pre_Tcas5.2_2_bng_contig.agp";
open (my $bng_contig_agp, ">", $bng_contig_agp_file) or die "can't open $bng_contig_agp_file:$!";
print $bng_contig_agp "##agp-version   2.0\n";
print $bng_contig_agp "# ORGANISM: Tribolium castaneum\n";
print $bng_contig_agp "# TAX_ID: 7070\n";
print $bng_contig_agp "# ASSEMBLY NAME: Tcas_5.2\n";
print $bng_contig_agp "# ASSEMBLY DATE: 11-November-2014\n";
print $bng_contig_agp "# GENOME CENTER: Bioinformatics Center at Kansas State University\n";
print $bng_contig_agp "# DESCRIPTION: AGP specifying the assembly of scaffolds from WGS contigs\n";
while (<$raw_contig_agp>)
{
    unless (/^#/)
    {
        chomp;
        my @gap = (split(/\t/))[0,1,2];
        my $start = join("\t",@gap);
        $start .= "\t";
        my $component_type = (split(/\t/))[4];
        if (($mapped_gaps{$start}) && ($component_type ne "W"))
        {
            my @columns= split(/\t/);
            $columns[4] = "$mapped_gaps{$start}";
            $columns[8] = "map";
            my $line = join("\t",@columns);
            print $bng_contig_agp "$line\n";
        }
        else
        {
            print $bng_contig_agp "$_\n";
        }
    }
}
###############################################################################
##############             Make raw ChLG AGP                 ##################
###############################################################################
my $chlg_agp_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/pre_Tcas5.2_2_chlg.agp";
open (my $chlg_agp, ">", $chlg_agp_file) or die "can't open $chlg_agp_file:$!";
print $chlg_agp "##agp-version   2.0\n";
print $chlg_agp "# ORGANISM: Tribolium castaneum\n";
print $chlg_agp "# TAX_ID: 7070\n";
print $chlg_agp "# ASSEMBLY NAME: Tcas_5.2\n";
print $chlg_agp "# ASSEMBLY DATE: 11-November-2014\n";
print $chlg_agp "# GENOME CENTER: Bioinformatics Center at Kansas State University\n";
print $chlg_agp "# DESCRIPTION: AGP specifying the assembly of ChLGs from WGS scaffolds\n";
my $current_chlg;
my $pos = 1;
my $agp_element=1;
my @current_chlg_agp;
my %seen;
my $placed = 1;
$chlg_in_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/ChLG_order.txt";
open ($chlg_in, "<", $chlg_in_file) or die "can't open $chlg_in_file: $!";
while (<$chlg_in>)
{
    chomp;
    if ($placed == 1)
    {
        if(/^#ChLGX/)
        {
            s/#//;
            $current_chlg = $_;
            $pos = 1;
            $agp_element=1;
        }
        elsif (/^#/)
        {
            pop(@current_chlg_agp);
            print $chlg_agp join('',@current_chlg_agp);
            s/#//;
            $current_chlg = $_;
            $pos = 1;
            $agp_element=1;
            @current_chlg_agp=();
        }
        else
        {
            $seen{$_} = 1; # log this sequence as added
            my $seq_obj = $db->get_Seq_by_id($_); #get fasta record
            my $length = $seq_obj->length;
            my $end_pos = $pos + $length - 1;
            unless( $current_chlg eq "ChLGY")
            {
                push (@current_chlg_agp, "$current_chlg\t$pos\t$end_pos\t$agp_element\tW\t$_\t1\t$length\t+\n");
                $pos = $end_pos + 1;
                ++$agp_element;
                $end_pos = $pos + 99;
                push (@current_chlg_agp, "$current_chlg\t$pos\t$end_pos\t$agp_element\tU\t100\tscaffold\tyes\tmap\n");
                $pos = $end_pos + 1;
                ++$agp_element;
            }
            else
            {
                push (@current_chlg_agp, "$current_chlg\t$pos\t$end_pos\t$agp_element\tW\t$_\t1\t$length\t?\n");
                $pos = $end_pos + 1;
                ++$agp_element;
                $end_pos = $pos + 99;
                push (@current_chlg_agp, "$current_chlg\t$pos\t$end_pos\t$agp_element\tU\t100\tcontig\tna\tmap\n");
                $pos = $end_pos + 1;
                ++$agp_element;
            }
        }
        if (/Unplaced/)
        {
            $placed = 0;
        }
    }
    else
    {
        $seen{$_} = 1; # log this sequence as added
        my $seq_obj = $db->get_Seq_by_id($_); #get fasta record
        my $length = $seq_obj->length;
        my $end_pos = $pos + $length - 1;
        print $chlg_agp "unplaced_seqid1\t$pos\t$end_pos\t1\tW\t$_\t1\t$length\t?\n";
    }
}

###############################################################################
##############      Add unplaced objects back to AGP         ##################
###############################################################################
$seq_in = Bio::SeqIO->new(-file => "<$fasta_in_file", -format => 'fasta');
my $unplaced_count = 2;
while (my $seq = $seq_in->next_seq)
{
    my $id = $seq->display_id();
    #    print "$id\n";
    
    unless ($seen{$id})
    {
#        $seq_out->write_seq($seq);
        my $seq_obj = $db->get_Seq_by_id($id); #get fasta record
        my $length = $seq_obj->length;
        print $chlg_agp "unplaced_seqid${unplaced_count}\t1\t$length\t1\tW\t$id\t1\t$length\t?\n";
        ++$unplaced_count;
    }
    
}
close $chlg_agp;
###############################################################################
##############         Set final names for contigs           ##################
###############################################################################
# open current contigs
my $contig_in_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/pre_Tcas5.2_2.fasta_contig.fasta";
open (my $contig_in, "<", $contig_in_file) or die "can't open $contig_in_file: $!";
# open final contigs
my $contig_out_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/Tcas5.2_contigs.fasta";
open (my $contig_out, ">", $contig_out_file) or die "can't open $contig_out_file: $!";
# open contig key
my $contig_keys_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/Tcas5.2_contigs_key.txt";
open (my $contig_keys, ">", $contig_keys_file) or die "can't open $contig_keys_file: $!";
print $contig_keys "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/pre_Tcas5.2_2.fasta_contig.fasta to /homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/Tcas5.2_contigs.fasta\n";
# replace contig headers
my $contig_counter = 1;
my %contig_key;
while (<$contig_in>)
{
    if (/^>/)
    {
        #>contig_seqid1 [organism=Tribolium castaneum] [strain=Georgia GA2] [country=USA: Kansas] [collection-date=Apr-2003]
        chomp;
        s/>//;
        my $new_header = "contig_seqid${contig_counter}";
        print $contig_out ">${new_header} [organism=Tribolium castaneum] [strain=Georgia GA2] [country=USA: Kansas] [collection-date=Apr-2003]\n";
        $contig_key{$_} = $new_header;
        print $contig_keys "$_\t$new_header\n";
        ++$contig_counter;
    }
    else
    {
        
#        print $contig_out "$_";
        chomp;
        $Text::Wrap::columns = 60;
        print $contig_out wrap('', '', $_) . "\n";
    }
}

###############################################################################
##############         Set final names for scaffolds         ##################
###############################################################################
# open current contigs
my $scaffold_in_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/pre_Tcas5.2_2.fasta";
open (my $scaffold_in, "<", $scaffold_in_file) or die "can't open $scaffold_in_file: $!";
# open final contigs
my $scaffold_out_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/Tcas5.2_scaffolds.fasta";
open (my $scaffold_out, ">", $scaffold_out_file) or die "can't open $scaffold_out_file: $!";
# open contig key
my $scaffold_keys_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/Tcas5.2_scaffolds_key.txt";
open (my $scaffold_keys, ">", $scaffold_keys_file) or die "can't open $scaffold_keys_file: $!";
print $scaffold_keys "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/pre_Tcas5.2_2.fasta to /homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/Tcas5.2_scaffolds.fasta\n";
# replace contig headers
my $scaffold_counter = 1;
my %scaffold_key;
while (<$scaffold_in>)
{
    if (/^>/)
    {
        #>scaffold_seqid1 [organism=Tribolium castaneum] [strain=Georgia GA2] [country=USA: Kansas] [collection-date=Apr-2003]

        chomp;
        s/>//;
        my $new_header = "scaffold_seqid${scaffold_counter}";
        print $scaffold_out ">${new_header} [organism=Tribolium castaneum] [strain=Georgia GA2] [country=USA: Kansas] [collection-date=Apr-2003]\n";
        $scaffold_key{$_} = $new_header;
        print $scaffold_keys "$_\t$new_header\n";
        ++$scaffold_counter;
    }
    else
    {
        print $scaffold_out "$_";
    }
}

###############################################################################
##########       Set final names for contigs to scaffold AGP          #########
###############################################################################
# open final contigs to scaffolds AGP
my $scaffolds_from_contigs_out_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/Tcas5.2_scaffolds_from_contigs.agp";
open (my $scaffolds_from_contigs_out, ">", $scaffolds_from_contigs_out_file) or die "can't open scaffolds_from_contigs_out_file: $!";
#open ($bng_agp, "<", $output_agp) or die "can't open $output_agp: $!";

close $bng_contig_agp;
#my $bng_contig_agp_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/pre_Tcas5.2_2_bng_contig.agp";
open ($bng_contig_agp, "<", $bng_contig_agp_file) or die "can't open $bng_contig_agp_file:$!";

while (<$bng_contig_agp>)
{
    unless(/^#/)
    {
        chomp;
        my @columns = (split(/\t/));
        my $component_type = (split(/\t/))[4];
        my $contig_id = (split(/\t/))[5];
        my $scaffold_id = (split(/\t/))[0];
        $columns[0]=$scaffold_key{$scaffold_id}; #swap scaffold id
        if ($component_type eq "W")
        {
            $columns[5]= $contig_key{$contig_id}; #swap contig id
        }
        my $corrected_agp = join("\t",@columns);
        print $scaffolds_from_contigs_out "$corrected_agp\n";
    }
    else
    {
        print $scaffolds_from_contigs_out "$_";
    }
}

###############################################################################
##########       Set final names for scaffolds to ChLG AGP          #########
###############################################################################
# open final contigs to scaffolds AGP
my $chlg_from_scaffolds_out_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/Tcas5.2_chlg_from_scaffolds.agp";
open (my $chlg_from_scaffolds_out, ">", $chlg_from_scaffolds_out_file) or die "can't open $chlg_from_scaffolds_out_file: $!";
# open intermediate contigs to scaffolds AGP
#my $chlg_agp_file = "/homes/bioinfo/bionano/Trib_cast_0002_september_2014/ncbi/pre_Tcas5.2_2_chlg.agp";
open (my $chlg_agp_in, "<", $chlg_agp_file) or die "can't open $chlg_agp_file:$!";

while (<$chlg_agp_in>)
{
    unless(/^#/)
    {
        chomp;
        my @columns = (split(/\t/));
        my $component_type = (split(/\t/))[4];
        my $scaffold_id = (split(/\t/))[5];
        if ($component_type eq "W")
        {
            $columns[5]= $scaffold_key{$scaffold_id}; #swap scaffold id
        }
        my $corrected_agp = join("\t",@columns);
        print $chlg_from_scaffolds_out "$corrected_agp\n";
    }
    else
    {
        print $chlg_from_scaffolds_out "$_";
    }
}

#verify that Y bin is correct 206
#awk '/Scaffold/{i++}i==206' /homes/bioinfo/bionano/Trib_cast_0002_september_2014/tcas4_scaffolds.fa | head

print "Done\n";

#export PERL5LIB=/usr/lib/perl5:/usr/lib/perl5/site_perl:/usr/lib/perl5/vendor_perl:/usr/lib64/perl5:/usr/lib64/perl5/site_perl:/usr/lib64/perl5/vendor_perl:/homes/bioinfo/perl5/lib/perl5:/homes/bioinfo/perl5/lib/perl5/x86_64-linux::/homes/bioinfo/bioinfo_software/perl/lib64/perl5:/homes/bioinfo/bioinfo_software/perl/lib64/perl5/site_perl:/homes/bioinfo/bioinfo_software/BioPerl/lib64/perl5/site_perl:/homes/bioinfo/bioinfo_software/BioPerl/lib64/perl5/site_perl/5.8.8/Bio:/homes/bioinfo/bioinfo_software/BioPerl/lib64/perl5/site_perl/5.8.8/x86_64-linux



