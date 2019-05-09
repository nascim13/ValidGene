use strict;
use warnings;
use Statistics::R;


# This script is used to validate incomplete DNA sequences. It is used on the output files of the <extract_incomplete_seq.pl> script. 
# Usage: validate_diagnostics_v8.pl -in <input.fasta> -fq <input.fastq> -ref_seq <reference.fasta> -refContig <referenceContig.fasta>

# Create a hash for the genetic code and reverse complement bases.
	
my %g = ('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');
my %reverse = ('C' => 'G', 'G' => 'C', 'A' => 'T', 'T' => 'A', 'c' => 'G', 'g' => 'C', 'a' => 'T', 't' => 'A');

# If not created yet, create directory for the output files.

if (! -e "Results") {

	system ("mkdir Results");
	system ("mkdir Results/PDF");
#	system ("mkdir Plots");
}

if (! -e "tmp") {

	system ("mkdir tmp");
	system ("mkdir tmp/Ref_blast");
	system ("mkdir tmp/Ref_seq");
	system ("mkdir tmp/bowtie2");
	system ("mkdir tmp/bowtie2/indexes");
	system ("mkdir tmp/bowtie2/sam_files");
	system ("mkdir tmp/bowtie2/bam_files");
	system ("mkdir tmp/bowtie2/bcf_files");
	system ("mkdir tmp/bowtie2/vcf_files");
	system ("mkdir tmp/bowtie2/depth_files");
	system ("mkdir tmp/bowtie2/depth_files/Plots");
}

#if (! -e "Results/Diagnostics.txt") {

#	open VALIDATION, ">Results/Diagnostics.txt";
#	open NT_SEQUENCES, ">Results/Validated_sequences.fasta";
#	open LOW_COV, ">Results/Low_coverage_sequences.fasta";
#}
#else {

#	open VALIDATION, ">>Results/Diagnostics.txt";
#       open NT_SEQUENCES, ">>Results/Validated_sequences.fasta";
#        open LOW_COV, ">>Results/Low_coverage_sequences.fasta";
#}

my $size = 0;
my $refSeq;

if (! $ARGV[5]) {

# Get a reference protein sequence
	system ("/home/renan/Software/ncbi-blast-2.7.1+/bin/blastx -query $ARGV[1].fasta -db nr -remote -outfmt 6 -max_target_seqs 1 -evalue 1e-05 -out tmp/Ref_blast/$ARGV[1].prot.blast");
	my $subjectProtein = (`cat tmp/Ref_blast/$ARGV[1].prot.blast | awk '{print \$2}'`);
	chomp($subjectProtein);
# Pull the reference protein sequence from NCBI
	system ("/home/renan/Software/edirect/efetch -db nucleotide -format fasta -id $subjectProtein > tmp/Ref_seq/$ARGV[1].prot.fasta");
# Get potential reference contigs
	system ("/home/renan/Software/ncbi-blast-2.7.1+/bin/blastn -query $ARGV[1].fasta -db nr -remote -outfmt 6 -evalue 1e-05 -out tmp/Ref_blast/$ARGV[1].initial.nt.blast");

# Calculate reference size
	if (-s "tmp/Ref_seq/$ARGV[1].prot.fasta") {

		open FASTA, "<tmp/Ref_seq/$ARGV[1].prot.fasta";
		my @fasta = <FASTA>;
		foreach my $line (@fasta) {
			chomp($line);

			if ($line !~ /^>/) {

				my @split = split //, $line;
				foreach my $aa (@split) {
					chomp($aa);

					$size++;
				}
			}
		}
	}
	$size = $size * 3 + 3;
}
	
my %hash;
my $output = 0;
my $marker = 0;
my $count = 0;
my $refStart;
my $refEnd;
my $refSubject;
my $refQuery;
my $fullSeq = 0;

if (-s "tmp/Ref_blast/$ARGV[1].initial.nt.blast") {
	
	open BLAST, "<tmp/Ref_blast/$ARGV[1].initial.nt.blast";
	my @blast = <BLAST>;
	foreach my $line (@blast) {
		chomp($line);

		last if $marker == 1;

		my @aux = split /\t/, $line;
		my $subjectContig = $aux[1];

# Pull reference contig
		system ("/home/renan/Software/edirect/efetch -db nucleotide -format fasta -id $subjectContig > tmp/Ref_seq/$ARGV[1].nt.fasta");

# BLAST reference protein against refenrece contig to get the correct coordinates for a full length gene
		system ("/home/renan/Software/ncbi-blast-2.7.1+/bin/makeblastdb -in tmp/Ref_seq/$ARGV[1].nt.fasta -dbtype nucl -out tmp/DB/$ARGV[1].nt.fasta");
		system ("/home/renan/Software/ncbi-blast-2.7.1+/bin/tblastn -query tmp/Ref_seq/$ARGV[1].prot.fasta -db tmp/DB/$ARGV[1].nt.fasta -outfmt 6 -max_hsps 1 -evalue 1e-05 -out tmp/Ref_blast/$ARGV[1].nt.$subjectContig.blast");

		open REFBLAST, "<tmp/Ref_blast/$ARGV[1].nt.$subjectContig.blast";
		my @refBlast = <REFBLAST>;
		foreach my $hit (@refBlast) {
			chomp($hit);

			my ($query, $subject, $identity, $length, $mismatch, $gapopen, $queryStart, $queryEnd, $subjectStart, $subjectEnd, $evalue, $score) = split /\t/, $hit;
			$refQuery = $query;
			$refSubject = $subject;
			$refStart = $subjectStart;
			$refEnd = $subjectEnd;
			$refSeq = ();

			if ($subjectStart < $subjectEnd) {

				if ($queryStart == 1 && $queryEnd*3+3 != $size) {  ######################## ParICMP4460.avrPto.2
					
					$refEnd = $subjectEnd + ($size-$subjectEnd) + ($subjectStart-1);
				}
				elsif ($queryStart != 1 && $queryEnd*3+3 == $size) {

					$refStart = $subjectStart - ($queryStart*3-3);
					$refEnd = $subjectEnd + 3;

				}
				elsif ($queryStart != 1 && $queryEnd*3+3 != $size) {

					$refStart = $subjectStart - ($queryStart*3);
					$refEnd = $subjectEnd + ($size-$queryEnd*3);
#print "HERE\n";
				}
				else {
					$refEnd = $subjectEnd + 3;

				}
			}
			else {

				if ($queryStart == 1 && $queryEnd*3+3 != $size) {

					$refEnd = $subjectEnd - ($size-($subjectStart-$subjectEnd+1));

				}
				elsif ($queryStart != 1 && $queryEnd*3+3 == $size) {

					$refStart = $subjectStart + ($queryStart*3-3);
					$refEnd = $subjectEnd - 3;
				}
				elsif ($queryStart != 1 && $queryEnd*3+3 != $size) {

					$refStart = $subjectStart + ($queryStart*3);
					$refEnd = $subjectEnd - ($size-$queryEnd*3);
				}
				else {

					$refEnd = $subjectEnd - 3;

				}
			}
		}
		open REF, "<tmp/Ref_seq/$ARGV[1].nt.fasta";
		my @reference = <REF>;
		foreach my $ref (@reference) {
			chomp($ref);

			if ($ref =~ /^>/) {

				$count = 0;
			}

			if ($ref !~ /^>/) {

				my @aux = split //, $ref;
				foreach my $base (@aux) {
					chomp($base);

					$count++;
					my $ucBase = uc($base);

					if ($refStart < $refEnd && $count >= $refStart && $count <= $refEnd) {

						$hash{$count} = $ucBase;
						$refSeq = $refSeq . $ucBase;
					}
					if ($refStart > $refEnd && $count <= $refStart && $count >= $refEnd) {
							
						$hash{$count} = $reverse{$ucBase};
						$refSeq = $refSeq . $reverse{$ucBase};

						if ($count eq $refStart) {

							$refSeq = reverse($refSeq);
						}
					}
				}
			}
		}

		print "$refQuery\t$refSubject\t$refStart\t$refEnd\t$size\t$refSeq\n";

		if (substr($refSeq, 0, 3) eq "ATG" || substr($refSeq, 0, 3) eq "GTG") {

			if ($g{substr($refSeq, -3)} eq "*" && $size == length($refSeq)) {

				$fullSeq = 1;
				last;
			}
		}
	}
}

else {
#        system ("/home/renan/Software/ncbi-blast-2.7.1+/bin/makeblastdb -in $ARGV[7] -dbtype nucl -out tmp/DB/$ARGV[7]");
#        system ("/home/renan/Software/ncbi-blast-2.7.1+/bin/tblastn -query $ARGV[5] -db tmp/DB/$ARGV[7] -outfmt 6 -evalue 1e-05 -out tmp/Ref_blast/$ARGV[1].nt.blast");
#	system ("cp $ARGV[7] tmp/Ref_seq/$ARGV[1].nt.fasta");

	my $results = (`cat tmp/Ref_blast/$ARGV[1].nt.blast | awk '{print \$1\"\t\"\$2\"\t\"\$9\"\t\"\$10}'`);
        chomp($results);

        my @split = split /\t/, $results;
        $refQuery = $split[0];
        $refSubject = $split[1];
        $refStart = $split[2];
        $refEnd = $split[3];
	$refSeq = ();

	if ($refStart < $refEnd) {

		$refEnd = $refEnd + 3;
	}
	else {

		$refEnd = $refEnd - 3;

	}

	open REF, "<$ARGV[7]";
        my @reference = <REF>;
        foreach my $ref (@reference) {
	        chomp($ref);

                if ($ref =~ /^>/) {

        	        $count = 0;
                }

                if ($ref !~ /^>/) {

                	my @aux = split //, $ref;
                        foreach my $base (@aux) {
                        	chomp($base);

                                $count++;
                                my $ucBase = uc($base);

                                if ($refStart < $refEnd && $count >= $refStart && $count <= $refEnd) {

                                	$hash{$count} = $ucBase;
                                        $refSeq = $refSeq . $ucBase;
                                }
                                if ($refStart > $refEnd && $count <= $refStart && $count >= $refEnd) {

                                	$hash{$count} = $reverse{$ucBase};
                                        $refSeq = $refSeq . $reverse{$ucBase};

                                        if ($count eq $refStart) {

                                        	$refSeq = reverse($refSeq);
                                        }
				}
			}
 		}
	}


	$size = length($refSeq);
}

#print "$refQuery\t$refStart\t$refEnd\t$size\t$refSeq\n";

if ($fullSeq) {
#if (substr($refSeq, 0, 3) eq "ATG" || substr($refSeq, 0, 3) eq "GTG" && $g{substr($refSeq, -3)} eq "*" && $size == length($refSeq)) {

	$output = 1;

	my $codonCount = 0;
	my $codon;
	my $refAASeq;

	my @split = split //, $refSeq;
	foreach my $base (@split) {
		chomp($base);

		$codonCount++;

		if ($codonCount <= 3) {

			$codon = $codon . $base;
		}
		if ($codonCount == 3) {

			$refAASeq = $refAASeq . $g{$codon};

			if ($g{$codon} eq "*") {

			      last;
			}
			$codon = ();
			$codonCount = 0;
		}
	}

	open RESULTS, ">Results/$ARGV[1].csv";

	print RESULTS "Query,$refQuery\nReference,$refSubject\nLocation in the reference,$refStart-$refEnd\nMapping visualization,samtools tview tmp/bowtie2/bam_files/$refQuery\_nt.sorted.bam tmp/bowtie2/indexes/nt.fasta\nStatus,";

	system ("bowtie2-build -q tmp/Ref_seq/$ARGV[1].nt.fasta tmp/bowtie2/indexes/$ARGV[1].nt.fasta");
	system ("bowtie2 -p 8 -k 10 --local -x tmp/bowtie2/indexes/$ARGV[1].nt.fasta $ARGV[3] > tmp/bowtie2/sam_files/$ARGV[1]\_nt.sam");


##			system ("bowtie2 -p 8 -k 10 --local -x tmp/bowtie2/indexes/$ref\_$gene_name\_$contig.fasta -f $sequencing_reads/$strain.fasta > tmp/bowtie2/sam_files/$strain\_$ref\_$gene_name\_$contig.sam");

	system ("samtools view -bS tmp/bowtie2/sam_files/$ARGV[1]\_nt.sam > tmp/bowtie2/bam_files/$ARGV[1]\_nt.bam");
	system ("samtools sort tmp/bowtie2/bam_files/$ARGV[1]\_nt.bam tmp/bowtie2/bam_files/$ARGV[1]\_nt.sorted");
	system ("samtools index tmp/bowtie2/bam_files/$ARGV[1]\_nt.sorted.bam");
	system ("samtools faidx tmp/Ref_seq/$ARGV[1].nt.fasta");
	system ("samtools mpileup -uf tmp/Ref_seq/$ARGV[1].nt.fasta tmp/bowtie2/bam_files/$ARGV[1]\_nt.sorted.bam > tmp/bowtie2/bcf_files/$ARGV[1]\_nt.bcf");
	system ("bcftools view -bvcg tmp/bowtie2/bcf_files/$ARGV[1]\_nt.bcf > tmp/bowtie2/bcf_files/$ARGV[1]\_nt\_2.bcf");
	system ("bcftools view tmp/bowtie2/bcf_files/$ARGV[1]\_nt\_2.bcf > tmp/bowtie2/vcf_files/$ARGV[1]\_nt.vcf");
	system ("samtools depth tmp/bowtie2/bam_files/$ARGV[1]\_nt.sorted.bam > tmp/bowtie2/depth_files/$ARGV[1]\_nt.depth");

	my %varCall;
	my %varDef;
	my %varRef;
	my $newSeq;
	my @vcfPos;
	my $covPercentage;

	open VCF, "<tmp/bowtie2/vcf_files/$ARGV[1]\_nt.vcf";
	my @vcf = <VCF>;
	foreach my $call (@vcf) {
		chomp($call);

		if ($call !~ /^#/) {

			my @split = split /\t/, $call;
			my $pos = $split[1];
			my $ref = $split[3];
			my $var = $split[4];
			my $qual = $split[5];

			if ($call =~ /.*DP4=(\d+),(\d+),(\d+),(\d+);MQ.*/) {

				my $supportRef = $1 + $2;
				my $supportAlt = $3 + $4;
				my $total = $1 + $2 + $3 + $4;
				$covPercentage = $supportAlt/$total;
			}
			if ($covPercentage >= 0.9 && $qual > 30) {

				if ($call !~ /INDEL/ && $var !~ /\,/) {
				
					$varDef{$pos} = "SNP";
					$varCall{$pos} = $var;
				
				}
				if ($call ~~ /INDEL/) {

					$varDef{$pos} = "INDEL";
					$varCall{$pos} = $var;
					$varRef{$pos} = $ref;
				}
			}
		}
	}

	my $coverageFilter = 0;
	my $avgContigCov = 0;
	my @lowCoverage;

	my $checkFile = (`cat tmp/bowtie2/depth_files/$ARGV[1]\_nt.depth | awk 'END { if (NR > 0) print \"1\"; else print \"0\"}'`);
	chomp($checkFile);

	if ($checkFile == 1) {

		$avgContigCov = (`cat tmp/bowtie2/depth_files/$ARGV[1]\_nt.depth | awk '{sum+=\$3} END {print sum/NR}'`);
		chomp($avgContigCov);

#				my $avgCovThird = $avgContigCov/3;
	}

	my %depth;	
# Create a line plot of the coverage depth per position in the effector sequence.

	open DEPTH_OUTPUT, ">tmp/bowtie2/depth_files/Plots/$ARGV[1]\_nt.txt";
	print DEPTH_OUTPUT "Position\tCoverage\tDefinition\n";

	open DEPTH, "<tmp/bowtie2/depth_files/$ARGV[1]\_nt.depth";
	my @depthOutput = <DEPTH>;
	foreach my $out (@depthOutput) {
		chomp($out);

		my @aux = split /\t/, $out;
		$depth{$aux[1]} = $aux[2];
	}

	if ($refStart > $refEnd) {

		for (my $i = $refStart; $i >= $refEnd; $i--) {
	
			if (defined $depth{$i}) {
			
				if ($depth{$i} >= 10) {

					print DEPTH_OUTPUT "$i\t$depth{$i}\tCoverage per position\n";
					print DEPTH_OUTPUT "$i\t$avgContigCov\tAverage contig coverage\n";
				}
				else {
					print DEPTH_OUTPUT "$i\t$depth{$i}\tCoverage per position\n";
					push(@lowCoverage, $i);
#							$coverage_filter = 1;
				}
			}
			
			else {
				print DEPTH_OUTPUT "$i\t0\tCoverage per position\n";
				print DEPTH_OUTPUT "$i\t$avgContigCov\tAverage contig coverage\n";
				$coverageFilter = 1;
				push (@lowCoverage, $i);

			}
		}
	}
	else {

		for (my $i = $refStart; $i <= $refEnd; $i++) {
	
			if (defined $depth{$i}) {

				if ($depth{$i} >= 10) {

					print DEPTH_OUTPUT "$i\t$depth{$i}\tCoverage per position\n";
					print DEPTH_OUTPUT "$i\t$avgContigCov\tAverage contig coverage\n";
				}
				else {
					
					print DEPTH_OUTPUT "$i\t$depth{$i}\tCoverage per position\n";
#							$coverage_filter = 1;
					push (@lowCoverage, $i);
				}
			}
			else {
					
				print DEPTH_OUTPUT "$i\t0\tCoverage per position\n";
				print DEPTH_OUTPUT "$i\t$avgContigCov\tAverage contig coverage\n";
				$coverageFilter = 1;
				push (@lowCoverage, $i);
			}
		}
	}

	my $R = Statistics::R -> new();
	$R -> start_sharedR;
	$R -> send ("library(ggplot2)");
	$R -> send ("data <- read.csv (file=\"tmp/bowtie2/depth_files/Plots/$ARGV[1]\_nt.txt\",head=TRUE,sep=\"\t\")");
	$R -> send ("ggplot (data = data, aes(x = Position, y = Coverage, fill = Definition)) + geom_line(aes(linetype=Definition, color=Definition)) + labs(x=\"Contig position\", y = \"Read coverage\") + scale_linetype_manual(values=c(\"solid\", \"solid\")) + scale_color_manual(values=c(\"#FF0000\", \"#0066CC\")) + theme_minimal() + theme(legend.title=element_blank())");
	$R -> send ("ggsave(\"Results/PDF/$ARGV[1].pdf\", width = 14)");

# Analyzing results
	my @indelCoordinates;

	if ($refStart < $refEnd) {

		my $subjectLength = $refEnd - $refStart;

		for (my $i = $refStart; $i <= $refEnd; $i++) {

			if ($varDef{$i}) {

				push (@vcfPos, $i);

				if ($varDef{$i} eq "SNP") {

					$newSeq = $newSeq . $varCall{$i};
				}

				if ($varDef{$i} eq "INDEL") {

					 my @splitVar = split //, $varCall{$i};
					 my @splitRef = split //, $varRef{$i};
					 my $varLength = scalar(@splitVar);
					 my $refLength = scalar(@splitRef);

					 if ($varCall{$i} !~ /\,/) {

						my $indelPos = $i - $refStart + 2;
						push (@indelCoordinates, $indelPos);
						$newSeq = $newSeq . $varCall{$i};
						$i = $i + $refLength -1;
					 }
					 else {

						 $newSeq = $newSeq . $hash{$i};
					 }	
				}

			}
			else {

				$newSeq = $newSeq . $hash{$i};
			}
		}
		my $codon;
		my $codonCount = 0;
		my $aaSeq;
		my $ntSeq;
		my @split = split //, $newSeq;
		foreach my $base (@split) {
			chomp($base);
			
			$codonCount++;
			$ntSeq = $ntSeq . $base;

			if ($codonCount <= 3) {

				$codon = $codon . $base;
			}
			if ($codonCount == 3) {
				
				$aaSeq = $aaSeq . $g{$codon};

				if ($g{$codon} eq "*") {

					last;
				}
				$codon = ();
				$codonCount = 0;
			}
		}
		my $querySize = length($ntSeq);

		if ($coverageFilter == 0) {

			if ($querySize == $size) {


				if (@vcfPos) {

					print RESULTS "SNP(s) on position(s) @vcfPos\nReference_nt_sequence,$refSeq\nReference_aa_sequence,$refAASeq\nValidated_nt_sequence,$ntSeq\nValidated_aa_sequence,$aaSeq\n";
#					print NT_SEQUENCES ">$ARGV[1].v\n$ntSeq\n";
#					print VALIDATION "$ARGV[1]\tSNP on position @vcfPos\n";
				}
				else {

					print RESULTS "Complete sequence\nReference_nt_sequence,$refSeq\nReference_aa_sequence,$refAASeq\nValidated_nt_sequence,$ntSeq\nValidated_aa_sequence,$aaSeq\n\n";
#					print VALIDATION "$ARGV[1]\tComplete sequence\n";
#					print NT_SEQUENCES ">$ARGV[1].v\n$ntSeq\n";

				}
			}
		
			else {

				if (@indelCoordinates) {

					print RESULTS "INDEL on position @indelCoordinates\nReference_nt_sequence,$refSeq\nReference_aa_sequence,$refAASeq\nValidated_nt_sequence,$ntSeq\nValidated_aa_sequence,$aaSeq\n";
#					print NT_SEQUENCES ">$ARGV[1].v\n$ntSeq\n";
#					print VALIDATION "$ARGV[1]\tINDEL on position @indelCoordinates\n";
				}
				else {

					print RESULTS "Early stop codon\nReference_nt_sequence,$refSeq\nReference_aa_sequence,$refAASeq\nValidated_nt_sequence,$ntSeq\nValidated_aa_sequence,$aaSeq\n";
#					print NT_SEQUENCES ">$ARGV[1].v\n$ntSeq\n";
#					print VALIDATION "$ARGV[1]\tEarly stop codon\n";

				}
			}
		}

		my $lowCovElements = scalar(@lowCoverage);
		my $countElements = 0;
		my $firstElement = 0;	
		my $elementMarker = 0;
		my $previousElement = 0;

		if ($coverageFilter == 1) {

			print RESULTS "Lack of coverage in the following positions: ";
#			print VALIDATION "$ARGV[1]\tLack of coverage in the following positions: ";
#			print LOW_COV ">$ARGV[1]\n";
		
			for (my $i = $refStart; $i <= $refEnd; $i++) {

				if ($i ~~ @lowCoverage) {

#					print LOW_COV "N";
				}
				else {
	
#					print LOW_COV "$hash{$i}";
				}
			}		

			if ($lowCovElements == 1) {

#				print VALIDATION "@lowCoverage";
				print RESULTS "@lowCoverage";
			}
			else {

				foreach my $an (@lowCoverage) {
					chomp($an);

					$countElements++;

					if ($elementMarker == 0) {

						$firstElement = $an;
						$elementMarker = 1;
					}
					my $diff = $an - $previousElement;

					if ($diff != 1 && $firstElement != $an && $lowCovElements != $countElements) {

#						print VALIDATION " $firstElement-$previousElement";
						print RESULTS " $firstElement-$previousElement";
						$firstElement = $an;
					}
					if ($diff == 1 && $firstElement != $an && $lowCovElements == $countElements) {

#						print VALIDATION " $firstElement-$previousElement";
						print RESULTS " $firstElement-$previousElement";
						$firstElement = $an;
					}
					$previousElement = $an;
				}
			}
			$coverageFilter = 0;
			@lowCoverage = ();
#			print VALIDATION "\n";
			print RESULTS "\n";
#			print LOW_COV "\n";
		}
	}
	
	if ($refStart > $refEnd) {

		my $subjectLength = $refStart - $refEnd;

		for (my $i = $refEnd; $i <= $refStart; $i++) {

			if ($varDef{$i}) {

				push (@vcfPos, $i);

				if ($varDef{$i} eq "SNP") {
					
					$newSeq = $newSeq . $reverse{$varCall{$i}};

				}
				if ($varDef{$i} eq "INDEL") {

					my $revVar = reverse($varCall{$i});
					my @splitVar = split //, $revVar;
					my @splitRef = split //, $hash{$i};
					my $varLength = scalar(@splitVar);
					my $refLength = scalar(@splitRef);

					if ($varCall{$i} !~ /\,/) {

						my $indelPos = $refStart - $i + 2;
						push (@indelCoordinates, $indelPos);
#								$newSeq =  $newSeq . $hash{$i};
							
						foreach my $an (@splitVar) {
							chomp($an);

							$newSeq = $newSeq . $reverse{$an};
						}
						$i = $i + $refLength;
					}
					else {

						$newSeq = $newSeq . $hash{$i};
					}		
				}
			
			}
			else {
				$newSeq = $newSeq . $hash{$i};
			}
		}

		my $codon = ();
		my $codonCount = 0;
		my $revSeq = reverse($newSeq);
		my $aaSeq;
		my $ntSeq;

		my @split = split //, $revSeq;
		foreach my $base (@split) {
			chomp($base);

			$codonCount++;	
			$ntSeq = $ntSeq . $base;

			if ($codonCount <= 3) {

				$codon = $codon . $base;
			}
			if ($codonCount == 3) {

				$aaSeq = $aaSeq . $g{$codon};

				if ($g{$codon} eq "*") {

					last;
				}
				$codon = ();
				$codonCount = 0;
			}
		}

		my $querySize = length($ntSeq);

		if ($coverageFilter == 0) {

			if ($querySize == $size) {

				if (@vcfPos) {

					print RESULTS "SNP(s) on position(s) @vcfPos\nReference_nt_sequence,$refSeq\nReference_aa_sequence,$refAASeq\nValidated_nt_sequence,$ntSeq\nValidated_aa_sequence,$aaSeq\n";
#					print NT_SEQUENCES ">$ARGV[1].v\n$ntSeq\n";
#					print VALIDATION "$ARGV[1]\tSNP on position @vcfPos\n";

				}
				else {

					print RESULTS "Complete sequence\nReference_nt_sequence,$refSeq\nReference_aa_sequence,$refAASeq\nValidated_nt_sequence,$ntSeq\nValidated_aa_sequence,$aaSeq\n";
#					print VALIDATION "$ARGV[1]\tComplete sequence\n";
#					print NT_SEQUENCES ">$ARGV[1].v\n$ntSeq\n";
				}
			}
			else {

				if (@indelCoordinates) {

					print RESULTS "INDEL on position @indelCoordinates\nReference_nt_sequence,$refSeq\nReference_aa_sequence,$refAASeq\nValidated_nt_sequence,$ntSeq\nValidated_aa_sequence,$aaSeq\n";
#					print NT_SEQUENCES ">$ARGV[1].v\n$ntSeq\n";
#					print VALIDATION "$ARGV[1]\tINDEL on position @indelCoordinates\n";
				}
				else {

					print RESULTS "Early stop codon\nReference_nt_sequence,$refSeq\nReference_aa_sequence,$refAASeq\nValidated_nt_sequence,$ntSeq\nValidated_aa_sequence,$aaSeq\n";
#					print NT_SEQUENCES ">$ARGV[1].v\n$ntSeq\n";
#					print VALIDATION "$ARGV[1]\tEarly stop codon\n";

				}
			}
		}
			
		my $firstElement = 0;
		my $elementMarker = 0;
		my $countElements = 0;
		my $previousElement = 0;
		my $lowCovElements = scalar(@lowCoverage);

		if ($coverageFilter == 1) {

			print RESULTS "Lack of coverage in the following positions: ";
#			print VALIDATION "$ARGV[1]\tLack of coverage in the following positions: ";
#			print LOW_COV ">$ARGV[1]\n";
		
			for (my $i = $refStart; $i >= $refEnd; $i--) {

				if ($i ~~ @lowCoverage) {

#					print LOW_COV "N";
				}
				else {
	
#					print LOW_COV "$hash{$i}";
				}
			}

			if ($lowCovElements == 1) {

#				print VALIDATION "@lowCoverage";
				print RESULTS "@lowCoverage";
			}
		
			else {

				foreach my $an (@lowCoverage) {
					chomp($an);

					$countElements++;

					if ($elementMarker == 0) {

						$firstElement = $an;
						$elementMarker = 1;
					}
					my $diff = $previousElement - $an;

					if ($diff != 1 && $firstElement != $an && $lowCovElements != $countElements) {

#						print VALIDATION " $firstElement-$previousElement";
						print RESULTS " $firstElement-$previousElement";
						$firstElement = $an;
					}
					if ($diff == 1 && $firstElement != $an && $lowCovElements == $countElements) {

#						print VALIDATION " $firstElement-$previousElement";
						print RESULTS " $firstElement-$previousElement";
						$firstElement = $an;
					}
					$previousElement = $an;
				}
			}	
			$coverageFilter = 0;
			@lowCoverage = ();
#			print VALIDATION "\n";
			print RESULTS "\n";
#			print LOW_COV "\n";
		}
	}
}

else {

	open RESULTS, ">Results/$ARGV[1].txt";
#	print RESULTS "Query:\n$ARGV[1]\n\nDiagnostics: Incomplete reference sequence\n";
#	print VALIDATION "$ARGV[1]\tIncomplete reference sequence\n";


}
#else {

#	open RESULTS, ">Results/$ARGV[1].txt";
#	print RESULTS "Query:\n$ARGV[1]\n\nSubject:\nNo reference hit for this sequence\n\nDiagnostics: There are no reference hits for this sequence\n";
#        print VALIDATION "$ARGV[1]\tThere are no reference hits for this sequence\n";
#}
