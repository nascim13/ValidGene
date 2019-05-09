<img src="DSG_doodle.png" width="400">

# ValidGene
Pipeline for the validation of incomplete nucleotide sequences.

## Requirements

* Software:<br>
`$ apt-get install ncbi-blast+`<br>
`$ apt-get install samtools`<br>
`$ apt-get install bowtie2`<br>

* Efetch:<br>
`$ cd ~/bin/bash`<br>
`$ perl -MNet::FTP -e \
    '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1);
    $ftp->login; $ftp->binary;
    $ftp->get("/entrez/entrezdirect/edirect.tar.gz");`<br>
`$ gunzip -c edirect.tar.gz | tar xf -`<br>
`$ rm edirect.tar.gz`<br>
`$ builtin exit`<br>
`$ export PATH=$PATH:$HOME/edirect >& /dev/null || setenv PATH "${PATH}:$HOME/edirect"`<br>
`$ ./edirect/setup.sh`<br>
`$ echo "export PATH=\$PATH:\$HOME/edirect" >> $HOME/.bash_profile`<br>


* Statistics R module for perl:<br>
`$ cpan install Statistics::R`<br>

* ggplot2 library for R:<br>
`$ R`<br>
`$ install.packages("ggplot2")`<br>

## Usage

`$ perl ValidGene.pl -in <fileName> -fq <fileName.fastq>`<br>

* fileName: name of the fasta file containing sequences to be analyzed (e.g fasta file = input.fasta; fileName = input).
* fileName.fastq: read file (.fastq extension).

## Test set

* This package comes with 3 incomplete sequences for testing.<br>
`$ perl ValidGene.pl -in seq1 -fq bmul.fastq`<br>
`$ perl ValidGene.pl -in seq2 -fq bmul.fastq`<br>
`$ perl ValidGene.pl -in seq3 -fq bmul.fastq`<br>

### Output files:<br>
* fileName.csv: Summarizes the results of the validation analysis.<br>
* fileName.pdf: A read coverage plot for the region containing the homolog sequence in the reference genome.<br>
