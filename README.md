<img src="DSG_doodle.png" width="400">

# ValidGene
Pipeline for the validation of incomplete nucleotide sequences

## Requirements

* Software:
`$ apt-get install ncbi-blast+`
`$ apt-get install samtools`
`$ apt-get install bowtie2`
* Efetch:
`$ cd ~/bin/bash
$ perl -MNet::FTP -e \
    '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1);
    $ftp->login; $ftp->binary;
    $ftp->get("/entrez/entrezdirect/edirect.tar.gz");'
$ gunzip -c edirect.tar.gz | tar xf -
$ rm edirect.tar.gz
$ builtin exit
$ export PATH=$PATH:$HOME/edirect >& /dev/null || setenv PATH "${PATH}:$HOME/edirect"
$ ./edirect/setup.sh
$ echo "export PATH=\$PATH:\$HOME/edirect" >> $HOME/.bash_profile`


* Statistics R module for perl:
`$ cpan install Statistics::R`

* ggplot2 library for R:
`$ R`
`$ install.packages("ggplot2")

##Usage

`$ perl validGene.pl -in <fileName> -fq <fileName.fastq>`

* fileName: name of the fasta file containing sequences to be analyzed (e.g fasta file = input.fasta; fileName = input)
* fileName.fastq: read file (.fastq extension)

# Test

* This package comes with 3 incomplete sequences for testing.

`$ perl ValidGene.pl -in a.partial -fq bmul.fastq`
`$ perl ValidGene.pl -in bc.partial -fq bmul.fastq`
`$ perl ValidGene.pl -in d.partial -fq bmul.fastq`

### Output files:
* fileName.csv: Summarizes the results of the validation analysis.
* fileName.pdf: A read coverage plot for the region containing the homolog sequence in the reference genome.
