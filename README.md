lskScripts
==========

A placeholder for all my random scripts.  Some scripts might eventually graduate and go to their own projects, so don't be surprised if anything leaves in the future.

HELP
==========
Run any script with -h or --help for these help menus.

USAGE
==========

    convertAlignment.pl
    Converts an alignment to another alignment format.
      /home/lkatz/bin/convertAlignment.pl -i alignmentfile -o outputPrefix [-f outputformat]
      -i the alignment file to input
      -o the output alignment file or directory
        Specify a directory by a trailing slash.
        If a directory, the output format will be determined by -f. Default: fasta.
      -f output format
        possible values are derived from bioperl.
        e.g. fasta, clustalw, phylip, xmfa
      -g input format to force the correct input format, but it's ok to let it guess
      -c to concatenate the alignment into one solid entry
      -l linkerSequence, to be used with -c, to separate entries.
        Default is no linker but a useful linker could be NNNNNNNNNN
      -r remove uninformative sites (Ns, gaps, and non-variable sites)
      -h this helpful menu
       at /home/lkatz/bin/convertAlignment.pl line 21.

    extractSequence.pl
    Usage: perl /home/lkatz/bin/extractSequence.pl -i inputGenomicFile -s start -e end
      -i input file: the file extension dictates the format
      -o outfile
      -s start coordinate. 1-based
      -e end coordinate. 1-based
      -c contig 
        put a negative sign in front of a contig to indicate revcom
        you may have multiple -c args
      -n name of organism (useful for genome browsers)
       at /home/lkatz/bin/extractSequence.pl line 21.

    fastacmd.pl
    /home/lkatz/bin/fastacmd.pl: main::main: 
    Finds a sequence in a fasta file and prints it  
    perl /home/lkatz/bin/fastacmd.pl -s search -d database
      -s search term in the defline.
        search can also be from stdin
      -d fasta file to search
        comma-separated databases if multiples
      -i case-insensitive

    fastqToFastaQual.pl
    Usage: /home/lkatz/bin/fastqToFastaQual.pl -i inputFastqFile [-n numCpus -q outputQualfile -f outputFastaFile] at /home/lkatz/bin/fastqToFastaQual.pl line 25.

    genomeDist.pl
    Finds the p-distance between two assemblies using mummer. With more genomes, it creates a table.
      Usage: /home/lkatz/bin/genomeDist.pl assembly.fasta assembly2.fasta [assembly3.fasta ...]
      -a to note averages. Switching the subject and query can reveal some artifacts in the algorithm.
      -q for minimal stdout
      -m method.  Can be mummer (default) or jaccard
        Mummer: uses mummer to discover SNPs and counts the total number
        Jaccard: (kmer method) counts 18-mers and calculates 1 - (intersection/union)
      -c minimum kmer coverage. Default: 2
      -k kmer length. Default: 18
       at /home/lkatz/bin/genomeDist.pl line 19.

    lasergeneToFna.pl
    Converts lasergene sequence files to a multifasta file
      Usage: lasergeneToFna.pl *.lasergene > file.fasta
       at /home/lkatz/bin/lasergeneToFna.pl line 9.

    lyve_splitgbk.pl
    Splits a genbank file from gendb back out into individual loci
      Usage: lyve_splitgbk.pl -i in -o out
      -i input genbank file
      -o output genbank file
        required: .gbk extension
      -d
        debug mode
      -l linker
        default: NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN
      -h
        this help menu
       at /home/lkatz/bin/lyve_splitgbk.pl line 23.

    phylovizFilePrep.pl
    Usage: /home/lkatz/bin/phylovizFilePrep.pl -i infile.fasta -o out.phyloviz
     where the first sequence is a reference sequence at /home/lkatz/bin/phylovizFilePrep.pl line 20.

    sortFastq.pl
    Allows for better compression of reads by sorting them. Uses smalt under the hood. CGP fast assembly for an assembler.
      Usage: /home/lkatz/bin/sortFastq.pl [-r reference.fasta] < reads.fastq | gzip -c > sorted.fastq.gz
      -r reference.fasta If not supplied, then the reads will be assembled to create one.
      -t tempdir/ Default: one will be created for you
      --numcpus How many cpus to use when mapping. Default: 1
      -p to describe paired-end
       at /home/lkatz/bin/sortFastq.pl line 19.
