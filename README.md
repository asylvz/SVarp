# SVarp
pangenome-based structural variant discovery


## Quick Start
	git clone https://git.hhu.de/liy34nag/svarp.git
	cd svarp
	make
	build/svarp -a xxx.gaf -g xxx.gfa --fasta xxx.fasta.gz --phase read_tags.tsv -i SAMPLE_NAME -o OUTPUT_FOLDER


## Requirements

* Only works in Linux.
* You need the fasta files of the reads bgzipped (i.e., bgzip NA12878.fasta)
* GAF file needs to be in unstable coordinate system. (more info here: https://github.com/lh3/gfatools/blob/master/doc/rGFA.md)
* For phased variants, you need to provide tags in ".tsv" format. E.g., whatshap haplotag (https://whatshap.readthedocs.io/). A sample command would be:
	````
	whatshap haplotag NA12878_longread.vcf.gz NA12878.bam -o NA12878_longread_phase --reference ref_genome.fasta --output-haplotag-list tags.tsv

## Dependencies

### wtdbg2 (https://github.com/ruanjue/wtdbg2)
	git clone https://github.com/ruanjue/wtdbg2
	cd wtdbg2 && make	
	(Copy wtdbg2.pl to your PATH)

### htslib (https://github.com/samtools/htslib)
	wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2
	mkdir htslib && tar -xvf htslib-1.17.tar.bz2 -C htslib --strip-components=1
	cd htslib && autoconf -i && ./configure --disable-lzma --disable-bz2 --disable-libcurl && make
	(htslib folder needs to reside inside SVarp's main folder)

### WFA2-lib (https://github.com/smarco/WFA2-lib)
	wget https://github.com/smarco/WFA2-lib/archive/refs/tags/v2.3.4.tar.gz --strip-components=1
	mkdir wfa && tar -xzf v2.3.4.tar.gz -C wfa
	cd wfa && make clean all
	(wfa folder needs to reside inside SVarp's main folder)

### zlib (http://www.zlib.net)
	Ubuntu: sudo apt-get install zlib1g-d

### GraphAligner (https://github.com/maickrau/GraphAligner)

### samtools (https://github.com/samtools/samtools)

    
	You can also use**make libs**to download htslib and WFA2-lib



