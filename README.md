# SVarp 
***Pangenome-based structural variant discovery***

The aim of SVarp is to discover haplotype resolved SVs on top of a pangenome graph reference using long sequencing reads. It outputs local assemblies of SV alleles, termed *svtigs*.

*Please open an issue for your questions or feel free to send me an e-mail (arda.soylev@hhu.de)*

## Quick Start
	git clone https://github.com/asylvz/SVarp.git
	cd svarp
	make
	build/svarp -a xxx.gaf -g xxx.gfa --fasta xxx.fasta.gz --phase read_tags.tsv -i SAMPLE_NAME -o OUTPUT_FOLDER


## Requirements

SVarp is developed and tested using Linux Ubuntu operating system

* You need the FASTA of the reads **bgzipped** (i.e., bgzip NA12878.fasta)
* GAF file needs to be in **unstable coordinate system**. (https://github.com/lh3/gfatools/blob/master/doc/rGFA.md)
* Pangenome ref. needs to be in GFA format. We tested SVarp with **HPRC Minigraph pangenome** (https://zenodo.org/record/6983934)
	- GAF alignments must be generated using the same GFA file that you use as input
* For **phased variants**, tags need to be generated in ".tsv" format using **whatshap haplotag** (https://whatshap.readthedocs.io/)
	- E.g., *whatshap haplotag NA12878.vcf.gz NA12878.bam -o NA12878_phase --reference ref_genome.fasta --output-haplotag-list tags.tsv*


## Dependencies

* **wtdbg2** 
 	- git clone https://github.com/ruanjue/wtdbg2 && cd wtdbg2 && make (copy *wtdbg2*, *wtpoa-cns* and *wtdbg2.pl* to your PATH)
* **zlib**
  	- sudo apt-get install zlib1g-d
* **GraphAligner**
  	- conda install -c bioconda graphaligner (https://github.com/maickrau/GraphAligner)
* **Samtools**
  	- https://github.com/samtools/samtools
* **Minimap2**
  	- git clone https://github.com/lh3/minimap2 && cd minimap2 && make
* **HTSlib and WFA2-lib**
  	- You can either use ***make libs*** or follow the steps below

### HTSlib (https://github.com/samtools/htslib)
	wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2
	mkdir htslib && tar -xvf htslib-1.17.tar.bz2 -C htslib --strip-components=1
	cd htslib && autoconf -i && ./configure --disable-lzma --disable-bz2 --disable-libcurl && make
	(htslib folder needs to reside inside SVarp's main folder)

### WFA2-lib (https://github.com/smarco/WFA2-lib)
	wget https://github.com/smarco/WFA2-lib/archive/refs/tags/v2.3.4.tar.gz --strip-components=1
	mkdir wfa && tar -xzf v2.3.4.tar.gz -C wfa
	cd wfa && make clean all
	(wfa folder needs to reside inside SVarp's main folder)

## All parameters

	Required arguments:
	--gaf (-a)                  : Alignment file in GAF format
	--graph (-g)                : Pangenome file in GFA format
	--fasta (-f)                : Fasta sequence file


	Optional arguments:
	--sample (-i)               : Sample name.
	--out (-o)                  : Output folder.
	--debug                     : Output multiple log files for debugging purpose.
	--skip-untagged             : Output only phased variants (~30% faster).
	--dist_threshold (-d)       : Distance threshold to merge SV breakpoints (default=100)
	--out (-o)                  : Output folder path
	--phase (-p)                : WhatsHap haplotag file in .tsv (https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-haplotag)
	--reads(-r)                 : Bgzipped FASTA file of reads for extensive mode (needed for WFA realignment)
	--sample (-i)               : Sample (Individual) name
	--support (-s)              : Minimum support for a cluster to be assembled (default=5 for diploid samples)
	--threads(-t)               : Number of threads for assembly and realignment (default:32)
	--help                      : Print this help menu

## Citation

*Soylev, A., Ebler, J., Pani, S., Rausch, T., Korbel, J., & Marschall, T. (2024). **SVarp: pangenome-based structural variant discovery**. bioRxiv, 2024-02.*
