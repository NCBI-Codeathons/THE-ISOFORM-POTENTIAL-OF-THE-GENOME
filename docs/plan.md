# Project 1 
## Find Alternative Starts
### Gervaise Henry

1. Search the human protein sequences for Ms
2. Can be done by taking the sequences files and spliting on Ms -- Create a files of the location of each M in each protein sequence
 - filter on 50aa
 - filter no Ms

# Project 2
## Alternative splicing
### Rong Lu
1. Take GTF file create a bed files with the chrom, start, stop, transcriptid_exon#
2. Use Bedfile and human fasta file to create a fasta file of every human transcript
  - https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html
3. Merge all exons per transcripid -- leaving out one exon at a time.  Ie missing 1, 2, 3, etc
4. Translate DNA to protein using R package (biostrings)

### Zahin Islam 
1. Pass through all species' fasta sequences into Blast for nucleotide alignment 
* All data files copied into new folder "fasta" 
* All files were unzipped from gunzip to fasta files and BlastX was run 
* output will be tabular file 

# Project 3
## Transcript conservation accross species
### Reed
1. Run BLASTX: various organism protein sequence against human genome
2. Parse BLASTX for alignments
3. Compare alignments to human genome against CDS sequences in GTF -- extract alignments not overlapping CDS
