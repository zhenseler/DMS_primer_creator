# DMS_primer_creator

PYTHON 2.X COMPATIBLE

Design primers for deep mutational scanning (Fowler 2014)

From an input gene sequence (containing at least 30 bp of flanking DNA) creates primers to be used for library construction for deep mutational scanning (DMS) experiments.  Primers can be used in conjuction with PALS mutagenesis (Kitzman 2015) or one-pot saturation mutagenesus (Wrenbeck 2016) library construction protocols.  Both protocols involve using primers which match input sequence perfectly, aside from one codon change.  

Has two modes: programmed or degenerate.
  - programmed: for each position, create individual primers for each codon mutation
  - degenerate: for each position, create one degenerate 'NNN' primer
  
Features
  - Primer creation
    - Makes sure tm is > 78oC
    - Extends to assure GC-clamp (if <6 bp away)
    - Creates primers for codon deletions as well
    
Arguments
		argv1 = input fasta filepath
		argv2 = mode (options: programmed or degenerate)
		argv3 = codon table filepath (codons must be in column one, no header)
		argv4 = position of first nucleotide to mutagenize
		argv5 = position of last nucleotide to mutagenize
		argv6 = sequence for 5\'-overhang
    argv7 = sequence for 3\'-overhang of the primer (do not reverse complement)
		argv8 = output filepath, ending with '.txt'
