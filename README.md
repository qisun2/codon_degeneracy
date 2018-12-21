# Codon Degeneracy

This script 0-fold 2-fold and 4-fold codon degeneracy sites from a genome. It requires a genome fasta file and a gff file.  

## Getting Started

### Prerequisites
Python3

### Installing

Download the python script

### Command line

python3 degeneracy myannotation.gff3 mygenome.fastas

The output are three files: 
deglist_0fold_sorted.bed
deglist_2fold_sorted.bed
deglist_4fold_sorted.bed

The position in the bed file conforms to bed file standard (0 based in the 2nd column, 1 based in the 3rd column), and it is the position of the mutable nucleotide not the codon start.
The Name column of the bed file includes the following field:
type, transcript ID, codon position in the CDS (0 based), codon position in the chromosome (0 based), mutation postion offset from codon start, codon

## Authors

* **Qi Sun** - *Initial work*

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments