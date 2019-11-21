from collections import defaultdict
from Bio import SeqIO
import shutil
import os
import sys
import argparse

class Orthologs:
    def __init__(self, ortho_file, seq_map_file):       
        """Store information on orthologs from another of strains """

        self.ortholog_map = defaultdict(lambda: defaultdict(list))
        with open(ortho_file) as infile:
            for line in infile:
                og, gene_id, spp = line.split()
                self.ortholog_map[og][spp].append(gene_id)

        self.sequence_map = {}
        with open(seq_map_file) as infile:
            for line in infile:
                spp, seq_file = line.split()
                self.sequence_map[spp] = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
        
        self.spp = self.sequence_map.keys()
        self.orth = self.ortholog_map.keys()
   
    def retrieve_seqs(self, ortho_id, spp_name = False):
        """Get sequences corresponding to a given ortho group.
        
        Arguments:
            ortho_id: the ID of a an ortho group present in this object.
            spp_name: boolean, if True, replace the ID (and name and description) of 
                      the sequence with the species name.
        Returns:
            A list of Biopython SeqRecords representing orthologs
        """
        gene_id_map = self.ortholog_map[ortho_id]
        seqs = []
        for sp in self.spp:
            for gid in gene_id_map[sp]:
                to_add = self.sequence_map[sp][gid]
                if spp_name:
                    to_add.id = sp
                    to_add.description = ""
                    to_add.name = ""
                seqs.append(to_add)
        return(seqs)

    def all_orthos(self, min_spp, include_paralogs = False):
        """Iterate over all ortholog groups
        
        Arguments:
            min_spp: minimmum number of species required for an
                     ortho group to be included.
            include_paralogs: boolean, if True, allow ortho groups 
                     that include paralogs (i.e. more than one gene
                     per species). If False skip such ortho groups
        Returns:
            This is a generator function, for use in for loops and simialr.
            The generator yields tuple with two elements
                1: the name of the ortho group being represented
                2: a list of SeqRecords representing sequences in this ortho group
        """
        for ortho_id in self.orth:
                gene_id_map = self.ortholog_map[ortho_id]
                if len(gene_id_map.keys()) < min_spp:
                    continue
                if include_paralogs:
                    for seq_list in gene_id_map.values():
                        if len(seq_list) > 1:
                            continue
                yield( (ortho_id,  self.retrieve_seqs(ortho_id, True) ))



"""
denote files for retrieve_seqs as variable
must dictate the name of the desired ortholog for retrieve.seqs and x=
"""
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = 'ortho_subset, sequence files, name of ortholog, minimum species')
	parser.add_argument('ortho_subset', help='tsv with ortholog information')
	parser.add_argument('sequence_path', help='tsv with sequence paths')
	parser.add_argument('min_species', help='minimum number of species', type=int)
	args = parser.parse_args()
	O = Orthologs(args.ortho_subset, args.sequence_path)
	O.retrieve_seqs(args.og_name)
	for ortho_name, seqs in O.all_orthos(args.min_species):
                out_name = "unaligned/{}.fna".foramt(ortho_name)
		SeqIO.write(seqs, out_name, "fasta")	


