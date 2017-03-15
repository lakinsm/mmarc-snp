#!/usr/bin/env python3

import sys
import glob
import os.path

reldir = '/'.join(os.path.realpath(__file__).split('/')[:-1])


def fastaParse(infile):
    with open(infile, 'r') as fastaFile:
        # Skip whitespace
        while True:
            line = fastaFile.readline()
            if line is "":
                return  # Empty file or premature end of file?
            if line[0] is ">":
                break
        while True:
            if line[0] is not ">":
                raise ValueError("Records in FASTA should begin with '>'")
            header = line[1:].rstrip()
            allLines = []
            line = fastaFile.readline()
            while True:
                if not line:
                    break
                if line[0] is ">":
                    break
                allLines.append(line.rstrip())
                line = fastaFile.readline()
            yield header, "".join(allLines).replace(" ", "").replace("\r", "")
            if not line:
                return  # Stop Iteration
        assert False, "Should not reach this line"


def load_snpsearch_metadata(infile):
    ret = {}
    with open(infile, 'r') as f:
        data = f.read().split('\n')
        for line in data:
            if line:
                entry = line.split(',')
                # name = class, mechanism, group, search_val, snp_name, wt_aa, mutant_aa, mmarc_codon_start,
                # mmarc_codon_end, gene
                ret.setdefault(entry[0], []).append(entry[1:10] + [entry[12]])
    return ret


def find_ungapped_position(seq, codon_start, codon_stop):
    gene_start, gene_stop = codon_start - 1, codon_stop - 1
    i = 0
    for i, nucleotide in enumerate(seq):
        if i < gene_stop and nucleotide == '-':
            gene_stop -= 1
        if i < gene_start and nucleotide == '-':
            gene_start -= 1
        if i == codon_stop:
            break
    if gene_stop - gene_start < 3:
        return False, False
    else:
        return gene_start, gene_stop


def traverse_MPAs(top_level_dir_path, snpsearch_metadata):
    sys.stdout.write('gene_header,model_name,0-based_start_index,0-based_stop_index,wild-type_amino_acid,' +
                     'mutant_amino_acid,symbolic_gene_name,class_annotation,mechanism_annotation,group_annotation,'+
                     'ungapped_sequence\n')
    for multiple_alignment_fasta in glob.glob(top_level_dir_path + '/*/*'):
        model_name = multiple_alignment_fasta.split('/')[-1].replace('.fasta', '')
        if model_name in snpsearch_metadata:
            if snpsearch_metadata[model_name][0][3] == '1' and snpsearch_metadata[model_name][0][-1]:  # if we need to search this model
                mpa_data = {header: seq for header, seq in fastaParse(multiple_alignment_fasta)}
                for header, seq in mpa_data.items():
                    for snp_entry in snpsearch_metadata[model_name]:
                        start, stop = find_ungapped_position(seq, int(snp_entry[7]), int(snp_entry[8]))
                        if start:
                            # header, modelname, 0-based start idx, 0-based stop idx, wild-type amino acid,
                            # mutant amino acid, gene, class, mechanism, group,ungapped seq
                            sys.stdout.write('{},{},{},{},{},{},{},{},{},{},{}\n'.format(
                                header,
                                model_name,
                                start,
                                stop,
                                snp_entry[5],
                                snp_entry[6],
                                snp_entry[9],
                                snp_entry[0],
                                snp_entry[1],
                                snp_entry[2],
                                seq.replace('-', '')
                            ))


if __name__ == '__main__':
    snpsearch_metadata = load_snpsearch_metadata(reldir + '/mmarc_snpsearch_metadata.csv')
    traverse_MPAs(reldir + '/models', snpsearch_metadata)

