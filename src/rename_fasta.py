# here we store functions to rename sequences in fasta files.

import os
import subprocess

class FastaIterator:
    '''This is an iterator for Fasta files.
    '''
    def __init__(self, directory, output_dir, mapping_file, silent=False):
        '''Initialize iterator by providing a directory.
        '''
        self.directory = directory
        self.output_dir = output_dir
        self.silent = silent  
        self.fasta_files = [file for file in os.listdir(directory) if file.endswith(".fasta")]
        self.mapping_file = mapping_file
   
    # neccessary methods    
    def prepend_taxo(self):
        '''Prepend taxonomy on each fasta header.
        '''

        # Populate dictionary with accession numbers and asscoiated taxonomy
        mapping_dict = {}

        with open(self.mapping_file, 'r') as file:
            for line in file:
                columns = line.split('\t')
                key = columns[0]
                value = columns[2].strip()

                # Populate the dictionary with the key and value
                mapping_dict[key] = value   
        
        for fasta_file in self.fasta_files:
            accession = os.path.splitext(fasta_file)[0]

            # Get the corresponding taxonomy from the dictionary
            if accession in mapping_dict:
                taxonomy = mapping_dict[accession]
            else:
                print(f"Accession number {accession} not found in the dictionary.")
                continue

            # Construct the output file name
            output_fasta_file = f"{accession}.taxo.fasta"

            # Execute the seqkit command to prepend the taxonomy to the header
            cmd = f"cat {os.path.join(self.directory, fasta_file)} | seqkit replace -p ^ -r '{taxonomy};' > {os.path.join(self.output_dir, output_fasta_file)}"
            subprocess.call(cmd, shell=True)

# end FastaIterator
