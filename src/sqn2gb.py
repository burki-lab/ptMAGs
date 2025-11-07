# here we store functions to convert sqn files 
# to genbank fasta files.

import csv
import os
import subprocess

class PlastomeSQNIterator:
    '''This is an iterator for Plastome SQN Files.
    '''
    def __init__(self, directory, output_dir, slurmlog_csv_file, silent=False):
        '''Initialize iterator by providing a directory.
        '''
        self.directory = directory
        self.output_dir = output_dir
        self.silent = silent  
        self.sqn_files = [file for file in os.listdir(directory) if file.endswith(".sqn")]
        self.slurmlog_csv_file = slurmlog_csv_file
   
    # neccessary methods    
    def run_asn2gb(self, restart_fails=True, force=False):
        '''Run asn2gb on all genomes.
        '''
        from utils_slurm import asn2gbSubmitter

        asn2gb = asn2gbSubmitter(self.slurmlog_csv_file, silent=self.silent)        
        
        for sqn_file in self.sqn_files:
            accession = sqn_file.split(os.path.extsep, 1)[0]
            output_file_name = f"{accession}.gb"
            output_gb_file = os.path.join(self.output_dir, output_file_name)

            # Check if the output file already exists
            if not force and os.path.exists(output_gb_file):
                if not self.silent:
                    print(f"Output file for {accession} already exists.")
                continue

            # Submit the job using asn2gbSubmitter
            script_kwargs = asn2gb.make_kwargs(accession, self.directory, self.output_dir)
            asn2gb.run(accession, self.directory, self.output_dir, force=force, restart_fails=restart_fails, **script_kwargs)

# end PlastomeSQNIterator
