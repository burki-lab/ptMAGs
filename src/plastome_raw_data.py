# here one finds a class that helps to iterate
# through all plastomes and it follow up analyses

from Bio import SeqIO
from Bio.SeqUtils import GC
import csv
import numpy as np
import os


from plastome_abstract_data import PlastomeAbstractFile
from plastome_abstract_data import PlastomeAbstractIterator
from plastome_name import PlastomeName

class PlastomeRawFile(PlastomeAbstractFile):
    '''This is a container class for raw plastome fasta files with contigs.
    '''
    def __init__(self, name, base_path, suffix=".fa", alias=""):
        '''Initialize abstract instance by name and base_path.
        '''
        super().__init__(name, base_path, suffix=suffix, alias=alias)
        return
    
    # abstract methods
    def summarize(self):
        '''Summarize file specific information.
        '''
        # lambda to make bin code from bin dict
        listize = lambda dct: [f"{i}:{j}" for i, j in dct.items()][0]

        # load data
        size = self.get_size()
        self.get_seq_summary()
        self.meta.get_data()
        
        if isinstance(self.meta.bin, dict):
            bin_name = [
                f"{i}:{listize(j)}" if isinstance(j, dict) else f"{i}:{j}" 
                for i, j in self.meta.bin.items()][0]
        else : bin_name = self.meta.bin
        # assemble the summary of file infos 
        summary_dict = {
            "id": self.id,
            "filename": self.file,
            "filesize": size,  # in bytes
            "region": self.meta.region,
            "bin": bin_name,
            "total_length": self.contigs_info["total_length"],
            "contig_count": self.contigs_info["contig_count"],
            "mean_gc": self.contigs_info["mean_gc"],
            "std_gc": self.contigs_info["std_gc"],
            "median_length": self.contigs_info["median_length"],
            "max_length": self.contigs_info["max_length"],
            "min_length": self.contigs_info["min_length"]
        }
        return summary_dict
    
    # additional methods
    def get_seq_summary(self):
        '''Summarize the sequences within the given file.
        '''
        # initialize dict for tracking sequence infor of contigs
        self.contig_sp_info = {
                "id": {},
                "length": {},
                "gc_content": {}
        }
        self.contigs_info = {}
        
        # collect info from contigs
        for i, seq in enumerate(SeqIO.parse(self.file, format="fasta")): 
            self.contig_sp_info["id"][i] = seq.id
            self.contig_sp_info["length"][i] = len(seq)
            self.contig_sp_info["gc_content"][i] = GC(seq.seq)/100
        
        # condense the sequence info to summary statistics
        self.contigs_info["contig_count"] = len(self.contig_sp_info["id"])
        self.contigs_info["total_length"] = sum(self.contig_sp_info["length"].values())
        
        self.contigs_info["mean_gc"] = np.mean(
            list(self.contig_sp_info["gc_content"].values()))
        self.contigs_info["std_gc"] = np.std(
            list(self.contig_sp_info["gc_content"].values()))

        self.contigs_info["median_length"] = np.median(
            list(self.contig_sp_info["length"].values()))
        self.contigs_info["max_length"] = np.max(
            list(self.contig_sp_info["length"].values()))
        self.contigs_info["min_length"] = np.min(
            list(self.contig_sp_info["length"].values()))
        
        return self.contigs_info
    
    def get_contig_seq(self, contig_name):
        '''Obtain sequence of a requested contig.
        '''
        for seq in SeqIO.parse(self.file, format="fasta"): 
            if contig_name == seq.id:
                return str(seq.seq)
        raise RuntimeError(f"Contig '{contig_name}' was not found in plastome {self.id}")
        
    def run_mfannot(self, mfannot_submitter, output_dir, force=False, **kwargs):
        '''Submit MFannot job for plastome.
        '''
        self.run_slurm_job(mfannot_submitter, output_dir, force=force, **kwargs)
        return
# end PlastomeRawFile


class PlastomeRawIterator(PlastomeAbstractIterator):
    '''This is an iterator for PlastomeRawFiles.
    '''
    def __init__(self, directory, suffix=".fa"):
        '''Initialize iterator by providing a directory.
        '''
        super().__init__(directory, file_container=PlastomeRawFile, suffix=suffix)
        return
   
    # neccessary methods
    def load_files(self):
        '''Select relevant files within directory that should be iterated through.
        '''
        super().load_files()
        return
    
    def store_information(self, csv_file_path, force=False):
        '''Store the class specific information for all plastomes.
        '''
        super().store_information(csv_file_path, force=force)
        return
    
    def run_mfannot(self, output_dir, slurmlog_csv_file, prot_coll_path, restart_fails=True, force=False):
        '''Run MFannot on all genomes.
        '''
        from utils_slurm import MFannotSubmitter
        mfs = MFannotSubmitter(slurmlog_csv_file, silent=self.silent)
        
        for i, plastome in enumerate(self.__iter__()):
            plastome.run_mfannot(mfs, output_dir, prot_coll_path=prot_coll_path, restart_fails=restart_fails, force=force)
# end PlastomeRawIterator