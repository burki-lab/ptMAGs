# this file includes (a) class(es) that allow to iterate through multiple gene files
# and perform some analyses on it.

from abc import ABC
from Bio import SeqIO
import os
import shutil

from __init__ import Scripts

class AbstractPathCarrier(ABC):
    '''This class simply provides the opportunity to check a file/dir.
    '''
    def check_path(self, path, path_type=""):
        '''Check if the path given exists.
        '''
        if not os.path.exists(path):
            raise FileNotFoundError(f"Provided {path_type} {path} does not exist.")
        return path
        
    def check_file(self, filepath):
        '''Check if the directory given exists.
        '''
        self.file = self.check_path(filepath, path_type="file")
        return
        
    def check_dir(self, directory):
        '''Check if the directory given exists.
        '''
        self.dir = self.check_path(directory, path_type="directory")
        return
    
    def run_slurm_job(self, slurm_submitter, output_dir, force=False, input_dir=None, **kwargs):
        '''Submit slurm job using SlurmSubmitter subclasses.
        '''
        from utils_slurm import SlurmSubmitter
        self_id = os.path.splitext(os.path.split(self.file)[-1])[0]
        if input_dir is None:
            self_dir = os.path.dirname(os.path.abspath(self.file))
        else:
            self_dir = input_dir
        
        slurm_submitter.silent = self.silent
        slurm_submitter.run(self_id, self_dir, output_dir, force=force, **kwargs)
        return
# end AbstractPathCarrier
    
    
class GeneIterator(AbstractPathCarrier):
    '''This class allows to iterate through a set of fasta files for genes.
    '''
    def __init__(self, gene_dir, gene_list=[], file_list=None, suffix="fasta", alias="",
                paralogue_seperator=None, silent=False):
        '''Initialize from directory and list of genes, or list of files.
        '''
        # label the instance with alias
        self.alias = alias
        
        # if there are paralogues of some genes, those can be defined in the sequence name
        # by a seperator that is often "@"
        self.par_sep = paralogue_seperator

        # load the files; either we load a selection of files, or we load genes from a directory
        self.check_dir(gene_dir)
        self.load_files(gene_list, file_list, suffix)
        
        self._unlocked = {
            "mafft": False,
            "trimal": False,
            "raxml": False,
            "iqtree": False
        }
        
        self.silent = silent
        return

    
    def load_files(self, gene_list, file_list, suffix):
        '''Select files that are relevant the instance.
        '''
        # given no file list, we load the genes from the directory
        if file_list is None:
            files = [fl for fl in map(lambda x : ".".join([x, suffix]), gene_list) if fl in os.listdir(self.dir)]
        # otherwise, we check the files for existence and use them
        else:
            files = [fl for fl in file_list if os.path.exists(fl)]
        
        # we wrap files into GeneWrappers for easy analysis
        self.genes = [
            GeneWrapper(fl, paralogue_seperator=self.par_sep) 
            for fl in map(lambda x : os.path.join(self.dir, x), files)]
        self.suffix = suffix
        return
    
    def unlock_pipeline(self):
        '''Unlock the run of the pipeline.
        '''
        self._unlocked = {key : True for key in self._unlocked.keys()}
        return
    
    def unlock_genes(self):
        '''Unlock the SLURM runs for all genes.
        '''
        list(map(lambda gene: gene.unlock(), self.genes))
        return
        
    def mafft(self, mafft_dir):
        '''Run mafft on all files.
        '''
        print("--- RUN MAFFT ---")
        
        if self._unlocked["mafft"] : self.unlock_genes()
        list(map(lambda gene: gene.mafft_align(mafft_dir), self.genes))
        self._unlocked["mafft"] = False
        return
    
    def set_mafft_dir(self, mafft_dir, suffix=None):
        '''Manually assign names of the mafft files.
        '''
        if suffix is None : suffix = gene.suffix.strip(".")
        list(map(
            lambda gene: gene.set_mafft_file(os.path.join(mafft_dir, f"{gene.gene}.{suffix}")),
            self.genes))
        self._unlocked["mafft"] = False
        return
    
    def trimal(self, trimal_dir, gt_dropoff=0.8, scr_key="TRIMAL_SINGLE_GENE_GAPFRACTION"):
        '''Run trimal on all files.
        '''
        print("--- RUN TRIMAL ---")
        
        if self._unlocked["trimal"] : self.unlock_genes()
        list(map(lambda gene: gene.trimal(trimal_dir, gt_dropoff=gt_dropoff, scr_key=scr_key), self.genes))
        self._unlocked["trimal"] = False
        return

    def set_trimal_dir(self, trimal_dir, suffix=None):
        '''Manually assign names of the trimal files.
        '''
        if suffix is None : suffix = gene.suffix.strip(".")

        list(map(
            lambda gene: gene.set_trimal_file(os.path.join(trimal_dir, f"{gene.gene}.{suffix}")),
            self.genes))
        self._unlocked["trimal"] = False
        return

    def raxml(self, raxml_dir, model="LG+G+F", force=False):
        '''Run raxml on all files.
        '''
        print("--- RUN RAXML-NG ---")
        
        if self._unlocked["raxml"] : self.unlock_genes()
        list(map(lambda gene: gene.tree_raxml(raxml_dir, model=model, force=force), self.genes))
        self._unlocked["raxml"] = False
        return
    
    def iqtree(self, iqtree_dir, force=False):
        '''Run iqtree on all files.
        '''
        print("--- RUN IQ-TREE ---")
        
        if self._unlocked["iqtree"] : self.unlock_genes()
        list(map(lambda gene: gene.tree_iqtree(iqtree_dir, force=force), self.genes))
        self._unlocked["iqtree"] = False
        return
    
    def run_mafft(self, output_dir, slurmlog_csv_file, amino_flag="--amino", force=False, **kwargs):
        '''Run single gene mafft on all genes.
        '''
        from utils_slurm import MAFFTSubmitter
        ms = MAFFTSubmitter(slurmlog_csv_file, silent=self.silent)
        
        list(map(lambda gene:
                 gene.run_mafft(ms, output_dir, amino_flag=amino_flag, force=force, **kwargs),
                 self.genes))
        return
    
    def run_trimal(self, output_dir, slurmlog_csv_file, alignment_dir=None, force=False, **kwargs):
        '''Run single gene trimal on all genes.
        '''
        from utils_slurm import TrimAlSubmitter
        ts = TrimAlSubmitter(slurmlog_csv_file, silent=self.silent)
        
        list(map(lambda gene : gene.run_trimal(ts, output_dir, force=force, alignment_dir=alignment_dir, **kwargs),
                 self.genes))
        return
    
        
    def run_bmge(self, output_dir, slurmlog_csv_file, alignment_dir=None, force=False, **kwargs):
        '''Run single gene BMGE trimming on all genes.
        '''
        from utils_slurm import BMGESubmitter
        bs = BMGESubmitter(slurmlog_csv_file, silent=self.silent)
        
        list(map(lambda gene : gene.run_bmge(bs, output_dir, force=force, alignment_dir=alignment_dir, **kwargs),
                 self.genes))
        return
             
    def run_siqtree(self, output_dir, slurmlog_csv_file, alignment_dir=None, force=False, **kwargs):
        '''Run single gene IQ tree on all genes.
        '''
        from utils_slurm import SinglePartitionIQTreeSubmitter
        siqs = SinglePartitionIQTreeSubmitter(slurmlog_csv_file, silent=self.silent)
        
        list(map(lambda gene : gene.run_siqtree(siqs, output_dir, force=force, alignment_dir=alignment_dir, **kwargs),
                 self.genes))
        return
    
    def run_divvier(self, output_dir, slurmlog_csv_file, force=False, **kwargs):
        '''Run Divvier trimming on all genes.
        '''
        from utils_slurm import DivvierSubmitter
        ds = DivvierSubmitter(slurmlog_csv_file, silent=self.silent)
        
        list(map(lambda gene:
                 gene.run_divvier(ds, output_dir, force=force, **kwargs),
                 self.genes))
        return
# end GeneIterator
    
class GeneWrapper(AbstractPathCarrier):
    '''This class wraps around a gene, organizes information about it and allows processing.
    '''
    def __init__(self, filepath, paralogue_seperator=None, silent=False):
        '''Initialize from fasta file name, that contains gene name.
        '''
        # store fasta file
        self.check_file(filepath)
        
        # extract gene namme and suffix from filename, as well as sequence type
        self.get_gene_info()
        self.infer_sequence_type()
        self.set_paralogue_splitter(paralogue_seperator=paralogue_seperator)
        
        # make a safety lock for slurm jobs
        self._lock = False
        
        # analyses performed ?
        self.maffted = False
        self.trimmed = False
        self.raxxed = False
        self.intelligent = False
        
        self.silent = silent
        return
        
    def unlock(self):
        '''Unlock for SLURM runs.
        '''
        self._lock = True
        return
        
    def get_gene_info(self):
        '''Load gene name and suffix.
        '''
        self.gene, self.suffix = os.path.splitext(
            os.path.split(self.file)[1])
        return
        
    def infer_sequence_type(self):
        '''Obtain the sequence type of the given file.
        '''
        from Bio.Seq import CodonTable
        
        first_sequence = SeqIO.parse(self.file, format="fasta").__next__()
        # by trying to translate we can infer, if sequence is nucleotide or not.
        try: first_sequence.translate()
        except CodonTable.TranslationError : self.seqtype = "amino_acid"
        else : self.seqtype = "nucleotide"
        self.length = len(first_sequence)
        return
    
    def get_sample_dict(self, force=False):
        '''Obtain map of names of all samples to paralogues of the sequence file.
        '''
        if hasattr(self, "sample_dict") and not force : return self.sample_dict
        # fast iteration though all seqs
        all_seqs = lambda : [sq.name for sq in SeqIO.parse(self.file, format="fasta")]
        # building keys and secondly adding entries
        self.sample_dict = {
            smp : [] for smp in set([self.splitter(sq) for sq in all_seqs()])}
        list(map(lambda x : self.sample_dict[self.splitter(x)].append(x), all_seqs()))
        
        return self.sample_dict
 
    def get_sample_list(self, force=False):
        '''Obtain the names of all samples of the sequence file.
        '''
        sample_list = sorted(self.get_sample_dict(force=force).keys())
        return sample_list
    
    def set_paralogue_splitter(self, paralogue_seperator=None):
        '''Define function that splits sequence names into sample names (paralogue free nomenclature).
        '''
        self.par_sep = paralogue_seperator
        # for known paralogue seperators those could be split of
        self.splitter = lambda x : x if self.par_sep is None else x.split(self.par_sep)[0]
        return
    
    
    # iqtree single gene
    def run_siqtree(self, siqtree_submitter, output_base, force=False, alignment_dir=None, **kwargs):
        '''Submit single gene IQ-Tree job for gene.
        '''
        assert alignment_dir is not None
        self.run_slurm_job(siqtree_submitter, output_base, force=force, input_dir=alignment_dir, **kwargs)
        return
    
    # mafft
    def run_mafft(self, mafft_submitter, output_dir, force=False, **kwargs):
        '''Submit MAFFT job for gene.
        '''
        self.run_slurm_job(mafft_submitter, output_dir, force=force, **kwargs)
        return


    def set_mafft_file(self, mafft_file):
        '''Set the file that contains mafft output
        '''
        self.mafft_file = mafft_file
        self.maffted = os.path.exists(mafft_file)
        return
       
    # trimal
    def run_trimal(self, trimal_submitter, output_dir, force=False, alignment_dir=None, **kwargs):
        '''Submit trimAl job for gene.
        '''
        assert alignment_dir is not None
        self.run_slurm_job(trimal_submitter, output_dir, force=force, input_dir=alignment_dir, **kwargs)
        return

    def set_trimal_file(self, trimal_file):
        '''Set the file that contains trimal output
        '''
        if self.maffted:
            self.trimal_file = trimal_file
            self.trimmed = os.path.exists(trimal_file)
            return
        raise RuntimeError((
            f"No maffted file provided. {self.file}"
            " Try self.mafft_align or manually set self.set_mafft_file"))
    
    # bmge  # trimming tool
    def run_bmge(self, bmge_submitter, output_dir, force=False, alignment_dir=None, **kwargs):
        '''Submit BMGE job for gene.
        '''
        assert alignment_dir is not None
        self.run_slurm_job(bmge_submitter, output_dir, force=force, input_dir=alignment_dir, **kwargs)
        return

    def set_bmge_file(self, trimmed_file):
        '''Set the file that contains trimmed output
        '''
        if self.maffted:
            self.trimmed_file = trimmed_file
            self.trimmed = os.path.exists(trimmed_file)
            return
        raise RuntimeError((
            f"No maffted file provided. {self.file}"
            " Try self.mafft_align or manually set self.set_mafft_file"))

    # raxml
    def tree_raxml(self, raxml_dir, model="LG+G+F", force=False):
        '''Build tree for an alignment using RAxML-NG.
        '''
        raxml_script = Scripts.DICT.value["RAXML_SINGLE_GENE"]
        
        # assign file
        if not os.path.exists(raxml_dir) : os.makedirs(raxml_dir)
        raxml_base = os.path.join(
            raxml_dir,
            "".join([self.gene, self.suffix]))
        
        if os.path.exists(raxml_base) and force : shutil.rmtree(raxml_base)
        self.set_raxml_base(raxml_base)
        
        # define script arguments
        args = [
            self.trimal_file,  # $1
            self.raxml_base,  # $2
            model  # $3
        ]
        
        self.submit_script(raxml_script, args=args)
        self.raxxed = True
        return
    
    def set_raxml_base(self, raxml_base):
        '''Set the basename of files that contains raxml output
        '''
        if self.trimmed:
            self.raxml_base = raxml_base
            self.raxxed = os.path.exists(raxml_base)
            return
        raise RuntimeError(
            ("No trimmed file provided."
             " Try self.trimal or manually set self.set_trimal_file"))

    # iqtree
    def tree_iqtree(self, iqtree_dir, slurmlog_csv_file=None, partitionfile_path=None, force=False, modelfinder=False, **kwargs):
        '''Build tree for an alignment using IQ-Tree.
        '''
        from utils_slurm import MultiPartitionIQTreeSubmitter
        
        assert slurmlog_csv_file is not None
        mpiqtree_submitter = MultiPartitionIQTreeSubmitter(slurmlog_csv_file, modelfinder=modelfinder, silent=self.silent)
        
        self.run_slurm_job(mpiqtree_submitter, iqtree_dir, partitionfile_path=partitionfile_path, force=force, **kwargs)
        return
    
    def set_iqtree_base(self, iq_base):
        '''Set the basename of files that contains raxml output
        '''
        if self.trimmed:
            self.iqtree_base = iq_base
            self.intelligent = os.path.exists(iq_base)
            return
        raise RuntimeError(
            ("No trimmed file provided."
             " Try self.trimal or manually set self.set_trimal_file"))

    # divvier  # trimming tool
    def run_divvier(self, divvier_submitter, output_dir, force=False, alignment_dir=None, **kwargs):
        '''Submit Divvier job for gene.
        '''
        assert alignment_dir is not None
        self.run_slurm_job(divvier_submitter, output_dir, force=force, input_dir=alignment_dir, **kwargs)
        return

    def set_divvier_file(self, trimmed_file):
        '''Set the file that contains trimmed output
        '''
        return self.set_trimal_file(trimmed_file)
# end GeneWrapper