# here one finds an abstract class that allows to submit 
# slurm jobs for plastome iterators

from abc import ABC
from abc import abstractmethod
import csv
from datetime import datetime
import numpy as np
import os
import pandas as pd
import regex as re
import shutil
import subprocess

from __init__ import Scripts

class SlurmSubmitter(ABC):
    '''This is an abstract class whose implementations should allow to submit a slurm script.
    '''
    def __init__(self, log_csv_file, silent=False):
        '''Initialize abstract instance by name and base_path.
        '''
        self.silent = silent
        
        # we load the implementation specific details
        self.set_default()
        self.set_header()

        # we define the script from config file
        self.set_script()
        
        # using pandas we load the csv file and keep it updated
        self.load_logfile(log_csv_file)
        return
    
    def set_script(self):
        '''Load the name of the slurm script for given SlurmSubmitter subclass.
        '''
        self.script = Scripts.DICT.value[self.name]
        return

    def set_header(self):
        '''Set header for log file.
        '''
        HEADER_COLS = [
            "job_id",
            "job_status",
            "time_submitted",
            "duration",
            "timelimit",
            "time_last_checked",
            "script_name",
            "script_file",
            "out_file",
            "err_file",
            "input_id",
            "n_cores",
            "software_output",
            "software_output_size"
        ]
        self.header = [*HEADER_COLS, *self.script_args]
        return    
    
    def load_logfile(self, log_csv_file):
        '''Load the logging csv file as pd.DataFrame into cache.
        '''
        self.log_file = log_csv_file
        
        # either a new log_file needs to be produced
        if not os.path.exists(log_csv_file):
            self.log_df = pd.DataFrame(columns=self.header)
            self.save_logfile()
        
        # or a good old one is checked and appended
        else:
            self.log_df = pd.read_csv(log_csv_file)
            self.log_df["job_id"] = self.log_df["job_id"].astype(int).astype(str)
            self.log_df["n_cores"] = self.log_df["n_cores"].astype(int).astype(str)
            self.log_df["software_output_size"] = self.log_df["software_output_size"].astype(int).astype(str)
            if sorted(self.header) != sorted(self.log_df.columns):
                raise IOError("Tracking csv file exists, but has wrong column names.")
        return
            
    def save_logfile(self):
        '''Save the updated logging DataFrame to its source file.
        '''
        self.log_df.to_csv(self.log_file, index=False)
        return

    def check(self, input_id, *script_args, restart_fails=True):
        '''Check if a job for given file id was already run and update in log file.
        '''
        jobs = self.log_df[self.log_df.loc[:,"input_id"] == input_id]
        
        # there can be no jobs ever run for a file id
        if jobs.empty:
            return False
        
        last_job_id = jobs.sort_values(
                by=["time_submitted"], ascending=False)[
                        "job_id"].values[0]
        
        return self.check_job_id(last_job_id, input_id, *script_args, restart_fails=restart_fails)
    
    def check_job_id(self, job_id, input_id, *script_args, restart_fails=True):
        '''Check and update job given its ID in the log file.
        '''
        from dateutil import parser

        REGEX_DICT = {
                "status": r"(?<=JobState\=)[^\ ]*",
                "submit_time": r"(?<=SubmitTime\=)[^\ ]*",
                "duration": r"(?<=RunTime\=)[^\ ]*",
                "timelimit": r"(?<=TimeLimit\=)[^\ ]*",
                "out_file": r"(?<=StdOut\=)[^\ ]*",
                "err_file": r"(?<=StdErr\=)[^\ ]*",
                "n_cores": r"(?<=NumCPUs\=)[^\ ]*"
                }
        SACCT_COLS = ["jobid", "start", "end", "state"]

        # for preexisting outputs, we do nothing.
        if job_id == "PREEXISTING_OUTPUT" : return True
    
        job = self.log_df[self.log_df["job_id"] == job_id]

        # we use scontrol for infos to document new job
        # otherwise we update using sacct, as this goes 
        # back far in time
        if job.empty:
            command_frags = ["scontrol", "show", "job", str(int(job_id))]
        # completed jobs are not updated
        elif job["job_status"].values[0] == "COMPLETED":
            return True
        elif restart_fails and (
            job["job_status"].values[0] in ["FAILED", "CANCELLED", "CANCELLED+", "TIMEOUT"]):
            return False
        else:
            command_frags = ["sacct", "-j", str(int(job_id)), "-n", "-o", ",".join(SACCT_COLS)]
        
        # catch the outputs
        srun = subprocess.Popen(
            command_frags, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sout, serr = srun.communicate()
        
        # for new job
        if job.empty:
            get_job_feature = lambda x : re.search(REGEX_DICT[x], sout.decode('utf-8')).group().strip()
            output_name, output_size = self.check_output(input_id,  # implementation specific output checker
                    **dict(zip(self.script_args, script_args)))
            row = [
                    str(int(job_id)),  # job_id
                    get_job_feature("status"),  # status
                    get_job_feature("submit_time"),  # time of submission
                    get_job_feature("duration"),  # duration of run
                    get_job_feature("timelimit"),  # timelimit
                    datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),  # time last checked
                    self.name,  # script name
                    self.script,  # script path
                    get_job_feature("out_file"),  # file for slurm output
                    get_job_feature("err_file"),  # file for error
                    input_id,  # id of input; neccessary!
                    str(int(get_job_feature("n_cores"))),  # number of cores/CPUs
                    output_name,  # (provisional) name of op file or dir
                    str(int(output_size)),  # size of op, if it exists
                    *script_args
                    ]
            self.log_df.loc[len(self.log_df)+1, :] = row
        
        # for simple status update
        else:
            sout_groups = [sg for sg in sout.decode('utf-8').split("\n")[0].split(" ") if sg!=""]
            souts = dict(zip(SACCT_COLS, sout_groups))
            if souts["end"] == "Unknown" and souts["start"] == "Unknown":
                self.log_df.loc[self.log_df["job_id"] == job_id, "duration"] = "00:00:00"
            else:
                if souts["end"] == "Unknown":
                    dur = datetime.now() - parser.parse(souts["start"])
                else:
                    dur = parser.parse(souts["end"]) - parser.parse(souts["start"])
                self.log_df.loc[self.log_df["job_id"] == job_id, "duration"] = (
                    f"{dur.days*24+dur.seconds//3600:02d}:{dur.seconds%3600//60:02d}:{dur.seconds%60:02d}")

            output_name, output_size = self.check_output(input_id,  # implementation specific output checker
                    **dict(zip(self.script_args, script_args)))
            
            self.log_df.loc[self.log_df["job_id"] == job_id, "job_status"] = souts["state"]
            self.log_df.loc[self.log_df["job_id"] == job_id, "time_last_checked"] = (
                datetime.now().strftime("%Y-%m-%dT%H:%M:%S"))
            self.log_df.loc[self.log_df["job_id"] == job_id, "software_output_size"] = output_size
        # write it down
        self.save_logfile()
        
        return True
    
    def cancel_running_jobs(self, input_id, *args):
        '''Cancel job for input ID if any is running.
        '''
        def cancel(job_id):
            '''Cancel a slurm job.
            '''
            command_frags = ["scancel", str(int(job_id))]
            srun = subprocess.Popen(
                command_frags, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            sout, serr = srun.communicate()
            if not self.silent:
                print(f"Job {job_id} was cancelled to run new job for {input_id}.") 
                
            self.check_job_id(job_id, input_id, *args)
            return
        
        # find running jobs ...
        running_jobs = self.log_df.loc[
            (self.log_df["input_id"] == input_id) &
            (self.log_df["job_status"] == "RUNNING"), ]
        
        # there can be no jobs ever run for a file id
        if running_jobs.empty : return
        
        # .. and cancel all of them.
        list(map(cancel, running_jobs["job_id"].unique().values))
        return
    
    def run(self, input_id, input_dir, output_dir, force=False, restart_fails=True, **kwargs):
        '''Submit a SLURM script (with list of args).
        '''        
        JOBID_REGEX = r"(?<=Submitted batch job )[0-9]*"
        
        # we design the input arguments for the slurm script
        script_kwargs = self.make_kwargs(input_id, input_dir, output_dir, **kwargs)
        args = [script_kwargs[k] for k in sorted(script_kwargs.keys())]

        # we avoid overwriting
        if not force:
            # .. of output
            if self.check_output(input_id, **script_kwargs)[1] > 0:
                if not self.silent:
                    print(f"Output for {input_id} already exists.")
                self.check(input_id, *args, restart_fails=restart_fails)
                return
            # .. or preexisting jobs
            elif self.check(input_id, *args, restart_fails=restart_fails):
                if not self.silent:
                    print(f"Jobs for {input_id} were already submitted. Check the outputs linked in logfile.") 
                return
        
        # we cancel running jobs if overwriting is forced.
        self.cancel_running_jobs(input_id, *args)
    
        # build command and run
        command_frags = ["sbatch", self.script, *args]
        srun = subprocess.Popen(
            command_frags, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sout, serr = srun.communicate()
            
        # report what was done
        if not self.silent:
            print(f"Command '{' '.join(command_frags)}' was submitted.\n")
            print("Response was:")
            print(f"[OUTPUT] {sout.decode('utf-8')}\n[ERROR] {serr.decode('utf-8')}\n--")
                
        # we extract the job_id and document the submission
        job_id_re = re.search(JOBID_REGEX, sout.decode('utf-8'))
        if job_id_re : job_id = job_id_re.group()
        else : raise RuntimeError(f"Did not manage to submit script for {input_id}.")
        
        self.check_job_id(job_id, input_id, *args, restart_fails=restart_fails)
        return
    
    def track_preexisting_output(self, input_id, **script_kwargs):
        '''Document pre exiting output in log file.
        '''
        # we write a custom row for preexisting output
        script_args = [script_kwargs[k] for k in sorted(script_kwargs.keys)]
        row = [
            "PREEXISTING_OUTPUT", "-", "2022-01-01T00:00:00", "-", "-", datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
            self.name, self.script, "-", "-", input_id, "-", output_name, str(output_size), *script_args]
        
        self.log_df.loc[len(self.log_df)+1, :] = row

        # write it down
        self.save_logfile()
        return
    
    # abstract methods
    @abstractmethod
    def set_default(self):
        '''Define name and argument names for slurm script.
        '''
        pass

    @abstractmethod
    def make_kwargs(self, input_id, input_dir, output_dir, **kwargs):
        '''Assemble script arguments.
        '''
        pass
    
    @abstractmethod
    def check_output(self, input_id, **script_kwargs):
        '''Check if output of file id already exisits.
        '''
        pass
# end SlurmSubmitter


class MFannotSubmitter(SlurmSubmitter):
    '''Class that submits MFannot jobs on slurm cluster.
    '''
    def set_default(self):
        '''Define name and argument names for slurm script.
        '''
        self.name = "MFANNOT_SINGLE_GENE"
        self.script_args = ["arg1_input_fasta", "arg2_output_name", "arg3_protein_collection_path"]
        return

    def make_kwargs(self, input_id, input_dir, output_dir, prot_coll_path=None, fasta_suffix="fa", masterfile_suffix="mf"):
        '''Assemble script arguments.
        '''
        # lazy kwarg check
        assert prot_coll_path is not None
        
        # make inputfile_path 
        input_fasta = os.path.join(input_dir, f"{input_id}.{fasta_suffix}")
        if not os.path.exists(input_fasta):
            raise FileNotFoundError(f"No {input_fasta}.")
        
        # make output file name
        output_name = os.path.join(output_dir, f"{input_id}.{masterfile_suffix}")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # check the peptide collection dir
        if not os.path.exists(prot_coll_path):
            raise FileNotFoundError(f"No path {prot_coll_path}.")
        
        # assemble all kwargs
        script_kwargs = dict(zip(
            self.script_args,
            [input_fasta, output_name, prot_coll_path]
        ))
        return script_kwargs
        
    def check_output(self, input_id, **kwargs):
        '''Return (provisional) output name and size.
        '''
        output_name = kwargs["arg2_output_name"]
        
        # by convention non existing files have -1 size
        if os.path.exists(output_name):
            output_size = os.stat(output_name).st_size
        else:
            output_size = -1
        return output_name, output_size
# end MFannotSubmitter



class asn2gbSubmitter(SlurmSubmitter):
    '''Class that submits asn2gb jobs on slurm cluster.
    '''
    def set_default(self):
        '''Define name and argument names for slurm script.
        '''
        self.name = "ASN_2_GB"
        self.script_args = ["arg1_input_sqn", "arg2_output_name"]
        return

    def make_kwargs(self, input_id, input_dir, output_dir, sqn_suffix="sqn", gb_suffix="gb", **kwargs):
        '''Assemble script arguments.
        '''
        
        # make inputfile_path 
        input_sqn = os.path.join(input_dir, f"{input_id}.mf.{sqn_suffix}")
        if not os.path.exists(input_sqn):
            raise FileNotFoundError(f"No {input_sqn}.")
        
        # make output file name
        output_name = os.path.join(output_dir, f"{input_id}.{gb_suffix}")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # assemble all kwargs
        script_kwargs = dict(zip(
            self.script_args,
            [input_sqn, output_name]
        ))
        return script_kwargs
        
    def check_output(self, input_id, **kwargs):
        '''Return (provisional) output name and size.
        '''
        output_name = kwargs["arg2_output_name"]
        
        # by convention non existing files have -1 size
        if os.path.exists(output_name):
            output_size = os.stat(output_name).st_size
        else:
            output_size = -1
        return output_name, output_size
# end asn2gbSubmitter


class SinglePartitionIQTreeSubmitter(SlurmSubmitter):
    '''Class that submits slurm jobs for single gene IQ tree analyses.
    '''
    def set_default(self):
        '''Define name and argument names for slurm script.
        '''
        self.name = "IQ_TREE_SINGLE_GENE"
        self.script_args = ["arg1_input_alignment", "arg2_output_base"]
        return
    
    def make_kwargs(self, input_id, input_dir, output_dir, fasta_suffix="fasta"):
        '''Assemble script arguments.
        '''
        # make inputfile_path 
        input_alignment = os.path.join(input_dir, f"{input_id}.{fasta_suffix}")
        if not os.path.exists(input_alignment):
            raise FileNotFoundError(f"No {input_alignment}.")
        
        # make output file name
        output_base = os.path.join(output_dir, f"{input_id}")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # assemble all kwargs
        script_kwargs = dict(zip(
            self.script_args,
            [input_alignment, output_base]
        ))
        return script_kwargs
    
    def check_output(self, input_id, **kwargs):
        '''Return (provisional) output name and size.
        '''
        output_base = kwargs["arg2_output_base"]
        output_suffices = ["treefile", "log", "ckp.gz", "iqtree", "bionj", "mldist"]
        
        tree_file = f"{input_id}.{output_suffices[0]}"
        # by convention non existing files have -1 size
        if all(map(os.path.exists, [f"{output_base}.{sfx}" for sfx in output_suffices])):
            output_size = os.stat(tree_file).st_size
        else:
            output_size = -1
        return tree_file, output_size
# end SinglePartitionIQTreeSubmitter

class MAFFTSubmitter(SlurmSubmitter):
    '''Class that submits slurm jobs for single gene MAFFT alignment.
    '''
    def set_default(self):
        '''Define name and argument names for slurm script.
        '''
        self.name = "MAFFT_SINGLE_GENE"
        self.script_args = ["arg1_input_fasta", "arg2_output_name", "arg3_amino_flag"]
        return
    
    def make_kwargs(self, input_id, input_dir, output_dir, amino_flag="--amino", fasta_suffix="fasta", aln_suffix="fasta"):
        '''Assemble script arguments.
        '''
        # lazy kwarg check
        assert amino_flag in ["--amino", ""]
        
        # make inputfile_path 
        input_fasta = os.path.join(input_dir, f"{input_id}.{fasta_suffix}")
        if not os.path.exists(input_fasta):
            raise FileNotFoundError(f"No {input_fasta}.")
        
        # make output file name
        output_name = os.path.join(output_dir, f"{input_id}.{aln_suffix}")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # assemble all kwargs
        script_kwargs = dict(zip(
            self.script_args,
            [input_fasta, output_name, amino_flag]
        ))
        return script_kwargs

    def check_output(self, input_id, **kwargs):
        '''Return (provisional) output name and size.
        '''
        output_name = kwargs["arg2_output_name"]
        
        # by convention non existing files have -1 size
        if os.path.exists(output_name):
            output_size = os.stat(output_name).st_size
        else:
            output_size = -1
        return output_name, output_size
# end MAFFTSubmitter

class TrimAlSubmitter(SlurmSubmitter):
    '''Class that submits slurm jobs for single gene alignment trimming using trimAl with a user-defined gap threshold.
    '''
    def set_default(self):
        '''Define name and argument names for slurm script.
        '''
        self.name = "TRIMAL_SINGLE_GENE_GAPFRACTION"
        self.script_args = ["arg1_input_fasta", "arg2_output_name"]
        return
    
    def make_kwargs(self, input_id, input_dir, output_dir, fasta_suffix="fasta", trimmed_suffix="fasta"):
        '''Assemble script arguments.
        '''
        # make inputfile_path 
        input_fasta = os.path.join(input_dir, f"{input_id}.{fasta_suffix}")
        if not os.path.exists(input_fasta):
            raise FileNotFoundError(f"No {input_fasta}.")
        
        # make output file name
        output_name = os.path.join(output_dir, f"{input_id}.{trimmed_suffix}")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # assemble all kwargs
        script_kwargs = dict(zip(
            self.script_args,
            [input_fasta, output_name]
        ))
        return script_kwargs

    def check_output(self, input_id, **kwargs):
        '''Return (provisional) output name and size.
        '''
        output_name = kwargs["arg2_output_name"]
        
        # by convention non existing files have -1 size
        if os.path.exists(output_name) and os.path.exists(f"{output_name}.html"):
            output_size = os.stat(output_name).st_size
        else:
            output_size = -1
        return output_name, output_size
# end TrimAlSubmitter


class MultiPartitionIQTreeSubmitter(SlurmSubmitter):
    '''Class that submits slurm jobs for single gene IQ tree analyses.
    '''
    def __init__(self, log_csv_file, modelfinder=False, silent=False):
        '''Initialize abstract instance by name and base_path.
        '''
        self.modelfinder = modelfinder
        super().__init__(log_csv_file, silent=silent)
        return
    
    def set_default(self):
        '''Define name and argument names for slurm script.
        '''
        if self.modelfinder : self.name = "IQ_TREE_BIG_MODELFINDER"
        else : self.name = "IQ_TREE_BIG"
        self.script_args = ["arg1_input_alignment", "arg2_partition_file", "arg3_output_base"]
        return
    
    def make_kwargs(self, input_id, input_dir, output_dir, partitionfile_path=None, fasta_suffix="fa", tree_base_name="out"):
        '''Assemble script arguments.
        '''
        # make inputfile_path 
        input_alignment = os.path.join(input_dir, f"{input_id}.{fasta_suffix}")
        if not os.path.exists(input_alignment):
            raise FileNotFoundError(f"No {input_alignment}.")
        
        # make output file name
        output_base = os.path.join(output_dir, input_id, tree_base_name)
        if not os.path.exists(os.path.join(output_dir, input_id)):
            os.makedirs(os.path.join(output_dir, input_id))
        
        # make partitionfile path
        if partitionfile_path is None:
            partitionfile_path = "_".join([os.path.splitext(input_alignment)[0], "partitions.nex"])
        assert os.path.exists(partitionfile_path)

        # assemble all kwargs
        script_kwargs = dict(zip(
            self.script_args,
            [input_alignment, partitionfile_path, output_base]
        ))
        return script_kwargs
    
    def check_output(self, input_id, **kwargs):
        '''Return (provisional) output name and size.
        '''
        output_base = kwargs["arg3_output_base"]
        output_suffices = ["treefile", "log", "ckp.gz", "iqtree", "bionj", "mldist"]
        
        tree_file = f"{output_base}.{output_suffices[0]}"
        # by convention non existing files have -1 size
        if all(map(os.path.exists, [f"{output_base}.{sfx}" for sfx in output_suffices])):
            output_size = os.stat(tree_file).st_size
        else:
            output_size = -1
        return tree_file, output_size
# end SinglePartitionIQTreeSubmitter


class DivvierSubmitter(SlurmSubmitter):
    '''Class that submits slurm jobs for single gene Divvier trimming.
    '''
    def set_default(self):
        '''Define name and argument names for slurm script.
        '''
        self.name = "DIVVIER_SINGLE_GENE"
        self.script_args = ["arg1_input_fasta"]
        return
    
    def make_kwargs(self, input_id, input_dir, output_dir, fasta_suffix="fa", aln_suffix="fa"):
        '''Assemble script arguments.
        '''
        # make inputfile_path 
        input_fasta = os.path.join(input_dir, f"{input_id}.{fasta_suffix}")
        if not os.path.exists(input_fasta):
            raise FileNotFoundError(f"No {input_fasta}.")
        
        # make output file name
        output_name = os.path.join(output_dir, f"{input_id}.{aln_suffix}")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        self.op = output_name
        
        # assemble all kwargs
        script_kwargs = dict(zip(
            self.script_args,
            [input_fasta]
        ))
        return script_kwargs

    def check_output(self, input_id, **kwargs):
        '''Return (provisional) output name and size.
        '''
        
        output_names = [f'{kwargs["arg1_input_fasta"]}.{sfx}' for sfx in 
                        ["PP", "partial.fas"]]
        
        # by convention non existing files have -1 size
        if os.path.exists(output_names[1]):
            output_size = os.stat(output_names[1]).st_size
        else:
            output_size = -1
        
        if output_size > -1 :
            for of, sfx in zip(output_names, [".PP", ""]):
                shutil.move(of, f"{self.op}{sfx}")
        if os.path.exists(self.op):
            output_size = os.stat(self.op).st_size
        
        return self.op, output_size
# end DivvierSubmitter


class BMGESubmitter(SlurmSubmitter):
    '''Class that submits slurm jobs for single gene BMGE trimming.
    '''
    def set_default(self):
        '''Define name and argument names for slurm script.
        '''
        self.name = "BMGE_SINGLE_GENE"
        self.script_args = ["arg1_input_fasta", "arg2_output_name"]
        return
    
    def make_kwargs(self, input_id, input_dir, output_dir, fasta_suffix="fasta", aln_suffix="fasta"):
        '''Assemble script arguments.
        '''
        # make inputfile_path 
        input_fasta = os.path.join(input_dir, f"{input_id}.{fasta_suffix}")
        if not os.path.exists(input_fasta):
            raise FileNotFoundError(f"No {input_fasta}.")
        
        # make output file name
        output_name = os.path.join(output_dir, f"{input_id}.{aln_suffix}")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # assemble all kwargs
        script_kwargs = dict(zip(
            self.script_args,
            [input_fasta, output_name]
        ))
        return script_kwargs

    def check_output(self, input_id, **kwargs):
        '''Return (provisional) output name and size.
        '''
        output_name = kwargs["arg2_output_name"]
        
        # by convention non existing files have -1 size
        if os.path.exists(output_name):
            output_size = os.stat(output_name).st_size
        else:
            output_size = -1
        return output_name, output_size
# end BMGESubmitter