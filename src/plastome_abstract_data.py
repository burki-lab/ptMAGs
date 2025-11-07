# here one finds a class that helps to iterate
# through all plastomes and it follow up analyses

import os
from abc import ABC
from abc import abstractmethod
import csv

from plastome_name import PlastomeName

class PlastomeAbstractFile(ABC):
    '''This is an abstract container class for plastome files.
    '''
    def __init__(self, name, base_path, suffix=None, alias="", silent=False):
        '''Initialize abstract instance by name and base_path.
        '''
        file_name = os.path.join(base_path, name)
        if os.path.exists(file_name):
            self.dir = base_path
            self.name = name
            self.file = file_name
            self.suffix = suffix
            self.make_id()
            self.meta = PlastomeName(name)
            self.alias = alias
            self.silent = silent
        else : ValueError(f"Requested file {file_name} does not exist.")
        return
    
    
    # overall methods
    def get_size(self):
        '''Obtain size of the file.
        '''
        return os.path.getsize(self.file)
    
    def make_id(self):
        '''Define the ID of a file given its name.
        '''
        base = [frg for frg in os.path.split(os.path.relpath(self.file, self.dir)) if frg!=""][0]
        self.id = os.path.splitext(base)[0]
        return
    
    def run_slurm_job(self, slurm_submitter, output_dir, force=False, **kwargs):
        '''Submit slurm job using SlurmSubmitter subclasses.
        '''
        from utils_slurm import SlurmSubmitter
        slurm_submitter.silent = self.silent
        slurm_submitter.run(self.id, self.dir, output_dir, force=force, **kwargs)
        return
        
    # abstract methods
    @abstractmethod
    def summarize(self):
        '''Summarize file specific information.
        '''
        return
# end PlastomeAbstractFile


class PlastomeAbstractIterator(ABC):
    '''This is an iterator for PlastomeAbstractFiles.
    '''
    def __init__(self, directory, file_container=PlastomeAbstractFile,
                 suffix=None, alias="", silent=False, **kwargs):
        '''Initialize abstract iterator by providing a directory.
        '''
        if not os.path.exists(directory):
            ValueError(f"Requested directory {directory} does not exist.")
        self.base_path = directory
        
        # some kinds of data organization are multiple files in a dir 
        # instead of just simple files.
        # thus we allow to list the directory in a "more recursive"
        # second depth level
        if "extended" in kwargs.keys() : self.extended = kwargs.pop("extended")
        else : self.extended = False
        
        self.suffix = suffix
        self.load_files()
        self.alias = alias
        
        self.container = file_container
        self.kwargs = kwargs
        
        self.silent = silent
        return
    
    
    # abstract methods
    @abstractmethod
    def load_files(self):
        '''Select relevant files within directory that should be iterated through.
        '''
        if self.extended : self.files = [
            os.path.join(b, f)  for b, dl, fl in os.walk(self.base_path) for f in fl]
        else : self.files = [os.path.join(self.base_path, f) for f in os.listdir(self.base_path)]
        
        if self.suffix is None : return
    
        self.files = sorted([os.path.relpath(fl, self.base_path) for fl in self.files if fl.endswith(self.suffix)])
        self.ids = [idn for idn, sf in map(os.path.splitext, self.files)]
        if self.extended : self.ids = [os.path.split(idn)[0] for idn in self.ids]
        return
    
    @abstractmethod
    def store_information(self, csv_file_path, force=False, buffered=True, **kwargs):
        '''Store the class specific information for all plastomes.
        '''
        self.check_csv_path(csv_file_path, force=force)
        
        with open(csv_file_path, "w") as cfile:
            header_index = 0
            i=0
            self.__iter__()
            
            # we can buffer mistaken iteration
            while True:
                try: plastome = self.__next__()
                except FileNotFoundError as err:
                    if buffered: 
                        header_index += 1
                        continue
                    else : raise err
                except StopIteration : break
                
                # then we summarize and write
                row = plastome.summarize(**kwargs)

                if i == header_index:
                    header = row.keys()
                    csv_writer = csv.DictWriter(cfile, fieldnames=header)
                    csv_writer.writeheader()
                csv_writer.writerow(row)
                i += 1
        return

    # iterator functions
    def __iter__(self):
        '''Iterate through all files.
        '''
        self.i = 0
        return self
    
    def __next__(self):
        '''Obtain next element in list as PlastomeAbstractFile.
        '''
        if self.i < len(self.files):
            self.i += 1
            current = self.__getitem__(self.i-1)
            return current
        else : raise StopIteration
        return
        
    def __getitem__(self, key):
        '''Obtain an element by given indexing key.
        '''
        # integer indexing
        if isinstance(key, int) and key >= 0 and key < len(self.files):
            return self.container(
                self.files[key], self.base_path, suffix=self.suffix, alias=self.alias, **self.kwargs)
        # id indexing
        elif isinstance(key, str) and key in self.ids:
            key_index = [i for i, el in enumerate(self.ids) if el==key][0]
            return self.__getitem__(key_index)
        
        raise IndexError(f"Wrong index: {key}")
    
    # additional methods
    def check_csv_path(self, csv_file_path, force=False):
        '''Check the path of a given csv file.
        '''
        if os.path.exists(csv_file_path) and not force:
            raise IOError("Information will not be stored, "
                          f"{csv_file_path} already exists.")
    
    def select_regions(self, regions):
        '''Subset the files by only occuring in given region(s).
        '''
        # small function that checks if plastome comes from any of given regions
        check_regions = lambda plname, regs: any(
            [plname.is_region(reg) for reg in regs])
        if isinstance(regions, str) : regions = [regions]
        self.files = [pl for pl in self.files if check_regions(pl, regions)]
        return 
# end PlastomeAbstractIterator
