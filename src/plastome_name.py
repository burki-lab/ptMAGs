class PlastomeName:
    '''PlastomeName helps to easily access simple Metadata of a pastome sequence by its name. 
    '''
    def __init__(self, name):
        '''Initialize PlastomeName by name.
        '''
        self.name = name
        return
        
    def get_data(self):
        '''Return georgraphic region a MAG was constructed from.
        '''
        if hasattr(self, "region"):
            return {self.region: self.bin}
        if self.name.startswith("CHL_"):
            splits = self.name.split("_")
            self.region = splits[1]
            if splits[5].startswith("c"):
                self.bin = {int(splits[3]): int(splits[4])}
            elif splits[6].startswith("c"):
                self.bin = {int(splits[3]): {
                    int(splits[4]): int(splits[5])}}
            else : self.bin = None
        else:
            self.region = "REFERENCE"
            self.bin = "REFERENCE"
        return {self.region: self.bin}
    
    def is_region(self, region):
        '''Check if a plastome was found in a specific region.
        '''
        if not hasattr(self, "region") : self.get_data()
        return self.region == region
# end PlastomeName