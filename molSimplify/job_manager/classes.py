import pickle
import os

def try_float(obj):
    # Converts an object to a floating point if possible
        try:
            floating_point = float(obj)
        except:
            floating_point = obj
        return floating_point
        
def strip_new_line(string):
        if string[-1] == '\n':
            return string[:-1]
        else:
            return string

class resub_history:
    # Class for saving information about previous resubmissions
    # should be initialized like:
    # resub = resub_history()
    # resub.read(outfile_path) do this step even if no pickel file exists already
    # update stored values
    # resub.save()
    
    def __init__(self,path = None):
        self.resub_number = 0
        self.status = 'Normal'
        self.needs_resub = False
        self.notes = []
        self.outfiles = []
        self.path = path
        
    def save(self):
        if self.path == None:
            raise Exception('The path for the resub_history pickel file is not specified!')
        with open(self.path,'wb') as handle:
            pickle.dump(self, handle, protocol = pickle.HIGHEST_PROTOCOL)
            
    def read(self,path):
        if path.endswith('.out'):
            path = path.rsplit('.',1)[0]+'.pickle'
        
        if os.path.isfile(path):
            with open(path,'rb') as handle:
                saved = pickle.load(handle)
            self.resub_number = saved.resub_number
            self.status = saved.status
            self.needs_resub = saved.needs_resub
            self.notes = saved.notes
            self.outfiles = saved.outfiles
        self.path = path

class textfile:
    ##Class for importing textfiles in a searchable way.
    def __init__(self,file_name=None):
        
        if file_name:
            raw_file = open(file_name,'r')
            self.lines = raw_file.readlines()
            raw_file.close()
            self.lines = [strip_new_line(i) for i in self.lines]
        
        else:
            self.lines = None
            
    def wordgrab(self,keywords,indices,last_line=False,first_line = False,min_value = False):
        ## takes two lists as an input
        # The first list is the keywords to look for
        # The second list is the indices to pull from the matching lines
        #  Returns a list of the resulting values. Numeric values are automatically converted to floats
        
        if type(keywords) != list:
            keywords = [keywords]
        if type(indices) != list:
            indices = [indices]
            
        results = dict()
        zipped_values = zip(keywords,indices)
        
        for line in self.lines:
            for keyword,index in zipped_values:
                if keyword in line:
                    
                    if type(index) == int:
                        matching_value = try_float(line.split()[index])
                    else:
                        matching_value = line.split()
                    
                    if keyword not in results.keys():
                        results[keyword] = [matching_value]
                    else:
                        results[keyword].append(matching_value)
        
        if (last_line and min_value) or (last_line and first_line) or (first_line and min_value):
            raise ValueError('Warning, incompatible options selected in text parsing')
            
        if last_line:
            for keyword in results.keys():
                results[keyword] = results[keyword][-1]
        if first_line:
            for keyword in results.keys():
                results[keyword] = results[keyword][0]
        if min_value:
            for keyword in results.keys():
                results[keyword] = min(results[keyword])
                
        results_to_return = []
        for key in keywords:
            if key in results.keys():
                results_to_return.append(results[key])
            else:
                if last_line or min_value or first_line:
                    results_to_return.append(None)
                else:
                    results_to_return.append([None])
        
        return results_to_return
