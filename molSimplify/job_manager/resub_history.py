import pickle
import os

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
