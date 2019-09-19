##Class for importing textfiles in a searchable way.
def try_float(obj):
    # Converts an object to a floating point if possible
        try:
            floating_point = float(obj)
        except:
            floating_point = obj
        return floating_point

class textfile:
    
    def __init__(self,file_name=None):
        
        if file_name:
            raw_file = open(file_name,'r')
            self.lines = raw_file.readlines()
            raw_file.close()
        
        else:
            self.lines = None
            
    def wordgrab(keywords,indices,last_line=False,first_line = False,min_value = False):
        ## takes two lists as an input
        # The first list is the keywords to look for
        # The second list is the indices to pull from the matching lines
        #  Returns a list of the resulting values. Numeric values are automatically converted to floats
        
        if type(keywords) != list:
            keywords = [keywords]
        if type(indices) != list:
            indices = [indices]
            
        results = dict()
        zipped_values == zip(keywords,indices)
        
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
            if key in results.keys()
                results_to_return.append(results[key])
            else:
                results[key] = None
            
        return results_to_return
