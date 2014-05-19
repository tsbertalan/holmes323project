'''
Some simple classes for handling the loading of parameter files.

Created on May 8, 2014
@author: tsbertalan
'''

def latexifyUnderscore(tex):
    assert tex[0] == tex[-1] == "$", "The string must begin and end in a dollar sign."  # TODO obviously not a real requirement
    assert len(tex) >= 2, "There must be 2 or more dollar signs in the string."
    if '_' in tex:
        tex = tex[1:-1].split('_')
        assert len(tex) == 2, "Only one or zero underscores are allowed."
        return '$' + tex[0] + '_{' + tex[1] + '}$'
    else:
        return tex


class Parameter(object):
    '''A trivial mixin class that adds parameter name fields.
    Meant to inherit from some numeric type like int or float.
    '''
    def __init__(self, *args, **kwargs):
        super(Parameter, self).__init__(*args, **kwargs)
        self._label = ""
        self._tex = None
    
    def slabel(self, label):
        self._label = label
    
    def stex(self, tex):
        self._tex = tex
    
    def glabel(self):
        return self._label
    
    def gtex(self):
        if self._tex is None:
            return latexifyUnderscore('$' + self.glabel() + '$')
        else:
            return self._tex
    
    def gval(self):  # note that there is not and should not be a sval()
        raise NotImplementedError('This method is virtual.')
#     def __repr__(self):
#         return self.label


class FloatParameter(float, Parameter):
    '''
    >>> m = FloatParameter(3.45)
    >>> m.slabel('em')
    >>> m.gval()
    3.45
    >>> m.gtex()
    '$em$'
    >>> m.stex('$m$')
    >>> m.gtex()
    '$m$'
    '''
    def gval(self):
        return float(self)


class IntParameter(int, Parameter):
    def gval(self):
        return int(self)


class Parameters(object):
    '''
    a parameters class:
    stackoverflow.com/questions/2597278/python-load-variables-in-a-dict-into-namespace
    '''
    def __init__(self, intParams=False):
        self.pnames = set()
        self.intParams = intParams
    
    def update(self, newVars):
        self.__dict__.update(newVars)
        
    def __getitem__(self, item):
        return self.__dict__[item]
    
    def addFromFile(self, filePath):
        '''File should have format something like this:
        
        a = 12
        b=64.0
        foo=-8
        
        with no blank lines. Later, we might split on hash marks to allow
        comments or descriptions after the paramter values.
        if the field intParams, then whether parameters are considered as floats or ints
        is determined by whether they contain a period or not. 
        '''
        lines = open(filePath).readlines()
        return self.addFromList(lines)
        
    def addFromList(self, lines):
        
        dictFromList = {}
        for line in lines:
            if len(line.strip()) > 0:  # ignore whitespace lines
                assert '=' in line, "line '%s' does not contain a '='" % line
                items = line.split('=')
                assert len(items) >= 2
                key, val = items[:2]
                key = key.strip()
                val = val.split('#')[0]  # ignore everything after a hashmark
                dictFromList[key] = val
                self.pnames.add(key)
        self.addFromDict(dictFromList)
             
    def addFromDict(self, updateDict):
        '''
        updateDict can contain either regular floats and ints, or the extended
        versions above. They'll get upgraded.
        ''' 
        upgradedDict = {}
        for key, val in updateDict.items():
            key = key.strip()
            if not self.intParams:  # TODO this should also accept the case that 
                ParamType = FloatParameter
            else:
                if isinstance(val, int):
                    ParamType = IntParameter
                elif isinstance(val, float):
                    ParamType = FloatParameter
                elif isinstance(val, str) and '.' not in val:
                    ParamType = IntParameter
                else:
                    ParamType = FloatParameter
                
            upgradedDict[key] = ParamType(val)
            upgradedDict[key].slabel(key)
        self.update(upgradedDict)
        
    def __repr__(self):  # this doesn't work for pretty-printing.
                         # Just call pp(self.__dict__) directly.
        return "Parameters:" + self.__dict__.__repr__()

    def __str__(self):
        return "Parameters:" + self.__dict__.__str__()


if __name__=="__main__":
    from doctest import testmod
    testmod()
    

