from __future__ import division
import sys
sys.path.append('../../')
from tools.userError import *
from copy import deepcopy

class compartment(object):
    """
    A class holding the information about a compartment (e.g., cytosol, mitochondria, etc) 

    Ali R. Zomorrodi - Segre Lab @ BU
    Last updated: 06-02-2016
    """

    def __init__(self, id, name = '', standard_id = '', name_aliases = [], model_id = '', notes = '', stdout_msgs = False): 
    
        # Gene id 
        self.id = id

        # na,e 
        if name == None or name == '':
            self.name = self.id
        else:
            self.name = name

        # Messages in the output
        self.stdout_msgs = stdout_msgs

        # Standard id: This is the standard compartment ids as detailed below
        # c: Cytosol (cytoplasm),   e: Extracellular,   g: Golgi,     m: Mitochondria
        # n: Nucleus,   p: Periplasm,    r: Endoplasmic reticulum,    x: Peroxisome
        if standard_id != '':
            self.standard_id = standard_id

        # If standard_id is not provided, try to find it out from id or name
        else:
            self.guess_standard_id() 

        # name aliases
        self.name_aliases = name_aliases

        # The model in which this compartment is used
        self.model_id = model_id

        # Notes and comments
        self.notes = notes

    def __setattr__(self,attr_name,attr_value):
        """
        Redefines funciton __setattr__
        INPUTS:
        -------
        attr_name: Attribute name
        attr_value: Attribute value
        """
        # id 
        if attr_name == 'id' and not isinstance(attr_value,str):
            raise TypeError("Invalid 'id' for compartment " + str(attr_value) + "! 'id' must be a string. A " + str(type(attr_value)) + " type object was entered instead")

        # Standard_id
        if attr_name == 'standard' and (attr_value is not None and not isinstance(attr_value,str)):
            raise TypeError("Invalid 'standard_id' for compartment " + self.id + "! 'standard_id' must be a string. A " + str(type(attr_value)) + " type object was entered instead")
        elif attr_name == 'standard_id' and attr_value not in ['c', 'p', 'e', 'g', 'm', 'n', 'p', 'r', 'x']:
            raise ValueError('Invalid value for standard_id. Allowed values are: {}'.format(['c', 'e', 'g', 'm', 'n', 'p', 'r', 'x']))   
 
        # Name
        if attr_name == 'name' and (attr_value is not None and not isinstance(attr_value,str)):
            raise TypeError("Invalid 'name' for compartment " + self.id + "! 'name' must be a string. A " + str(type(attr_value)) + " type object was entered instead")
    
        # Model id
        if attr_name == 'model_id' and (attr_value is not None and not isinstance(attr_value,str)):
            raise TypeError("Invalid 'model_id' for compartment " + self.id + "! 'model_id' must be a string. A " + str(type(attr_value)) + " type object was entered instead")

        # stdout_msgs
        if attr_name == 'stdout_msgs' and not isinstance(attr_value, bool):
            raise TypeError("Invalid stdout_msgs! True or False expected, a variable of type {} was entered instead".format(type(stdout_msgs)))
    
        self.__dict__[attr_name] = attr_value


    def guess_standard_id(self):
        """
        Guesses the standard id of the compartment
        """
        standard_id = ''
        if self.id.lower() in ['c0','[c]', '_c'] or self.name.lower() in ['cytosol', 'cytoplasm']: 
            standard_id = 'c'
        elif self.id.lower() in ['p0','[p]', '_p'] or self.name.lower() == 'periplasm': 
            standard_id = 'p'
        elif self.id.lower() in ['e0','[e]', '_e'] or self.name.lower() in ['extracellular', 'extra cellular', 'extra-cellular'] or 'extra' in self.name.lower(): 
            standard_id = 'e'
        elif self.id.lower() in ['g0','[g]', '_g'] or self.name.lower() == 'golgi': 
            standard_id = 'g'
        elif self.id.lower() in ['m0','[m]', '_m'] or self.name.lower() == 'mitochondria': 
            standard_id = 'm'
        elif self.id.lower() in ['n0','[n]', '_n'] or self.name.lower() == 'nucleus': 
            standard_id = 'm'
        elif self.id.lower() in ['r0','[r]', '_r'] or self.name.lower() in ['endoplasmic reticulum', 'Endoplasmic_reticulum'] or 'endoplasmic' in self.name.lower() or 'reticulum' in self.name.lower(): 
            standard_id = 'm'
        elif self.id.lower() in ['x0','[x]', '_x'] or self.name.lower() == 'peroxisome': 
            standard_id = 'x'

        if standard_id == '':
            print "WARNING! (compartment.py): No standard_id could be assigned to compartment with id = {} , name = {}".format(self.id, self.name)
        else:
            self.standard_id = standard_id
            if self.stdout_msgs:
                print "{} was assigned as standard_id to compartment with id = {} , name = {}".format(standard_id, self.id, self.name)
       

