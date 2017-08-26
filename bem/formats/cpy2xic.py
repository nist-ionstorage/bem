import math     # Math functions
import re       # Pattern matching
import string   # String functions
import sys      # need the exit() function
import optparse # Command line parsing

from TrapGeometry import * # Jason: all the actual code has been moved here 

#################################################
## Parse command line
#################################################
def standaloneDriver():
    usage="usage: %prog cpy_file xic_file"
    parser = optparse.OptionParser(usage=usage)
    (options,args) = parser.parse_args()

    # Get the file names
    if len(args)!=2:
        parser.error("Needs two arguments.")

    file_cpy = args[0]
    file_xic = args[1]
    
    # Load .CPY file
    trapGeometry=TrapGeometry()
    trapGeometry.readCPY(file_cpy,re.compile('(.*)') ) 
    
    # Output geometry to xic
    trapGeometry.savexic(file_xic)
    

if __name__ == "__main__":
    standaloneDriver()
