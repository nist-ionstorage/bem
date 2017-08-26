#! python
#No comments before the klafbang, please

import sys
import os, os.path, shutil
import re
import string
import optparse
import codecs #needed to write xml (2-byte unicode)
import math
import copy
import operator
import random
    
class WiringDB:
    def __init__(self):
    	self.columnnames={} #Map: column name -> index
    	self.tuples=[]      #List of lists.
        pass

    def readFile(self, filename):
        columnnames={}
        tuples=[]
        hasheader=False
        numcols=0
        inputfile=file(filename,'r')
        for line in inputfile:
            if line.startswith('#') or len(line.strip())==0 :
                continue
            #We have a content line:
            v=map(string.strip, line.split(','))
            if not hasheader:
                hasheader=True
                numcols=len(v)
                for (idx,name) in enumerate(v):
                    self.columnnames[name.lower()]=idx
            else:
                buf=['']*numcols
                buf[0:len(v)]=v
                self.tuples.append(buf)
        inputfile.close()

    def getColumn(self, columnname):
        idx=self.columnnames[columnname.lower()] #Exception if unknown
        return map( lambda v : v[idx], self.tuples )
#End of class WiringDB    

def makeStyXML(enames, outputfilename):
    s=u''
    s=s+"""<?xml version="1.0" encoding="utf-16"?>""" + '\n'
    s=s+"<Export>"+'\n'
    for ename in enames:
	print 'Adding electrode: %s\n' %ename
        s=s+"<Colors Version=\"11\">"
        s=s+"<Style InternalName=\"1:"
        s=s+ename
        ars=''
        for i in range(3):
            ars=ars+str(round(random.random(),6))+' '
        ars=ars+'1'
        s=s+'\" EditFlag=\"0\" Ambient=\"'+ars+'\" Diffuse=\"0.560784 0.521569 0.537255 1\" Specular=\"0.780392 0.780392 0.780392 1\" Emissive=\"0 0 0 1\" Shininess=\"4\" '
        s=s+""" IsMasked="0" Mask="0 0 0 1" Opacity="100" Scale="1" Rotation="0" TextureFile="AlloySand2.bmp" TexturePath="1" ShowBackFaces="0" ReflectionImageFile="" ReflectionImagePath="0" Refraction="1" """
        s=s+""" BumpMapImageFile="" BumpMapImagePath="0" BumpMapAmount="1" BumpMapScale="1" BumpMapRotation="0"> """
        s=s+ """<Name Language="1033" Description="">"""
        s=s+ename 
        s=s+"</Name>"
        s=s+"</Style>"
        s=s+"</Colors>"+'\n'
    s=s+"</Export>"+'\n'
    #write to file -- XML requires 2 byte unicode
    vFile = codecs.open(outputfilename,"w","utf-16")
    vFile.write(s)
    vFile.close()

#############################
# MAIN 
#############################

parser=optparse.OptionParser(usage=
"""%prog [options] file1.wdb [file2.wdb [...]]
A script for manipulating wiring databases.
Currently, only function is to generate a .styxml file
based on the NAME column.

.wdb Format:
Ignored lines: Empty lines and lines starting with #
Content lines: <v1> [,<v2> [,...]]
First content line must be column names.
Remaining content lines are db tuples.
- Missing fields are set to empty string.
- Extra fields are dropped

Column names should be alphanummeric. Case insensitive.
Suggested convention:
ELECTRODE_NAME : Electrode name. Form should be DC*, RF, GND
CHIP_PIN       : Chip carrier pin if applicable
CABLE_PIN      : Format: Ann, Bnn, ...
DAC_CHANNEL    : Format: Devn/n, ...

The following preextensions are suggested:
_trap.wdb   : trap config:    ELECTRODE_NAME [CHIP_PIN]
_vacum.wdb  : vacum config:   CHIP_PIN CABLE_PIN
_wiring.wdb : mux box config: CABLE_PIN DAC_CHANNEL
""")

parser.add_option("-x","--writeStyXML", dest="styxml", action ="store_true",
                  help="Generate .styxml file based on ELECTRODE_NAME channel")
parser.add_option("--StyXMLFile", dest="styxmlfilename", type="string",
                  help="Name of .styxml. Default: <input1>.styxml")
parser.add_option("--StyXMLPrefix", dest="styxmlprefix", type="string",
                  help="Electrode name prefix for .styxml. Default: %default")
parser.set_defaults(styxml=False)
parser.set_defaults(styxmlprefix="TRAPELECTRODE_")
(options,args)=parser.parse_args()

if len(args)<1:
    parser.error('Incorrect number of arguments')
# Read the datebase:
inputdb=WiringDB()
inputdb.readFile(args[0])

if not options.styxmlfilename :
    [options.styxmlfilename,inputfileext]=os.path.splitext(args[0])
    options.styxmlfilename+='.styxml'

if options.styxml :
    enames=inputdb.getColumn('ELECTRODE_NAME')
    enames=map(lambda s: options.styxmlprefix+s,enames)
    makeStyXML(enames,options.styxmlfilename)


