# -*- coding: utf-8 -*-
"""
configfile.py - Human-readable text configuration file library 
Copyright 2010  Luke Campagnola
Distributed under MIT/X11 license. See license.txt for more infomation.

Used for reading and writing dictionary objects to a python-like configuration
file format. Data structures may be nested and contain any data type as long
as it can be converted to/from a string using repr and eval.
"""

import re, os, sys, datetime
import numpy
from collections import OrderedDict
# from .pgcollections import OrderedDict
# from . import units
# from python2_3 import asUnicode, basestring
# from .Qt import QtCore
# from .Point import Point
# from .colormap import ColorMap
# GLOBAL_PATH = None # so not thread safe.

# -*- coding: utf-8 -*-
## Very simple unit support:
##  - creates variable names like 'mV' and 'kHz'
##  - the value assigned to the variable corresponds to the scale prefix
##    (mV = 0.001)
##  - the actual units are purely cosmetic for making code clearer:
##  
##    x = 20*pA    is identical to    x = 20*1e-12

## No unicode variable names (μ,Ω) allowed until python 3

SI_PREFIXES = 'yzafpnum kMGTPEZY'
UNITS = 'm,s,g,W,J,V,A,F,T,Hz,Ohm,S,N,C,px,b,B'.split(',')
allUnits = {}

def addUnit(p, n):
    g = globals()
    v = 1000**n
    for u in UNITS:
        g[p+u] = v
        allUnits[p+u] = v
    
for p in SI_PREFIXES:
    if p ==  ' ':
        p = ''
        n = 0
    elif p == 'u':
        n = -2
    else:
        n = SI_PREFIXES.index(p) - 8

    addUnit(p, n)

cm = 0.01

def evalUnits(unitStr):
    """
    Evaluate a unit string into ([numerators,...], [denominators,...])
    Examples:
        N m/s^2   =>  ([N, m], [s, s])
        A*s / V   =>  ([A, s], [V,])
    """
    pass
    
def formatUnits(units):
    """
    Format a unit specification ([numerators,...], [denominators,...])
    into a string (this is the inverse of evalUnits)
    """
    pass
    
def simplify(units):
    """
    Cancel units that appear in both numerator and denominator, then attempt to replace 
    groups of units with single units where possible (ie, J/s => W)
    """
    pass
    

def ColorMap(p, c):
    return p, list(c)

def Point(p):
    return list(p)






class ParseError(Exception):
    def __init__(self, message, lineNum, line, fileName=None):
        self.lineNum = lineNum
        self.line = line
        #self.message = message
        self.fileName = fileName
        Exception.__init__(self, message)
        
    def __str__(self):
        if self.fileName is None:
            msg = "Error parsing string at line %d:\n" % self.lineNum
        else:
            msg = "Error parsing config file '%s' at line %d:\n" % (self.fileName, self.lineNum)
        msg += "%s\n%s" % (self.line, self.message)
        return msg
        #raise Exception()
        

def writeConfigFile(data, fname):
    s = genString(data)
    fd = open(fname, 'w')
    fd.write(s)
    fd.close()
    
def readConfigFile(fname):
    #cwd = os.getcwd()

    if not os.path.exists(fname):
        fname = fname2        
    try:
        #os.chdir(newDir)  ## bad.
        fd = open(fname)
        s = fd.read() # .encode("latin1")
        s = str(s)
        fd.close()
        s = s.replace("\r\n", "\n")
        s = s.replace("\r", "\n")
        data = parseString(s)[1]
    except ParseError:
        sys.exc_info()[1].fileName = fname
        raise
    except:
        print("Error while reading config file %s:"% fname)
        raise
    #finally:
        #os.chdir(cwd)
    return data

def appendConfigFile(data, fname):
    s = genString(data)
    fd = open(fname, 'a')
    fd.write(s)
    fd.close()


def genString(data, indent=''):
    s = ''
    for k in data:
        sk = str(k)
        if len(sk) == 0:
            print(data)
            raise Exception('blank dict keys not allowed (see data above)')
        if sk[0] == ' ' or ':' in sk:
            print(data)
            raise Exception('dict keys must not contain ":" or start with spaces [offending key is "%s"]' % sk)
        if isinstance(data[k], dict):
            s += indent + sk + ':\n'
            s += genString(data[k], indent + '    ')
        else:
            s += indent + sk + ': ' + repr(data[k]) + '\n'
    return s
    
def parseString(lines, start=0):
    
    data = OrderedDict()
    # lines = str(lines)
    if isinstance(lines, str):
        lines = lines.split('\n')
        lines = [l for l in lines if re.search(r'\S', l) and not re.match(r'\s*#', l)]  ## remove empty lines
    indent = measureIndent(lines[start])
    ln = start - 1
    
    try:
        while True:
            ln += 1
            #print ln
            if ln >= len(lines):
                break
            
            l = lines[ln]
            
            ## Skip blank lines or lines starting with #
            if re.match(r'\s*#', l) or not re.search(r'\S', l):
                continue
            
            ## Measure line indentation, make sure it is correct for this level
            lineInd = measureIndent(l)
            if lineInd < indent:
                ln -= 1
                break
            if lineInd > indent:
                #print lineInd, indent
                raise ParseError('Indentation is incorrect. Expected %d, got %d' % (indent, lineInd), ln+1, l)
            
            ## set up local variables to use for eval
            local =allUnits.copy()
            local['OrderedDict'] = OrderedDict
            local['readConfigFile'] = readConfigFile
            local['Point'] = Point
            # local['QtCore'] = QtCore
            local['ColorMap'] = ColorMap
            local['datetime'] = datetime
            # Needed for reconstructing numpy arrays
            local['array'] = numpy.array
            for dtype in ['int8', 'uint8', 
                          'int16', 'uint16', 'float16',
                          'int32', 'uint32', 'float32',
                          'int64', 'uint64', 'float64']:
                local[dtype] = getattr(numpy, dtype)            
            
            if ':' not in l:
                raise ParseError('Missing colon', ln+1, l)
            
            (k, p, v) = l.partition(':')
            k = k.strip()
            v = v.strip()

            if len(k) < 1:
                raise ParseError('Missing name preceding colon', ln+1, l)
            if k[0] == '(' and k[-1] == ')':  ## If the key looks like a tuple, try evaluating it.
                try:
                    k1 = eval(k, local)
                    if type(k1) is tuple:
                        k = k1
                except:
                    pass
            if re.search(r'\S', v) and v[0] != '#':  ## eval the value
                if (v[0] == '[' or v[0] == '(') and (v[-1] == ']' or v[-1] == ')'):
                     v = v.replace('L', '')
                try:
                    val = eval(v, local)
                except:
                    ex = sys.exc_info()[1]
                    raise ParseError("Error evaluating expression '%s': [%s: %s]" % (v, ex.__class__.__name__, str(ex)), (ln+1), l)
            else:
                if ln+1 >= len(lines) or measureIndent(lines[ln+1]) <= indent:
                    #print "blank dict"
                    val = {}
                else:
                    #print "Going deeper..", ln+1
                    (ln, val) = parseString(lines, start=ln+1)
            data[k] = val
        #print k, repr(val)
    except ParseError:
        raise
    except:
        ex = sys.exc_info()[1]
        raise ParseError("%s: %s" % (ex.__class__.__name__, str(ex)), ln+1, l)
    #print "Returning shallower..", ln+1
    return (ln, data)
    
def measureIndent(s:str):
    n = 0
    while n < len(s) and s[n] == ' ':
        n += 1
    return n
    
def test():

    fn = "/Users/pbmanis/Desktop/2018.09.27_000/ImageSequence_000/.index"
    data = readConfigFile(fn)
    print(data)
    fn = "/Users/pbmanis/Desktop/Python/mrk-nf107/data_for_testing/CCIV/.index"
    data = readConfigFile(fn)
    print(data)
    
if __name__ == '__main__':
    import tempfile
    fn = tempfile.mktemp()
    tf = open(fn, 'w')
    cf = """
key: 'value'
key2:              ##comment
                   ##comment
    key21: 'value' ## comment
                   ##comment
    key22: [1,2,3]
    key23: 234  #comment
    """
    tf.write(cf)
    tf.close()
    print("=== Test:===")
    num = 1
    for line in cf.split('\n'):
        print("%02d   %s" % (num, line))
        num += 1
    print(cf)
    print("============")
    data = readConfigFile(fn)
    # print(data)
    os.remove(fn)

# now test some real files

