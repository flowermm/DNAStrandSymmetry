'''
RMfileReader.py

Created on Dec 26, 2011

@author: Kristin Helling
'''

import sys
import re
import string
import bz2
import subprocess
import pickle
import glob
####################

# Two parameters for the program:
# 1. Align file (.bz2 format)
# 2. Out file



# Useful regular expressions

symbols = "UATGCKMRYSWBVHDXN"
r_coords = re.compile("(\d+)\s([" + symbols + symbols.lower() + "\-]+)\s(\d+)")

name = re.compile("^C?\s*(\S*(?:[^\s\d]|(?:\d\s)))\s*\d")  # REMEMBER: This can add an extra space to the name

float_re = "(\d+\.\d+)"
top_line_re = re.compile("(\d+)\s+" + float_re + "\s+" + float_re + "\s+" +
                         float_re + "\s+(chr(\d+))\s+(\d+)\s+(\d+)\s+\((\d+)\)\s+(C\s+|\+\s+|\s*)" #(C?|\+?)\s+"
                         + "(\([A-Z]+\)n|[A-Z]+_rich|[A-Z|a-z|0-9|-]+|[A-Z|a-z|0-9]+_|[A-Z|a-z|0-9|-]+_[A-Z|a-z|0-9]+)"
                         + "#(\w+/\w+_type|Low_complexity|Simple_repeat|\w+/\w+)\s+" 
                         + "\(?(\d+)\)?" + "\s+(\d+)\s+\(?(\d+)\)?\s+(\d+)\s*")
complement_re = re.compile("^C?\s+(\(|[A-Z]+)")
first_lineAL = re.compile("^\s+[chr]")
comp_line = re.compile("^C?\s+(\(|[A-Z]+)")
rm_prnth = re.compile("^\(")
out_header = re.compile("^[A-Z|a-z]")
rm_endsName = re.compile("(_3end|_5end)$")

class AlignFile:

    """ This class takes the file pointer from a standard .align file and 
    retrieves all of the inforrmation from each data entry. Each entry usually 
    begins with a whitespace separated line containing statistical information
    about each repeat and its family along with its class name. Several lines will then 
    follow that include the repeat and its complement DNA. The last two 
    lines of each entry will end with data about the transition/tranversion
    ratio and the initial gap rate along with the average gap size, 
    respectfully."""
    
    def __init__(self, fp):

        """ This method initializes all of the statistical and string data of 
        each .align entry. This includes the Smith-Waterman score, query 
        sequence name, beginning position in the query sequence, ending 
        position in the query sequence, the 'left' position in the query (length), the 
        determination of the query's complement (either on the positive or 
        complement string), the repeat's class and family, the repeat's 
        complement 'left' (length), beginning, and end positions, the percent deletion, 
        insertion, and division to the hundredths position, both the query and 
        its repeat (transposable element) sequence, the transition/transversion ratio, initial gap 
        rate, average gap size, and another integer denoted 'info' that 
        currently has an unknown reason for being in the entry besides a 'LineageID'. It then calls 
        the readAlignFile on the file pointer that has been passed to this
        class as a parameter. """

        self.SWscore = 0 #Smith-Waterman score
        self.querySeq = "" #Query sequence name
        self.posQueryBeg = 0 #Beginning position in modern query sequence
        self.posQueryEnd = 0 #Ending position in modern query sequence
        self.leftQuery = 0 #Number of bases in query sequence past the end position of the matched repeat
        self.posOrC = "" #Match aligns with either the complement (C) or positive (+) of query sequence
        self.matchRepeat = "" #Repeat name
        self.repClassFam = "" #Repeat family/class
        self.leftR = 0 #Remaining bases in (complement of) repeated ancestor sequence prior to the beginning of the match
        self.percDiv = 0.0 #Percent of substitutions in the matching region
        self.percDel = 0.0 #Percent deletions in query sequence according to the repeat sequence
        self.percIns = 0.0 #Percent insertions in query sequence according to the repeat sequence       
        self.posRepBeg = 0 #Starting position in the repeat sequence
        self.posRepEnd = 0 #Ending position in the repeat sequence      
        self.info = 0 #Lineage identifier
        self.seq = "" #Modern sequence
        self.CorPseq = "" #Ancestral sequence
        self.transition_tranversion = 0.0 #Ratio of transitions to transversions
        self.Gap_init_rate = 0.0 #Rate for gaps occurring
        self.avg_gap_size = 0.0 #Average gap size in match
        
        self.readAlignFile(fp)


    def readAlignFile(self, fp):

        """This method is used for reading through and parsing the data in the 
        .fa.align file. It groups the first line into a regular expression (if 
        it exists for the entry) and retrieves the percent division, percent 
        deletion, and percent insertion,(percDiv, percDel, and percIns, 
        respectively, all floats). It creates the left-closed, right open 
        interval based on if there is a 'C' or a '+'/empty string in the file for
        the data found in the entry (posRepBeg and posRepEnd, both ints). It also takes
        in the number that is last within the first line of the repeat data, but it is 
        unclear of what that number specifies (info, int). It also creates 
        and sets the string and values for the query sequence, its complement repeat 
        (seq and CorPseq, str), and the values for the rate of transition versus
        transversion (transition_transversion, float), the initial gap rate 
        (Gap_init_rate, float), and the average gap size (avg_gap_size, float), 
        and then breaks to finish the Align object."""

        chrAL = ""
        CchrAL = ""     
        firstLineFlag = 0
        atypStartCountQ = 0
        atypStartCountC = 0
        
        for line in fp:
            if sys.version_info[0] == 3 and type(line) == bytes:
                line = line.decode()

            #parses the data for the top line and sets the needed information
            if re.search(top_line_re, line):
                startLine = re.search(top_line_re, line)
                
                self.SWscore = int(startLine.group(1))
                self.percDiv = float(startLine.group(2))
                self.percDel = float(startLine.group(3))
                self.percIns = float(startLine.group(4))
                
                self.querySeq = startLine.group(5)
                self.posQueryBeg = int(startLine.group(7))-1
                self.posQueryEnd = int(startLine.group(8))
                self.leftQuery = int(startLine.group(9))
                
                self.matchRepeat = startLine.group(11)
                if re.search(rm_endsName, self.matchRepeat):
                    self.matchRepeat = self.matchRepeat[:-5]
                self.repClassFam = startLine.group(12)
                
                
                #if "+"
                if startLine.group(10) == "" or startLine.group(10) == "+":
                    self.posOrC = "+"
                    self.leftR = startLine.group(15)
                    self.posRepBeg = int(startLine.group(13))-1
                    self.posRepEnd = int(startLine.group(14))
                # if "C"
                else:
                    self.posOrC = startLine.group(10)
                    self.leftR = startLine.group(13)
                    self.posRepBeg = int(startLine.group(15))-1
                    self.posRepEnd = int(startLine.group(14))
                    
                #unsure of what this number correlates to besides LineageID
                self.info = startLine.group(16)
                
                firstLineFlag = 1
        
            #looks for the query sequence if the typical first line is missing
            elif re.match(first_lineAL,line) and firstLineFlag == 0:

                if r_coords.search(line):
                    r = r_coords.search(line)

                if name.search(line):
                    SeqName = name.search(line)
                                                
                chrAL = chrAL + r.group(2)
                
                if atypStartCountQ == 0:
                    self.querySeq = SeqName.group(1).rstrip()
                    self.posQueryBeg = int(r.group(1))-1
                    self.posQueryEnd = int(r.group(3))
                    atypStartCountQ = atypStartCountQ + 1
                else:
                    self.posQueryEnd = int(r.group(3))
                        
        
            #looks for the complement string if the typical first line is missing
            elif complement_re.search(line) and firstLineFlag == 0:
                group2 = line.split()
                if group2[0] == "C":
                    self.posOrC = group2[0]
                    
                    if r_coords.search(line):
                        r = r_coords.search(line)
                    
                    if atypStartCountC == 0 and name.search(line):
                         
                        SeqInfo = name.search(line)
                        getNames = SeqInfo.group(1).rstrip()
                        
                        if not getNames.find("#") == -1 and not getNames.find("#") == len(SeqInfo.group(1).rstrip())-1:
                            
                            # this is if the first "#" splits the two names, and there are more than one "#" in the name
                            twoNames  = getNames.split("#",1)
                            self.matchRepeat = twoNames[0]
                            self.repClassFam = twoNames[1]
                            
                        elif getNames.find("#") == -1 and not getNames.find("-") == -1:
                            
                            # this is if the first "-" splits the two names, and there are more than one "-"
                            twoNames  = getNames.split("-",1)
                            self.matchRepeat = twoNames[0]
                            self.repClassFam = getNames                                             
                            
                        elif not getNames.find("#") == -1 and getNames.find("#") == len(SeqInfo.group(1).rstrip())-1:
                            
                            # this is if the last character is a "#", and can be split like the one above
                            newNames = getNames.rstrip("#")
                            twoNames  = newNames.split("-",1)
                            self.matchRepeat = twoNames[0]
                            self.repClassFam = newNames                                             
                            
                        else:
                            
                            #this is a catch all for my data, in case the name has been cut off in another odd format
                            self.matchRepeat = getNames
                            self.repClassFam = "Null"
                        
                        if re.search(rm_endsName, self.matchRepeat):
                            self.matchRepeat = self.matchRepeat[:-5]
                        self.posRepEnd = int(r.group(1))-1
                        self.posRepBeg = int(r.group(3))
                        atypStartCountC = atypStartCountC + 1
                        
                    else:
                        self.posRepBeg = int(r.group(3))
                        
                    CchrAL = CchrAL + r.group(2)
                    
                    
                else:
                    #same scenario that first line is missing, but for "+" instance
                    self.posOrC = "+"
                    
                    if r_coords.search(line):
                        r = r_coords.search(line)   
                        
                    if atypStartCountC == 0 and name.search(line):
                        
                        SeqInfo = name.search(line)
                        getNames = SeqInfo.group(1).rstrip()                                            
                        
                        if not getNames.find("#") == -1 and not getNames.find("#") == len(SeqInfo.group(1).rstrip())-1:
                                                                            
                        # this is if the first "#" splits the two names, and there are more than one "#"
                            twoNames  = getNames.split("#",1)
                            self.matchRepeat = twoNames[0]
                            self.repClassFam = twoNames[1]                                      
                        
                        elif getNames.find("#") == -1 and not getNames.find("-") == -1:
                            
                            # this is if the first "-" splits the two names, and there are more than one "-"
                            twoNames  = getNames.split("-",1)
                            self.matchRepeat = twoNames[0]
                            self.repClassFam = getNames                                                 
                            
                        elif not getNames.find("#") == -1 and getNames.find("#") == len(SeqInfo.group(1).rstrip())-1:
                            
                            # this is if the last character is a "#", and can be split like the one above
                            newNames = getNames.rstrip("#")
                            twoNames  = newNames.split("-",1)
                            self.matchRepeat = twoNames[0]
                            self.repClassFam = newNames                                                 
                            
                        else:
                            
                            #this is a catch all for my data, in case the name has been cut off in another odd format
                            self.matchRepeat = getNames
                            self.repClassFam = "Null"                                           
                        
                        if re.search(rm_endsName, self.matchRepeat):
                            self.matchRepeat = self.matchRepeat[:-5]
                        self.posRepBeg = int(r.group(1))-1
                        self.posRepEnd = int(r.group(3))
                        atypStartCountC = atypStartCountC + 1
                    else:   
                        self.posRepEnd = int(r.group(3))
                
                    CchrAL = CchrAL + r.group(2)            
                    
            #looks for the query sequence if the first line is a typical entry
            elif re.match(first_lineAL,line) and firstLineFlag == 1:
                if r_coords.search(line):
                    r = r_coords.search(line)   
                    chrAL = chrAL + r.group(2)
                else:
                    print("ERROR IN LINE FORMAT")
        
        
            #looks for the complement string if the first line is a typical entry
            elif re.search(comp_line,line) and firstLineFlag == 1:
                if r_coords.search(line):
                    r = r_coords.search(line)
                    CchrAL = CchrAL + r.group(2)                        
                else:
                    print("ERROR IN LINE FORMAT")
                
            #looks for the transition/transversion ratio
            elif line.find("Transitions") != -1:
                trans_info = line.split()
                self.transition_transversion = float(trans_info[4])
        
            #sets the values for the two repeat strings and the gap data
            elif line.find("Gap") != -1:
                gap_info = line.split()
                self.seq = chrAL
                self.CorPseq = CchrAL
        
                self.Gap_init_rate = float(gap_info[3])
                self.avg_gap_size = float(gap_info[11])         
                
                break    

    def __str__(self):

        """This method creates a string that will be used as a unique identifier
        when creating the Repeat objects later in this library. It creates a string
        that holds the beginning and end of the query sequences (ints) and the 
        the repeat's class name and family (strings), all separated by the '#'
        character. """

        if self.repClassFam != "":
            string = str(self.posQueryBeg) + "#" + str(self.posQueryEnd) + "#" + \
                    self.matchRepeat + "#" + self.repClassFam
        else:
            string = str(self.posQueryBeg) + "#" + str(self.posQueryEnd) + "#" + \
                    self.matchRepeat + "#"      
        return string

class OutFile:

    """This class creates an object specifically for the .out file. Because each
    line of the .out file is an entry, the class takes in one line from the 
    file as it is passed in and creates one object per each line (after the 
    header lines of the file have been removed)."""

    def __init__(self, outStr):

        """The __init__ method constructs the .out entry by calling 
        readOutFile() on the outStr (string) that has been passed in as a 
        parameter. """  

        self.readOutFile(outStr)

    def readOutFile(self, outStr):
        
        """This method is used for parsing the .fa.out string. It splits
        the string passed in to an array, and assigns values to the 
        Smith-Waterman value (SWscore, int), percent division (percDiv, float),
        percent deletion (percDel, float), percent insertion (percIns, float),
        the query sequence name (querySeq, str), where in the query the repeat 
        begins, ends, and 'left' (posQueryBeg, posQueryEnd, and leftQuery, all 
        ints respectively), whether the repeat lies on the complement or the 
        positive (posOrC, str), the matching repeat name (matchRepeat, str), its 
        class/family name (repClassFam, str), position in repeat 'left' (leftR, 
        int), the beginning and end position in the repeat's complement 
        (posRepBeg and posRepBeg, ints), and the ID number in the file (ID, 
        int) that is an identifier for groupings with locations on sequence/repeat
        class and family. This method also creates the left-open, right-closed interval 
        based on the existence of parenthesis in the repeat's complement 
        beginning and end positions.
        
        Calls the castToInt method to remove the parenthesis and cast the number 
        to an int."""

        lineArr = outStr.split()
        
        # Set SWscore
        self.SWscore = int(lineArr[0])
        
        # Set percDiv
        self.percDiv = float(lineArr[1])
        
        # Set percDel
        self.percDel = float(lineArr[2])

        # Set percIns
        self.percIns = float(lineArr[3])        
        
        # Set querySeq
        self.querySeq = lineArr[4]              
        
        # Set posQueryBeg
        self.posQueryBeg = int(lineArr[5])-1
        
        # Set posQueryEnd
        self.posQueryEnd = int(lineArr[6])

        # Set leftQuery
        self.leftQuery = self.castToInt(lineArr[7])
                
        # Set posOrC
        self.posOrC = lineArr[8]
        
        # Set matchRepeat
        self.matchRepeat = lineArr[9]
        
        # Set repClassFam
        self.repClassFam = lineArr[10]
        
        if lineArr[11].find("(") == -1:
            # Set posRepBeg
            self.posRepBeg = int(lineArr[11])-1
                    
            # Set posRepEnd
            self.posRepEnd = int(lineArr[12])
            
            # Set leftR
            self.leftR = self.castToInt(lineArr[13])                        
                            
        else:
            # Set posRepBeg
            self.posRepBeg = int(lineArr[13])-1
    
            # Set posRepEnd
            self.posRepEnd = int(lineArr[12])
            
            # Set leftR
            self.leftR = self.castToInt(lineArr[11])                        
                                                                                    
        # Set ID - identifies groups with similar locations on sequence and repeat family
        self.ID = int(lineArr[14])
    

    def castToInt(self,data):
        
        """This method takes the 'data' and returns an integer that removes the
        parathesis surrounding it and casts it to the int type."""
        
        if re.search(rm_prnth,data):
            intData = int(data[1:len(data)-1])
            return int(intData)
        else:
            return int(data)    

    def __str__(self):
        
        """This method creates a string that will be used as a unique identifier
        when creating the Repeat objects later in this tool and overrides the 
        str method. It creates a string that holds the beginning and end of the 
        query sequences (ints) and the the repeat's class name and family 
        (strings), all separated by the '#' character. """
        
        string = str(self.posQueryBeg) + "#" + str(self.posQueryEnd) + "#" + \
                        self.matchRepeat + "#" + self.repClassFam
        return string

class Repeat:

    """This object takes an AlignFile and an OutFile object and takes data from 
    each class in order to retain the most accurate data of the two objects, 
    along with combining the two to get a complete data set for the repeat
    sequence. """  
    
    def __init__(self, alignObj = None, outObj = None):

        """This is the __init__ takes in two arguements - the object created 
        from the .align file and the object from the .out file. If the Align object
        is not given, then it sets all attributes to None. Otherwise, it calls the 
        method setRepObj to begin setting the data found from both files in 
        order to create the complete Repeat data object. """    
        
        if(alignObj == None):
            self.SWscore = None #Smith-Waterman score
            self.querySeq = None #Query sequence name
            self.posQueryBeg = None #Beginning position in modern query sequence
            self.posQueryEnd = None #Ending position in modern query sequence
            self.leftQuery = None #Number of bases in query sequence past the end position of the matched repeat
            self.posOrC = None #Match aligns with either the complement (C) or positive (+) of query sequence
            self.matchRepeat = None #Repeat name
            self.repClassFam = None #Repeat family/class
            self.leftR = None #Remaining bases in (complement of) repeated ancestor sequence prior to the beginning of the match
            self.ID = None #Identifies groups with similar locations on sequence and repeat family
            self.percDiv = None #Percent of substitutions in the matching region
            self.percDel = None #Percent deletions in query sequence according to the repeat sequence
            self.percIns = None #Percent insertions in query sequence according to the repeat sequence  
            self.posRepBeg = None #Starting position in the repeat sequence
            self.posRepEnd = None #Ending position in the repeat sequence
            self.info = None #LineageID
            self.seq = None #Modern sequence
            self.CorPseq = None #Ancestral sequence
            self.transition_tranversion = None #Ratio of transitions to transversions
            self.Gap_init_rate = None #Rate for gaps occurring
            self.avg_gap_size = None #Average gap size in match             
        else:
            self.setRepObj(alignObj, outObj)
     # Inspectors / Mutators
    def ancestor_seq(self, seq = None):
        """Set/Return aligned ancestor sequence."""
        if seq != None:
            self.CorPseq = seq
        return self.CorPseq

    def ancestor_coords(self, begin = None, end = None):
        """Set/Return alignment corrds w.r.t. sequences."""
        if begin != None:
            self.posRepBeg, self.posRepEnd = begin, end
        return (self.posRepBeg, self.posRepEnd)

    def modern_seq(self, seq = None):
        """Set/Return aligned modern sequences."""
        if seq != None:
            self.seq = seq
        return self.seq

    def modern_chr(self, seq_name = None):
        """Set/Return name of chromsome / sequence contiaining modern instance"""
        if seq_name != None:
            self.querySeq = seq_name
        return self.querySeq

    def modern_coords(self, begin = None, end = None):
        """Set/Return coordinates of modern sequence on genome."""
        if begin != None:
            self.posQueryBeg, self.posQueryEnd = begin,end
        return (self.posQueryBeg, self.posQueryEnd)

    def start(self):
        """Return the start position on the modern strand"""
        return self.posQueryBeg

    def finish(self):
        """Return the finish position on the modern strand"""
        return self.posQueryEnd

    def modern_strand(self, strand = None):
        """Return strand."""
        if strand != None:
            self.posOrC = strand
        return self.posOrC

    def rep_name(self, rep_name = None):
        if rep_name != None:
            self.repClassFam = rep_name
        return self.repClassFam
    
    def class_name(self, class_name = None):
        if class_name != None:
            self.matchRepeat = class_name
        return self.matchRepeat


    # Helpers   
    def setRepObj(self, alignObj, outObj):

        """Sets all of the attributes for the Repeat object from the combined
        .align and .out objects. From the .out file, this includes: the Smith-
        Waterman score (SWscore, int), the query sequence name (querySeq, str),
        where in the query the repeat begins, ends, and 'left' (posQueryBeg, 
        posQueryEnd, and leftQuery, all ints respectively), whether the repeat 
        lies on the complement or the positive (posOrC, str), the matching 
        repeat name (matchRepeat, str), its class/family name (repClassFam, 
        str), position in repeat 'left' (leftR, int), and the ID number in the 
        .out file (ID, int). As for the .align file, it retrieves and sets the 
        following: it retrieves the percent division, percent deletion, and 
        percent insertion,(percDiv, percDel, and percIns, respectively, all 
        floats), the repeat's complement beginning and end positions (posRepBeg 
        and posRepEnd, ints), the number that is placed last within the first 
        line of the repeat data, but it is unclear of what that number 
        represents (info, int), the repeat's sequence and its complement (seq 
        and CorPseq, str), the value for the rate of transition versus 
        transversion (transition_transversion, float), the initial gap rate (
        Gap_init_rate, float), and, finally, the average gap size (avg_gap_size,
        float) for the Repeat object."""
        
        #from .out file
        self.SWscore = outObj.SWscore
        self.querySeq = outObj.querySeq
        #self.posQueryBeg = outObj.posQueryBeg
        #self.posQueryEnd = outObj.posQueryEnd
        #self.leftQuery = outObj.leftQuery
        self.posOrC = outObj.posOrC
        self.matchRepeat = outObj.matchRepeat
        self.repClassFam = outObj.repClassFam
        self.leftR = outObj.leftR
        self.ID = outObj.ID
        
        #from .align file
        self.percDiv = alignObj.percDiv
        self.percDel = alignObj.percDel
        self.percIns = alignObj.percIns
        self.posQueryBeg = alignObj.posQueryBeg
        self.posQueryEnd = alignObj.posQueryEnd
        self.leftQuery = alignObj.leftQuery
        self.posRepBeg = alignObj.posRepBeg
        self.posRepEnd = alignObj.posRepEnd
        self.info = alignObj.info
        self.seq = alignObj.seq
        self.CorPseq = alignObj.CorPseq
        self.transition_tranversion = alignObj.transition_tranversion
        self.Gap_init_rate = alignObj.Gap_init_rate
        self.avg_gap_size = alignObj.avg_gap_size
    
    def __str__(self):

        """This method overloads the str operator and prints out a string that 
        is human readable about the data for each Repeat object. It prints out
        the query sequence name, beginning and end positions for the query 
        sequence, where the repeat sequence is located (positive or complement
        strand), the name and its family/class, and its beginning and ending
        positions. """
        
        classString = "Information on Repeat object:\n Query Sequence: %s\n " + \
                "Beginning Position for Query Sequence: %i\n End Position for Query Sequence: %i\n " + \
                "Repeat Sequence complement located on: %s\n Repeat name: %s\n " + \
                "Repeat Class/Family: %s\n Beginning Position for Repeat Sequence: " + \
                "%i\n End Position for Repeat Sequence: %i\n" 
        return classString % (self.querySeq, self.posQueryBeg, self.posQueryEnd,
                                                  self.posOrC, self.matchRepeat, self.repClassFam, 
                                                  self.posRepBeg, self.posRepEnd)



def OutFileGenerator(outFile, excludeSet = {}):
    fp = open(outFile)
    
    # skip over header
    fp.readline()
    fp.readline()

    for line in fp:
        if line.rstrip():
            O = OutFile(line)
            if not excludeSet or O.repClassFam not in excludeSet:
                yield O

def openBZ2(file):
    if sys.version_info[0] == 2:
        return bz2.BZ2File(file)
    else:
        return bz2.open(file)

def AlignFileGenerator(alignFile, excludeSet = {}):
    fp = openBZ2(alignFile) if alignFile.endswith(".bz2") else open(alignFile)

    while 1:
        A = AlignFile(fp)
        if str(A) == "0#0##":
            break
        if not excludeSet or A.repClassFam not in excludeSet:
            yield A


def createObjArray(alignFilePointer, outFilePtr1):
   
    """Takes the .align and the .out file pointers and creates the designated 
    objects for each file. Removes the two header lines from the .out file 
    before calling the OutFile class on the .out line. Returns the OutFile and 
    AlignFile object arrays. """
    
    alignOut = []
    
    #arrays created to hold all .align and .out objects 
    alignObjs = []
    outObjs = []    
     
    #reads through .align file to create each Align object
    while 1:
        A = AlignFile(alignFilePointer)
        if str(A) == "0#0##":
            break
        else:
            alignObjs.append(A)
            
    #skips over header of .out object
    outFilePtrF = outFilePtr1.readline().strip()
    
    #skips over any whitespace lines found in .out file
    while re.match(out_header,outFilePtrF):
        outFilePtrF = outFilePtr1.readline()

    #reads through .out file to create each Out object
    for outFilePtrF in outFilePtr1:
        B = OutFile(outFilePtrF)
        outObjs.append(B)

    alignOut.append(outObjs)
    alignOut.append(alignObjs)
    
    return alignOut

def createRepObjsOld(outObjs, alignObjs):
    
    """This is a standalone method that takes in two arrays - one of .out 
    objects and one of align objects. It then takes both of those arrays and 
    creates a Repeat class object for each of the matching strings found in the 
    two arrays. Takes the unused entries in the .align file and creates a 
    separate align-object only array. Will return both arrays."""  
    
    repObjs = []
    unusedAlignObjs = []
    repDict = {}

    #using dictionary to create Repeat objects
    
    #creates Repeat dictionary based on out objects (will always have fewer entries
    #than .align file)
    for i in range(len(outObjs)):
        repDict[str(outObjs[i])] = outObjs[i]
        
    #creates Repeat objects based on the matching strings of the .align and .out
    #objects
    for j in range(len(alignObjs)):
#        if repDict.has_key(str(alignObjs[j])):
        if str(alignObjs[j]) in repDict:
            C = Repeat(alignObjs[j], repDict[str(alignObjs[j])])
            repObjs.append(C)    
        else:
            unusedAlignObjs.append(alignObjs[j])

    #array that returns an array of created Repeat objects
    #and unused Align objects
    RepAndUnalignObjs = []
    RepAndUnalignObjs.append(repObjs)
    RepAndUnalignObjs.append(unusedAlignObjs)
    
    return RepAndUnalignObjs

def createRepObjs(alignFile, outFile, excludeFile = None):
    excludeSet = {line.rstrip() for line in open(excludeFile)} if excludeFile else {}
    alignObjs = list(AlignFileGenerator(alignFile, excludeSet))
    outGen = OutFileGenerator(outFile, excludeSet)

    repObjs = []
    unusedObjs = []
    
    j = 0
    for o in outGen:
        while j < len(alignObjs) and alignObjs[j].posQueryBeg < o.posQueryBeg: j+=1
        k = j
        while k < len(alignObjs) and alignObjs[k].posQueryEnd <= o.posQueryEnd: k+=1
        if k > j:
            A = max([(r.SWscore,j-k-i,r) for i,r in enumerate(alignObjs[j:k])])[-1]  # Use the longest repeat
            repObjs.append(Repeat(A, o))
            j = k
        if j == len(alignObjs):
            break

    return repObjs
        

def pickleRepeats(repeatList, fileName):
    """
    Creates a pickled file of the repeat list.
    Will compress the file into .bz2 output if:
    1) compressed - True
    2) fileName has a .bz2" at the end.
    """
    fp = open(fileName, "wb")
    for obj in repeatList:
        pickle.dump(obj, fp)
    fp.close()


def create_prmFile(alignFile, outFile, prmFile, excludeFile = None):
    pickleRepeats(createRepObjs(alignFile, outFile, excludeFile), prmFile) 
    
    
def prm_generator(chr, dir = "."):
    fileName = "%s/%s.%d.prm" % (dir, chr, sys.version_info[0])
    fp = open(fileName, 'rb')

    while True:
        try:
            yield pickle.load(fp)
        except EOFError:
            break

def unpickleRepeats(chr, dir = "."):
    """Will unpickle the .bz2 version if:
    1) The name ends with .bz2.
    2) No non-compressed version is present.
    """
    return list(prm_generator(chr, dir))
                

if __name__ == "__main__":
    excludeFile = "hg18.rpt_exclude.txt"
    for file in glob.glob("chr*.fa.out"):
        if re.match("chr\d+\.fa*", file):
            file = file.rstrip(".fa.out")
            print("Converting: %s" % (file))
            create_prmFile(file + ".fa.align.bz2", file + ".fa.out", file + ".%d.prm" % (sys.version_info[0]), excludeFile = excludeFile)


