import pickle 
from RMfileReader import *
from numpy import *
from read_files import *

"""makes c++ global count file - must include chromosomes to include in global count,
what families are to be included, the directory with the prmfiles to read from, and an
optional output file"""
def global_counts(chromosomes,includedfams,directory,fileToUse = None):
    Cglobal = {}
    bases = {'A':0,'C':1,'G':2,'T':3}
    for c in chromosomes:
        print(c)
        if(directory=="Karro"):
            List = unpickleRepeats(str(c),"/shared/karrojegrp/cache/human/hg18/seq/rmsk")
        else:
            List = unpickleRepeats(str(c),str(directory))
        for r in List:
            ancestor = r.ancestor_seq()
            modern = r.modern_seq()
            special = r.rep_name()
            family = r.class_name()
            if family in includedfams:
                M = zeros((4,4))
                for c1,c2 in zip(ancestor,modern):
                    if(c1 in bases and c2 in bases):
                        M[bases[c1],bases[c2]] += 1
                if(M.sum() >=40):
                    if family not in Cglobal: 
                        Cglobal[family] = zeros((4,4))
                    Cglobal[family] += M
    if(fileToUse == None):
        fileToUse = open("globalFamilies.ccm","w")
    fileToUse.write("OneSection ")
    for entry in Cglobal:
        fileToUse.write(str(entry)+" ")
        for i in range(4):
            for j in range(4):
                fileToUse.write(str(Cglobal[entry][i,j])+" ")
    fileToUse.write("EndSection")

"""Makes a c++ count file based on a whole chromosome"""
def globalChromosome(includedfams,chromosome,directory):
    fileToUse = open(str(chromosome[0]) + "Whole.ccm","w")
    global_counts(chromosome,includedfams,directory,fileToUse)

"""Makes a c++ count file for a given chromosome from Martin data in a given directory"""
def makeMartinFileReadable(includedfams,chromosome,directory):
    Cgamma = {}
    for partNo,start,finish,R in read_mat(str(directory)+"hg18_gene.mat",target_chr = str(chromosome)):
        family = R.class_name
        access = str(start) + '-' + str(finish)
        if(family in includedfams and R.M.sum() >= 40):
            if Cgamma.get(access) == None:
                Cgamma[access] = {}
            if Cgamma[access].get(family) == None:
                Cgamma[access][family] = matrix(R.M)
            else:
                Cgamma[access][family] += matrix(R.M)
    #write to file
    file = open(str(chromosome)+"Martin.ccm","w")
    for section in Cgamma:
        file.write(str(section)+" ")
        for family in Cgamma[section]:
            file.write(str(family)+" ")
            for i in range(4):
                for j in range(4):
                    file.write(str(Cgamma[section][family][i,j])+" ")
        file.write("EndSection"+" ")


"""Makes a c++ local count file for a given chromosome with given partitions from a prmfile in a given directory"""        
def count_by_chrm(chromosome,includedfams,partitions,directory):
    bases = {'A':0,'C':1,'G':2,'T':3}
    tracker = 0
    if(directory=="Karro"):
        List = unpickleRepeats(str(chromosome),"/shared/karrojegrp/cache/human/hg18/seq/rmsk")
    else:
        List = unpickleRepeats(str(chromosome),str(directory))
    Cgamma = {}
    for r in List:
        partNotFound = True
        strand = r.modern_strand()
        ancestor = r.ancestor_seq()
        modern = r.modern_seq()
        start,finish = r.modern_coords()
        special = r.rep_name()
        family = r.class_name()
        if(family in includedfams and tracker < len(partitions)):
            M = zeros((4,4)) #Create count matrix
            for c1,c2 in zip(ancestor,modern): #Fill count matrix:
                if(c1 in bases and c2 in bases):
                   M[bases[c1],bases[c2]] += 1 
            #add to a partition
            while(partNotFound and tracker < len(partitions)):
                if(int(start) >= int(partitions[tracker][0]) and int(finish) <= int(partitions[tracker][1])):
                    access = str(chromosome) + ':' + str(partitions[tracker][0]) + '-' + str(partitions[tracker][1]).strip()
                    if(M.sum() >= 40):  #If instance has more than 40 bases add to overall count of family
                        if access not in Cgamma:
                            Cgamma[access] = {}
                        if family not in Cgamma[access]:
                            Cgamma[access][family] = M
                        else:
                            Cgamma[access][family] += M
                    partNotFound = False
                elif(int(start) >= int(partitions[tracker][1])):
                    tracker+=1 
                else:
                    partNotFound = False
    #write file for c++
    file = open(str(chromosome)+".ccm","w")
    for section in Cgamma:
        file.write(str(section)+" ")
        for family in Cgamma[section]:
            file.write(str(family)+" ")
            for i in range(4):
                for j in range(4):
                    file.write(str(Cgamma[section][family][i,j])+" ")
        file.write("EndSection"+" ")
                        
                   
                        
"""creates the partitions for the local count files based on input partition file"""                        
def partition_setup(part_file):
    parts = []
    f = open(part_file,"r")
    for line in f:
        line = line.split("\t")
        parts.append((line[1],line[2]))
    return parts           

"""main function to create files based on input commands, see help"""
def main(todo,directory,chromosome=None):
    chromosomes = ["chr22","chr21","chr20","chr19","chr18","chr17","chr16","chr15","chr14","chr13","chr12","chr11","chr10","chr9","chr8","chr7","chr6","chr5","chr4","chr3","chr2","chr1"]
    martinFamilies = []
    f1 = open("martin_repeats.txt","r")
    for line in f1:
        martinFamilies.append(line.strip())
    f1.close()
    if(todo=="local"):
        if(chromosome != None):
            part = partition_setup(str(chromosome)+".txt")
            count_by_chrm(chromosome,martinFamilies,part,directory)
        else:
            for c in chromosomes:
                print(c)
                part = partition_setup(c+".txt")
                count_by_chrm(c,martinFamilies,part,directory)
    elif(todo=="globalChrom" and chromosome != None):
        chromosomes = [str(chromosome)]
        globalChromosome(martinFamilies,chromosomes,directory)
    elif(todo=="globalChrom" and chromosome == None):
        for c in chromosomes:
            globalChromosome(martinFamilies,[c],directory)
    elif(todo=="global"):
        global_counts(chromosomes,martinFamilies,directory)
    elif(todo=="Martin" and chromosome != None):
        makeMartinFileReadable(martinFamilies,chromosome,directory)
    else:
        print("Incorrect parameters, use help for more information.")

if __name__ == "__main__":
    if(len(sys.argv) == 3):
        main(sys.argv[1],sys.argv[2])
    elif(len(sys.argv)==4):
        main(sys.argv[1],sys.argv[2],sys.argv[3])
    elif(len(sys.argv)==2 and sys.argv[1] == "help"):
        print("This program takes 2-3 parameters to create readable \n" +
        "count files for the c++ program to find asymmetric regions. \n"+
        "Assumes prm files and partition files are ordered by coordinates. \n" +
        "Parameter 1 - Please supply a task: \n"+
        "\"global\" - creates global family count file using all 22 chromosomes \n"+
        "\"globalChrom\" - creates overall chromosome count file \n" +
        "\"local\" - creates local count files \n"+
        "\"Martin\" - creates readable files from Martin data (must supply chromosome number).  Note that this has multiple python file dependencies.\n\n"+
        "Parameter 2 -  Supply a directory to find prm files in: \n"+
        "\"Karro\" - uses prm files in this directory \n"+
        "~\ - or specify other directory. \n"+
        "Note: there needs to be the correct repeat masker files from RMFileReader (prms) or hg18_gene.mat file in the given directory.\n\n" +
        "Parameter 3 (optional except for the Martin task) - Specify what chromosome (chr) you want to create \n"+ 
        "a file for if using the \"globalChrom\", \"local\" or \"Martin\" tasks -  \n"+
        "note, if not specified, files will be created for all chromosomes. \n\n"+
        "Also note, this program uses partition files for the local option - \n" +
        "these files simply have to be in the directory the program is running from, and are of the title chromosome.txt \n"+
        "(Example, \"chr22.txt\" ) which contains a chromosome number, the beginning, and end position of the section separated by tabs on each line \n" +
        "(Example, chr22\t1234\t4567). \n\n")
    else:
        print("Incorrect parameters, please supply at least a task and a directory, and (optionally for local), a chromosome --- or use parameter \"help\"")
