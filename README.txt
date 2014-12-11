DNAStrandSymmetry
=================

Masters Project for finding the type of symmetry of a given region, and for finding regions that may switch their symmetry status.
Authors: Michelle Flowers, and Dr. John Karro

================
Firstly, our current algorithm for finding asymmetric regions:
1. Make a “local count matrix” of each family in each partition leaving out instances with less than 40 bases – use makeReadableFile.py for this, symmetryTester.cpp will do the rest.
2.  For each partition, take each count matrix and create the asymmetric model P matrix (sum rows and divide each index by sum of their row).
3.  Calculate “local ds” for each family in each partition by the logdet method.
4.  Take the log of P, aka find the “Rt matrix” of each family in each section.
5.  Calculate the “q matrix” for each family in each section - use the local ds for this (Rt/d)
6.  Now for all the families in one section, calculate the “average Q matrix”.
7.  With this average Q matrix, and global distances (“dg”) of the families (found using the global count matrices and steps 2 and 3), create a “Phat matrix” for each family - that is, expm(Q*dg) for each family in the section.
8.  Find the log-likelihood, “logL” - multiply the local count matrix of the family with the log of the pHat matrix (element wise log!!) found for the family to get a new matrix.  Now sum each element of the new matrix to get a number per family.  For the section, the log-likelihood is the sum of all of these numbers.
9.  Calculate the BIC - which is, -2*(logL for section) + (f+k)*log(16*k) where k is the number of families used, and f is 5 or 11 depending on your model.
10. Repeat with symmetric model P – see code for symmetric model P calculation.
11.  Which ever has the lowest score is the preferred model.

================
Now, to get the files to use with the c++ code of this algorithm use the makeReadableFiles.py program (every file created has the extension “.ccm” - Note, the directory ccmFiles has the global counts "globalFamilies.ccm" file and the local counts per chromosome (partitioned by genes) .ccm files).  This program does a lot of preprocessing on the "prmfiles" made by Dr. John Karro and as such needs access to these files (which are too large to be provided here, hence the few given .ccm files - if others are needed, please let me know) and a few of his python scripts (RMfileReader.py, RptMatrix.py, read_file.py, and matrix_fun.py).  It also needs access to a list of acceptable families, which is currently hardcoded as the “martin_repeats.txt” file (so modify this file if you want to change allowable overall families, otherwise leave it be), and needs access to partition files for making local count files (more on this later).  All extras needed (except the prmfiles) are stored in the directory called "makeFilesNeededExtras". This program, in short, creates family count matrices (note, this merges the same families per section, and does not store the information by each repeat instance) which count how the bases in a given ancestor sequence of a family transition to the bases in its modern sequence. Here are the commands for the program (again, you must have the prmfiles to actually use these):

To get a global count file (i.e. combining all families in every chromosome) use the command:
python makeReadableFile.py global directory
Where directory can be “Karro” if using Miami University’s RedHawk cluster, or can be specified otherwise to where your prmfiles are located.

To get local count files (i.e. partition a chromosome into sections that have separate family counts based on a partition file - used to find symmetry status of those particular sections), use the command:
python makeReadableFile.py local directory chromosome
Where chromosome should be formatted as “chr#” (example, “chr22”).  It should be noted that this is optional, and if excluded this command will create local files for each chromosome.  Note that the partition files have to be in the directory that the program is running from, and are of the title chromosome.txt (Example, "chr22.txt”) which contains a chromosome number, the beginning, and end position of the section separated by tabs on each line (Example, chr22	1234	4567).

This program can also be used to generate global count files for specific chromosomes only using the “globalChrom” option, as well as used to generate readable Martin count matrices from the “hg18_gene.mat” file (also not included due to size) using the “Martin” option.  Use the “help” option to see more on how these are called:
python makeReadableFile.py help

=====================

symmetryTester.cpp – A program implementing the symmetry finding algorithm described above.
To use, compile the symmetryTester.cpp file (it is important to include the Armadillo package):
g++ -std=c++11 symmetryTester.cpp –o main –larmadillo

To just find asymmetric sections, use the commad:
./main chr#.ccm global.ccm
Where “chr#.ccm” is the file that contains the counts of the families in each section made with the program makeReadableFile.py as discussed above (note, you can combine all local files into a giant file for this to use – “cat chr*.ccm all.ccm” – or use global files here to check the entire genome/specific chromosomes for overall symmetry type).  And where “global.ccm” represents the global count file you want to use for the global distances discussed in the algorithm – global distances that are calculated in the program will be stored in the file “GlobalDs.txt”.  The program will print to the command line how many asymmetric sections were found.  The asymmetric sections that are found will be stored in the file “AsymmetricSections.txt”

To find sections that switch symmetry, use the command:
./main chr#.ccm global.ccm –l
This defaults the limits to use 30% of the youngest/newest families, and 30% of the oldest families, after removing 10% at each end (removal in case of outliers).  The limits are based on the ages/distances found using the global file.

To specify what percent of families to use and how many to discard as outliers, use the command:
./main chr#.ccm global.ccm -1 l1 l2
Where l1 is the percent of new and old families to use in decimal form (example: .4) and l2 is the percent to consider outliers (example: .05).

The program will print to the command line how many asymmetric sections were found (overall), how many were found using just new families, and how many were found using just old families.  Particular section names will be stored in “AsymmetricSections.txt”, “AsymmetricSectionsNew.txt”, and “AsymmetricSectionsOld.txt” respectively, with the age/distance limits listed first for the New and Old files.  Finally, the program will compare the new and old output.  It will print to the command line how many sections have switched, and how many in what direction (when looking at old to new).  The file “EasyReadCompare.txt” lists the sections that have changed symmetry in an easy to read format – it will list the new outcome, the old outcome, and then the overall outcome.  Note that each outcome will contain the number of bases used in the calculations, followed by the number of families used in the calculations, and then the symmetry result.  “LineCompareFile.txt” has the same information but contains it all on a single line per section, starting with section name, then the new information, the old information next, and the overall information last.  The “AllFile.txt” contains all sections in the same format as “LineCompareFile.txt”, regardless if that section has switched symmetry.
