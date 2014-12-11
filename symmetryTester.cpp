#include <iostream>
#include <stdio.h>
#include <armadillo>
#include <map>
#include <math.h>
#include <fstream>
#include <string>

//Helpful typedefs for main Map types
typedef std::map<std::string,arma::mat> MAP;
typedef std::map<std::string,float> DMAP;


/**
Read in the next section to work with or read in global counts
 - the file being read has the counts of each family (values stored by rows) by section, 
ending the section with "EndSection" - global files only have one section
@param myfile The file to be read that has the counts for each family in each section
-various files can be created from makeReadableFiles.py.
@return countDictionary A MAP mapping a family to a count matrix for a 
particular section - note the section name is not stored in this MAP, please read it 
in prior to calling function!
*/
MAP dealWithPartitions(std::ifstream& myfile){
  MAP countDictionary;
  std::string familyName;
  double a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p;
  if(!myfile.eof()){
    myfile >> familyName;
    while(familyName!="EndSection"){
      myfile >> a >> b >> c >> d >> e >> f >> g >> h >> i >>
	j >> k >> l >> m >> n >> o >> p;
      arma::mat filledMat;
      filledMat << a << b << c << d << arma::endr
		<< e << f << g << h << arma::endr
		<< i << j << k << l << arma::endr
		<< m << n << o << p << arma::endr;
      countDictionary.emplace(familyName,filledMat);
      myfile >> familyName;
    }
  }
  return countDictionary;
}

/**
Modify the input matrix to be the P matrix for the asymmetric model - 
sum each row, and divide each index by its row value.
@param M The Armadillo mat matrix to calculate P from (the count matrix)- 
note this matrix will be modified to reflect the newly calculated P matrix.
*/
void calculate_P(arma::mat& M){
  arma::colvec rowSums = arma::sum(M,1);
  if(rowSums(0)!=0)
    M.row(0) = M.row(0)/rowSums(0);
  if(rowSums(1)!=0)
    M.row(1) = M.row(1)/rowSums(1);
  if(rowSums(2)!=0)
    M.row(2) = M.row(2)/rowSums(2);
  if(rowSums(3)!=0)
    M.row(3) = M.row(3)/rowSums(3);
  
}

/**
Modify the input matrix to be the P matrix for the symmetric model (PSym).
@param M The Armadillo mat matrix to calculate PSym from (the count matrix) - 
note this matrix will be modified to reflect the newly calculated PSym matrix.
 */
void calculate_PSym(arma::mat& M){
  double lambda1=1,lambda2=1;
  arma::colvec rowSums = arma::sum(M,1);
  arma::mat rowCalcs;
  //Ensure we aren't dividing by zero, otherwise, leave row zero - deal with later
  if(rowSums(0)!=0 || rowSums(3)!=0) 
    lambda1 = rowSums(0) + rowSums(3);
  if(rowSums(1)!=0 || rowSums(2)!=0)
    lambda2 = rowSums(1)+rowSums(2);
  rowCalcs = {(M(0,0)+M(3,3))/lambda1,(M(0,1)+M(3,2))/lambda1,
	      (M(0,2)+M(3,1))/lambda1,(M(0,3)+M(3,0))/lambda1};
  M.row(0) = rowCalcs;
  rowCalcs = {(M(1,0)+M(2,3))/lambda2,(M(1,1)+M(2,2))/lambda2,
	      (M(1,2)+M(2,1))/lambda2,(M(1,3)+M(2,0))/lambda2};
  M.row(1) = rowCalcs;
  rowCalcs = {M(1,3),M(1,2),M(1,1),M(1,0)};
  M.row(2) = rowCalcs;
  rowCalcs = {M(0,3),M(0,2),M(0,1),M(0,0)};
  M.row(3) = rowCalcs;
}

/**
Modify the input matrix to be the log (true) or exp (false) of the input matrix
--useful for calculating the Rt matrices and pHats for a section and its families.
@param M An Armadillo mat matrix to be used in the calculations.
@param log A boolean for log function (true) or exp function (false).
- note, modifies input, no return!
*/
void calculate_logORexp(arma::mat& M,bool log){
  arma::cx_mat eigvec, vInverse, Aprime, A;
  arma::cx_vec eigval;
  arma::eig_gen(eigval,eigvec,M);
  vInverse = arma::inv(eigvec); 
  A = vInverse*M*eigvec;
  if(log)
    Aprime = arma::diagmat(arma::log(A));
  else
    Aprime = arma::diagmat(arma::exp(A));
  A = eigvec*Aprime*vInverse;
  M = arma::conv_to<arma::mat>::from(A);
}

/**
Modify the input matrix to be the q matrix given a distance d and an Armadillo matrix 
representing the Rt matrix (log(P)).
@param M An Armadillo mat matrix representing the Rt matrix of a family in an section.
@param d A double representing distance for a particular family in a section.
- note, modifies input, no return!
*/
void calculate_q(arma::mat& M, double d){
  //divides M by d - if d is not finite, makes M zero matrix
  if(std::isfinite(d)){
    M = M/d;
  }else{
    M.zeros();
  }
}

/**
Calculates and returns a distance d given a matrix representing the P matrix
@param M The Armadillo mat matrix representing the P matrix of a family in a section.
@return double A value representing the distance for the given M.
*/
double calculate_d(arma::mat& M){
  //Uses the logdet method to find distance
  double val,sign;
  arma::log_det(val,sign,M);
  return .25*val*-1;
}

/** 
Finds and returns the average matrix from a given matrix list
--useful for finding avgerage Q matrix.  Also allows for filtering based
on age/distance of the family - will only include families that fall into the 
limits into the average.
@param matList A MAP - when finding average Q, this has all families and their q matrices for a given section.
@param dGlobalList A DMAP containing the global distances for all families.
@param lowerlimit A double representing the lower limit on global distance.
@param upperlimit A double representing the upper limit on global distance.
@param applyLimits A bool representing whether or not to apply limits to our averaging.
@return avgM An Armadillo mat matrix containing the average value matrix from matList.
*/
arma::mat calculate_avg(MAP& matList,DMAP& dGlobalList,double lowerlimit,double upperlimit,bool applyLimits){
  //add up all matrices and divide by number of matrices in list
  MAP::iterator it;
  arma::mat avgM(4,4);
  avgM.zeros();
  for(it=matList.begin();it!=matList.end();++it){
    if(it->second.is_finite()){
      if(applyLimits && dGlobalList[it->first] > lowerlimit && dGlobalList[it->first] < upperlimit){
	avgM += it->second;
      }else if(!applyLimits){
	avgM += it->second;
      }
    }
  }
  return avgM/matList.size();
}

/** 
Calculates p hats - the new approximate p matrices - for the given
d values and the avgerage Q (rate) matrix of
the section. Returns a new MAP object with these approximations 
--limits can be applied to only include certain p hats in the new MAP object - useful
for looking at old and new age/distance groups and their impact on symmetry.
--if using just average D (an alternative algorithm), see extra code at bottom.
@param avgQ An Armadillo mat matrix representing the average Q of a section.
@param dList A DMAP containing all the families in a section and their local distances.
@param dGlobalList A DMAP containing all families and their global distances.
@param global A bool, true to use global distances in calculations, false for local ones (true in Karro algorithm).
@param lowerlimit A double representing the lower limit on global distance.
@param upperlimit A double representing the upper limit on global distance.
@param applyLimits A bool representing whether or not to apply limits to our averaging.
@return pHats A MAP containing families and their respective pHat matrix values.
*/
MAP calculate_pHats(arma::mat& avgQ,DMAP& dList,DMAP& dGlobalList,bool global, double lowerlimit, double upperlimit, bool applyLimits){
  DMAP::iterator it;
  arma::mat tmp;
  MAP pHats;
  //for each family in section, use the average Q and unique family d to calculate a new matrix
  for(it=dList.begin();it!=dList.end();++it){
    if(it->second!= 0 && std::isfinite(it->second)){
      if(global){
	  tmp = avgQ*dGlobalList[it->first];
      }else{
	  tmp = avgQ*it->second;
      }
      //take the exp of the new matrix
      calculate_logORexp(tmp,false);
      //don't keep families out of the limits if there are any limits
      if(applyLimits && dGlobalList[it->first] > lowerlimit && dGlobalList[it->first] < upperlimit)
	pHats.emplace(it->first,tmp);
      else if(!applyLimits)
	pHats.emplace(it->first,tmp);
    }
  }
  return pHats;
}

/** 
Calculates the log-likelihood of a given list of family counts and
and pHats for a section.
--if using average pHat of a section (alternative algorithm), see bottom for alternative code.
@param countList A MAP containing the original count matrices for the families in a section.
@param pHatList A MAP containing the pHat matrices for the families in a section.
@return logL The log-likelihood of a given section.
*/ 
double calculate_logL(MAP& countList,MAP& pHatList){
  //for each family's pHat matrix in a section, 
  //take its element-wise log and multiply by its original count matrix
  //and sum this matrix for its logL - sum all logL in a section.
  MAP::iterator it;
  double logL=0;
  for(it=pHatList.begin();it!=pHatList.end();++it){
    if(it->second.is_finite()){
      logL+= arma::sum(arma::sum(countList[it->first]%arma::log(it->second)));
    }
  }
  return logL;
}

/**
Calculate the BIC given the log-likelihood of a section, the number of free variables
in the model, and the number of families in the section.
@param logL A double representing the log-likelihood of the section.
@param f A double representing the number of free variables used in the model (11 for asymmetric,
5 for symmetric).
@param k A double representing the number of families in the section.
@return double The BIC score for the above input.
*/
double calculate_BIC(double logL, double f, double k){
  //BIC formula according to Mugal et al.
  double BIC=0;
  if(logL != 0)//In case of limits excluding all families in a section
    BIC = -2*logL + (f+k)*std::log(16*k);
  return BIC;
}

/**
Method to return the q-balance of a q matrix.
@param Q The Armadillo mat matrix to find the q-balance for.
@return balance The double representing q-balance.
 */
double findQBalance(arma::mat& Q){
  double balance = std::abs(Q(0,1)-Q(3,2)) + std::abs(Q(0,2)-Q(3,1)) 
    + std::abs(Q(0,3)-Q(3,0)) + std::abs(Q(1,0)-Q(2,3)) + std::abs(Q(1,2)-Q(2,1))
    + std::abs(Q(1,3)-Q(2,0)); //+ std::abs(Q(0,0)-Q(3,3)) + std::abs(Q(1,1)-Q(2,2));
  return balance;
}


/**
Work loop to find symmetry of the a section. - Calculates p according to model specification,
then finds local d, the Rt matrices, and q matrices.  Next averages q matrices to get the 
average Q of the section, and uses that Q and global ds of each family in the section 
to calculate the pHats.  The pHats and original counts of the families are then used to find the
log-likelihood, which is used to find the BIC score of the section - the BIC score is returned.
--limits can also be applied in order to include only families with certain global age/distance ranges in 
the calculations - useful for looking at old and new age/distance groups and their impact on symmetry.
@param countList A MAP containing families and their counts for a section. 
@param dGlobalList A DMAP containing all families and their global distances.
@param asymm A bool indicating if the model is asymmetric (true) or symmetric (false).
@param lowerlimit A double representing the lower limit on global distance.
@param upperlimit A double representing the upper limit on global distance.
@param applyLimits A bool representing whether or not to apply limits to our averaging.
@return double The BIC score for the section under the model indicated in the input.
*/
double calculations_loop(MAP& countList,DMAP& dGlobalList,bool asymm,
			 double lowerlimit, double upperlimit, bool applyLimits,double& qb,int& numFamiliesUsed){
  MAP qList,pHatsList;
  DMAP dList;
  double BICVal=0,logL;
  MAP::iterator it;
  arma::mat avgQ;
  //save the countList, we need this later for log likelihood
  qList = countList;
  //for each family in section, do calculations:
  for(it = qList.begin();it != qList.end();++it){
    if(asymm){
      calculate_P(it->second);
    }else{
      calculate_PSym(it->second);
    }
    dList.emplace(it->first,calculate_d(it->second));
    calculate_logORexp(it->second,true);
    calculate_q(it->second,dList[it->first]); //or dGlobalList[it->first] for alternative algorithm
  }
  //Calculate the average Q for the section:
  avgQ = calculate_avg(qList,dGlobalList,lowerlimit,upperlimit,applyLimits);
  qb = findQBalance(avgQ);
  //Calculate pHats for families in section:
  pHatsList = calculate_pHats(avgQ,dList,dGlobalList,true,lowerlimit,upperlimit,applyLimits);
  numFamiliesUsed = pHatsList.size();
  //Alternative method below use: arma::mat avgP = calculate_pHat(avgQ,dList,dGlobalList,true);
  //Calculate logL:
  logL = calculate_logL(countList,pHatsList);
  //Alternative method below use: logL = calculate_logL2(countList,avgP);
  //BIC calculations:
  if(asymm){
    BICVal = calculate_BIC(logL,11,pHatsList.size());
  }else{
    BICVal = calculate_BIC(logL,5,pHatsList.size());
  }
  return BICVal;
}

/**
The main method to determine section symmetry type.  Reads in the file responsible to have sections counts
and calls the calculations loop on each section - once as an asymmetric model and again as a symmetric model.
@param myfile An ifstream object that is the file with the count information to be read for each section.
@param outfile An ostream object that asymmetric results are written to.
@param dGlobalList A DMAP containing all families and their global distances.
@param lowerlimit A double representing the lower limit on global distance.
@param upperlimit A double representing the upper limit on global distance.
@param applyLimits A bool representing whether or not to apply limits to our averaging.
@return results A std::map<std::string,std::string> containing each section and type of symmetry + number of bases.
*/  
std::map<std::string,std::string> findRegionType(std::ifstream& myfile,std::ostream& outfile,DMAP& dGlobalList,
					  double lowerlimit,double upperlimit,bool applyLimits){
  std::map<std::string,std::string> results;
  int count = 0;
  double qbA=0,qbS=0;
  int numFamUsedA=0,numFamUsedS=0;
  if(myfile.is_open()){
    //Go through this section by section
    while(!myfile.eof()){
      MAP countDictionary;
      MAP::iterator it;
      double BICValA,BICValS;
      std::string sectionName;
      myfile >> sectionName;
      //std::cout << sectionName << std::endl;
      if(!myfile.eof()){
	//read in section:
	countDictionary = dealWithPartitions(myfile);
	//calculate under asymmetric model:
	BICValA = calculations_loop(countDictionary,dGlobalList,true,lowerlimit,upperlimit,applyLimits,qbA,numFamUsedA);
	//calculate under symmetric model:
	BICValS = calculations_loop(countDictionary,dGlobalList,false,lowerlimit,upperlimit,applyLimits,qbS,numFamUsedS);
	//find total number of bases used in calculations:
	arma::mat allC(4,4);
	allC.zeros();
	for(it=countDictionary.begin();it!=countDictionary.end();++it){
	  if(applyLimits && dGlobalList[it->first] > lowerlimit && dGlobalList[it->first] < upperlimit)
	    allC += it->second;
	  else if(!applyLimits)
	    allC += it->second;
	}
	int bases = sum(sum(allC));
	//bases in section, followed by q-balance of S, followed by q-balance of A,
	//followed by symmetry type found
	std::string model = std::to_string(bases) + " " + std::to_string(numFamUsedA) + " " + std::to_string(qbA) 
	  + " ";
	if(BICValA==0 && BICValS==0){//in case limits cause no families
	  results.emplace(sectionName,"N");
	}else if(BICValA < BICValS){//asymmetric case
	  model += "A";
	  results.emplace(sectionName, model);
	  outfile << sectionName << "\n";
	  count++;
	}else{//symmetric case
	  model += "S";
	  results.emplace(sectionName,model);
	}
      }
    }
  }else{
    std::cout << "FILE NOT FOUND" << std::endl;
    return results;
  }
  if(!applyLimits)
    std::cout << "There are " << count << " asymmetric regions using all families." << std::endl;
  else
    std::cout << "There are " << count << " asymmetric regions using only families between the distances: "
	      << lowerlimit << ", " << upperlimit << std::endl;
  return results;
}

/**
A method to roughly find limits based on user input.  The percentToTake indicates the percent of total families to 
include from the front and end of the given global distances list (which will be sorted by ages/distances), 
and the cutOffPercent indicates how much to skip in the front and end (it's the starting index for the percentToTake).
--useful to calculate symmetry based on age of the families: Is symmetry different when we only use old or new families?
@param dGlobalList A DMAP containing of all families and their global distances.
@param percentToTake A double the user specifies to indicate the percent of families to include in the limits.
@param cutOffPercent A double the user specifies to indicate how much of the front/end to skip.
@return limits A std::vector<double> that stores the age/distance limits of the families to use.
*/
std::vector<double> findLimits(DMAP& dGlobalList,double percentToTake, double cutOffPercent){
  int index = 0;
  std::vector<double> valList;
  std::vector<double> limits (4,0);
  for(DMAP::iterator it = dGlobalList.begin(); it != dGlobalList.end(); ++it)
    valList.push_back(it->second);
  std::sort(valList.begin(),valList.end());
  index = valList.size()*cutOffPercent;
  limits[0] = valList[index];
  limits[1] = valList[index + (int)(valList.size()*percentToTake)];
  limits[2] = valList[valList.size()-index-(int)(valList.size()*percentToTake)];
  limits[3] = valList[valList.size()-index-1];
 
  return limits;
}

/**
Main function -can be used to asymmetric regions and can apply limits of the ages/distances of 
the families used in the calculations.
Must supply command line with a file to run calculations on and a global count file.
Optional: -l and percent of families to take/outlier cutoff.  
See README for more info on files returned and command line.
 */
int main(int argc, char *argv[]){
  MAP globalDictionary;
  DMAP dGlobalList;
  int countcompare = 0,countSChange=0,countAChange=0;
  std::vector<double> limits (4,0);
  std::string random;
  std::map<std::string,std::string> results,results1, results2;
  MAP::iterator it1;
  DMAP::iterator d1;
  std::ifstream myfile(argv[1]);
  std::ifstream myfile2(argv[2]);
  std::ofstream outfile("AsymmetricSections.txt");//file for asymmetric findings
  std::ofstream outfile2("GlobalDs.txt");//readable global d file
  //deal with global file and get global distances:
  if(myfile2.is_open()){
    myfile2 >> random;
    globalDictionary = dealWithPartitions(myfile2);
    for(it1 = globalDictionary.begin();it1 != globalDictionary.end();++it1){
      calculate_PSym(it1->second);
      dGlobalList.emplace(it1->first,calculate_d(it1->second));
    }
    for(d1 = dGlobalList.begin();d1!=dGlobalList.end();++d1)
        outfile2 << d1->first << " " << d1->second << "\n";
  }
  results = findRegionType(myfile,outfile,dGlobalList,limits[0],limits[0],false);
  myfile.close();
  outfile.close();
  //if limits are used (use -l):
  if(argc > 3){
    std::ifstream myfile2(argv[1]);
    std::ofstream outfile3("AsymmetricSectionsNew.txt");
    std::ofstream outfile4("AsymmetricSectionsOld.txt");
    std::ofstream outfile5("EasyReadCompare.txt");
    std::ofstream karrofile("LineCompareFile.txt");
    std::ofstream karrofile2("AllFile.txt");
    if(argc > 4)
      limits = findLimits(dGlobalList,atof(argv[4]),atof(argv[5]));
    else
      limits = findLimits(dGlobalList,.3,.1);
    outfile3 << "d limits: " << limits[0] << ", " << limits[1] << "\n";
    outfile4 << "d limits: " << limits[2] << ", " << limits[3] << "\n";
    results1 = findRegionType(myfile2,outfile3,dGlobalList,limits[0],limits[1],true);
    myfile2.close();
    std::ifstream myfile3(argv[1]);
    results2 = findRegionType(myfile3,outfile4,dGlobalList,limits[2],limits[3],true);
    outfile3.close();
    outfile4.close();
    //compare the new, old, and all family results!
    for(std::map<std::string,std::string>::iterator it2 = results.begin(); it2 != results.end(); it2++){
      if(results1[it2->first].back() != results2[it2->first].back()){
	outfile5 << "Using 'new' families " << it2->first << " is: " << results1[it2->first] << "\n";
	outfile5 << "Using 'old' families " << it2->first << " is: " << results2[it2->first] << "\n";
	outfile5 << "Using all families " << it2->first << " is: " << it2->second << "\n\n";
	karrofile << it2->first << " " << results1[it2->first] << " " << results2[it2->first] << " " << it2->second << "\n\n";
	countcompare++;
	if(results2[it2->first].back()=='S')
	  countSChange++;
	if(results2[it2->first].back()=='A')
	  countAChange++;
      }
      karrofile2 << it2->first << " " << results1[it2->first] << " " << results2[it2->first] << " " << it2->second << "\n\n";
    }
    std::cout << "There are " << countcompare << " regions that change symmetry when comparing new/old family sets."<< std::endl;
    std::cout << "There are " << countSChange << " regions that change from symmetric to asymmetric over time (old->new)" << std::endl;
    std::cout << "There are " << countAChange << " regions that change from asymmetric to symmetric over time (old->new)" << std::endl;
    outfile5.close();
    karrofile.close();
    
  }
  return 0;
}



/**
Alternative to calculating pHats and logL above- this method creates a single
pHat matrix from the average of a specified d list and an average Q instead of multiple 
pHats from individual family ds and an average Q.  The following method calculates the 
logL from the single pHat using the total base count of the section instead of using 
multiple pHats and each family's count.  Note: limits not available for these methods.


arma::mat calculate_pHat(arma::mat& avgQ,DMAP& dList,DMAP& dGlobalList,bool global){
  DMAP::iterator it;
  double sum=0,davg=0,count=0;
  for(it=dList.begin();it!=dList.end();++it){
    if(it->second!= 0 && std::isfinite(it->second)){
      if(global){
        count++;
	sum+=dGlobalList[it->first];
      }else{
        count++;
	sum += it->second; 
      }
    }
  }
  davg = sum/count;
  arma::mat tmp = avgQ*davg;
  calculate_logORexp(tmp,false);
  return tmp;
}


double calculate_logL2(MAP& countList, arma::mat phat){
  MAP::iterator it;
  double logL;
  arma::mat allC(4,4);
  allC.zeros();
  for(it=countList.begin();it!=countList.end();++it)
    allC += it->second;
  logL = arma::sum(arma::sum(allC%arma::log(phat)));
  return logL;
}

*/
