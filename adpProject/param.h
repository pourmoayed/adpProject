#ifndef PARAM_HPP
#define PARAM_HPP

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
#include <iostream>
#include <armadillo>
#include <vector>
using namespace std;
using namespace arma;

/** 
* Class: Initialize the parameters: 
* @param pigs Number of pigs per pen
* @param tMax maximum length of growing period
* @param pNum Number of pens in section
* @param secNum Number of sections in herd
* @param simNum Number of simulation runs in AVI algorithm
* @param epochNum Number of decision epochs in AVI algorithm
* @param truckNum Maximum number of available trucks
* @param k1Feed Feed conversion rate related to metabolic weight 
* @param k2Feed Feed conversion rate related to growth
* @param numSampleDGLSM Number of samples in DGLM model
* @param tStartMarketing Start time of marketing decisions 
* @param minMarketingSize Minimum number of marketed pigs 
* @param cleaningPeriod Cleaning period at pen after terminiation (days) 
* @param marketingLength Time period between marketing decisions until sending to abattoir(days)
* @param avgGRate Average growth rate at herd
* @param avgLeanP Average leanness percentage at herd
* @param avgInsWeight Average insertion weight into the pens when starting a production process
* @param avgInsSd Average standard deviation in the first week of production
* @param rndIniSd Initial standard deviation values when starting a production process(randomely selected)
* @param avgInsG Average growh rate at herd per week (kg)
* @param adpPost Boolean variable showing the ADP is running under post-decision states (true) or only pre-decision states     
* @param comparePolicy Boolean variable used after for numerical experiments. True: comparing different marketing policies. False: An example for marketing decisions at herd  
* @param meanWeights A vector containing the mean weight of pigs during growing period.    
* @param sdWeights A vector containing standard deviations of weights during growing period   
* @param observedWeight A sample of weights during growing period used in DLM model   
* @param observedFeed A sample of feed intake values during growing period used in DLM model   
* @param sdMeasure Measurement error of weighing 
* @param feedPrice Unit price of feed (FEsv)  
* @param porkPrice Maximum unit price per kg carcass
* @param pigletPrice Unit price of piglet 
* @param truckCost Fixed cost of truck 
* @param truckCapacity Capacity of one truck 
* @param discount Discount value showing time value of money used in the AVI algorithm 
* @param sL Possible discrete values for the deviation of weights from average weight during growing period 
* @param sG Possible discrete values for the average growth of pigs in pens 
* @param sSd Possible discrete values for the deviation of standard deviations of weights from the average stadard deviations during growing period 
* @param sL0 Initial index value related to sL  
* @param sG0 Initial index value related to sG 
* @param sSd0 Initial index value related to sSd 
* @param iniL Initial values of sL for 5 sections   
* @param iniG Initial values of sG for 5 sections   
* @param iniSd Initial values of sSd for 5 sections   
* @param dL A matrix inculding possible discrete intervals for sL  
* @param dSd A matrix inculding possible discrete intervals for sSd
* @param dG A matrix inculding possible discrete intervals for sG
* @param m0 Initial posterior mean for DLM model
* @param C0 Initial posterior covariance matrix for DLM model
* @param W Covariance matrix related to system error in DLM model
* @param V Covariance matrix related to observation error in DLM model
* @param FF Design matrix of DLM model (observation equation)
* @param GG FF Design matrix of DLM model (system equation) 
* @param L Covariance matrix for calculating transition probabilities in DLM 
*
* @author Reza Pourmoayed (rpourmoayed@econ.au.dk)
*/
class PARAM
{ 
	/** A structure to return three cubes from function \code(DLMkalman()) used to store the required information in DLM model.  
     */  
	struct myList {
	  arma::cube g;
	  arma::cube f;
	  arma::cube l;
   };

public:
	/** Constructor. Store the recquired parameters .
	* 
     */
	PARAM(){
	  pigs = 18;
      tMax = 15;
      penNum = 20;
      secNum = 3;
	  simNum = 300;
	  epochNum = 120;
	  truckNum = 6;
	  k1Feed = 0.044;
	  k2Feed = 1.4;
	  numSampleDGLSM = 35;
      tStartMarketing = 9; //we start from t=0!
      minMarketingSize = 0;
      cleaningPeriod = 4;
      marketingLength = 3; //1;
      avgGRate = 6;
      avgLeanP = 61;
      avgInsWeight = 30;
      avgInsSd = 2.2;
	  rndIniSd << 1.8 << 2.8 << 3.8 << 4.8 <<endr; 
	  avgInsG = 6;
	  adpPost = true;
	  comparePolicy = true;
	  //meanWeights << 28.00<<  35.13<<  42.27<<  49.40<<  56.53<<  63.67<<  70.80<<  77.93<<  85.07<<  92.20<<  99.33<< 106.47<< 113.60<< 120.73<< 127.87<< 135.00<< 143 << endr;	  
	  meanWeights << 26.4 <<  33.4 <<  40.4 <<   47.4 <<   54.4 <<   61.4 <<   68.4 <<   75.4 <<  82.4 <<   89.4<<   96.4 <<  103.4 <<  110.4 <<  117.4 <<  124.4 <<  131.4 << endr;
	  sdWeights <<  6 << 6.5 <<  7 <<  7.5 <<  8 << 8.5 << 9 <<  9.5 << 10 << 10.5 << 11 << 11.5 << 12 << 12.5 << 13 << 13.5 <<  endr;
      observedWeight << 28.91122 <<  34.73722 <<  40.55647 <<  46.43588 <<  52.60592 <<  58.90283 <<  64.73567 <<  70.88710 <<  76.19884 <<  82.30708 <<  88.41322 <<  94.68169 <<  100.68169 <<  106.68169 <<  112.68169 << endr;
	  observedFeed << 11.17050 <<  11.57909 <<  11.91425 <<  12.31422 <<  12.76750 <<  13.18425 <<  13.57245 <<  13.94777 <<  14.36087 <<  14.76617 <<  15.10727 <<   15.56838 <<  15.96838 <<  16.57838 <<  17.06838  << endr;
	  sdMeasure = 1;
	  feedPrice = 1.8;
      porkPrice = 10.3;
      pigletPrice = 375;
	  truckCost = 400;
	  truckCapacity = 205;
	  discount = (double)1/(double)(1+0.05); 
	  sL << -20 << -18 << -16 << -14 << -12 << -10 << -8 << -6 << -4 << -2 << 0 << 2 << 4 << 6 << 8 << 10 << 12 << 14 << 16 << 18 << 20 << endr;
	  //sSd << -0.8*6  <<  -0.8*5 <<  -0.8*4 << -0.8*3 << -0.8*2  << -0.8*1 << 0 << 0.8*1 << 0.8*2  <<  0.8*3 <<  0.8*4 << 0.8*5 << 0.8*6  << endr;
	  sG << 4.2<< 5 << 5.8 << 6.6 << 7.4 << 8.2 << endr;
	  //sL << -10 << -8 << -6 << -4 << -2 << 0 << 2 << 4 << 6 << 8 << 10 << endr;
	  sSd << -4 << -3 << -2 << -1 << 0 << 1 << 2 << 3 << 4 << endr;
	  //sG << 5 << 5.5 << 6 << 6.5 << 7 << 7.5 << 8 << 8.5 << 9 << endr;
	  sL0 = 10;
	  sG0 = 3;
	  sSd0 = 1;
	  iniL << 0 << 0 << 0 << 0 << 0 << endr;
	  iniG << 6 << 6 << 6 << 6 << 6 << endr;
	  iniSd << -3 << -3 << -3 << -3 << -3 <<endr; 
	  coefL << 1 << 1 << 1 << 1 << 1 << endr;
	  coefG << 1 << 1 << 1 << 1 << 1 << endr;
	  coefSd << 1 << 1 << 1 << 1 << 1 <<endr;
	  dL = discritize(sL,-1000,1000);
	  dSd = discritize(sSd, -1000 , std::numeric_limits<float>::infinity() );
	  dG = discritize(sG,-1000, 1000);

	  m0 << 26.4 << endr 
	      << 5.8 << endr;

	  
	  C0 << 5 << 0.9 << endr
		  << 0.9 << 2 << endr;	  
	  
	  //C0 << 4.26 << 0.9 << endr      
		 // << 0.9 << 0.53 << endr;

	  W << 2.1 << -0.124 <<  endr
		 << -0.124 << 0.112 << endr;

	  V << (double)1/(double)pigs << 0.027 << endr
		 << 0.027 << 0.012 << endr;

	  FF = DLMkalman().f;
	  GG = DLMkalman().g;
	  L = DLMkalman().l;

	}
	// methods

	//variables
  int pigs;
  int tMax, tStartMarketing, sL0, sG0, sSd0;
  int penNum, secNum, simNum, epochNum, truckNum, numSampleDGLSM;
  int minMarketingSize;
  double cleaningPeriod;
  double marketingLength;
  double avgGRate;
  double k1Feed;
  double k2Feed; 
  double avgLeanP;
  double avgInsWeight;
  double avgInsSd;
  double avgInsG;
  double feedPrice;
  double porkPrice;
  double pigletPrice;
  double truckCost;
  double sdMeasure;
  double discount;
  int truckCapacity;
  bool adpPost;
  bool comparePolicy;

  vec sdWeights;
  vec meanWeights;
  vec observedWeight;
  vec observedFeed;
  vec rndIniSd;
  vec iniL;
  vec iniG;
  vec iniSd;
  vec coefL;
  vec coefG;
  vec coefSd;  
  vec sL;
  vec sSd;
  vec sG;
  mat dL;
  mat dSd;
  mat dG;

  mat m0;
  mat C0;
  mat W;
  mat V;  
  cube FF;
  cube GG;
  cube L;
 
private:
	//methods 
	/** Build design matrices of the SSM and use the Kalman filter to find the covariance matrices 
	*/
   myList DLMkalman(){
     vec timeVectorGG(tMax);
	 cube FF(2,2,tMax);
	 cube GG(2,2,tMax);

	 int t;
	 myList mt;
 
   cube m(2,1,tMax); // means of posterior - preallocate an empty list of length tMax
   cube C(2,2,tMax); // covariance matrices of posterior
   cube L(2,2,tMax); // covariance matrices of transition probabilities of latent variables for one ahead forecast.
   cube Q(2,2,tMax); // covariance matrices of transition probabilities of observable variables.
   cube A(2,2,tMax); // a Matrix defined in the steps of DLM.

   mat a(2,1);
   mat R(2,2);
   mat f(2,1);  
   cube y(2,1,tMax);	 
   
	 for(t=0; t<tMax;t++){
      if( (t==0) || (t==1) || (t==1) ){
         timeVectorGG[t] = sqrt(1.8);  
      }else{
         timeVectorGG[t] = (double)(t)/(double)(t-1);
      }
    }

    for(int i=0; i<tMax; i++){
		FF.slice(i) << 1 << k1Feed << endr
		            << 0 << k2Feed << endr;
		GG.slice(i) << 1 << 1 << endr
	             	<< 0 << 1 << endr;
		y.slice(i) << observedWeight[i]  << endr
			        << observedFeed[i] << endr;
		}

  
   for(t=0; t<tMax; t++){
      //Prior
      if(t==0){
		  a = GG.slice(t) * m0;
          R = GG.slice(t) * C0 * GG.slice(t).t() + W;
      }
      else{
         a = GG.slice(t) * m.slice(t-1);
         R = GG.slice(t) * C.slice(t-1) * GG.slice(t).t() + W;
      }
      // One step forcast
	  FF.slice(t)(0,1) = 7*k1Feed/(pow(a(0,0),0.25)); 
	  f = FF.slice(t).t() * a;
      Q.slice(t) =  FF.slice(t).t() * R * FF.slice(t) + V;
      ////Posterior (we see Y_t here)
      A.slice(t) = R * FF.slice(t) * inv(Q.slice(t));
      m.slice(t) = a + A.slice(t) * (y.slice(t)-f);
      C.slice(t) = R - A.slice(t) *  Q.slice(t) * A.slice(t).t(); 
      //Compute the covariance matrix for transition probabilities
      L.slice(t) = R - C.slice(t);
   }

   mt.f = FF;
   mt.g = GG;
   mt.l = L;

   return(mt);
}

   /** find the matrix including the discritized values of state variables given the vector of center points
   * @param givenVec Vector of center points for a state variable
   * @param lower Lower bound of possible values
   * @param upper Upper bound of possible values
   *
   * @return A matrix of discritized intervals and center points 
	*/
   mat discritize(vec givenVec, double lower, double upper ){
	   mat discrite(givenVec.n_elem,3);
	   for(int j=0;j<givenVec.size(); j++)
		   discrite(j,0) = givenVec[j];
		   for(int j=0;j<givenVec.size(); j++){
			   if(j==0){
				   discrite(j,1) = lower; //-infinty
				   discrite(j,2) = discrite(j,0) + ( discrite(j+1,0) -  discrite(j,0) )/2;
			   }
			   if(j==(givenVec.size() -1) )
				   discrite(j,2) = upper; //infinity
			   if(j!=0)
				   discrite(j,1) = discrite(j-1,2);
			   if(j!=(givenVec.size() -1 ) )
				   discrite(j,2) = discrite(j,0) + ( discrite(j+1,0) -  discrite(j,0)  )/2;
		   }
		   return(discrite);
   }



   /**  Descritize the state space of the HMDP based on the number of the stage: Subfunction : Find the centerpoints for the live weight
    * and stansard deviation  
    * 
    * @param Week The current stage (week) in the production system
	* @param numPoints Number of intervals for discritization
    * @param length The length of the related intervals in the current styage
	* @param givenInfo Information for the possible values for a state variable.
    *
    * @return A vector included the center points to discritize the state space in the current stage       
    * @author Reza Pourmoayed \email{rpourmoayed@econ.au.dk}
   */

   arma::vec pointSW(int week, int numPoints, int length, arma::vec givenInfo){  
   int x;
   arma::vec centerPoints;
   x = 2*numPoints+1; 
   for(int i=0; i<numPoints; i++){
      centerPoints[i]= givenInfo[week] - length*( numPoints - i );
      centerPoints[x-(i+1)]= givenInfo[week] + length*( numPoints - i );
   }
   centerPoints[numPoints]= givenInfo[week];
   return(centerPoints);
}

//  C0 << 5 << 0.9 << endr
//		  << 0.9 << 2 << endr;

   //variables:

};

#endif