#ifndef SAMPLEPATH_H
#define SAMPLEPATH_H

#include <armadillo>  
#include <iostream>
#include <vector>
#include <random>
#include <boost/math/distributions.hpp>
using namespace boost::math;

#include "param.h"
#include "cumNorm.h"
using namespace std;

/**
* This class uses simulation methods to simulate the next states and compute the reward coeficiants of the objective function.
*
* @author Reza Pourmoayed (rpourmoayed@econ.au.dk)
*/
class samplePath
{

public:  // methods

    /** Constructor of the class sample path
     *
     */
	samplePath(): p() {}

//methods
	
  /** simulate a value from the posterior distribution of the GSSM.
   *
   * @param t Index time.
   * @param m Posterior at time t
   *
   * @return A vector of simulated posterior for time t+1
   * @author Reza Pourmoayed (rpourmoayed@econ.au.dk)
   */
  arma::vec simGSSM(int t, arma::vec m){
	arma::vec mean;
    arma::mat sigma(2,2);
    arma:: mat transferMat(2,2);
	if(t==0){
	  mean = p.m0;
	  sigma = p.C0;
    }else{
	  mean  = p.GG.slice(t) * m; // G_{t+1}m_t
	  sigma = p.L.slice(t); // L_{t+1}  
    }
	//if choleski decompositio possibe: 
	transferMat = arma::chol(sigma);
	// if  spectral decomposition possible
	//arma::vec eigval;
    //arma::mat eigvec;
    //arma::eig_sym(eigval, eigvec, sigma);
	//transferMat =  eigvec * sqrt( abs(diagmat(eigval) ) );
	return( mean + (arma::mean( arma::randn<arma::mat>(1,sigma.n_cols) )  * transferMat).t() ); 
  }

     /** Transition probability related to GSSM, i.e. bivariate cumulative normal 
	 */  
   double logPrGSSM(int t, arma::vec lower, arma::vec upper, arma::vec mt) { 

	   arma::vec mean;
       arma::mat sigma(2,2);
       arma:: mat transferMat(2,2);
       if(t==0){
    	  mean = p.m0;
    	  sigma = p.C0;
       }else{
    	  mean  = p.GG.slice(t) * mt; // G_{t+1}m_t
    	  sigma = p.L.slice(t); // L_{t+1}  
       } 
      return log( pNorm2D_arma(lower, upper, mean, sigma) * pow(10.0,10) );             
   }

  
  /** simulate a value from the posterior distribution of the nGSSM (i.e. simulate from beta prime distribution).
   *
   * @param t Index time.
   * @param sd posterior at time t
   *
   * @return a simulated posterior for time t+1
   * @author Reza Pourmoayed (rpourmoayed@econ.au.dk)
   */
  double simNGSSM(int t, double sd){ // note: t refers to tN that is the time in the next week!
	  double G, a, alpha, gamma, var, oShape, iShape, randBeta, randBetaPrime, randUniform;
	  int RandIndexSd;
	  oShape = (double) (p.numSampleDGLSM-1)/(2);   //shape parameter of observation distribution
	  iShape = (double) (p.numSampleDGLSM-1)/(2);//(double) (numSample-3)/(numSample-5); //("shape parameter of prior at t=1"): c_1 in the paper

	  if( t<p.tStartMarketing-1 ){  
         G = 1.1;
      }else{
         G = pow( (double) (t)/(t-1),2);
      }       


	  if(t==0){
		  RandIndexSd = rand() % (p.rndIniSd.n_elem-1); 
		  var = pow(p.rndIniSd[RandIndexSd],2);
		  a = G * var; // location
	  }else{
		  var = pow(sd,2);
		  a = (double) (G * var * ( iShape + oShape*(t-1) ) ) / ( iShape + oShape*(t) ) ; // location
	  }
	       
      alpha = oShape; // shape 1
      gamma = iShape + oShape*(t-1) + 1; // shape 2
      //beta =1 ;  // Weibul parameter

	
	  //generate a random sample from beta distribution using the properties of gamma distribution
	  
	  /*boost::random::mt19937 rng; 
	  uniform_01_distribution<> dis();
	  randFromUnif = uniform_01(rng); */
	  randUniform = ( (double) rand() / (double) (RAND_MAX) );
	  beta_distribution<> dist(alpha, gamma);
	  randBeta = quantile(dist, randUniform );
	  if(randBeta>0.99) randBeta = 0.99; 

	  //typedef std::mt19937 Ga;    
	  //typedef std::gamma_distribution<> D;
	  //Ga g;  	    
	  //D x(alpha, 1);	
	  //D y(gamma, 1);
	  //randBeta = (double)x(g) / (double)(x(g) + y(g));

	  //find a sample from beta prime distribution using a random value from beta distribution
	  randBetaPrime = (double)a / (double)(1-randBeta);
    
 	return( sqrt(randBetaPrime) ); 
  }


     /** Transition probability related to nGSSM, i.e. inverse gamma cdf. 
    *
    * @param t Week (time t).
    * @param n Number of pigs at time t.
    * @param lower Lower bound on variance at time t+1.
    * @param upper Upper bound on variance at time t+1.
    * @param var Estimate of variance at time t.
    */
   double logPrNGSSM(int t, double lower, double upper, double var) { //SOLVED[Reza] : Based on the formulatiom for this probability in the paper, I changed "n" to "nf" (nf is the sample size). 

      double logProbSd, xUpper, xLower, a, s, alpha, gamma, beta, G; 
	  double oShape = (double) (p.numSampleDGLSM-1)/(2);   //shape parameter of observation distribution
	  double iShape = (double) (p.numSampleDGLSM-1)/(2);//(double) (numSample-3)/(numSample-5); //("shape parameter of prior at t=1"): c_1 in the paper

	  if( t<p.tStartMarketing-1 ){ // 
         G = 1.1;
      }else{
         G = pow( (double) (t+1)/(t),2);
      }          
            
	  //if(t==0){
		 // a = G * var; // location
	  //}else{
		 // a = (double) (G * var * ( iShape + oShape*t ) ) / ( iShape + oShape*(t+1) ) ; // location
	  //}

	  a = (double) (G * var * ( iShape + oShape*t ) ) / ( iShape + oShape*(t+1) ) ; // location
      s = a; // scale
      alpha = oShape; // shape 1
      gamma = iShape + oShape*t + 1; // shape 2
      beta =1 ;  // Weibul parameter
      xUpper = (double) (1)/( 1 + pow ( (double)(upper-a)/(s),-beta ) );
	  if(xUpper<0) xUpper = 0;
      xLower = (double) (1)/( 1 + pow ( (double)(lower-a)/(s),-beta ) ); 
	  if(xLower<0) xLower = 0;

	  beta_distribution<> dist(alpha, gamma);  	  
	  logProbSd = log( (cdf(dist,xUpper) -  cdf(dist,xLower) ) * pow(10.0,10) ); // ?
      if( exp(logProbSd)<0 ) cout<<"error_minus"<<endl;
     
      return (logProbSd); 
      }


  /** Simulate the pen and find the revenue of selling and cost of feeding.
   *
   * @param time Week number
   * @param weighDev Posterior mean related to weight deviation from average weight at herd
   * @param weightG Posterior mean related to growth
   * @param weightSdDev Posterior mean related to deviations of sd from average sd at herd
   *
   * @author Reza Pourmoayed (rpourmoayed@econ.au.dk)
   */
 void simulateRew( int time, double weighDev, double weightG, double weightSdDev ){
    
	int samples=1000;
	double weightMean = p.meanWeights[time] + weighDev;
	double weightSd = p.sdWeights[time] +  weightSdDev;
	double dailyGrowth = (double)weightG/(double)7;

	unsigned int pigNums =  p.pigs;

    arma::mat leanT(samples,pigNums);  // leanness
    arma::mat lWT(samples,pigNums);    // live weight
	arma::mat oWT(samples,pigNums);    // live weight
    arma::mat sortOWT(samples,pigNums); 
	arma::mat sortLWT(samples,pigNums);
	arma::mat sumFeed(samples,pigNums);  // feed sums for 7 days
    arma::mat sumFeed3(samples,pigNums);  // feed sums for 3 days
	arma::mat rewT(samples,pigNums);   // reward
	arma::mat cWT(samples,pigNums);    // carcass weight
    arma::mat sortCWT(samples,pigNums);    // carcass weight of sorted pigs
	arma::mat dailyGT(samples,pigNums);
	arma::rowvec rev(pigNums,arma::fill::zeros);   // avg. reward of culling the k'th pig
    arma::rowvec cost(pigNums,arma::fill::zeros);   // avg. cost of feeding the k'th pig
	arma::rowvec lWeight(pigNums,arma::fill::zeros);   // avg. cost of feeding the k'th pig
	
    //  simulation 
    oWT = arma::randn<arma::mat>(samples,pigNums);   // use std. norm dist arma function
	oWT = oWT*sqrt( pow(weightSd,2) + pow(p.sdMeasure,2) ) + weightMean; // transform to general norm dist
	if (weightMean<50) oWT = abs(oWT);
	sortOWT = arma::sort(oWT,"ascend",1);
	sortLWT = sortOWT - arma::randn<arma::mat>(samples,pigNums)*p.sdMeasure;
	sumFeed.fill(0);//zeros();  // initialize
    sumFeed3.fill(0);  // feed sums for 7 days
	dailyGT.fill(dailyGrowth);

	for (int d=1;d<8;d++) { // days in week
        sumFeed = sumFeed + 1.549*dailyGT + 0.044*pow(abs(sortLWT),0.75);
        if (d==p.marketingLength) {
          // carcass weight
          sortCWT = 0.84*sortLWT-5.89 + arma::randn<arma::mat>(samples,pigNums)*1.4;
		  sortCWT.elem( arma::find(sortCWT < 0) ).zeros();
          // leanness
          leanT = (-30*(dailyGT-0.8571429))/4 + 61;   // assume constant growth during the week
          sumFeed3 = sumFeed;
        }
        sortLWT = sortLWT + dailyGT;
       // sortLWT.elem( arma::find(sortLWT < 0) ).zeros();  // set negative to zero - hack to not get negative weights for very rare samples
      }
	for(unsigned int i=0; i<leanT.n_elem; ++i) {
        leanT[i] = priceLeanness(leanT[i]);
		rewT[i] = priceCarcass(sortCWT[i]);
      }
	 rewT = sortCWT % (leanT + rewT) - sumFeed3*p.feedPrice;
     rev = arma::mean(rewT);
     cost = arma::mean(sumFeed*p.feedPrice);
	 lWeight = arma::mean(sortLWT);
	 cost.insert_cols(0,arma::colvec(1,arma::fill::zeros) );
	 rev.insert_cols(0,arma::colvec(1,arma::fill::zeros) );
	 lWeight.insert_cols(0,arma::colvec(1,arma::fill::zeros) );
	 revCoef = rev;
     costCoef = cost;
	 sortLW = lWeight;
  }

  /** Simulate the pen and find the revenue of selling and cost of feeding for an averge weighted pig.
   *
   * @param time Week number
   * @param weighDev Posterior mean related to weight deviation from average weight at herd
   * @param weightG Posterior mean related to growth
   * @param weightSdDev Posterior mean related to deviations of sd from average sd at herd
   * 
   * @author Reza Pourmoayed (rpourmoayed@econ.au.dk)
   */
 void simulateRewAvg( int time, double weighDev, double weightG, double weightSdDev ){
    
	int samples=1000;
	double weightMean = p.meanWeights[time] + weighDev;
	double weightSd = p.sdWeights[time] +  weightSdDev;
	double dailyGrowth = (double)weightG/(double)7;

	unsigned int pigNums =  p.pigs;

    arma::vec leanT(samples);  // leanness
    arma::vec lWT(samples);    // live weight
	arma::vec oWT(samples);    // live weight
    arma::vec sortOWT(samples); 
	arma::vec sortLWT(samples);
	arma::vec sumFeed(samples);  // feed sums for 7 days
    arma::vec sumFeed3(samples);  // feed sums for 3 days
	arma::vec rewT(samples);   // reward
	arma::vec cWT(samples);    // carcass weight
    arma::vec sortCWT(samples);    // carcass weight of sorted pigs
	arma::vec dailyGT(samples);
	double avgRev; // avg. reward of culling
    double avgCost; // avg. cost of feeding 
	
    //  simulation 
    oWT = arma::randn<arma::mat>(samples);   // use std. norm dist arma function
	oWT = oWT*sqrt( pow(weightSd,2) + pow(p.sdMeasure,2) ) + weightMean; // transform to general norm dist
	if (weightMean<50) oWT = abs(oWT);
	lWT = oWT - arma::randn<arma::mat>(samples)*p.sdMeasure;
	sumFeed.fill(0);//zeros();  // initialize
    sumFeed3.fill(0);  // feed sums for 7 days
	dailyGT.fill(dailyGrowth);

	for (int d=1;d<8;d++) { // days in week
        sumFeed = sumFeed + 1.549*dailyGT + 0.044*pow(abs(lWT),0.75);
        if (d==p.marketingLength) {
          // carcass weight
          cWT = 0.84*lWT-5.89 + arma::randn<arma::mat>(samples)*1.4;
		  cWT.elem( arma::find(cWT < 0) ).zeros();
          // leanness
          leanT = (-30*(dailyGT-0.8571429))/4 + 61;   // assume constant growth during the week
          sumFeed3 = sumFeed;
        }
        lWT = lWT + dailyGT;
       // sortLWT.elem( arma::find(sortLWT < 0) ).zeros();  // set negative to zero - hack to not get negative weights for very rare samples
      }
	for(unsigned int i=0; i<leanT.n_elem; ++i) {
        leanT[i] = priceLeanness(leanT[i]);
		rewT[i] = priceCarcass(cWT[i]);
      }
	 rewT = cWT % (leanT + rewT) - sumFeed3*p.feedPrice;
     avgRev = arma::mean(rewT);
     avgCost = arma::mean(sumFeed*p.feedPrice);
	 avgRevCoef = avgRev;
     avgCostCoef = avgCost;
  }


//variables:
  arma::rowvec revCoef;
  arma::rowvec costCoef;
  arma::rowvec sortLW;
  double avgRevCoef;
  double avgCostCoef;

private:
//methods

  /** Find selling price of an individual pig using market price p.
   */
  double priceCarcass(double sWeight){
   //Interval [0,60) :
   if(0<=sWeight && sWeight<60) return( ( 10.3 - 0.1*(70-60) - 0.2*(60-sWeight) ) );
   //Interval [60,70) :  
   if(60<=sWeight && sWeight<70) return( ( 10.3 - 0.1*(70-sWeight) ) );
   //Interval [70,86) : 
   if(70<=sWeight && sWeight<86) return( 10.3);
   //Interval [86,95) :  
   if(86<=sWeight && sWeight<95) return( ( 10.3 - 0.1*(sWeight-86) ) );
   //Interval [95,100) :
   if(95<=sWeight && sWeight<100) return(9.3);
   //Interval [100,Inf) : 
   if(100<=sWeight) return(9.1);
   return(0);
}

  /** Find  leannes bonus of an individual pig.
   */
  double priceLeanness(const double & lean) {
    return(0.1*(lean - 61));
  }

private:
//variables:
  PARAM p;
};

#endif
