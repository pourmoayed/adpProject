
#include <ilcplex/ilocplex.h>
#include<iostream>
using namespace std; 
#include <armadillo>
using namespace arma;
#include <random>
//#include <boost/math/distributions.hpp>
//using namespace boost::math; 
#include <fstream>      // fstream
#include <vector>
#include <string>
#include <algorithm>    // copy
#include <iterator>     // ostream_operator
#include <boost/tokenizer.hpp>
using namespace boost;
#include "param.h"
#include "samplePath.h"
#include "adp.h"
#include "max.h"
#include "cumNorm.h"

void main()
{    
  PARAM p;

  //Estimate the slope values of the approximate value function (Note: results for  Experiment 1 in the paper can be obtained by running ap.valueADP() in one single pen)
  //srand( time(0) );
  //ADP ap;
  //ap.valueADP();

  //Find the results for marketing decisions at herd level (Experiment 2 in the paper)
  if(!p.comparePolicy){
  srand( 1234567890 ); // rand for making plot at herd level.
  ADP ap;
  ap.optActions();
  }
  
  //Find the results for comparing marketing policies (Experiment 3 in the paper) 
  if(p.comparePolicy){
  srand( 12345678 ); // rand for comparing policies.
  ADP ap;
  ap.rewCoefficients();
  ap.readFinalSlopes();
  int dim = 100;
  arma::vec rewVec(dim);
  arma::vec feedVec(dim);
  arma::vec avgWeightVec(dim);
  arma::vec cycleVec(dim);
  arma::vec cycleLengthVec(dim);
  arma::vec truckVec(dim);
  arma::vec marketVec(dim);
  arma::vec truckUtilizeVec(dim);

  //made samples for srand;
  arma::vec randNum(dim);
  for(int j=0; j<dim; j++){
  randNum[j] = rand() % (99999999999);
  }
 

  for(int i=0; i<dim; i++){
  srand( randNum[i] );
  ap.optActions();
  rewVec[i] = ap.RewardWeek;
  feedVec[i] = ap.FeedWeek;
  avgWeightVec[i] = ap.avgMarkWeight;
  cycleVec[i] = ap.cycleNumber;
  cycleLengthVec[i] = ap.cycleLength;
  truckVec[i] = ap.sumTruckNum;
  marketVec[i] = ap.sumMarket;
  truckUtilizeVec[i] = (double) (ap.sumMarket) / (double) (ap.sumTruckNum*p.truckCapacity);
  cout<< " simNum " << i << endl;
  }

  cout<<" meanReward "<< arma::mean(rewVec) <<endl;
  cout<<" sdReward "<< sqrt(arma::var(rewVec)) <<endl;
  cout<<" meanFeed "<< arma::mean(feedVec) <<endl;
  cout<<" sdFeed "<< sqrt(arma::var(feedVec)) <<endl;
  cout<<" meanavgWeight "<< arma::mean(avgWeightVec) <<endl;
  cout<<" sdavgWeight "<< sqrt(arma::var(avgWeightVec)) <<endl;
  cout<<" meancycle "<< arma::mean(cycleVec) <<endl;
  cout<<" sdcycle "<< sqrt(arma::var(cycleVec)) <<endl;
  cout<<" meancycleLength "<< arma::mean(cycleLengthVec) <<endl;
  cout<<" sdcycleLength "<< sqrt(arma::var(cycleLengthVec)) <<endl;
  cout<<" meantruckNum "<< arma::mean(truckVec) <<endl;
  cout<<" sdtruckNum "<< sqrt(arma::var(truckVec)) <<endl;
  cout<<" meanmrketNum "<< arma::mean(marketVec) <<endl;
  cout<<" sdmarketkNum "<< sqrt(arma::var(marketVec)) <<endl;
  cout<<" meantruckUtilize "<< arma::mean(truckUtilizeVec) <<endl;
  cout<<" sdmtruckUtilize "<< sqrt(arma::var(truckUtilizeVec)) <<endl;

  //Write the results in two csv files:
	ofstream  policyInfo;
	policyInfo.open("C:\\Academic\\policyInfo.csv", ios::trunc);

	policyInfo << "meanReward" << ";" << "sdReward" << ";" << "meanFeed" << ";" << "sdFeed" <<  ";" << "meanavgWeight" << ";" << "sdavgWeight" << ";" << "meancycle" << ";" << "sdcycle" 
		<< ";" << "meancycleLength" << ";" << "sdcycleLength" << ";" << "meantruckNum" << ";" << "sdtruckNum" << ";" << "meanmrketNum" << ";" << "sdmrketNum" <<";" << "meantruckUtilize" << ";" << "sdtruckUtilize" << endl;

	policyInfo << arma::mean(rewVec) << ";" << sqrt(arma::var(rewVec)) << ";" << arma::mean(feedVec) << ";" << sqrt(arma::var(feedVec)) << ";" << arma::mean(avgWeightVec) << ";" 
		<< sqrt(arma::var(avgWeightVec)) << ";" << arma::mean(cycleVec) << ";" << sqrt(arma::var(cycleVec)) << ";" << arma::mean(cycleLengthVec) << ";" << sqrt(arma::var(cycleLengthVec)) << ";" 
		<< arma::mean(truckVec) << ";" << sqrt(arma::var(truckVec)) << ";" << arma::mean(marketVec) << ";" << sqrt(arma::var(marketVec)) << ";"
		<< arma::mean(truckUtilizeVec) << ";" << sqrt(arma::var(truckUtilizeVec)) << endl;

  }

  system("PAUSE");

}
