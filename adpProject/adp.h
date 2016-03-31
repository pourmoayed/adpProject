#ifndef ADP_HPP
#define ADP_HPP


#include <armadillo>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "param.h"
#include "max.h"
#include "samplePath.h"
#include "cumNorm.h"
using namespace std;

// ===================================================

/**
* Class for ADP algorithm
*
* This class includes all the functions and variables that we need to implement the ADP.
*
* @author Reza Pourmoayed (rpourmoayed@econ.au.dk)
*/
class ADP
{
	  
public:  // methods

    /** Constructor. For the input arguments see file adp.cpp     
	*/
	ADP(); 

   /** Run the ADP using the parameters defined in \code{setParameters} by te value iteration method.
    *
    *  @return return the approximate value functions for all states.
    *  @author Reza Pourmoayed (rpourmoayed@econ.au.dk)
    *  @export
    */
   void valueADP();

   /** Run the ADP using the parameters defined in \code{setParameters} by te policy iteration method. This function can be completed in future.
    *
    *  @return return the policy including the best action for each state.
    *  @author Reza Pourmoayed (rpourmoayed@econ.au.dk)
    *  @export
    */
   void policyADP();


  /** Initialize slopes (\var(slope[t][iL][iG][iSd]) )
   */
  void initializeValue();

  /** Initialize state values (\var(state[i][j]) ) for the first decision epoch randomly
   */
  void initializeState();

  /** Find the index of a given state based on the discritization matrix .
 * 
 * @param st A value in a given discretized matrix
 * @param dis A matrix containing discretized values of a continuous state variable.  
 *
 * @return The index number of a specefic interval in dis that includes st. 
 */
  int findIndex(double st, arma::mat dis);

  /** Calculate the coeficiants of reward and cost for the kth pig and store them in arrays \var(revC[t][iL][iG][iSd]) and \var(costC[t][iL][iG][iSd]) 
   */
  void rewCoefficients();

  /** Upate the slopes of the value function in each iteration of the algorithm using pre-decision state method 
   */
  void updaeSlope(vector < vector <arma::vec> > & stateVec, vector < vector < vector < vector<double> > > >  slopeVec, double objective, double alphaUpdate, arma::mat expSlope);

  /** Upate the slopes of the value function in each iteration of the algorithm using post-decision state method
   */
  void updaeSlopePost(vector < vector <arma::vec> > & stateVec, vector < vector < vector < vector<double> > > >  slopeVec, double objective, double alphaUpdate, arma::mat expSlope);

   /** Transition function to define the values of states in the next decision epoch. This function is impelemented in each iteration of the algorithm.  
   */
  void transState();

   /** Compute the maximum diff between observed slopes considered as the stop criterion.   
   */
  void stopCriteria();


     /** Calculate the log transformed probabilities Pr(m_t+1 | m_t, t).
    * 
    *  Values are stored in the 7-dim vector \var(prM[t][iRS][iSWt][iSGt][iSW][iSG][iR]).
    */
   void CalcTransPrM();
   
   
     /** Calculate the log transformed probabilities Pr(var_t+1 | var_t, t, n).
    * 
    *  Values are stored in the 3-dim vector \var(prSd[t][iSSdt][iSSd]).
    */
   void CalcTransPrSd();


   /** Read the initial slope values from the related csv files used in value iteration algorithm.
    */
   void readIniSlopes();

    /** Read the final slope values from the related csv files used in the maximization problem to find best actions.
    */
   void readFinalSlopes();


   /** Store the final slope values in the related csv files
    */
   void storeSlopes();


   /** find the best actions using the given slope values by solving a maximization problem. 
    */
   void optActions(); 


 // variables
   int sizeSL;
   int sizeSG;
   int sizeSSd;
   double maxDiffSlope;
   arma::mat expSlope;

   double RewardWeek, FeedWeek, avgMarkWeight, cycleNumber, cycleLength;  //used for comparing policies
   int sumTruckNum, sumMarket; //used for comparing policies

 
   vector < vector <arma::vec> > state; // state[i][j]  state vector.
   vector < vector < vector < vector<double> > > > slope; // slope[t][iL][iG][iSsd] slope of the state components in the approximate value function.
   vector < vector < vector < vector<arma::vec> > > > slopeVal; // slope[t][iL][iG][iSsd] history of slopes calculated using approximate value iteration algorithm.
   vector < vector < vector < vector<arma::rowvec> > > > revC; // revC[t][iL][iG][iSsd]  coeficiants related to revenue of selling for an unselected pen.
   vector < vector < vector < vector<arma::rowvec> > > > costC; // costC[t][iL][iG][iSsd]  coeficiants related to cost of feeding for an unselected pen.
   vector < vector < vector < vector<arma::rowvec> > > > lWeight; // lWeight[t][iL][iG][iSsd]  live weight of the sorted pigs.

   vector< vector< vector<double> > > prSd; // prSd[t][iSSdt][iSSd] log probability (var_t+1|var_t, t, n)
   vector < vector < vector < vector < vector<double> > > > > prM; 
   private:
   
   int iL0, iG0, iSd0;
   arma::mat iLPost;
   arma::mat iGPost; 
   arma::mat iSdPost; 
   arma::mat tPost;
   // define an object of the required classes:
   samplePath sample;
   double alphaDummy;
   PARAM p;
   MAX opt;
};
#endif
