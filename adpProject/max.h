#ifndef MAX_H
#define MAX_H

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
#include <iostream>
using namespace std;
#include <vector>
#include <armadillo>
#include "param.h"

// ===================================================

/**
* In this class maximisation problems are solved to find the best marketing decisions used in the AVI algorithm and the final results.
*
* For each iteration of ADP, a maximization problem is solved and the optimal actions are reported
*
* @author Reza pourmoayed
*/
class MAX
{
public:  // methods
	
    /** Constructor. Store the parameters.
     *  
     */
   MAX(): p()  {
   termSec = arma::vec(p.secNum,arma::fill::zeros);
   numPigPen = arma::mat(p.secNum,p.penNum,arma::fill::zeros);
   }

//mhetods 

  /** optimisation problem in each iteration of ADP algorithm to find the optimal actions.
   *
   * @param rev An array of coeficiants related to the unit revenie of selling the pigs.
   * @param cost An array of coeficiants related to the unit cost of keeping the pigs.
   * @param staeVec An multi dimensinal array related to the state structure of the model.  
   * @param slopeVec A vector related to the slope vector of value function
   * @param expSlope A vector related to the expected slope values used when the AVI is run based on only pre-decision states.
   *
   * The optimal values related to the objective function and the number of pigs for selling are stored in the global 
   * variables numPigPen, termSec and objective.  
   * @author Reza Pourmoayed (rpourmoayed@econ.au.dk)
   */
   void maxCplex( vector < vector < vector < vector<arma::rowvec> > > > & rev, vector < vector < vector < vector<arma::rowvec> > > > & cost,
	   vector < vector <arma::vec> > & stateVec, vector < vector < vector < vector<double> > > > & slopeVec, arma::mat & expSlope);

   /** Find the slope of the value function by solving a set of linear models(currently this function is not used in the AVI algorithm).
   *
   * @param rev An array of coeficiants related to the unit revenie of selling the pigs.
   * @param cost An array of coeficiants related to the unit cost of keeping the pigs.
   * @param staeVec An multi dimensinal array related to the state structure of the model.  
   * @param slopeVec An multi dimensinal array related to the slopes of value function
   * @param alphaUpdate Step size needed to update the slope values
   * @param expSlope A vector related to the expected slope values used when the AVI is run based on only pre-decision states.
   *
   * The value of slope parameters are changed in "slopeVec". Note that we have reference for "slopeVec".   
   * @author Reza Pourmoayed (rpourmoayed@econ.au.dk)
   */
  void slopesUpdate(vector < vector < vector < vector<arma::rowvec> > > > & rev, vector < vector < vector < vector<arma::rowvec> > > > & cost,
	 vector < vector <arma::vec> > & stateVec, vector < vector < vector < vector<double> > > > &  slopeVec, double & alphaUpdate, arma::mat & expSlope);

 

  /** Linear relaxation of the optimization problem only for updating the slope parametes using dual values(currently this function is not used in the AVI algorithm).
   *
   * @param rev An array of coeficiants related to the unit revenie of selling the pigs.
   * @param cost An array of coeficiants related to the unit cost of keeping the pigs.
   * @param staeVec An multi dimensinal array related to the state structure of the model.  
   * @param slopeVec A vector related to the slope vector of value function
   * @param trucksNum Number of trucks
   * @param alphaUpdate step size value in the AVI algorithm
   * @param expSlope A vector related to the expected slope values used when the AVI is run based on only pre-decision states
   * @param numPigs Number of pigs at the herd
   *
   * The optimal values related to the objective function is stored in the global 
   * variable objective.  
   * @author Reza Pourmoayed (rpourmoayed@econ.au.dk)
   */
  void maxCplexUpdate( vector < vector < vector < vector<arma::rowvec> > > > & rev, vector < vector < vector < vector<arma::rowvec> > > > & cost,
	   vector < vector <arma::vec> > & stateVec, vector < vector < vector < vector<double> > > > & slopeVec, int & trucksNum, double & alphaUpdate, arma::mat & expSlope, arma::mat & numPigs);
    
  
  /** optimization problem for updating the slope parametes using the differnece between the objective functions.
   *
   * @param rev An array of coeficiants related to the unit revenie of selling the pigs.
   * @param cost An array of coeficiants related to the unit cost of keeping the pigs.
   * @param staeVec An multi dimensinal array related to the state structure of the model.  
   * @param slopeVec A vector related to the slope vector of value function
   *
   * The optimal values related to the objective function is stored in the global 
   * variable objective.  
   * @author Reza Pourmoayed (rpourmoayed@econ.au.dk)
   */
  void maxCplexRelax( vector < vector < vector < vector<arma::rowvec> > > > & rev, vector < vector < vector < vector<arma::rowvec> > > > & cost,
	   vector < vector <arma::vec> > & stateVec, vector < vector < vector < vector<double> > > > & slopeVec, arma::mat & expSlope);



   /** optimisation problem in each iteration of ADP algorithm to find the optimal actions (used for myopic policies).
   *
   * @param rev An array of coeficiants related to the unit revenie of selling the pigs.
   * @param cost An array of coeficiants related to the unit cost of keeping the pigs.
   * @param lW live weight of sorted pigs needed to evaluate the myopic policies.
   * @param staeVec An multi dimensinal array related to the state structure of the model.
   * @param expSlope A vector related to the expected slope values used when the AVI is run based on only pre-decision states
   * @param slopeVec A vector related to the slope vector of value function
   * @param ths Threshold weight values used for evaluating threshold policies.
   *
   * The optimal values related to the objective function and the number of pigs for selling are stored in the global 
   * variables numPigPen, termSec and objective.  
   * @author Reza Pourmoayed (rpourmoayed@econ.au.dk)
   */
   void maxMyopic( vector < vector < vector < vector<arma::rowvec> > > > & rev, vector < vector < vector < vector<arma::rowvec> > > > & cost,
	   vector < vector < vector < vector<arma::rowvec> > > > & lW, vector < vector <arma::vec> > & stateVec, vector < vector < vector < vector<double> > > > & slopeVec, arma::mat & expSlope, double & ths);

     
//   variables:
   arma::mat numPigPen;
   arma::mat dualPigPen;
   arma::vec termSec;
   double objective;
   double reward;
   double objectiveRelax;
   int truckNeed;
   IloEnv   envRelax;


private:
	//variables:
	PARAM p;	
};

#endif