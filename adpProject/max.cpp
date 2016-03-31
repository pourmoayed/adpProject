

#include "max.h"
typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<NumVarMatrix>   NumVar3Matrix;
typedef IloArray<IloNumArray> NumMatrix;


void MAX::maxCplex( vector < vector < vector < vector<arma::rowvec> > > > & rev, vector < vector < vector < vector<arma::rowvec> > > > & cost,
	   vector < vector <arma::vec> > & stateVec, vector < vector < vector < vector<double> > > > &  slopeVec, arma::mat & expSlope){

	   IloEnv   env;
	   IloInt i, j, k, t, iL, iG, iSd;
	   IloModel model(env);
	   IloNumArray sumSec(env,p.secNum);
     
      //define dicision variables 
      NumVar3Matrix x(env,p.secNum);
      for(i=0; i<p.secNum; i++){
    	  x[i] = NumVarMatrix(env,p.penNum);
    	  for(j=0; j<p.penNum; j++){
     		  x[i][j] = IloNumVarArray(env,stateVec[i][j][4]+1);
    		  for(k=0; k<=stateVec[i][j][4]; k++){
    			  x[i][j][k] = IloNumVar(env,0,1,ILOINT); //binary the for number of remaining pigs in each pen
	    	  }
       	  }
       }
      IloNumVar y(env,0,p.truckNum,ILOINT); // integer for the number of trucks
      IloNumVarArray g(env,p.secNum, 0, 1, ILOINT); // binary for terminating decisions


	  IloExpr obj1(env);
      IloExpr obj2(env);
      IloExpr obj3(env);
      IloExpr obj4(env);
	  IloExpr conExp1(env);


	  for(i=0; i<p.secNum; i++)
	  {
		  obj3 += g[i]*p.pigletPrice*p.pigs*p.penNum; //piglet cos
		  IloExpr conExp2(env);
		  sumSec[i]=0;
		  for(j=0; j<p.penNum; j++)
		  {			  
			  iL = stateVec[i][j][0]; iG = stateVec[i][j][1]; iSd = stateVec[i][j][2]; t = stateVec[i][j][3];
			  
			  if(p.adpPost){
				  if( stateVec[i][j][3] == (p.tMax-1) ) obj4 += g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; else 
				  obj4 += (stateVec[i][j][4] - IloSum(x[i][j]) ) * slopeVec[t+1][iL][iG][iSd]  +    g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; //value function
			  }
			  if(!p.adpPost){
				  if( stateVec[i][j][3] == (p.tMax-1) ) obj4 += g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; else 
				  obj4 += (stateVec[i][j][4] - IloSum(x[i][j]) ) * expSlope(i,j)  +    g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; //value function			  			  
				  //cout<<" expSlope(i,j) Max " << expSlope(i,j) << endl;
			  }
			  			  
			  if( stateVec[i][j][3]==(p.tMax-1) ) model.add(g[i]==1);
			  model.add( x[i][j][0]==0 );
			  sumSec[i] = sumSec[i] + stateVec[i][j][4];
			  for(k=0; k<=stateVec[i][j][4]; k++)
			  {
				  obj1 += ( x[i][j][k] * rev[t][iL][iG][iSd][k] ) - ( (1-x[i][j][k]) * cost[t][iL][iG][iSd][k] ); //net profit
				  conExp1 += x[i][j][k];
				  if( (stateVec[i][j][3]>=0) & (stateVec[i][j][3]<p.tStartMarketing-1) ) model.add(x[i][j][k]==0);
				  if( (k>0) ) model.add( x[i][j][k-1] <= x[i][j][k] );
				  conExp2 = x[i][j][k] + conExp2;				  
			  }
		  }
		  model.add( sumSec[i]  - conExp2 <= 1000000000*(1-g[i]) ); // constrints imply an empty pen after termination 
		  model.add( sumSec[i]  - conExp2 >= 1-g[i] ); // constrints imply an empty pen after termination
		  conExp2.end();
	  }

	  obj2 = y*p.truckCost; //transport cost
	  model.add( IloMaximize(env, obj1 - obj2 - obj3 + p.discount*obj4 ) );
	  model.add(conExp1 <= p.truckCapacity * y );
	  conExp1.end();

	  // Solve the model
	  IloCplex cplex(env);
      cplex.extract(model);
	  cplex.setOut(env.getNullStream());
      cplex.solve();
	  IloCplex clearModel();
	  //cplex.out() << "Solution status: " << cplex.getStatus() << endl;	  

	  termSec.zeros();
	  numPigPen.zeros();
	  truckNeed = 0;

	  truckNeed = cplex.getValue(y);
	  objective = cplex.getObjValue(); 
	  reward = cplex.getValue(obj1 - obj2 - obj3);
	  for(i=0; i<p.secNum; i++){
		  termSec[i] =  cplex.getValue(g[i]);
		  for(j=0; j<p.penNum; j++){
    			 numPigPen(i,j) = cplex.getValue(IloSum(x[i][j]));
		  }
	  }

	  obj1.end();
	  obj2.end();
      obj3.end();
      obj4.end();

	  env.end();   
   }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void MAX::maxCplexRelax( vector < vector < vector < vector<arma::rowvec> > > > & rev, vector < vector < vector < vector<arma::rowvec> > > > & cost,
	   vector < vector <arma::vec> > & stateVec, vector < vector < vector < vector<double> > > > &  slopeVec, arma::mat & expSlope){

	   
	   //IloModel model(envRelax);
	   //modelRelax(model, rev, cost, stateVec, slopeVec);

	   IloInt i, j, k, t, iL, iG, iSd;
	   IloModel model(envRelax);
	   IloNumArray sumSec(envRelax,p.secNum);
     
      //define dicision variables 
      NumVar3Matrix x(envRelax,p.secNum);
      for(i=0; i<p.secNum; i++){
    	  x[i] = NumVarMatrix(envRelax,p.penNum);
    	  for(j=0; j<p.penNum; j++){
     		  x[i][j] = IloNumVarArray(envRelax,stateVec[i][j][4]+1);
    		  for(k=0; k<=stateVec[i][j][4]; k++){
    			  x[i][j][k] = IloNumVar(envRelax,0,1,ILOINT); //binary, float: the for number of remaining pigs in each pen
	    	  }
       	  }
       }
      IloNumVar y(envRelax,0,p.truckNum,ILOINT); // integer, float: for the number of trucks
      IloNumVarArray g(envRelax,p.secNum, 0, 1,ILOINT); // binary, float: for terminating decisions


	  IloExpr obj1(envRelax);
      IloExpr obj2(envRelax);
      IloExpr obj3(envRelax);
      IloExpr obj4(envRelax);
	  IloExpr conExp1(envRelax);


	  for(i=0; i<p.secNum; i++)
	  {
		  obj3 += g[i]*p.pigletPrice*p.pigs*p.penNum; //piglet cos
		  IloExpr conExp2(envRelax);
		  sumSec[i]=0;
		  for(j=0; j<p.penNum; j++)
		  {			  
			  iL = stateVec[i][j][0]; iG = stateVec[i][j][1]; iSd = stateVec[i][j][2]; t = stateVec[i][j][3];

			  if(p.adpPost){
				  if( stateVec[i][j][3] == (p.tMax-1) ) obj4 += g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; else  
				  obj4 += (stateVec[i][j][4] - IloSum(x[i][j]) ) * slopeVec[t+1][iL][iG][iSd]  +    g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; //value function 
			  }
			  if(!p.adpPost){
				  if( stateVec[i][j][3] == (p.tMax-1) ) obj4 += g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; else 
				  obj4 += (stateVec[i][j][4] - IloSum(x[i][j]) ) * expSlope(i,j)  +    g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; //value function			  			  
			  }
			  			  
			  if( stateVec[i][j][3]==(p.tMax-1) ) model.add(g[i]==1);
			  model.add( x[i][j][0]==0 );
			  sumSec[i] = sumSec[i] + stateVec[i][j][4];
			  for(k=0; k<=stateVec[i][j][4]; k++)
			  {
                  //obj1 += ( x[i][j][k] * rev[t][iL][iG][iSd][k] ) - ( (1-x[i][j][k]) * cost[t][iL][iG][iSd][k] ); //net profit
                  obj1 += ( x[i][j][k] * arma::mean(rev[t][iL][iG][iSd]) ) - ( (1-x[i][j][k]) * arma::mean(cost[t][iL][iG][iSd]) ); //net profit
				  conExp1 += x[i][j][k];
				  if( (stateVec[i][j][3]>=0) & (stateVec[i][j][3]<p.tStartMarketing-1) ) model.add(x[i][j][k]==0);
				  if( (k>0) ) model.add( x[i][j][k-1] <= x[i][j][k] );
				  conExp2 = x[i][j][k] + conExp2;				  
			  }
		  }
		  model.add( sumSec[i]  - conExp2 <= 1000000000*(1-g[i]) ); // constrints imply an empty pen after termination 
		  model.add( sumSec[i]  - conExp2 >= 1-g[i] ); // constrints imply an empty pen after termination
		  conExp2.end();
	  }

	  obj2 = y*p.truckCost; //transport cost
	  model.add( IloMaximize(envRelax, obj1 - obj2 - obj3 + p.discount*obj4 ) );
	  model.add(conExp1 <= p.truckCapacity * y );
	  conExp1.end();

	  obj1.end();
	  obj2.end();
      obj3.end();
      obj4.end();

	  // Solve the model
	  IloCplex cplex(envRelax);
      cplex.extract(model);
	  cplex.setOut(envRelax.getNullStream());
      cplex.solve();
	  IloCplex clearModel();
	  //cplex.out() << "Solution status: " << cplex.getStatus() << endl;	  

	  objectiveRelax = cplex.getObjValue(); 

	   //envRelax.end();       
   }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void MAX::maxMyopic( vector < vector < vector < vector<arma::rowvec> > > > & rev, vector < vector < vector < vector<arma::rowvec> > > > & cost,
	   vector < vector < vector < vector<arma::rowvec> > > > & lW, vector < vector <arma::vec> > & stateVec, vector < vector < vector < vector<double> > > > &  slopeVec, arma::mat & expSlope, double & ths){

	   IloEnv   env;
	   IloInt i, j, k, t, iL, iG, iSd;
	   IloModel model(env);
	   IloNumArray sumSec(env,p.secNum);
	   int sumHerd=0;
     
      //define dicision variables 
      NumVar3Matrix x(env,p.secNum);
      for(i=0; i<p.secNum; i++){
    	  x[i] = NumVarMatrix(env,p.penNum);
    	  for(j=0; j<p.penNum; j++){
     		  x[i][j] = IloNumVarArray(env,stateVec[i][j][4]+1);
    		  for(k=0; k<=stateVec[i][j][4]; k++){
    			  x[i][j][k] = IloNumVar(env,0,1,ILOINT); //binary the for number of remaining pigs in each pen
	    	  }
       	  }
       }
      IloNumVar y(env,0,p.truckNum,ILOINT); // integer for the number of trucks
      IloNumVarArray g(env,p.secNum, 0, 1, ILOINT); // binary for terminating decisions

	  //Calculate total number of pigs at herd that can be marketed 
	  for(int ii=0; ii<p.secNum; ii++){
		  for(int jj=0; jj<p.penNum; jj++){
			  if( stateVec[ii][jj][3]>=(p.tStartMarketing-1) ) sumHerd += stateVec[ii][jj][4];
		  }
	  }


	  IloExpr obj1(env);
      IloExpr obj2(env);
      IloExpr obj3(env);
      IloExpr obj4(env);
	  IloExpr conExp1(env);


	  for(i=0; i<p.secNum; i++)
	  {
		  obj3 += g[i]*p.pigletPrice*p.pigs*p.penNum; //piglet cos
		  IloExpr conExp2(env);
		  sumSec[i]=0;
		  for(j=0; j<p.penNum; j++)
		  {			  
			  iL = stateVec[i][j][0]; iG = stateVec[i][j][1]; iSd = stateVec[i][j][2]; t = stateVec[i][j][3];
			  
			  if(p.adpPost){
				  if( stateVec[i][j][3] == (p.tMax-1) ) obj4 += g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; else 
				  obj4 += (stateVec[i][j][4] - IloSum(x[i][j]) ) * slopeVec[t+1][iL][iG][iSd]  +    g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; //value function
			  }
			  if(!p.adpPost){
				  if( stateVec[i][j][3] == (p.tMax-1) ) obj4 += g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; else 
				  obj4 += (stateVec[i][j][4] - IloSum(x[i][j]) ) * expSlope(i,j)  +    g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; //value function			  			  
				  //cout<<" expSlope(i,j) Max " << expSlope(i,j) << endl;
			  }
			  
			  if( stateVec[i][j][3]==(p.tMax-1) ) model.add(g[i]==1);			  
			  model.add( x[i][j][0]==0 );
			  sumSec[i] = sumSec[i] + stateVec[i][j][4];

			  //constraints for policy (all in all out)
			  //if( stateVec[i][j][3]== 14 ) model.add(g[i]==1);

			  for(k=0; k<=stateVec[i][j][4]; k++)
			  {
				  obj1 += ( x[i][j][k] * rev[t][iL][iG][iSd][k] ) - ( (1-x[i][j][k]) * cost[t][iL][iG][iSd][k] ); //net profit
				  conExp1 += x[i][j][k];
				  if( (stateVec[i][j][3]>=0) & (stateVec[i][j][3]<p.tStartMarketing-1) ) model.add(x[i][j][k]==0);				  
				  if( (k>0) ) model.add( x[i][j][k-1] <= x[i][j][k] );				  
				  //constraints for policy (threshould weight)
     			  //if( (stateVec[i][j][3]>=p.tStartMarketing-1) & (stateVec[i][j][3]!=(p.tMax-1)) & (lW[t][iL][iG][iSd][k]>ths)  ) model.add(x[i][j][k]==1);
	    		  //if( (stateVec[i][j][3]>=p.tStartMarketing-1) & (stateVec[i][j][3]!=(p.tMax-1)) & (lW[t][iL][iG][iSd][k]<=ths)  ) model.add(x[i][j][k]==0);

				  //constraints for policy (all in all out)
				  //if( (stateVec[i][j][3]>=p.tStartMarketing-1) & (stateVec[i][j][3]!=(p.tMax-1)) & (stateVec[i][j][3]==(p.tMax-1))  ) model.add(x[i][j][k]==1);
	    		  //if( (stateVec[i][j][3]>=p.tStartMarketing-1) & (stateVec[i][j][3]!=(p.tMax-1)) & (stateVec[i][j][3]!= 14 ) ) model.add(x[i][j][k]==0);

				  conExp2 = x[i][j][k] + conExp2;				  
			  }
		  }

		  //constraints for myopic policy:
		  //if(sumSec[i]<=1) model.add(g[i]==1);

		  model.add( sumSec[i]  - conExp2 <= 1000000000*(1-g[i]) ); // constrints imply an empty pen after termination 
		  model.add( sumSec[i]  - conExp2 >= 1-g[i] ); // constrints imply an empty pen after termination
		  conExp2.end();
	  }

	  obj2 = y*p.truckCost; //transport cost
	  model.add( IloMaximize(env, obj1 - obj2 - obj3  + p.discount*obj4 ) ); // -100000000*(p.truckCapacity * y - conExp1)   ) );   //  -10000000000000000000
	  //model.add(conExp1 <= p.truckCapacity * y );
	  
	  //constraints for new myopic policy (FTC):
	  int tSec1, tSec2, tSec3;
	  tSec1 = stateVec[0][0][3]; tSec2 = stateVec[1][0][3]; tSec3 = stateVec[2][0][3];

	  if( ( ( tSec1==(p.tMax-1) || tSec1<(p.tStartMarketing-1) ) & ( tSec2==(p.tMax-1) || tSec2<(p.tStartMarketing-1) ) & ( tSec3==(p.tMax-1) || tSec3<(p.tStartMarketing-1) ) ) 
		  || ( ( (double)sumHerd / (double)p.truckCapacity )< 1 ) ) {model.add(conExp1 == sumHerd); if(sumHerd==0) model.add(y==0); else model.add(y==1);}  else {model.add(conExp1 == p.truckCapacity * y ); }

	  //if( ( (double)sumHerd / (double)p.truckCapacity )>=1 ) model.add(conExp1 == p.truckCapacity * y );
	  //cout<<(double)sumHerd / (double)p.truckCapacity<<endl;
	  //if( 0<=( (double)sumHerd / (double)p.truckCapacity )< 1 ) model.add(conExp1 <= p.truckCapacity * y );
	  
	  conExp1.end();

	  // Solve the model
	  IloCplex cplex(env);
      cplex.extract(model);
	  cplex.setOut(env.getNullStream());
      cplex.solve();
	  IloCplex clearModel();
	  //cplex.out() << "Solution status: " << cplex.getStatus() << endl;	  

	  termSec.zeros();
	  numPigPen.zeros();
	  truckNeed = 0;

	  truckNeed = cplex.getValue(y);
	  cplex.out() << "cplex.getValue(y) " << cplex.getValue(y) << endl;
	  objective = cplex.getObjValue(); 
	  reward = cplex.getValue(obj1 - obj2 - obj3);
	  for(i=0; i<p.secNum; i++){
		  termSec[i] =  cplex.getValue(g[i]);
		  for(j=0; j<p.penNum; j++){
    			 numPigPen(i,j) = cplex.getValue(IloSum(x[i][j]));
		  }
	  }

	  obj1.end();
	  obj2.end();
      obj3.end();
      obj4.end();

	  env.end();   
   }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






	   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   



















   
   
   
   
   
   void  MAX::slopesUpdate(vector < vector < vector < vector<arma::rowvec> > > > & rev, vector < vector < vector < vector<arma::rowvec> > > > & cost,
	 vector < vector <arma::vec> > & stateVec, vector < vector < vector < vector<double> > > > & slopeVec, double & alphaUpdate, arma::mat & expSlope){

	   IloEnv   envR;
	   IloInt i, j, iR, jR, k, t, iL, iG, iSd;
	   IloModel model(envR);
	   IloNumArray sumSec(envR,p.secNum);
	   int rNew; // new value of number of pig in a specific pen- for finding the deviation
	   double extraValue; // extra value of objective function by adding a pig to a specific prn, for finding the deviation
	   double deviation; // deviation of the objective values.
	   double objectiveMain; // value of objective function of the main model
     
      //define dicision variables 
      NumVar3Matrix x(envR,p.secNum);
      for(i=0; i<p.secNum; i++){
    	  x[i] = NumVarMatrix(envR,p.penNum);
    	  for(j=0; j<p.penNum; j++){
     		  x[i][j] = IloNumVarArray(envR,p.pigs+1);
    		  for(k=0; k<=p.pigs; k++){
    			  x[i][j][k] = IloNumVar(envR,0,1,ILOFLOAT); //binary the for number of remaining pigs in each pen
	    	  }
       	  }
       }

	  NumVar3Matrix z(envR,p.secNum);
      for(i=0; i<p.secNum; i++){
    	  z[i] = NumVarMatrix(envR,p.penNum);
    	  for(j=0; j<p.penNum; j++){
     		  z[i][j] = IloNumVarArray(envR,p.pigs+1);
    		  for(k=0; k<=p.pigs; k++){
    			  z[i][j][k] = IloNumVar(envR,0,1,ILOFLOAT); //binary the for number of remaining pigs in each pen
	    	  }
       	  }
       }

      IloNumVar y(envR,0,p.truckNum,ILOFLOAT); // integer for the number of trucks
      IloNumVarArray g(envR,p.secNum, 0, 1, ILOFLOAT); // binary for terminating decisions

	  IloExpr obj1(envR);
      IloExpr obj2(envR);
      IloExpr obj3(envR);
      IloExpr obj4(envR);
	  IloExpr conExp1(envR);
	  IloExprArray conExp2(envR,p.secNum);
	  for (int i = 0; i < p.secNum; i++)
       conExp2[i] = IloExpr(envR);



	  for(i=0; i<p.secNum; i++)
	  {
		  obj3 += g[i]*p.pigletPrice*p.pigs*p.penNum; //piglet cost
		  //IloExpr conExp2(envR);
		  sumSec[i]=0;
		  for(j=0; j<p.penNum; j++)
		  {	
			  iL = stateVec[i][j][0]; iG = stateVec[i][j][1]; iSd = stateVec[i][j][2]; t = stateVec[i][j][3];
			  
			  if(p.adpPost){
				  if( stateVec[i][j][3] == (p.tMax-1) ) obj4 += g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; else 
				  obj4 += (stateVec[i][j][4] - IloSum(x[i][j]) ) * slopeVec[t+1][iL][iG][iSd]  +    g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; //value function
			  }
			  if(!p.adpPost){
				  if( stateVec[i][j][3] == (p.tMax-1) ) obj4 += g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; else 
				  obj4 += (stateVec[i][j][4] - IloSum(x[i][j]) ) * expSlope(i,j)  +    g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; //value function			  			  
			  }
			  			  
			  if( stateVec[i][j][3]==(p.tMax-1) ) model.add(g[i]==1);
			  model.add( x[i][j][0]==0 ); model.add( z[i][j][0]==0 );
			  model.add( IloSum(x[i][j]) + IloSum(z[i][j]) == stateVec[i][j][4] );
			  sumSec[i] = sumSec[i] + stateVec[i][j][4];
			  for(k=0; k<=p.pigs; k++)
			  {
				  obj1 += ( x[i][j][k] * rev[t][iL][iG][iSd][k] ) - ( (z[i][j][k]) * cost[t][iL][iG][iSd][k] ); //net profit
				  conExp1 += x[i][j][k];
				  if( (stateVec[i][j][3]>=0) & (stateVec[i][j][3]<p.tStartMarketing-1) ) model.add(x[i][j][k]==0); 
				  if( (k>0) & ( k <= stateVec[i][j][4]) ) model.add( x[i][j][k-1] - x[i][j][k] <= 0  );
				  if(  k > stateVec[i][j][4] ) {model.add( x[i][j][k] == 0); model.add( z[i][j][k] == 0) ;} 
				  conExp2[i] = x[i][j][k] + conExp2[i] ;				  
			  }
		  }
		  model.add( sumSec[i]  - conExp2[i] <= 1000000000*(1-g[i]) ); // constrints imply an empty pen after termination 
		  model.add( sumSec[i]  - conExp2[i] >= 1-g[i] ); // constrints imply an empty pen after termination
		  //conExp2.end();
	  }

	  model.add(conExp1 <= p.truckCapacity * y );	  
	  obj2 = y*p.truckCost; //transport cost
	  model.add( IloMaximize(envR, obj1 - obj2 - obj3 + p.discount*obj4 ) );
	  //conExp1.end();
	  //conExp2.end();
	  //obj1.end();
	  //obj2.end();
   //   obj3.end();
   //   obj4.end();

	  // Solve the main model
	  IloCplex cplex(envR);
      cplex.extract(model);
	  cplex.setOut(envR.getNullStream());
      cplex.solve();
	  objectiveMain = cplex.getObjValue(); 


	  //Change the number of pigs in pen jR of section iR
	  for(iR=0; iR<p.secNum; iR++){
		  for(jR=0; jR<p.penNum; jR++){
			  if( (stateVec[iR][jR][4]!=p.pigs) & (stateVec[iR][jR][3]!=p.tMax-1)   ){
			  rNew = stateVec[iR][jR][4] + 1;
			  extraValue = 0;

			  //modify the model
			  model.remove( x[iR][jR][rNew] == 0); model.remove( z[iR][jR][rNew] == 0 ); 
			  model.add(x[iR][jR][rNew-1] - x[iR][jR][rNew] <= 0);
			   
			  model.remove( IloSum(x[iR][jR]) + IloSum(z[iR][jR]) == stateVec[iR][jR][4] );
			  model.add( IloSum(x[iR][jR]) + IloSum(z[iR][jR]) == stateVec[iR][jR][4] + 1  );

			  model.remove( sumSec[iR]  - conExp2[iR] <= 1000000000*(1-g[iR]) ); 
		      model.remove( sumSec[iR]  - conExp2[iR] >= 1-g[iR] );
			  //cout<< "rNew" << rNew <<endl;
			  //cout<< "sumSec[iR]" << sumSec[iR]<<endl;
			  sumSec[iR] = sumSec[iR] + 1;
			  //cout<< "sumSec[iR]" << sumSec[iR]<<endl;
			  model.add( sumSec[iR]  - conExp2[iR] <= 1000000000*(1-g[iR]) ); 
		      model.add( sumSec[iR]  - conExp2[iR] >= 1-g[iR] );

			  iL = stateVec[iR][jR][0]; iG = stateVec[iR][jR][1]; iSd = stateVec[iR][jR][2]; t = stateVec[iR][jR][3];
			  if( stateVec[iR][jR][3] != (p.tMax-1) ) extraValue = slopeVec[t+1][iL][iG][iSd]; //obj4+=slopeVec[t+1][iL][iG][iSd];

			  //solve new model 
               cplex.extract(model);
			   cplex.setOut(envR.getNullStream());
			   cplex.solve();
			   //update slopes
			   deviation = cplex.getObjValue() + extraValue - objectiveMain;
			   slopeVec[t][iL][iG][iSd] = alphaUpdate * slopeVec[t][iL][iG][iSd] + (1-alphaUpdate) * deviation;
			   //cout<<" slopeVec[t][iL][iG][iSd] "<<slopeVec[t][iL][iG][iSd]<<endl;

			   //change the model to the orginal one
			   model.add( x[iR][jR][rNew] == 0); model.add( z[iR][jR][rNew] == 0 ); 
			   model.remove(x[iR][jR][rNew-1] - x[iR][jR][rNew] <= 0);
			   
			   model.remove( IloSum(x[iR][jR]) + IloSum(z[iR][jR]) == stateVec[iR][jR][4] + 1 );
			   model.add( IloSum(x[iR][jR]) + IloSum(z[iR][jR]) == stateVec[iR][jR][4] );

			   model.remove( sumSec[iR]   - conExp2[iR] <= 1000000000*(1-g[iR]) ); 
		       model.remove( sumSec[iR]  - conExp2[iR] >= 1-g[iR] );
			   //cout<< "sumSec[iR]" << sumSec[iR]<<endl;
			   sumSec[iR] = sumSec[iR] - 1;
			   //cout<< "sumSec[iR]" << sumSec[iR]<<endl;
			   model.add( sumSec[iR]  - conExp2[iR] <= 1000000000*(1-g[iR]) ); 
		       model.add( sumSec[iR]  - conExp2[iR] >= 1-g[iR] );			   
			  
			  }

		  }
	  }

envR.end();

 }

 //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

   void MAX::maxCplexUpdate( vector < vector < vector < vector<arma::rowvec> > > > & rev, vector < vector < vector < vector<arma::rowvec> > > > & cost,
	   vector < vector <arma::vec> > & stateVec, vector < vector < vector < vector<double> > > > &  slopeVec, int & trucksNum, double & alphaUpdate, arma::mat & expSlope, arma::mat & numPigs){

	   IloEnv   env;
	   IloInt i, j, k, t, iL, iG, iSd;
	   IloModel model(env);
	   IloNumArray sumSec(env,p.secNum);
     
      //define dicision variables 
      NumVar3Matrix x(env,p.secNum);
      for(i=0; i<p.secNum; i++){
    	  x[i] = NumVarMatrix(env,p.penNum);
    	  for(j=0; j<p.penNum; j++){
     		  x[i][j] = IloNumVarArray(env,p.pigs+1);
    		  for(k=0; k<=p.pigs; k++){
    			  x[i][j][k] = IloNumVar(env,0,1,ILOFLOAT); //binary the for number of remaining pigs in each pen
	    	  }
       	  }
       }

	  //define dicision variables 
      NumVar3Matrix z(env,p.secNum);
      for(i=0; i<p.secNum; i++){
    	  z[i] = NumVarMatrix(env,p.penNum);
    	  for(j=0; j<p.penNum; j++){
     		  z[i][j] = IloNumVarArray(env,p.pigs+1);
    		  for(k=0; k<=p.pigs; k++){
    			  z[i][j][k] = IloNumVar(env,0,1,ILOFLOAT); //binary the for number of remaining pigs in each pen
	    	  }
       	  }
       }

      //IloNumVar y(env,0,p.truckNum,ILOFLOAT); // integer for the number of trucks
      IloNumVarArray g(env,p.secNum, 0, 1, ILOFLOAT); // binary for terminating decisions


	  IloExpr obj1(env);
      IloExpr obj2(env);
      IloExpr obj3(env);
      IloExpr obj4(env);
	  IloExpr conExp1(env);


	  for(i=0; i<p.secNum; i++)
	  {
		  obj3 += g[i]*p.pigletPrice*p.pigs*p.penNum; //piglet cos
		  IloExpr conExp2(env);
		  sumSec[i]=0;
		  for(j=0; j<p.penNum; j++)
		  {			  
			  iL = stateVec[i][j][0]; iG = stateVec[i][j][1]; iSd = stateVec[i][j][2]; t = stateVec[i][j][3];
			  
			  if(p.adpPost){
				  if( stateVec[i][j][3] == (p.tMax-1) ) obj4 += g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; else 
				  obj4 += (stateVec[i][j][4] - IloSum(x[i][j]) ) * slopeVec[t][iL][iG][iSd]  +    g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; //value function
			  }
			  if(!p.adpPost){
				  if( stateVec[i][j][3] == (p.tMax-1) ) obj4 += g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; else 
				  obj4 += (stateVec[i][j][4] - IloSum(x[i][j]) ) * expSlope(i,j)  +    g[i] * p.pigs * slopeVec[0][p.sL0][p.sG0][p.sSd0]; //value function			  			  
			  }
			  
			  model.add( IloSum(x[i][j]) == numPigs(i,j) );

			  if( stateVec[i][j][3]==(p.tMax-1) ) model.add(g[i]==1);
			  model.add( x[i][j][0]==0 ); model.add( z[i][j][0]==0 );
			  sumSec[i] = sumSec[i] + stateVec[i][j][4];
			  for(k=0; k<=p.pigs; k++)
			  {
				  obj1 += ( x[i][j][k] * rev[t][iL][iG][iSd][k] ) - ( (z[i][j][k]) * cost[t][iL][iG][iSd][k] ); //net profit
				  conExp1 += x[i][j][k];
				  if( (stateVec[i][j][3]>=0) & (stateVec[i][j][3]<p.tStartMarketing-1) ){ model.add(x[i][j][k]==0); }
				  if( (k>0) & ( k <= stateVec[i][j][4]) ) model.add( x[i][j][k-1] - x[i][j][k] <= 0  );
				  if(  k > stateVec[i][j][4] ) {model.add( x[i][j][k] == 0); model.add( z[i][j][k] == 0) ;}
				  conExp2 = x[i][j][k] + conExp2;				  
			  }
		  }
		  model.add( sumSec[i]  - conExp2 <= 1000000000*(1-g[i]) ); // constrints imply an empty pen after termination 
		  model.add( sumSec[i]  - conExp2 >= 1-g[i] ); // constrints imply an empty pen after termination
		  conExp2.end();
	  }

	  //obj2 = y*p.truckCost; //transport cost
	  model.add( IloMaximize(env, obj1 - trucksNum*p.truckCost  - obj3 + p.discount*obj4 ) );
	  //model.add(conExp1 <= p.truckCapacity * y );
	  model.add(conExp1 <= p.truckCapacity *trucksNum  );
	  conExp1.end();

	  obj1.end();
	  obj2.end();
      obj3.end();
      obj4.end();


	  //We add extra constraints to find relevant dual variables
	  IloRangeArray c(env);
	  for(i=0; i<p.secNum; i++){
		  for(j=0; j<p.penNum; j++){
			  c.add( IloSum(x[i][j]) + IloSum(z[i][j]) == stateVec[i][j][4] );
		  }
	  }
	  model.add(c);


	  // Solve the model
	  IloCplex cplex(env);
      cplex.extract(model);
	  cplex.setOut(env.getNullStream());
      cplex.solve();
	  IloCplex clearModel();
	  //cplex.out() << "Solution status: " << cplex.getStatus() << endl;	  

	  //objective = cplex.getObjValue(); 

	  //cout<<"step1"<<endl;
	   //Extract dual from linear relaxation
	  double deviation;
	  dualPigPen.zeros();
	  IloNumArray vals(env);
	  cplex.getDuals(vals, c);
	  int xx=0;
	  int term;
	  for(i=0; i<p.secNum; i++){
		  term = cplex.getValue(g[i]);
		  for(j=0; j<p.penNum; j++){
			  iL = stateVec[i][j][0]; iG = stateVec[i][j][1]; iSd = stateVec[i][j][2]; t = stateVec[i][j][3];
			  if((!p.adpPost) & (t>=p.tStartMarketing-1) ){
			  //cout<<"step2"<<endl;
			  //dualPigPen(i,j) = vals[xx];
			  if(term!=1) deviation = vals[xx] + expSlope(i,j); else deviation = vals[xx]; 
			  slopeVec[t][iL][iG][iSd] = (1-alphaUpdate) * slopeVec[t][iL][iG][iSd] + alphaUpdate * deviation;
			  //cout<< "slopeVec " << slopeVec[t][iL][iG][iSd]<<endl; 
			  //cout<< " iL: " << iL << " iG: " << iG << " iSd: " << iSd << " t: " << t <<  " slope1 " <<slopeVec[t][iL][iG][iSd]<<endl;
			  env.out() << "duals: " <<i <<j << " : " << vals[xx]  << endl;
			  //slopeVec[t+1][iL][iG][iSd]
			  //cout<<"step3"<<endl;
			  }
			  xx=xx+1;
		  }
	  }

	     
	  env.end();

   }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   

 