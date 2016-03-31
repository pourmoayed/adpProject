#include "adp.h"

// ===================================================

ADP::ADP() : sample(), p(), opt() {
	sizeSL = p.sL.size();
	sizeSG = p.sG.size();
    sizeSSd = p.sSd.size();
	
	iL0 = p.sL0; //findIndex( p.avgInsWeight - p.meanWeights(0), p.dL );
	iG0 = p.sG0; //findIndex( p.avgInsG, p.dG );
	iSd0 = p.sSd0; //findIndex( p.avgInsSd - p.sdWeights(0), p.dSd );
	
	iLPost = arma::mat(p.secNum,p.penNum,arma::fill::zeros);
	iGPost = arma::mat(p.secNum,p.penNum,arma::fill::zeros);
	iSdPost = arma::mat(p.secNum,p.penNum,arma::fill::zeros);
	tPost = arma::mat(p.secNum,p.penNum,arma::fill::zeros);
	expSlope = arma::mat(p.secNum,p.penNum,arma::fill::zeros);


   revC = vector < vector < vector < vector<arma::rowvec> > > > (p.tMax,
           vector < vector < vector<arma::rowvec> > > (sizeSL,
           vector < vector<arma::rowvec> > (sizeSG,
           vector <arma::rowvec>(sizeSSd) ) ) );  // revC[t][iL][iG][iSd]

  costC = vector < vector < vector < vector<arma::rowvec> > > > (p.tMax,
           vector < vector < vector<arma::rowvec> > > (sizeSL,
           vector < vector<arma::rowvec> > (sizeSG,
           vector <arma::rowvec>(sizeSSd) ) ) );  // costC[t][iL][iG][iSd]

  lWeight = vector < vector < vector < vector<arma::rowvec> > > > (p.tMax,
           vector < vector < vector<arma::rowvec> > > (sizeSL,
           vector < vector<arma::rowvec> > (sizeSG,
           vector <arma::rowvec>(sizeSSd) ) ) );  // sWeight[t][iL][iG][iSd]


   state = vector < vector <arma::vec> >(p.secNum,
           vector <arma::vec>(p.penNum) ) ;// state[i][j]<iL,iG,iSd,t,r>


   slope = vector < vector < vector < vector<double> > > > (p.tMax,
           vector < vector < vector<double> > > (sizeSL,
           vector < vector<double> > (sizeSG,
           vector <double>(sizeSSd) ) ) );  // slope[t][iL][iG][iSd]


   slopeVal = vector < vector < vector < vector<arma::vec> > > > (p.tMax,
           vector < vector < vector<arma::vec> > > (sizeSL,
           vector < vector<arma::vec> > (sizeSG,
           vector <arma::vec>(sizeSSd) ) ) );  // slopeVal[t][iL][iG][iSd]<arma:vec>


   prSd = vector< vector< vector<double> > > (p.tMax-1, 
          vector< vector<double> >(sizeSSd, 
          vector<double>(sizeSSd) ) ); //prSd[t][iSdt][iSd]


   prM = vector < vector < vector < vector < vector<double> > > > > (p.tMax-1, 
      vector < vector < vector < vector<double> > > > (sizeSL, 
      vector < vector < vector<double> > > (sizeSG, 
      vector< vector<double> >(sizeSL,
      vector <double>(sizeSG) ) ) ) );  // prM[t][iLt][iGt][iL][iG]
}

// ===================================================

void ADP::valueADP() {
  if(!p.adpPost){ CalcTransPrSd(); CalcTransPrM(); } // || p.adpPost 
  rewCoefficients();  
  initializeValue();
  //readIniSlopes(); //for reading initial slope values from a csv file
  int h = 1;
  int n;
  double alpha; // = 0.6; //0.6; // 0.6; //= 0.1;
  while( (h<p.simNum) || (maxDiffSlope<5) ){ //stopping criteria
    initializeState();
	 alpha = (double)100/(double)(100+h-1);
	//alpha = (double)alpha/(double)(1+alpha-0.9);
	cout<<"iterNum"<<h<<endl;
	for(n=1; n<p.epochNum; n++){  //decision epochs.
		cout<<"epochNum"<<n<<endl;
		opt.maxCplex(revC, costC, state, slope, expSlope);
		//opt.maxCplexUpdate(revC, costC, state, slope, opt.truckNeed, alpha, expSlope, opt.numPigPen);
		if( (p.adpPost) & (n>0) ) updaeSlopePost(state, slope, opt.objective, alpha, expSlope);
		if( (!p.adpPost) ) updaeSlope(state, slope, opt.objective, alpha, expSlope);
		//opt.slopesUpdate(revC, costC, state, slope, alpha);
		transState();
    }
	stopCriteria();
	h+=1;
  }
  storeSlopes();
}
// ===================================================

void ADP::policyADP() {
 
}

// ===================================================
void ADP::optActions() {
	int i, j, iL, iG, iSd, t, n, numPigs;
	int sumSec, sumTotal, numMarkPen;
	int cycleNumbers = 0;
	double sumW, sumG, sumSd, sumReward;
	double sumFeed = 0;
	double sumMarkWeight = 0;
	int numWeeks = 120;
	double threshold = 110;
	if(!p.adpPost){ CalcTransPrSd(); CalcTransPrM(); } // || p.adpPost 
    
	arma::vec totalCull(numWeeks);
	arma::vec trucks(numWeeks);

 //find the best actions for a herd of three sections	
 if( !p.comparePolicy ){

	//initializee the states in the first week with different conditions:
	
	vector< vector< vector<arma::vec> > > numRemainPen;
	vector< vector< vector<arma::rowvec> > > numRemainPenWeight;
	vector< vector<arma::vec> > numRemainSec;

	numRemainPen = vector< vector< vector<arma::vec> > > (numWeeks, 
		vector< vector<arma::vec> >(p.secNum, 
		vector<arma::vec>(p.penNum) ) );

	numRemainPenWeight = vector< vector< vector<arma::rowvec> > > (numWeeks, 
		vector< vector<arma::rowvec> >(p.secNum, 
		vector<arma::rowvec>(p.penNum) ) );

    numRemainSec = vector< vector<arma::vec> > (numWeeks,  
		vector<arma::vec>(p.secNum) ); 
	
	rewCoefficients();	
	readFinalSlopes();

    for(i=0; i<p.secNum;i++){
		t = 0; 
		for(j=0; j<p.penNum; j++){
			state[i][j][0] =  p.sL0; //findIndex(p.iniL[i], p.dL); //findIndex(  meanNext[0]*p.coefL[i] - p.meanWeights(t), p.dL ); // findIndex(p.iniL[i], p.dL);
			state[i][j][1] =  findIndex(p.iniG[i], p.dG); //p.sG0; // // findIndex( meanNext[1]*p.coefG[i], p.dG ); //findIndex(p.iniG[i], p.dG);
			state[i][j][2] =  p.sSd0; //findIndex(p.iniSd[i], p.dSd); //findIndex(  sample.simNGSSM(t,sdnGSSM)*p.coefSd[i] - p.sdWeights(t), p.dSd ); // findIndex(p.iniSd[i], p.dSd);
			state[i][j][3] = t;
			state[i][j][4] = p.pigs;
		}
	}

	for(n=0; n<numWeeks; n++){
		opt.maxCplex(revC, costC, state, slope, expSlope);		
		//opt.maxMyopic(revC, costC, lWeight, state, slope, expSlope, threshold);
		//cout<<" reward "<<opt.reward<<endl;
		sumTotal=0;
		for(i=0; i<p.secNum; i++){
			t = state[i][0][3];
			sumSec = 0; sumW = 0; sumG=0; sumSd=0;
			for(int x=0; x<p.penNum; x++){
				sumSec += state[i][x][4];  sumW += p.dL(state[i][x][0],0) + p.meanWeights(t); sumG += p.dG(state[i][x][1],0); sumSd += p.dSd(state[i][x][2],0) + p.sdWeights(t);  
			}
			numRemainSec[n][i][0] = t;
			numRemainSec[n][i][1] = sumSec; //- (int)sum( opt.numPigPen.row(i) );
			numRemainSec[n][i][2] = sumW/p.penNum;
			numRemainSec[n][i][3] = sumG/p.penNum;
			numRemainSec[n][i][4] = sumSd/p.penNum;
			numRemainSec[n][i][5] = (int)sum( opt.numPigPen.row(i) );
			sumTotal += numRemainSec[n][i][5];

			if( n==0 & t==0) {numRemainSec[n][i][6] = 0; cycleNumbers =0; }
			if( t!=0) {numRemainSec[n][i][6] = numRemainSec[n-1][i][6] ; cycleNumbers = numRemainSec[n-1][i][6]; }
			if( n!=0 & t==0) {numRemainSec[n][i][6] = numRemainSec[n-1][i][6] + 1;  cycleNumbers = numRemainSec[n-1][i][6] +1; } 		

			for(j=0; j<p.penNum; j++){
				iL = state[i][j][0]; iG = state[i][j][1]; iSd = state[i][j][2]; numPigs = state[i][j][4];
				numMarkPen = (int)opt.numPigPen(i,j);
				numRemainPen[n][i][j][0] = t;
				numRemainPen[n][i][j][1] = numPigs; // - opt.numPigPen(i,j);
				numRemainPen[n][i][j][2] =  p.dL(iL,0) + p.meanWeights(t);   
				numRemainPen[n][i][j][3] = p.dG(iG,0); 
				numRemainPen[n][i][j][4] = p.dSd(iSd,0) + p.sdWeights(t);
				numRemainPen[n][i][j][5] = numMarkPen;
				numRemainPenWeight[n][i][j] = lWeight[t][iL][iG][iSd]; //.subvec(1,numPigs)


				//calculate the total feed consumption and average weight at slaughter
				sumFeed += (double)sum( costC[t][iL][iG][iSd].subvec(0, state[i][j][4] - numMarkPen) ) / (double)p.feedPrice;
				if( numMarkPen !=0 ) sumMarkWeight += sum( lWeight[t][iL][iG][iSd].subvec(state[i][j][4] - numMarkPen + 1 ,state[i][j][4]) );
				 
			}
		}
		totalCull[n] = sumTotal;
		trucks[n] = opt.truckNeed;
		transState();
	}

	//Write the results in two csv files:
	ofstream  numPigsPens;
	ofstream  numPigsSecs;
    numPigsPens.open("C:\\Academic\\numPigsPens.csv", ios::trunc);
    numPigsSecs.open("C:\\Academic\\numPigsSecs.csv", ios::trunc);
    numPigsPens << "epoch" << ";" << "t" << ";" << "section" << ";" << "pen" <<  ";" << "numPigs" << ";" << "weight" << ";" << "growth" << ";" << "sd" << ";" << "opt" << ";" << "cycleNum" << ";" << "pigWeight" <<  endl;
    numPigsSecs << "epoch" << ";" << "t" << ";" << "section" << ";" << "numPigs" << ";" <<  "weight" << ";" << "growth" << ";" << "sd" << ";" << "opt" << ";" << "totalCull" << ";" << "cycleNum" << ";" << "truckNum" << endl;
  
	for(n=0; n<numWeeks; n++){
		for(i=0; i<p.secNum; i++){
			numPigsSecs << n << ";" << numRemainSec[n][i][0] << ";" << i << ";" << numRemainSec[n][i][1] << ";" << numRemainSec[n][i][2] << ";" << numRemainSec[n][i][3] << ";" << numRemainSec[n][i][4] << ";" << numRemainSec[n][i][5] << ";" << totalCull[n] << ";" << numRemainSec[n][i][6] << ";" << trucks[n] << endl;  
			for(j=0; j<p.penNum; j++){
				for(int k=1; k<=p.pigs; k++ ){
					numPigsPens << n << ";" << numRemainPen[n][i][j][0] << ";" << i << ";" << j << ";" << numRemainPen[n][i][j][1] << ";" << numRemainPen[n][i][j][2] << ";" << numRemainPen[n][i][j][3] << ";" << numRemainPen[n][i][j][4] << ";" << numRemainPen[n][i][j][5] <<  ";" << numRemainSec[n][i][6] << ";" << numRemainPenWeight[n][i][j][k] <<  endl;  
				}
			}
		}
	}
	numPigsPens.close();
	numPigsSecs.close();
 }
		

 
 // compare policies

 if(p.comparePolicy){ 

	//initializee the states in the first week with different conditions: 
    for(i=0; i<p.secNum;i++){
		t = 8; // ? 
		for(j=0; j<p.penNum; j++){
			state[i][j][0] =  p.sL0; //findIndex(p.iniL[i], p.dL); //findIndex(  meanNext[0]*p.coefL[i] - p.meanWeights(t), p.dL ); // findIndex(p.iniL[i], p.dL);
			state[i][j][1] =  findIndex(p.iniG[i], p.dG); //p.sG0; // // findIndex( meanNext[1]*p.coefG[i], p.dG ); //findIndex(p.iniG[i], p.dG);
			state[i][j][2] =  p.sSd0; //findIndex(p.iniSd[i], p.dSd); //findIndex(  sample.simNGSSM(t,sdnGSSM)*p.coefSd[i] - p.sdWeights(t), p.dSd ); // findIndex(p.iniSd[i], p.dSd);
			state[i][j][3] = t;
			state[i][j][4] = p.pigs;
		}
	}

	
	double sumRewardPen, sumRewardSec, sumRewardEpoch;
	double cycleLgth = 0; 
	sumRewardEpoch = 0;
	for(n=0; n<numWeeks; n++){
		//opt.maxCplex(revC, costC, state, slope, expSlope);		
		opt.maxMyopic(revC, costC, lWeight, state, slope, expSlope, threshold);
		//cout<<" reward "<<opt.reward<<endl;
		sumTotal=0;
		sumRewardSec=0;
		for(i=0; i<p.secNum; i++){						
			t =  state[i][0][3];
			sumTotal += (int)sum( opt.numPigPen.row(i) );
			if( n==0 ) { cycleNumbers = 0; cycleLgth = 0; }
			if( (n!=0) & (opt.termSec[i]==1) ) { cycleNumbers = cycleNumbers + 1; cycleLgth = cycleLgth + (t+1); } 		
			sumRewardPen=0;
			for(j=0; j<p.penNum; j++){
				iL = state[i][j][0]; iG = state[i][j][1]; iSd = state[i][j][2]; numPigs = state[i][j][4]; t =  state[i][j][3];
				numMarkPen = (int)opt.numPigPen(i,j);
				if( numMarkPen!=0 ) sumRewardPen += sum( revC[t][iL][iG][iSd].subvec(numPigs - numMarkPen + 1 ,numPigs) ) - sum( costC[t][iL][iG][iSd].subvec(0, numPigs - numMarkPen) );
				if( numMarkPen==0 ) sumRewardPen += -sum( costC[t][iL][iG][iSd].subvec(0, numPigs - numMarkPen) );
	
				//calculate the total feed consumption and average weight at slaughter
				sumFeed += (double)sum( costC[t][iL][iG][iSd].subvec(0, state[i][j][4] - numMarkPen) ) / (double)p.feedPrice;
				if( numMarkPen !=0 ) sumMarkWeight += sum( lWeight[t][iL][iG][iSd].subvec(state[i][j][4] - numMarkPen + 1 ,state[i][j][4]) );
				 
			}
			sumRewardSec += sumRewardPen - opt.termSec[i]*p.pigletPrice*p.pigs*p.penNum;
		}		
		sumRewardEpoch = sumRewardSec - opt.truckNeed*p.truckCost;
		if(n==0) sumReward = sumRewardEpoch; else sumReward += pow(p.discount,n)*sumRewardEpoch;   

		totalCull[n] = sumTotal;
		trucks[n] = opt.truckNeed;
		transState();
	}
	cycleLgth = (double)cycleLgth / (double)cycleNumbers;
	//cout<<" sumReward "<< sumReward <<endl;
	//cout<<" sumFeed " << sumFeed << endl;
	//cout<<" RewardWeek "<< sumReward/(double)numWeeks <<endl;
	//cout<<" FeedWeek " << sumFeed/(double)numWeeks << endl;
	//cout<<" sumMarkWeight " << sumMarkWeight << endl;
	//cout<<" avgMarkWeight " << (double)sumMarkWeight/(double)sum(totalCull) << endl;

	RewardWeek = sumReward/(double)numWeeks;
	FeedWeek = sumFeed/(double)numWeeks;
	avgMarkWeight = (double)sumMarkWeight/(double)sum(totalCull);
	cycleNumber = cycleNumbers;
	cycleLength = cycleLgth;
	sumTruckNum = sum(trucks);
	sumMarket = sum(totalCull);
 }
}
// ===================================================


void ADP::initializeValue(){
	int iL, iG, iSd, t; 
	for(t=0; t<p.tMax; t++){
		for(iL=0; iL<sizeSL; iL++ ){
			for(iG=0; iG<sizeSG; iG++ ){
				for(iSd=0; iSd<sizeSSd; iSd++ ){
					slope[t][iL][iG][iSd] = 450; //revC[t][iL][iG][iSd][14]; // It is better we start from the slopes found in the feed paper
				}
			}
		}
	}
	//slope[0][iL0][iG0][iSd0] = 413;
}
// ===================================================

void ADP::initializeState(){
	//srand( time(0) );
  int i, j, iL, iG, iSd, t, iPigs;
  for(i=0; i<p.secNum;i++){
	  t = rand() % (p.tMax-2); //rand() % (p.tMax-1); 
     for(j=0; j<p.penNum; j++){
		 iL = rand() % (sizeSL-1);
		 iG = rand() % (sizeSG-1);		 
    	 if( (t>=0) & (t<p.tStartMarketing) ){
        	  iPigs = p.pigs;
			  iSd = rand() % ( (int)(sizeSSd-1)/(int)2 );
		 }else{
		     iPigs = rand() % p.pigs;
		     iSd = rand() % (sizeSSd-1);
		 }
		 if(t==0){
		   	iL = p.sL0;
		    iG = p.sG0;
			iSd = p.sSd0;
		 }
	  
	  state[i][j][0] = iL; // iL
	  state[i][j][1] = iG; // iG
      state[i][j][2] = iSd; //iSd
      state[i][j][3] = t;//t; //t 
	  state[i][j][4] = p.pigs; //iPigs; //r
    }
  }
}
 
// ===================================================
void ADP::rewCoefficients(){
	int iL, iG, iSd, t; 
	for(t=0; t<p.tMax; t++){
		cout<< " t "<<t<<endl;
		for(iL=0; iL<sizeSL; iL++ ){
			for(iG=0; iG<sizeSG; iG++ ){
				for(iSd=0; iSd<sizeSSd; iSd++ ){
					sample.simulateRew( t, p.dL(iL,0), p.dG(iG,0), p.dSd(iSd,0) );
					revC[t][iL][iG][iSd] = sample.revCoef;
					costC[t][iL][iG][iSd] = sample.costCoef;
					lWeight[t][iL][iG][iSd] = sample.sortLW;
				}
			}
		}
	}
}

// ===================================================

void ADP::updaeSlope( vector < vector <arma::vec> > & stateVec, vector < vector < vector < vector<double> > > > slopeVec, double object, double alphaUpdate, arma::mat expSlope){
	int i, j, iL, iG, iSd, t;
	int repeatNum = 1;
	unsigned int sizeVec;
	double deviation;
	double objecRelax;
	MAX optRelax; 
	optRelax.maxCplexRelax(revC, costC, stateVec, slopeVec, expSlope);
	objecRelax = optRelax.objectiveRelax;
	for(i=0; i<p.secNum;i++){
		for(j=0; j<p.penNum; j++){
			iL = stateVec[i][j][0]; iG = stateVec[i][j][1]; iSd = stateVec[i][j][2]; t = stateVec[i][j][3];
			slopeVal[t][iL][iG][iSd][0] = 0;
			if( slopeVal[t][iL][iG][iSd].n_elem >1 ) repeatNum = arma::max( arma::hist(slopeVal[t][iL][iG][iSd], arma::unique(slopeVal[t][iL][iG][iSd]) ) ); //compute number of repeats for the same values of slopes
			
	
			//if( t>=p.tStartMarketing-1 ){ //( t>=p.tStartMarketing-1 )// (opt.termSec[i]==1) & (t==p.tStartMarketing-1)
			//	cout<< " i " << i << " j " << j <<" iLt " << p.dL(iL,0) +  p.meanWeights(t)  << " iGt " << p.dG(iG,0) << " iSdt " <<  p.dSd(iSd,0) + p.sdWeights(t) << " t " << t<< " numPigst " << stateVec[i][j][4] << " slope0 " << slope[0][iL0][iG0][iSd0] << " slopeexp " << expSlope(i,j) << endl;
			//	//cout<< " rev "  << arma::mean(revC[t][iL][iG][iSd]) << endl;
			//	//cout<< " rev "  << arma::mean(costC[t][iL][iG][iSd]) << endl;
			//}
    		
			if( (stateVec[i][j][4]!=p.pigs) & (repeatNum <10) ){   
				stateVec[i][j][4] = stateVec[i][j][4] + 1; 
         		optRelax.maxCplexRelax(revC, costC, stateVec, slopeVec, expSlope);
	        	deviation = optRelax.objectiveRelax - objecRelax;
				slope[t][iL][iG][iSd] =  (1-alphaUpdate) * slope[t][iL][iG][iSd] + alphaUpdate *  deviation;
				sizeVec = slopeVal[t][iL][iG][iSd].n_elem;
				slopeVal[t][iL][iG][iSd].resize(sizeVec + 1);
				slopeVal[t][iL][iG][iSd][sizeVec] = slope[t][iL][iG][iSd];
				stateVec[i][j][4] = stateVec[i][j][4] - 1;
		   }
			if( (stateVec[i][j][4]==p.pigs) & (t>=p.tStartMarketing-1) & (repeatNum <10)  ){  
				stateVec[i][j][4] = stateVec[i][j][4] - 1; 
         		optRelax.maxCplexRelax(revC, costC, stateVec, slopeVec, expSlope);
	        	deviation = objecRelax - optRelax.objectiveRelax;
				slope[t][iL][iG][iSd] = (1-alphaUpdate) * slope[t][iL][iG][iSd] +  alphaUpdate * deviation;
				sizeVec = slopeVal[t][iL][iG][iSd].n_elem;
				slopeVal[t][iL][iG][iSd].resize(sizeVec + 1);
				slopeVal[t][iL][iG][iSd][sizeVec] = slope[t][iL][iG][iSd];
				stateVec[i][j][4] = stateVec[i][j][4] + 1; 
		   }
			if(opt.termSec[i]==1){				
				//deviation = (double) ( stateVec[i][j][4] * slope[t][iL][iG][iSd]  -  ( sum( revC[t][iL][iG][iSd].subvec(0,stateVec[i][j][4]) ) - p.pigletPrice*p.pigs ) ) / (double)(p.pigs*p.discount);
				slopeVal[0][iL0][iG0][iSd0][0] = 0;
				deviation = (double) ( stateVec[i][j][4] * slope[t][iL][iG][iSd]  -  ( stateVec[i][j][4] * arma::mean(revC[t][iL][iG][iSd]) - p.pigletPrice*p.pigs ) ) / (double)(p.pigs*p.discount);
				slope[0][iL0][iG0][iSd0] = (1-alphaUpdate) * slope[0][iL0][iG0][iSd0] +  alphaUpdate * deviation;
				sizeVec = slopeVal[0][iL0][iG0][iSd0].n_elem;
				slopeVal[0][iL0][iG0][iSd0].resize(sizeVec + 1);
				slopeVal[0][iL0][iG0][iSd0][sizeVec] = slope[0][iL0][iG0][iSd0];				
			}
	  	}
	}
	optRelax.envRelax.end();
}


// ===================================================

void ADP::updaeSlopePost( vector < vector <arma::vec> > & stateVec, vector < vector < vector < vector<double> > > > slopeVec, double object, double alphaUpdate, arma::mat expSlope){
	int i, j, iL, iG, iSd, t, iLP, iGP, iSdP, tP;
	int repeatNum = 1;
	unsigned int sizeVec;
	double deviation;
	double objecRelax, sumTermSlope, sumRewardTerm;
	MAX optRelax; 
	optRelax.maxCplexRelax(revC, costC, stateVec, slopeVec, expSlope);
	objecRelax = optRelax.objectiveRelax;
	for(i=0; i<p.secNum;i++){
		sumTermSlope = 0, sumRewardTerm=0; 
		for(j=0; j<p.penNum; j++){			
			iLP = iLPost(i,j); iGP = iGPost(i,j); iSdP = iSdPost(i,j); tP = tPost(i,j);
			iL = stateVec[i][j][0]; iG = stateVec[i][j][1]; iSd = stateVec[i][j][2]; t = stateVec[i][j][3];
			slopeVal[tP][iLP][iGP][iSdP][0] = 0;
			if( slopeVal[tP][iLP][iGP][iSdP].n_elem >1 ) repeatNum = arma::max( arma::hist(slopeVal[tP][iLP][iGP][iSdP], arma::unique(slopeVal[tP][iLP][iGP][iSdP]) ) ); //compute number of repeats for the same values of slopes

			//if( t>=p.tStartMarketing-1 ){ //( (opt.termSec[i]==1) & (t==p.tStartMarketing-1) ){
			//	cout<< " iLt " << p.dL(iL,0) +  p.meanWeights(t)  << " iGt " << p.dG(iG,0) << " iSdt " <<  p.dSd(iSd,0) + p.sdWeights(t) << " t " << t<< " numPigst " << stateVec[i][j][4] << " slope0 " << slope[0][iL0][iG0][iSd0] << " slopeexp " << expSlope(i,j)  << " slopePost "<<  slope[tP][iLP][iGP][iSdP] << endl;
			//	//cout<< " rev "  << revC[t][iL][iG][iSd] << endl;
			//	//cout<< " rev "  << costC[t][iL][iG][iSd] << endl;
			//}

			if( (stateVec[i][j][4]!=p.pigs) &  (repeatNum <10) ){   
				stateVec[i][j][4] = stateVec[i][j][4] + 1; 
         		optRelax.maxCplexRelax(revC, costC, stateVec, slopeVec, expSlope);
	        	deviation = optRelax.objectiveRelax - objecRelax;
				slope[tP][iLP][iGP][iSdP] =  (1-alphaUpdate) * slope[tP][iLP][iGP][iSdP] +  alphaUpdate * deviation;
				sizeVec = slopeVal[tP][iLP][iGP][iSdP].n_elem;
				slopeVal[tP][iLP][iGP][iSdP].resize(sizeVec + 1);
				slopeVal[tP][iLP][iGP][iSdP][sizeVec] = slope[tP][iLP][iGP][iSdP];
				stateVec[i][j][4] = stateVec[i][j][4] - 1;
		   }
			if( (stateVec[i][j][4]==p.pigs)  & (t>=p.tStartMarketing-1) & (repeatNum <10) ){  
				stateVec[i][j][4] = stateVec[i][j][4] - 1; 
         		optRelax.maxCplexRelax(revC, costC, stateVec, slopeVec, expSlope);
	        	deviation = objecRelax - optRelax.objectiveRelax;
				slope[tP][iLP][iGP][iSdP] =  (1-alphaUpdate) * slope[tP][iLP][iGP][iSdP] +  alphaUpdate * deviation;
				sizeVec = slopeVal[tP][iLP][iGP][iSdP].n_elem;
				slopeVal[tP][iLP][iGP][iSdP].resize(sizeVec + 1);
				slopeVal[tP][iLP][iGP][iSdP][sizeVec] = slope[tP][iLP][iGP][iSdP];
				stateVec[i][j][4] = stateVec[i][j][4] + 1; 
		   }
			if( (opt.termSec[i]==1)  ) { sumTermSlope += stateVec[i][j][4] * slope[t][iL][iG][iSd];  sumRewardTerm += sum( revC[t][iL][iG][iSd].subvec(0,stateVec[i][j][4]) ) - p.pigletPrice*p.pigs; }   
	  	}
		if( (opt.termSec[i]==1)  ){
			deviation = (double) ( sumTermSlope  -  sumRewardTerm ) / (double)( p.pigs * p.penNum * p.discount);
			slopeVal[0][iL0][iG0][iSd0][0] = 0;
			slope[0][iL0][iG0][iSd0] = (1-alphaUpdate) * slope[0][iL0][iG0][iSd0] +  alphaUpdate * deviation;
			sizeVec = slopeVal[0][iL0][iG0][iSd0].n_elem;
			slopeVal[0][iL0][iG0][iSd0].resize(sizeVec + 1);
			slopeVal[0][iL0][iG0][iSd0][sizeVec] = slope[0][iL0][iG0][iSd0];				
			}
	}
	optRelax.envRelax.end();
}

// ===================================================

void ADP::transState(){
	int i, j, iL, iG, iSd, t, iLN, iGN, iSdN, numPigs, tN, numPigsOpt;
	arma:: vec meanGSSM(2);
	arma:: vec meanNext(2);
	double sdnGSSM;
	//srand( time(0) );
	iLPost.zeros(); iGPost.zeros(); iSdPost.zeros(); tPost.zeros(); expSlope.zeros();
	for(i=0; i<p.secNum;i++){
		for(j=0; j<p.penNum; j++){			
			iL = state[i][j][0]; iG = state[i][j][1]; iSd = state[i][j][2]; t = state[i][j][3]; numPigs = state[i][j][4];	
			meanGSSM[0] = p.dL(iL,0) + p.meanWeights(t);
			meanGSSM[1] = p.dG(iG,0);
			sdnGSSM = p.dSd(iSd,0) + p.sdWeights(t);
			if(opt.termSec[i]==1){
				tN=0;
				//meanNext = sample.simGSSM(tN, meanGSSM );
				state[i][j][0] = p.sL0; //findIndex(p.iniL[i], p.dL); //findIndex(  meanNext[0] - p.meanWeights(tN), p.dL );  // findIndex(  meanNext[0]*p.coefL[i] - p.meanWeights(tN), p.dL ); //
				state[i][j][1] = p.sG0; //findIndex(p.iniG[i], p.dG); // p.sG0; // //findIndex( meanNext[1], p.dG ); //   // findIndex( meanNext[1]*p.coefG[i], p.dG );  //
				state[i][j][2] = p.sSd0; //findIndex(p.iniSd[i], p.dSd); //findIndex(  sample.simNGSSM(tN,sdnGSSM) - p.sdWeights(tN), p.dSd );   //  findIndex(  sample.simNGSSM(tN,sdnGSSM)*p.coefSd[i] - p.sdWeights(tN), p.dSd ); // 
				state[i][j][3] = 0;
				state[i][j][4] = p.pigs;
			}else{
				tN=t+1;
				numPigsOpt = opt.numPigPen(i,j);
				meanNext = sample.simGSSM(tN, meanGSSM);
				state[i][j][0] = findIndex( meanNext[0] - p.meanWeights(tN), p.dL );
				state[i][j][1] = findIndex( meanNext[1], p.dG );
				state[i][j][2] = findIndex( sample.simNGSSM(tN,sdnGSSM) - p.sdWeights(tN), p.dSd );
				state[i][j][3] = t+1;
				state[i][j][4] = numPigs - numPigsOpt;
				//if(t==0) cout<< " i " << i << " j " << j << " w " <<   meanNext[0] << " g " << meanNext[1] << endl << " sd " << sample.simNGSSM(t,sdnGSSM) << endl; 
				//if we use ADP with post decisions:
				if(p.adpPost)
				iLPost(i,j) = iL; iGPost(i,j) = iG; iSdPost(i,j) = iSd; tPost(i,j) = t+1; 

				//if we use ordinary ADP based on pre-decision state:
				if(  (!p.adpPost ) & (state[i][j][3]>=p.tStartMarketing-1) & ( state[i][j][3]<p.tMax-1 )  ){ // 
					iL = state[i][j][0]; iG = state[i][j][1]; iSd = state[i][j][2]; t = state[i][j][3];
					for (iLN=0; iLN<sizeSL; iLN++){
						for (iGN=0; iGN<sizeSG; iGN++){
							for (iSdN=0; iSdN<sizeSSd; iSdN++){
								expSlope(i,j) += prSd[t][iSd][iSdN] * prM[t][iL][iG][iLN][iGN]  * slope[t+1][iLN][iGN][iSdN];
							}
						}
					}

					expSlope(i,j) = (double) expSlope(i,j)/ (double)pow(10.0,20);
				}

			}
		}
	}
}

// ===================================================

void ADP::stopCriteria(){
	int iL, iG, iSd, t;
	arma:: vec slopeDiff;
	unsigned int sizeVec, sizeVec1;
	for(t=0; t<p.tMax; t++){
		for(iL=0; iL<sizeSL; iL++ ){
			for(iG=0; iG<sizeSG; iG++ ){
				for(iSd=0; iSd<sizeSSd; iSd++ ){
					sizeVec = slopeVal[t][iL][iG][iSd].n_elem;
					if(sizeVec>1){
						sizeVec1 = slopeDiff.n_elem;
				        slopeDiff.resize(sizeVec1 + 1);
						slopeDiff[sizeVec1] = (slopeVal[t][iL][iG][iSd][sizeVec-1] - slopeVal[t][iL][iG][iSd][sizeVec-2]);
					}  
				}
			}
		}
	}
	maxDiffSlope =  arma::max(slopeDiff);
}


// ===================================================

int ADP::findIndex(double st, arma::mat dis){
  for(int i=0; i<dis.n_rows; i++){
    if( ( st>=dis(i,1) ) & ( st<dis(i,2) )  )
      return(i);
    }
	  cout<<"error in index: "<< st << " dis " << dis << endl;
      return(-1);
  }

// ===================================================
void ADP::CalcTransPrSd() {  // calc values for prSd[t][iSdt][iSd]
   //cpuTime.Reset(0); cpuTime.StartTime(0);
   double upperSd, lowerSd, meanSd;
	for(int t = (p.tStartMarketing-1); t<(p.tMax-1); t++){
      for (int iSdt=0; iSdt<sizeSSd; iSdt++) {
		  meanSd = p.dSd(iSdt,0) + p.sdWeights(t);
         for (int iSd=0; iSd<sizeSSd; iSd++) {
			 if (iSd==0) lowerSd = 1; else lowerSd = p.dSd(iSd,1) + p.sdWeights(t+1);
			 upperSd = p.dSd(iSd,2) + p.sdWeights(t+1);
			 prSd[t][iSdt][iSd] = exp( sample.logPrNGSSM(t, pow(lowerSd,2), pow(upperSd,2), pow(meanSd,2) ) );   				
         }
      }   
   } 
   //Rcout << "Time for calculating prSd: " << cpuTime.StopAndGetTotalTimeDiff(0) << endl;
} 

// ===================================================

void ADP::CalcTransPrM() {   // calc values prM[t][iLt][iGt][iL][iG]
   //cpuTime.Reset(0); cpuTime.StartTime(0);
   int t, iLt, iGt, iL, iG;
   arma::vec lower(2), upper(2), mt(2);
    for(t = (p.tStartMarketing-1); t<(p.tMax-1); t++) {
             for (iLt=0; iLt<sizeSL; iLt++) {  
               for (iGt=0; iGt<sizeSG; iGt++) {                   
                      mt[0] = p.dL(iLt,0) + p.meanWeights(t); mt[1] = p.dG(iGt,0);	  
                      for (iL=0; iL<sizeSL; iL++) {   
                        for (iG=0; iG<sizeSG; iG++) {  
                         lower[0] = p.dL(iL,1) + p.meanWeights(t+1); lower[1] = p.dG(iG,1);
						 upper[0] = p.dL(iL,2) + p.meanWeights(t+1); upper[1] = p.dG(iG,2);
						 prM[t][iLt][iGt][iL][iG] = exp( sample.logPrGSSM(t, lower, upper, mt) );    // use log transf to avoid underflow       
						}
					  }
		    	   }
    	        }
           	}	   
    //Rcout << "Time for calculating prM: " << cpuTime.StopAndGetTotalTimeDiff(0) << endl;
}

// ===================================================

void ADP::readIniSlopes(){
  int iL, iG, iSd, t;
  double slopeGiven;
  ifstream  myFile;
  myFile.open("C:\\Academic\\slopeHMDP.csv");
  vector<vector<double> > values;
  vector<double> valueline;
  string line;
  while (getline(myFile, line))
    {
        stringstream linestream(line);
        string item;
        while (getline(linestream, item, ';'))
        {
			valueline.push_back(atof(item.c_str()));
        }
		values.push_back(valueline);
		valueline.clear();
    }    
  myFile.close();

  for(int i=1; i<(values.size()-1); i++){
	  t = values[i][0];
	  iL = values[i][1];
	  iG = values[i][2];
	  iSd = values[i][3];
	  slopeGiven = values[i][4];
	  slope[t][iL][iG][iSd] = slopeGiven;
  }
}

// ===================================================

void ADP::readFinalSlopes(){
  int iL, iG, iSd, t;
  double slopeGiven;
  ifstream  myFile;
  myFile.open("C:\\Academic\\dataADP.csv");
  vector<vector<double> > values;
  vector<double> valueline;
  string line;
  while (getline(myFile, line))
    {
        stringstream linestream(line);
        string item;
        while (getline(linestream, item, ';'))
        {
			valueline.push_back(atof(item.c_str()));
        }
		values.push_back(valueline);
		valueline.clear();
    }    
  myFile.close();

  for(int i=1; i<values.size(); i++){
	  t = values[i][0];
	  iL = values[i][1];
	  iG = values[i][2];
	  iSd = values[i][3];
	  slopeGiven = values[i][4];
	  slope[t][iL][iG][iSd] = slopeGiven;
  }
}

// ===================================================

void ADP::storeSlopes(){

//Store the resalts in the csv files:
  ofstream  myFile;
  ofstream  myFileVal;
  myFile.open("C:\\Academic\\dataADP.csv", ios::trunc);
  myFileVal.open("C:\\Academic\\dataValADP.csv", ios::trunc);
  myFile << "t" << ";" << "iL" << ";" << "iG" << ";" << "iSd" << ";" << "slope" << endl;
  myFileVal << "t" << ";" << "iL" << ";" << "iG" << ";" << "iSd" <<  endl;
  
  int iL, iG, iSd, t; 
	for(t=0; t<p.tMax; t++){
		for(iL=0; iL<sizeSL; iL++ ){
			for(iG=0; iG<sizeSG; iG++ ){
				for(iSd=0; iSd<sizeSSd; iSd++ ){
					if(slopeVal[t][iL][iG][iSd].n_elem > 5) 
						myFile << t << ";" << iL << ";" << iG << ";" << iSd << ";" << slope[t][iL][iG][iSd] << endl;
					if(slopeVal[t][iL][iG][iSd].n_elem > 5){
						myFileVal << t << ";" << iL << ";" << iG << ";" << iSd << endl;
						myFileVal << slopeVal[t][iL][iG][iSd]<< endl;
					} 
				}
			}
		}
	}

	myFile.close();
	myFileVal.close();

}

// ===================================================
