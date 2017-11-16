#include <Rcpp.h>
#include <cmath>
using namespace std;
using namespace Rcpp;

//Predefine a function for transition probability for later use
double transProb(bool eq, double t1, double t2, double theta, int n) {
    double res=(1-exp((t1-t2)*theta))/n;
    if (!eq) {
        return res;
    }
    else {
        res = 1-(n-1)*res;
        return res;
    }
}

//' Viterbi algorithm implementated in Rcpp
//'
//' This function implements Viterbi algorithm on a continuous-time	HMM, given a vector of timestamps, parameter theta for transition probability function, and a matrix of probability for observed outcome given each state at each timestamp
//'
//' @param ts size m NumericVector of timestamp associated with each observed outcome
//' @param theta a double type parameter for the calculation of transition probability
//' @param obs a NumericMatrix with # of columns the same as the length of ts, and # of rows the same as # of possible states, each element is the probability for observed outcome given each state at each timestamp
//' @return a vector of length m (the length of ts), representing the Viterbi path of maximum likelihood
//' @examples
//' obs <- matrix(c(0.88,0.10,0.88,0.10,0.02,0.30,0.02,0.30,0.10,0.60),2,5)
//' theta <- log(2)
//' ts <- c(1,2,3,4,5)
//' ctmcViterbi(ts,theta,obs)
//' @references BIOS 615 lecture notes
//' @export
// [[Rcpp::export]]
IntegerVector ctmcViterbi(NumericVector ts, double theta, NumericMatrix obs) {
  int m = (int)ts.size();
  int n = (int)obs.nrow();
  if ( obs.ncol() != m )
    stop("The input matrix does not conform to the other parameters");

  IntegerVector viterbiPath(m);
  NumericMatrix delta(n,m);
  IntegerMatrix phi(n,m);
  //build a 2*m matrix to save the top two most likely state at each timestamp
  //the purpose of this matrix is to reduce time complexity to mn
  IntegerMatrix TopStates(2,m);

  double max1=-INFINITY,max2=-INFINITY;
  double pSameState = transProb(true,0,ts[0],theta,n);
  double pDiffState = transProb(false,0,ts[0],theta,n);
  double initProb = pSameState>pDiffState ? pSameState:pDiffState;
  for (int i=0;i<n;i++) {
    delta(i,0) = log(initProb*obs(i,0));
    //fill in topstates matrix for the convenience of calculation in next timestamp
    if (delta(i,0)>max1) {
        max2=max1;
        max1=delta(i,0);
        TopStates(1,0) = TopStates(0,0);
        TopStates(0,0) = i;
    }
    else if (delta(i,0)>max2) {
        max2=delta(i,0);
        TopStates(1,0)=i;
    }
  }

  int prevMaxS;
  for (int j=1;j<m;j++) {
    max1=-INFINITY;
    max2=-INFINITY;
    pSameState = transProb(true,ts[j-1],ts[j],theta,n);
    pDiffState = transProb(false,ts[j-1],ts[j],theta,n);
    for (int i=0;i<n;i++) {
        //fill in delta and phi matrix
        //prevMaxS is the state with the largest delta value in previous timestamp and different from current state
        prevMaxS = i==TopStates(0,j-1)?TopStates(1,j-1):TopStates(0,j-1);
        if (delta(i,j-1)+log(pSameState) > delta(prevMaxS,j-1)+log(pDiffState)) {
            delta(i,j)= delta(i,j-1)+log(pSameState*obs(i,j));
            phi(i,j)=i;
        }
        else {

            delta(i,j)=delta(prevMaxS,j-1)+log(pDiffState*obs(i,j));
            phi(i,j)=prevMaxS;
        }
        //fill in topstates matrix for the convenience of calculation in next timestamp
        if (delta(i,j)>max1) {
            max2=max1;
            max1=delta(i,j);
            TopStates(1,j) = TopStates(0,j);
            TopStates(0,j) = i;
        }
        else if (delta(i,j)>max2) {
            max2=delta(i,j);
            TopStates(1,j)=i;
        }
    }
  }

  double v=-INFINITY;
  for (int i=0;i<n;i++) {
    if (delta(i,m-1)>v) {
        v=delta(i,m-1);
        viterbiPath[m-1]=i;
    }
  }

  for (int j=m-1;j>0;j--) {
    viterbiPath[j-1]=phi(viterbiPath[j],j);
  }

  return viterbiPath;
}

//' Forward Backward algorithm implementated in Rcpp
//'
//' This function implements Forward Backward algorithm on a continuous-time HMM, given a vector of timestamps, parameter theta for transition probability function, and a matrix of probability for observed outcome given each state at each timestamp
//'
//' @param ts size m NumericVector of timestamp associated with each observed outcome
//' @param theta a double type parameter for the calculation of transition probability
//' @param obs a NumericMatrix with # of columns the same as the length of ts, and # of rows the same as # of possible states, each element is the probability for observed outcome given each state at each timestamp
//' @return a m*n NumeriMatrix, where m is the length of ts and n is the # of possible state. The matrix represents the conditional probability of each state at each timestamp given the observed outcome
//' @examples
//' obs <- matrix(c(0.88,0.10,0.88,0.10,0.02,0.30,0.02,0.30,0.10,0.60),2,5)
//' theta <- log(2)
//' ts <- c(1,2.95,3,4,5)
//' ctmcForwardBackward(ts,theta,obs)
//' @references BIOS 615 lecture notes
//' @export
// [[Rcpp::export]]
NumericMatrix ctmcForwardBackward(NumericVector ts, double theta, NumericMatrix obs) {
  int m = (int)ts.size();
  int n = (int)obs.nrow();
  if ( obs.ncol() != m )
    stop("The input matrix does not conform to the other parameters");

  NumericMatrix condProb(n,m);
  NumericMatrix alpha(n,m);
  NumericMatrix beta(n,m);

  //Start forward algorithm
  //Initiate first column in alpha
  double pSame = transProb(true,0,ts[0],theta,n);
  double pDiff = transProb(false,0,ts[0],theta,n);
  for (int i=0;i<n;i++) {
    alpha(i,0)=(pSame/n+pDiff*(n-1)/n)*obs(i,0);
  }

  double col_sum;
  for (int j=1;j<m;j++) {
    //calculate the sum of column j-1 in alpha to reduce the time complexity within each j loop from O(n^2) to O(n)
    //also for the purpose of calculating conditional probability of each column, which would address precision issue
    //the sum idea was came from myself originally, the conditional probability idea was found online, from lecture notes of other universities
    col_sum=0;
    for (int i=0;i<n;i++) {
        col_sum+=alpha(i,j-1);
    }
    //fill in column j of alpha using forward algorithm
    pSame = transProb(true,ts[j-1],ts[j],theta,n);
    pDiff = transProb(false,ts[j-1],ts[j],theta,n);
    for (int i=0;i<n;i++) {
        alpha(i,j)=(pSame*alpha(i,j-1)+pDiff*(col_sum-alpha(i,j-1)))/col_sum;
        alpha(i,j)*=obs(i,j);
    }
  }
  //End forward algorithm
  //Start backward algorithm
  for (int i=0;i<n;i++) {
    beta(i,m-1)=1;
  }

  for (int j=m-2;j>=0;j--) {
    //calculate the column sum of beta(i,j+1)*obs(i,j+1) to reduce the time complexity within each j loop from O(n^2) to O(n)
    //also for the purpose of calculating conditional probability of each column, which would address precision issue
    //the sum idea was came from myself originally, the conditional probability idea was found online, from lecture notes of other universities
    col_sum=0;
    for (int i=0;i<n;i++) {
        col_sum+=beta(i,j+1)*obs(i,j+1);
    }
    //fill in column j of beta using backward algorithm
    pSame = transProb(true,ts[j],ts[j+1],theta,n);
    pDiff = transProb(false,ts[j],ts[j+1],theta,n);
    for (int i=0;i<n;i++) {
        beta(i,j)=(pSame*beta(i,j+1)*obs(i,j+1)+pDiff*(col_sum-beta(i,j+1)*obs(i,j+1)))/col_sum;
    }
  }
  //End backward algorithm

  //Fill in condProb using alpha and beta
  for (int j=0;j<m;j++) {
    col_sum=0;
    for (int i=0;i<n;i++) {
        col_sum+=alpha(i,j)*beta(i,j);
    }
    for (int i=0;i<n;i++) {
        condProb(i,j)=alpha(i,j)*beta(i,j)/col_sum;
    }
  }
  return condProb;
}
