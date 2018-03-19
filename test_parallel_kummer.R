library(rbenchmark)
library(Rcpp)

Rcpp::sourceCpp(
  code = '#include <Rcpp.h>
          #include <cmath>
          #include <iostream>
          using namespace Rcpp;

        // [[Rcpp::export]]
        NumericVector Kummer(NumericVector z, NumericVector a, NumericVector b, double relTol) {
        double sumando;
        int n=z.size(),i,j;
        NumericVector suma(n);

        for (i=0;i<n;++i) {
          if (z[i]<80) {
            for (j=0,sumando=1,suma[i]=1;fabs(sumando/suma[i])>relTol;++j){
              //Rcout << j;
              sumando*=(z[i] / (j+1)) * ((a[i] + j) / (b[i] + j));
              suma[i]+=sumando;
            }
          }
          else {
            for (j=0,sumando=1,suma[i]=1;fabs(sumando/suma[i])>relTol;++j){
              //Rcout << j;
              sumando*=((1-a[i]+j) / (j+1)) * ((b[i] - a[i] + j) / (z[i]));
              suma[i]+=sumando;
            }
          suma[i]*=exp(z[i]+(a[i]-b[i])*log(z[i])+lgamma(b[i])-lgamma(a[i]));
          }
        }
        return(suma);
      }
  ')

Rcpp::sourceCpp(
  code = '
  #include "Rcpp.h"
  #include "RcppParallel.h"
  #include <cmath>
  using namespace RcppParallel;

// [[Rcpp::depends(RcppParallel)]]
  struct wKummer : public Worker {

  //input
    const RVector<double> z;
    const RVector<double> a;
    const RVector<double> b;
    const double relTol;

  // output
    RVector<double> suma;

  // initialize from Rcpp input and output
    wKummer(Rcpp::NumericVector z, Rcpp::NumericVector a, Rcpp::NumericVector b, double relTol, Rcpp::NumericVector suma)
       : z(z), a(a), b(b), relTol(relTol), suma(suma) {}

  // function call operator that work for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; i++) {
          suma[i] = 1;
          double sumando = 1;
          //double dif = a[i] - b[i];
          //double dif2 = 1- a[i];

        if (z[i]<80) {
          for (long j=0;fabs(sumando/suma[i])>relTol;++j){
            sumando*=(z[i] / (j+1)) * ((a[i] + j) / (b[i] + j));
            suma[i]+=sumando;
           }
         }
         else {
           for (long j=0;fabs(sumando/suma[i])>relTol;++j){
             sumando*=((1-a[i]+j) / (j+1)) * ((b[i] - a[i] + j) / (z[i]));
             suma[i]+=sumando;
           }
           suma[i]*=exp(z[i]+(a[i]-b[i])*log(z[i])+lgamma(b[i])-lgamma(a[i]));
        }
      }
    }
};

// [[Rcpp::export]]
  Rcpp::NumericVector pKummer(Rcpp::NumericVector z, Rcpp::NumericVector a, Rcpp::NumericVector b, double relTol) {
     int n=z.size();
     Rcpp::NumericVector suma(n);

     wKummer wk(z, a, b, relTol, suma);

     // call it with parallelFor
     parallelFor(0, z.length(), wk, 2000);

     return(suma);
  }

  ')


skummer <- function(x, a, b, relTol = 1e-6){

  output <- Kummer(x, a, b, relTol)
  output[is.infinite(output)] <- .Machine$double.xmax
  return(output)
}

pkummer <- function(x, a, b, relTol = 1e-6){

  output <- pKummer(x, a, b, relTol)
  output[is.infinite(output)] <- .Machine$double.xmax
  return(output)
}


z<-rnorm(100000, 18, 5)
a<-rnorm(100000, 50, 10)
b<-rnorm(100000, 10, 5)

#stopifnot(identical(skummer(z, a, b, 1e-6), pkummer(z, a, b, 1e-6)))
stopifnot(identical(pestim::kummer(z, a, b, 1e-6), pestim::pkummer(z, a, b, 1e-6)))

#res <- benchmark(skummer(z, a, b, 1e-6), pkummer(z, a, b, 1e-6), order="relative")
res <- benchmark(pestim::kummer(z, a, b, 1e-6, nThreads = 1), pestim::kummer(z, a, b, 1e-6, nThreads = 8), order="relative")
res[, 1:4]


