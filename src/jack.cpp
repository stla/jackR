#include "jack.h"
#include "symbolicQspray.h"

using namespace SYMBOLICQSPRAY;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
typedef Qspray<int> Zspray;
typedef std::unordered_map<std::pair<int, int>, Zspray, pairHasher> Zij;

Zspray sch(Partition lambda, Zij S, int m, int k, Partition nu) {
  const int nusize = nu.size();
  if(nusize == 0 || nu[0] == 0 || m == 0) {
    return Zspray(1);
  }
  if(nusize > m && nu[m] > 0) {
    return Zspray(0);
  }
  if(m == 1){
    return Qlone<int>(1).power(nu[0]);
  }
  int N = _N(lambda, nu);
  std::pair<int, int> Nm = std::make_pair(N, m);
  if(auto search = S.find(Nm); search != S.end()) {
    return S[Nm];
  }
  Zspray s = sch(lambda, S, m-1, 1, nu);
  int i = k;
  while(nusize >= i && nu[i-1] > 0) {
    if(nusize == i || nu[i-1] > nu[i]) {
      Partition _nu(nu);
      _nu[i-1] = nu[i-1] - 1;
      if(nu[i-1] > 1) {
        s += Qlone<int>(m) * sch(lambda, S, m, i, _nu);
      } else {
        s += Qlone<int>(m) * sch(lambda, S, m-1, 1, _nu);
      }
    }
    i++;
  }
  if(k == 1) {
    S[Nm] = s;
  }
  return s;
}

Zspray SchurPol(int n, Partition lambda) {
  Zij S;
  return sch(lambda, S, n, 1, lambda);
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// [[Rcpp::export]]
Rcpp::List SchurPolRcpp(int n, Rcpp::IntegerVector lambda) {
  Partition lambdaP(lambda.begin(), lambda.end());
  Zspray S = SchurPol(n, lambdaP);
  Polynomial<int> P = S.get();
  Polynomial<gmpq> Q;
  for(auto it = P.begin(); it != P.end(); it++) {
    Q[it->first] = gmpq(it->second);
  }
  return returnQspray(Qspray<gmpq>(Q));
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename T>
using Qij = std::unordered_map<std::pair<int, int>, Qspray<T>, pairHasher>;

template <typename T>
Qspray<T> jac(
  Partition lambda, Qij<T> S, T alpha,
  int m, int k, Partition mu, Partition nu, T beta
) {
  const int nusize = nu.size();
  if(nusize == 0 || nu[0] == 0 || m == 0) {
    return Qspray<T>(T(1));
  }
  if(nusize > m && nu[m] > 0) {
    return Qspray<T>(T(0));
  }
  T oneT(1);
  if(m == 1) {
    T al(0);
    T prod(1);
    for(int i = 1; i < nu[0]; i++) {
      al += alpha;
      prod *= (al + oneT);
    }
    return Qspray<T>(T(prod)) * (Qlone<T>(1).power(nu[0]));
  }
  int N = _N(lambda, nu);
  std::pair<int, int> Nm = std::make_pair(N, m);
  if(k == 0) {
    if(auto search = S.find(Nm); search != S.end()) {
      return S[Nm];
    }
  }
  Qspray<T> s = 
    jac(lambda, S, alpha, m-1, 0, nu, nu, oneT) * Qspray<T>(beta) *
      (Qlone<T>(m).power(weight(mu) - weight(nu)));
  int i = k > 1 ? k : 1;
  while(nusize >= i && nu[i-1] > 0) {
    if(nusize == i || nu[i-1] > nu[i]) {
      Partition _nu(nu);
      _nu[i-1] = nu[i-1] - 1;
      T gamma = beta * _betaratio<T>(mu, nu, i, alpha);
      if(nu[i-1] > 1) {
        s += jac(lambda, S, alpha, m, i, mu, _nu, gamma);
      } else {
        s += jac(lambda, S, alpha, m-1, 0, _nu, _nu, oneT) * 
              Qspray<T>(gamma) * (Qlone<T>(m).power(weight(mu) - weight(_nu)));
      }
    }
    i++;
  }
  if(k == 0) {
    S[Nm] = s;
  }
  return s;
}

template <typename T>
Qspray<T> JackPol(int n, Partition lambda, T alpha) {
  Qij<T> S;
  return jac(lambda, S, alpha, n, 0, lambda, lambda, T(1));
}

SymbolicQspray JackSymPol(int n, Partition lambda) {
  Qij<RatioOfQsprays<gmpq>> S;
  RatioOfQsprays<gmpq> alpha(Qlone<gmpq>(1), Qspray<gmpq>(gmpq(1)));
  return jac(lambda, S, alpha, n, 0, lambda, lambda, RatioOfQsprays<gmpq>(1));
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// [[Rcpp::export]]
Rcpp::List JackPolRcpp(int n, Rcpp::IntegerVector lambda, std::string alpha) {
  Partition lambdaP(lambda.begin(), lambda.end());
  gmpq alphaQ(alpha);
  Qspray<gmpq> P = JackPol<gmpq>(n, lambdaP, alphaQ);
  return returnQspray(P);
}

// [[Rcpp::export]]
Rcpp::List JackSymPolRcpp(int n, Rcpp::IntegerVector lambda) {
  Partition lambdaP(lambda.begin(), lambda.end());
  SymbolicQspray P = JackSymPol(n, lambdaP);
  return returnSymbolicQspray(P);
}

