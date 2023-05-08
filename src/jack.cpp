#include "jack.h"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
typedef std::unordered_map<std::pair<int, int>, Zpoly, pairHasher> Zij;

Zpoly sch(Partition lambda, Zij S, int m, int k, Partition nu) {
  const int nusize = nu.size();
  if(nusize == 0 || nu[0] == 0 || m == 0) {
    return unitPoly<int>();
  }
  if(nusize > m && nu[m] > 0) {
    return zeroPoly<int>();
  }
  if(m == 1){
    return polyPow<int>(lonePoly<int>(1), nu[0]);
  }
  int N = _N(lambda, nu);
  std::pair<int, int> Nm = std::make_pair(N, m);
  if(auto search = S.find(Nm); search != S.end()) {
    return S[Nm];
  }
  Zpoly s = sch(lambda, S, m-1, 1, nu);
  int i = k;
  while(nusize >= i && nu[i-1] > 0) {
    if(nusize == i || nu[i-1] > nu[i]) {
      Partition _nu(nu);
      _nu[i-1] = nu[i-1] - 1;
      if(nu[i-1] > 1) {
        s = polyAdd<int>(
          s, polyMult<int>(lonePoly<int>(m), sch(lambda, S, m, i, _nu))
        );
      } else {
        s = polyAdd<int>(
          s, polyMult<int>(lonePoly<int>(m), sch(lambda, S, m-1, 1, _nu))
        );
      }
    }
    i++;
  }
  if(k == 1) {
    S[Nm] = s;
  }
  return s;
}

Zpoly SchurPol(int n, Partition lambda) {
  Zij S;
  return sch(lambda, S, n, 1, lambda);
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// [[Rcpp::export]]
Rcpp::List SchurPolRcpp(int n, Rcpp::IntegerVector lambda) {
  Partition lambdaP(lambda.begin(), lambda.end());
  Zpoly P = SchurPol(n, lambdaP);
  int nterms = P.size();
  Rcpp::List Exponents(nterms);
  Rcpp::IntegerVector Coeffs(nterms);
  int i = 0;
  for(auto it = P.begin(); it != P.end(); it++) {
    Powers pows = it->first;
    Rcpp::IntegerVector expnts(pows.begin(), pows.end());
    Exponents(i) = expnts;
    Coeffs(i) = it->second;
    i++;
  }
  return Rcpp::List::create(
    Rcpp::Named("exponents") = Exponents,
    Rcpp::Named("coeffs")    = Coeffs
  );
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
typedef std::unordered_map<std::pair<int, int>, Qpoly, pairHasher> Qij;

Qpoly jac(
  Partition lambda, Qij S, gmpq alpha,
  int m, int k, Partition mu, Partition nu, gmpq beta
) {
  const int nusize = nu.size();
  if(nusize == 0 || nu[0] == 0 || m == 0) {
    return unitPoly<gmpq>();
  }
  if(nusize > m && nu[m] > 0) {
    return zeroPoly<gmpq>();
  }
  gmpq oneq(1, 1);
  if(m == 1) {
    gmpq al(0, 1);
    gmpq prod(1, 1);
    for(int i = 1; i < nu[0]; i++) {
      al += alpha;
      prod *= (al + oneq);
    }
    return polyMult<gmpq>(
      constantPoly<gmpq>(prod), polyPow<gmpq>(lonePoly<gmpq>(1), nu[0])
    );
  }
  int N = _N(lambda, nu);
  std::pair<int, int> Nm = std::make_pair(N, m);
  if(k == 0) {
    if(auto search = S.find(Nm); search != S.end()) {
      return S[Nm];
    }
  }
  Qpoly s = polyMult<gmpq>(
    jac(lambda, S, alpha, m-1, 0, nu, nu, oneq),
    polyMult<gmpq>(
      constantPoly<gmpq>(beta),
      polyPow<gmpq>(lonePoly<gmpq>(m), weight(mu) - weight(nu))
    )
  );
  int i = k > 1 ? k : 1;
  while(nusize >= i && nu[i-1] > 0) {
    if(nusize == i || nu[i-1] > nu[i]) {
      Partition _nu(nu);
      _nu[i-1] = nu[i-1] - 1;
      gmpq gamma = beta * _betaratio(mu, nu, i, alpha);
      if(nu[i-1] > 1) {
        s = polyAdd<gmpq>(
          s, jac(lambda, S, alpha, m, i, mu, _nu, gamma)
        );
      } else {
        s = polyAdd<gmpq>(
          s,
          polyMult<gmpq>(
            jac(lambda, S, alpha, m-1, 0, _nu, _nu, oneq),
            polyMult<gmpq>(
              constantPoly<gmpq>(gamma),
              polyPow<gmpq>(lonePoly<gmpq>(m), weight(mu) - weight(_nu))
            )
          )
        );
      }
    }
    i++;
  }
  if(k == 0) {
    S[Nm] = s;
  }
  return s;
}

Qpoly JackPol(int n, Partition lambda, gmpq alpha) {
  Qij S;
  gmpq oneq(1, 1);
  return jac(lambda, S, alpha, n, 0, lambda, lambda, oneq);
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// [[Rcpp::export]]
Rcpp::List JackPolRcpp(int n, Rcpp::IntegerVector lambda, std::string alpha) {
  Partition lambdaP(lambda.begin(), lambda.end());
  gmpq alphaQ(alpha);
  Qpoly P = JackPol(n, lambdaP, alphaQ);
  int nterms = P.size();
  Rcpp::List Exponents(nterms);
  Rcpp::CharacterVector Coeffs(nterms);
  int i = 0;
  for(auto it = P.begin(); it != P.end(); it++) {
    Powers pows = it->first;
    Rcpp::IntegerVector expnts(pows.begin(), pows.end());
    Exponents(i) = expnts;
    Coeffs(i) = q2str(it->second);
    i++;
  }
  return Rcpp::List::create(
    Rcpp::Named("exponents") = Exponents,
    Rcpp::Named("coeffs")    = Coeffs
  );
}
