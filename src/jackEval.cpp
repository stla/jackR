#include "jack.h"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
typedef std::unordered_map<std::pair<int, int>, gmpq, pairHasher> Gij;

gmpq jacEval(
    std::vector<gmpq> x, Partition lambda, Gij S, gmpq alpha,
    int m, int k, Partition mu, Partition nu, gmpq beta
) {
  const int nusize = nu.size();
  if(nusize == 0 || nu[0] == 0 || m == 0) {
    return gmpq(1);
  }
  if(nusize > m && nu[m] > 0) {
    return gmpq(0);
  }
  gmpq oneq(1, 1);
  if(m == 1) {
    gmpq al(0, 1);
    gmpq prod(1, 1);
    for(int i = 1; i < nu[0]; i++) {
      al += alpha;
      prod *= (al + oneq);
    }
    return gmpqpow(x[0], nu[0]) * prod;
  }
  int N = _N(lambda, nu);
  std::pair<int, int> Nm = std::make_pair(N, m);
  if(k == 0) {
    if(auto search = S.find(Nm); search != S.end()) {
      return S[Nm];
    }
  }
  gmpq s = jacEval(x, lambda, S, alpha, m-1, 0, nu, nu, oneq) * beta *
    gmpqpow(x[m-1], weight(mu) - weight(nu));
  int i = k > 1 ? k : 1;
  while(nusize >= i && nu[i-1] > 0) {
    if(nusize == i || nu[i-1] > nu[i]) {
      Partition _nu(nu);
      _nu[i-1] = nu[i-1] - 1;
      gmpq gamma = beta * _betaratio(mu, nu, i, alpha);
      if(nu[i-1] > 1) {
        s = s + jacEval(x, lambda, S, alpha, m, i, mu, _nu, gamma);
      } else {
        s = s + jacEval(x, lambda, S, alpha, m-1, 0, _nu, _nu, oneq) *
          gamma * gmpqpow(x[m-1], weight(mu) - weight(_nu));
      }
    }
    i++;
  }
  if(k == 0) {
    S[Nm] = s;
  }
  return s;
}

gmpq JackEval(std::vector<gmpq> x, Partition lambda, gmpq alpha) {
  Gij S;
  gmpq oneq(1, 1);
  return jacEval(x, lambda, S, alpha, x.size(), 0, lambda, lambda, oneq);
}

// [[Rcpp::export]]
std::string JackEvalRcpp(
  Rcpp::StringVector x, Rcpp::IntegerVector lambda, std::string alpha
) {
  int n = x.size();
  std::vector<gmpq> xQ;
  xQ.reserve(n);
  for(int i = 0; i < n; i++) {
    xQ.emplace_back(gmpq(Rcpp::as<std::string>(x[i])));
  }
  Partition lambdaP(lambda.begin(), lambda.end());
  gmpq alphaQ(alpha);
  gmpq result = JackEval(xQ, lambdaP, alphaQ);
  return q2str(result);
}


// [[Rcpp::export]]
void test() {
  std::vector<gmpq> x = {gmpq(2), gmpq(3), gmpq(4), gmpq(5)};
  std::vector<int> lambda = {3, 1};
  gmpq alpha(5, 2);
  Rcpp::Rcout << JackEval(x, lambda, alpha);
}
