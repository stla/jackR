#include <Rcpp.h>
#include <boost/multiprecision/gmp.hpp>
typedef std::vector<int>                    Partition;
typedef std::vector<signed int>             Powers;
typedef boost::multiprecision::mpq_rational gmpq;
typedef boost::multiprecision::mpz_int      gmpi;


class vecHasher {
public:
  size_t operator()(const Powers& exponents) const {
    // thanks to Steffan Hooper for advice
    std::size_t seed = 0;
    for(auto& i : exponents) {
      seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

template <typename CoeffT>
using Poly = std::unordered_map<Powers, CoeffT, vecHasher>;

typedef Poly<int>  Zpoly;
typedef Poly<gmpq> Qpoly;


class pairHasher {
public:
  size_t operator()(const std::pair<int, int>& ij) const {
    size_t seed = 0;
    seed ^= ij.first + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= ij.second + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
int _N(Partition lambda, Partition mu) {
  size_t n = lambda.size();
  int out = mu[n-1];
  int product = 1;
  for(size_t i = n-1; i > 0; i--) {
    product *= lambda[i] + 1;
    out += mu[i-1] * product;
  }
  return out;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
void simplifyPowers(Powers& pows) {
  int n = pows.size();
  if(n == 0) {
    return;
  }
  n--;
  Powers::iterator it = pows.end();
  bool zero = pows[n] == 0;
  while(zero && n >= 0) {
    it--;
    n--;
    zero = pows[n] == 0;
  }
  if(n == -1) {
    pows = {};
  } else {
    pows.erase(it, pows.end());
  }
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename CoeffT>
Poly<CoeffT> polyAdd(Poly<CoeffT> P1, Poly<CoeffT> P2) {
  Powers pows;
  for(auto it = P2.begin(); it != P2.end(); ++it) {
    pows = it->first;
    P1[pows] += P2[pows];
    if(P1[pows] == 0) {
      P1.erase(pows);
    }
  }
  return P1;
}

template Zpoly polyAdd<int>(Zpoly, Zpoly);

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
Powers growPowers(Powers pows, int m, int n) {
  Powers gpows;
  gpows.reserve(n);
  for(int i = 0; i < m; i++) {
    gpows.emplace_back(pows[i]);
  }
  for(int i = m; i < n; i++) {
    gpows.emplace_back(0);
  }
  return gpows;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename CoeffT>
Poly<CoeffT> polyMult(const Poly<CoeffT> P1, const Poly<CoeffT> P2) {

  Poly<CoeffT> Pout;
  Powers powssum;
  int i;

  for(auto it1 = P1.begin(); it1 != P1.end(); ++it1) {
    const CoeffT c1 = it1->second;
    if(c1 != 0) {
      Powers pows1 = it1->first;
      int n1 = pows1.size();
      for(auto it2 = P2.begin(); it2 != P2.end(); ++it2) {
        const CoeffT c2 = it2->second;
        if(c2 != 0) {
          Powers pows2 = it2->first;
          int n2 = pows2.size();
          powssum.clear();
          if(n1 < n2) {
            Powers gpows = growPowers(pows1, n1, n2);
            powssum.reserve(n2);
            for(i = 0; i < n2; i++) {
              powssum.emplace_back(gpows[i] + pows2[i]);
            }
          } else if(n1 > n2) {
            Powers gpows = growPowers(pows2, n2, n1);
            powssum.reserve(n1);
            for(i = 0; i < n1; i++) {
              powssum.emplace_back(pows1[i] + gpows[i]);
            }
          } else {
            powssum.reserve(n1);
            for(i = 0; i < n1; i++) {
              powssum.emplace_back(pows1[i] + pows2[i]);
            }
          }
          Pout[powssum] += c1 * c2;
        }
      }
    }
  }

  return Pout;
}

template Zpoly polyMult<int>(const Zpoly, const Zpoly);


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename CoeffT>
Poly<CoeffT> zeroPoly() {
  Poly<CoeffT> out;
  return out;
}

template Zpoly zeroPoly<int>();


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename CoeffT>
Poly<CoeffT> unitPoly() {
  Poly<CoeffT> out;
  Powers pows(0);
  CoeffT one(1);
  out[pows] = one;
  return out;
}

template Zpoly unitPoly<int>();


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename CoeffT>
Poly<CoeffT> lonePoly(const int n) {
  Poly<CoeffT> out;
  Powers pows(n, 0);
  pows[n-1] = 1;
  CoeffT one(1);
  out[pows] = one;
  return out;
}

template Zpoly lonePoly<int>(const int);


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename CoeffT>
Poly<CoeffT> polyPow(const Poly<CoeffT> P, int n) {
  Poly<CoeffT> out;
  if(n >= 1) {
    if(n == 1) {
      out = P;
    } else {
      out = unitPoly<CoeffT>();
      for(; n > 0; n--) {
        out = polyMult(P, out);
      }
    }
  } else {
    out = unitPoly<CoeffT>();
  }
  return out;
}

template Zpoly polyPow<int>(const Zpoly, int);


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
  if(S.contains(Nm)) {
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
gmpq _betaratio(Partition kappa, Partition mu, int k, gmpq alpha) {
  std::vector<gmpq> mu_q;
  std::vector<gmpq> kappa_q;
  std::vector<gmpq> s;
  mu_q.reserve(k); kappa_q.reserve(k); s.reserve(k);
  for(int i = 0; i < k; i++) {
    mu_q.emplace_back(gmpq(mu[i], 1));
    kappa_q.emplace_back(gmpq(kappa[i], 1));
    s.emplace_back(gmpq(i+1, 1));
  }
  gmpq oneq(1, 1);
  gmpq t = s[k-1] - alpha * mu_q[k-1];
  std::vector<gmpq> u;
  u.reserve(k);
  for(int i = 0; i < k; i++) {
    u.emplace_back(t + oneq - s[i] + alpha * kappa_q[i]);
  }
  std::vector<gmpq> v;
  v.reserve(k-1);
  for(int i = 0; i < k-1; i++) {
    v.emplace_back(t - s[i] + alpha * mu_q[i]);
  }
  int musize = mu.size();
  int muk = mu[k-1];
  std::vector<gmpq> w;
  w.reserve(muk-1);
  gmpq al(0, 1);
  for(int i = 1; i < muk; i++) {
    int j = 0;
    while(j < musize && mu[j] >= i) {
      j++;
    }
    al += alpha;
    w.emplace_back(gmpq(j, 1) - t - al);
  }
  gmpq prod1(1, 1);
  gmpq prod2(1, 1);
  gmpq prod3(1, 1);
  for(int i = 0; i < k; i++) {
    prod1 *= (u[i] / (u[i] + alpha - oneq));
  }
  for(int i = 0; i < k-1; i++) {
    prod2 *= (oneq + alpha / v[i]);
  }
  for(int i = 0; i < muk-1; i++) {
    prod3 *= (oneq + alpha / w[i]);
  }
  return alpha * prod1 * prod2 * prod3;
}

// [[Rcpp::export]]
void test() {
  Partition kappa = {4, 2, 2};
  Partition mu    = {3, 2, 1};
  int k = 3;
  gmpq alpha(2, 3);
  Rcpp::Rcout << _betaratio(kappa, mu, k, alpha);
}
