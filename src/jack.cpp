#include <Rcpp.h>
typedef std::vector<int>        Partition;
typedef std::vector<signed int> Powers;

class Hasher {
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

typedef std::unordered_map<Powers, int, Hasher> Zpoly;


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
std::unordered_map<Powers, CoeffT, Hasher> polyAdd(
    std::unordered_map<Powers, CoeffT, Hasher> P1,
    std::unordered_map<Powers, CoeffT, Hasher> P2
) {
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
std::unordered_map<Powers, CoeffT, Hasher> polyMult(
  const std::unordered_map<Powers, CoeffT, Hasher> P1,
  const std::unordered_map<Powers, CoeffT, Hasher> P2
) {

  std::unordered_map<Powers, CoeffT, Hasher> Pout;
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
std::unordered_map<Powers, CoeffT, Hasher> zeroPoly() {
  std::unordered_map<Powers, CoeffT, Hasher> out;
  return out;
}

template Zpoly zeroPoly<int>();


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename CoeffT>
std::unordered_map<Powers, CoeffT, Hasher> unitPoly() {
  std::unordered_map<Powers, CoeffT, Hasher> out;
  Powers pows(0);
  CoeffT one(1);
  out[pows] = one;
  return out;
}

template Zpoly unitPoly<int>();


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
template <typename CoeffT>
std::unordered_map<Powers, CoeffT, Hasher> lonePoly(const int n) {
  std::unordered_map<Powers, CoeffT, Hasher> out;
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
std::unordered_map<Powers, CoeffT, Hasher> polyPow(
  const std::unordered_map<Powers, CoeffT, Hasher> P,
  int n
) {
  std::unordered_map<Powers, CoeffT, Hasher> out;
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
  while(nusize >= i && nu[i] > 0) {
    if(nusize == i || nu[i] > nu[i+1]) {
      Partition _nu(nu);
      _nu[i] = nu[i] - 1;
      if(nu[i] > 1) {
        s = polyAdd<int>(
          s, polyMult<int>(lonePoly<int>(m), sch(lambda, S, m, i, _nu))
        );
      } else {
        s = polyAdd<int>(
          s, polyMult<int>(lonePoly<int>(m), sch(lambda, S, m-1, 0, _nu))
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

Zpoly SchurPol(int n, Partition lambda) {
  Zij S;
  return sch(lambda, S, n, 0, lambda);
}

// [[Rcpp::export]]
Rcpp::List SchurPolRcpp(int n, std::vector<int> lambda) {
  Zpoly P = SchurPol(n, lambda);
  int nterms = P.size();
  Rcpp::List Exponents(nterms);
  Rcpp::IntegerVector Coeffs(nterms);
  int i = 0;
  for(auto it = P.begin(); it != P.end(); it++) {
    Powers pows = it->first;
    Rcpp::IntegerVector expnts(pows.begin(), pows.end());
    Exponents(i) = expnts;
    Coeffs(i) = it->second;
  }
  return Rcpp::List::create(
    Rcpp::Named("exponents") = Exponents,
    Rcpp::Named("coeffs")    = Coeffs
  );
}

