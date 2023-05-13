#include "jack.h"

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
template <typename numT>
numT _betaratio(Partition kappa, Partition mu, int k, numT alpha) {
  std::vector<numT> muT;
  std::vector<numT> kappaT;
  std::vector<numT> s;
  muT.reserve(k); kappaT.reserve(k); s.reserve(k);
  for(int i = 0; i < k; i++) {
    muT.emplace_back(numT(mu[i]));
    kappaT.emplace_back(numT(kappa[i]));
    s.emplace_back(numT(i+1));
  }
  numT oneT(1);
  numT t = s[k-1] - alpha * muT[k-1];
  std::vector<numT> u;
  u.reserve(k);
  for(int i = 0; i < k; i++) {
    u.emplace_back(t + oneT - s[i] + alpha * kappaT[i]);
  }
  std::vector<numT> v;
  v.reserve(k-1);
  for(int i = 0; i < k-1; i++) {
    v.emplace_back(t - s[i] + alpha * muT[i]);
  }
  int musize = mu.size();
  int muk = mu[k-1];
  std::vector<numT> w;
  w.reserve(muk-1);
  numT al(0);
  for(int i = 1; i < muk; i++) {
    int j = 0;
    while(j < musize && mu[j] >= i) {
      j++;
    }
    al += alpha;
    w.emplace_back(numT(j) - t - al);
  }
  numT prod1(1);
  numT prod2(1);
  numT prod3(1);
  for(int i = 0; i < k; i++) {
    prod1 *= (u[i] / (u[i] + alpha - oneT));
  }
  for(int i = 0; i < k-1; i++) {
    prod2 *= (oneT + alpha / v[i]);
  }
  for(int i = 0; i < muk-1; i++) {
    prod3 *= (oneT + alpha / w[i]);
  }
  return alpha * prod1 * prod2 * prod3;
}

template gmpq _betaratio<gmpq>(Partition, Partition, int, gmpq);
template double _betaratio<double>(Partition, Partition, int, double);


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
int weight(Partition mu) {
  int w = 0;
  int musize = mu.size();
  for(int i = 0; i < musize; i++) {
    w += mu[i];
  }
  return w;
}
