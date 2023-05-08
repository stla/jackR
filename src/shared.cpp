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
