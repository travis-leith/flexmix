// Core numerical kernels for flexmix, parallelised via RcppParallel.
// RcppParallel uses tbb or tinythread so parallelism is non-optional.

#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

// ---------------------------------------------------------------------------
// log_row_sums worker
// ---------------------------------------------------------------------------
struct LogRowSumsWorker : public Worker {
  const RMatrix<double> m;
  RVector<double>       out;
  
  LogRowSumsWorker(const NumericMatrix m_, NumericVector out_)
    : m(m_), out(out_) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    const int k = m.ncol();
    for (std::size_t i = begin; i < end; ++i) {
      double mx = m(i, 0);
      for (int j = 1; j < k; ++j)
        if (m(i, j) > mx) mx = m(i, j);
        double s = 0.0;
        for (int j = 0; j < k; ++j)
          s += std::exp(m(i, j) - mx);
        out[i] = mx + std::log(s);
    }
  }
};

// [[Rcpp::export]]
NumericVector cpp_log_row_sums(NumericMatrix m) {
  NumericVector out(m.nrow());
  LogRowSumsWorker worker(m, out);
  parallelFor(0, m.nrow(), worker);
  return out;
}

// ---------------------------------------------------------------------------
// E-step worker
// Fuses: add log(prior), exp, normalise, NaN repair, log_row_sums
// ---------------------------------------------------------------------------
struct EStepWorker : public Worker {
  const RMatrix<double> postunscaled_log;
  const RMatrix<double> prior;
  const bool            prior_per_row;
  RMatrix<double>       logpost;
  RMatrix<double>       postuns;
  RMatrix<double>       posts;
  RVector<double>       lrs;
  
  EStepWorker(const NumericMatrix& pul,
              const NumericMatrix& pr,
              bool ppr,
              NumericMatrix& lp,
              NumericMatrix& pu,
              NumericMatrix& ps,
              NumericVector& l)
    : postunscaled_log(pul), prior(pr), prior_per_row(ppr),
      logpost(lp), postuns(pu), posts(ps), lrs(l) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    const int k = postunscaled_log.ncol();
    for (std::size_t i = begin; i < end; ++i) {
      double mx = R_NegInf;
      for (int j = 0; j < k; ++j) {
        const double pj = prior_per_row ? prior(i,j) : prior(0,j);
        const double v  = postunscaled_log(i,j) + std::log(pj);
        logpost(i,j) = v;
        if (v > mx) mx = v;
      }
      double s = 0.0;
      for (int j = 0; j < k; ++j)
        s += std::exp(logpost(i,j) - mx);
      const double lr = mx + std::log(s);
      lrs[i] = lr;
      
      bool any_nan = false;
      for (int j = 0; j < k; ++j) {
        postuns(i,j) = std::exp(logpost(i,j));
        posts(i,j)   = std::exp(logpost(i,j) - lr);
        if (!std::isfinite(posts(i,j))) any_nan = true;
      }
      if (any_nan) {
        for (int j = 0; j < k; ++j) {
          posts(i,j)  = prior_per_row ? prior(i,j) : prior(0,j);
          postuns(i,j) = std::numeric_limits<double>::min();
        }
      }
    }
  }
};

// [[Rcpp::export]]
List cpp_estep(NumericMatrix postunscaled_log,
               NumericMatrix prior,
               IntegerVector groupfirst_idx,   // 0-based
               NumericVector weights) {         // length 0 = unweighted
  const int n = postunscaled_log.nrow();
  const int k = postunscaled_log.ncol();
  const bool prior_per_row = (prior.nrow() == n);
  
  NumericMatrix logpost(n, k);
  NumericMatrix postuns(n, k);
  NumericMatrix posts(n, k);
  NumericVector lrs(n);
  
  EStepWorker worker(postunscaled_log, prior, prior_per_row,
                     logpost, postuns, posts, lrs);
  parallelFor(0, n, worker);
  
  // llh sum is small (groupfirst rows only) — serial is fine
  const bool weighted = (weights.size() == n);
  double llh = 0.0;
  const int gn = groupfirst_idx.size();
  for (int t = 0; t < gn; ++t) {
    const int i = groupfirst_idx[t];
    llh += weighted ? lrs[i] * weights[i] : lrs[i];
  }
  
  return List::create(
    _["logpost"]      = logpost,
    _["postunscaled"] = postuns,
    _["postscaled"]   = posts,
    _["llh"]          = llh
  );
}

// ---------------------------------------------------------------------------
// Hard / CEM assignment worker
// ---------------------------------------------------------------------------
struct HardAssignWorker : public Worker {
  const RMatrix<double> p;
  const RVector<double> weights;
  const bool            weighted;
  RMatrix<double>       z;
  
  HardAssignWorker(const NumericMatrix& p_,
                   const NumericVector& w_,
                   bool wt,
                   NumericMatrix& z_)
    : p(p_), weights(w_), weighted(wt), z(z_) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    const int k = p.ncol();
    for (std::size_t i = begin; i < end; ++i) {
      int best = 0;
      double bv = p(i, 0);
      for (int j = 1; j < k; ++j)
        if (p(i,j) > bv) { bv = p(i,j); best = j; }
        z(i, best) = weighted ? weights[i] : 1.0;
    }
  }
};

// [[Rcpp::export]]
NumericMatrix cpp_hard_assign(NumericMatrix p, NumericVector weights) {
  const bool weighted = ((int)weights.size() == p.nrow());
  NumericMatrix z(p.nrow(), p.ncol());
  HardAssignWorker worker(p, weights, weighted, z);
  parallelFor(0, p.nrow(), worker);
  return z;
}

// ---------------------------------------------------------------------------
// Group posteriors worker
// Two-phase: serial reduce into sums, parallel broadcast back.
// The reduce must be serial (race condition on sums); broadcast is parallel.
// ---------------------------------------------------------------------------
struct GroupBroadcastWorker : public Worker {
  const RMatrix<double> sums;
  const RVector<int>    group;
  RMatrix<double>       out;
  
  GroupBroadcastWorker(const NumericMatrix& s,
                       const IntegerVector& g,
                       NumericMatrix& o)
    : sums(s), group(g), out(o) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    const int k = out.ncol();
    for (std::size_t i = begin; i < end; ++i) {
      const int g = group[i] - 1;
      for (int j = 0; j < k; ++j)
        out(i,j) = sums(g,j);
    }
  }
};

// [[Rcpp::export]]
NumericMatrix cpp_group_posteriors(NumericMatrix x,
                                   IntegerVector group,
                                   int n_groups) {
  const int n = x.nrow();
  const int k = x.ncol();
  
  // serial reduce (no races)
  NumericMatrix sums(n_groups, k);
  for (int i = 0; i < n; ++i) {
    const int g = group[i] - 1;
    for (int j = 0; j < k; ++j)
      sums(g,j) += x(i,j);
  }
  
  // parallel broadcast
  NumericMatrix out(n, k);
  GroupBroadcastWorker worker(sums, group, out);
  parallelFor(0, n, worker);
  return out;
}