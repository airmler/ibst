#ifndef IBST_KMESH_H
#define IBST_KMESH_H

#include <cassert>
#include <functional>
#include <numeric>
#include <vector>

namespace ibst {
using std::accumulate;
using std::multiplies;
class Kmesh {
public:
  int64_t Nk;
  std::vector<int64_t> mesh;
  std::vector<std::vector<int64_t>> kPoints;
  Kmesh(std::vector<int64_t> _mesh) {
    mesh = _mesh;
    Nk = accumulate(mesh.begin(), mesh.end(), 1L, multiplies<int64_t>());
    for (int64_t i(0); i < Nk; i++) kPoints.push_back(idxToK(i));
  }

  void print() {
    for (int64_t i(0); i < Nk; i++)
      printf("%ld : %ld %ld %ld\n",
             i,
             kPoints[i][0],
             kPoints[i][1],
             kPoints[i][2]);
  }
  int64_t backfold(const int64_t i, const int64_t N) { return (i + 4 * N) % N; }

  std::vector<int64_t> backfold(const std::vector<int64_t> in) {
    return {backfold(in[0], mesh[0]),
            backfold(in[1], mesh[1]),
            backfold(in[2], mesh[2])};
  }

  std::vector<int64_t> idxToK(const int64_t i) {
    return {i % mesh[0], (i / mesh[0]) % mesh[1], i / mesh[0] / mesh[1]};
  }

  int64_t kToIdx(const std::vector<int64_t> kPoint) {
    return backfold(kPoint[0], mesh[0]) + backfold(kPoint[1], mesh[1]) * mesh[0]
         + backfold(kPoint[2], mesh[2]) * mesh[0] * mesh[1];
  }

  int64_t getMinusIdx(const int64_t i) {
    std::vector<int64_t> k(idxToK(i));
    std::vector<int64_t> d({-k[0], -k[1], -k[2]});
    return kToIdx(backfold(d));
  }

  int64_t idxMinusIdx(const int64_t i, const int64_t j) {
    std::vector<int64_t> ki(idxToK(i));
    std::vector<int64_t> kj(idxToK(j));
    std::vector<int64_t> d({kj[0] - ki[0], kj[1] - ki[1], kj[2] - ki[2]});
    return kToIdx(backfold(d));
  }

  int64_t getForthIdx(const int64_t k, const int64_t i, const int64_t j) {
    // kk + ko = ki + kj ==>  ko = ki + kj - kk
    std::vector<int64_t> ki(idxToK(i));
    std::vector<int64_t> kj(idxToK(j));
    std::vector<int64_t> kk(idxToK(k));
    std::vector<int64_t> d(
        {ki[0] + kj[0] - kk[0], ki[1] + kj[1] - kk[1], ki[2] + kj[2] - kk[2]});
    return kToIdx(backfold(d));
  }

  std::vector<std::vector<int64_t>> getDefaultNzc(const int64_t order) {
    if (order == 0)
      return {{}};
    if (order == 1) {
      std::vector<std::vector<int64_t>> out(Nk, std::vector<int64_t>(2));
      for (int64_t o(0); o < Nk; o++) out[o] = {o, o};
      return out;
    }
    if (order == 2) {
      std::vector<std::vector<int64_t>> out(Nk, std::vector<int64_t>(3));
      for (int64_t o(0); o < Nk; o++) out[o] = {o, o, o};
      return out;
    }
    if (order == 3) {
      std::vector<std::vector<int64_t>> out(Nk * Nk, std::vector<int64_t>(4));
      int64_t it(0);
      for (int64_t m(0); m < Nk; m++)
        for (int64_t n(0); n < Nk; n++)
          out[n + m * Nk] = {idxMinusIdx(n, m), n, m, it++};
      return out;
    }
    if (order == 4) {
      std::vector<std::vector<int64_t>> out(Nk * Nk * Nk,
                                            std::vector<int64_t>(5));
      int64_t it(0);
      for (int64_t m(0); m < Nk; m++)
        for (int64_t n(0); n < Nk; n++)
          for (int64_t o(0); o < Nk; o++)
            out[o + n * Nk + m * Nk * Nk] = {m,
                                             n,
                                             o,
                                             getForthIdx(o, n, m),
                                             it++};
      return out;
    }
    assert(0);
    return {};
  }
};

} // namespace ibst

#endif /* IBST_KMESH_H */
