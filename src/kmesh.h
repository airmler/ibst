#ifndef IBST_KMESH_H
#define IBST_KMESH_H

#include <numeric>
#include <vector>
#include <cassert>
#include <functional>

namespace ibst{
using std::accumulate;
using std::multiplies;
class Kmesh {
 public:
  size_t Nk;
  std::vector<size_t> mesh;
  std::vector< std::vector<size_t> > kPoints;
  Kmesh(std::vector<size_t> _mesh) {
    mesh = _mesh;
    Nk = accumulate(mesh.begin(), mesh.end(), 1UL, multiplies<size_t>());
    for (size_t i(0); i < Nk; i++) kPoints.push_back(idxToK(i));
  }

  void print(){
    for (size_t i(0); i < Nk; i++)
      printf("%ld : %ld %ld %ld\n", i, kPoints[i][0], kPoints[i][1], kPoints[i][2]);
  }
  size_t backfold(const size_t i, const size_t N) { return (i+4*N)%N; }

  std::vector<size_t> backfold(const std::vector<size_t> in){
    return { backfold(in[0], mesh[0])
           , backfold(in[1], mesh[1])
           , backfold(in[2], mesh[2])
           };
  }

  std::vector<size_t> idxToK(const size_t i){
    return {  i%mesh[0]
           , (i/mesh[0])%mesh[1]
           ,  i/mesh[0]/mesh[1]
           };
  }

  size_t kToIdx(const std::vector<size_t> kPoint){
    return backfold(kPoint[0], mesh[0])
         + backfold(kPoint[1], mesh[1])*mesh[0]
         + backfold(kPoint[2], mesh[2])*mesh[0]*mesh[1];
  }


  size_t getMinusIdx(const size_t i){
    std::vector<size_t> k(idxToK(i));
    std::vector<size_t> d({-k[0], -k[1], -k[2]});
    return kToIdx(backfold(d));
  }

  size_t idxMinusIdx(const size_t i, const size_t j){
    std::vector<size_t> ki(idxToK(i));
    std::vector<size_t> kj(idxToK(j));
    std::vector<size_t> d( { kj[0] - ki[0]
                      , kj[1] - ki[1]
                      , kj[2] - ki[2]
                      } );
    return kToIdx(backfold(d));
  }

  size_t getForthIdx(const size_t k, const size_t i, const size_t j){
    // kk + ko = ki + kj ==>  ko = ki + kj - kk
    std::vector<size_t> ki(idxToK(i));
    std::vector<size_t> kj(idxToK(j));
    std::vector<size_t> kk(idxToK(k));
    std::vector<size_t> d( { ki[0] + kj[0] - kk[0]
                        , ki[1] + kj[1] - kk[1]
                        , ki[2] + kj[2] - kk[2]
                        } );
    return kToIdx(backfold(d));
  }

  std::vector< std::vector<size_t> > getDefaultNzc(const size_t order){
    if ( order == 0) return {{}};
    if ( order == 1) {
      std::vector< std::vector<size_t> > out(Nk, std::vector<size_t>(2));
      for (size_t o(0); o < Nk; o++) out[o] = {o, o};
      return out;
    }
    if ( order == 2) {
      std::vector< std::vector<size_t> > out(Nk, std::vector<size_t>(3));
      for (size_t o(0); o < Nk; o++) out[o] = {o, o, o};
      return out;
    }
    if ( order == 3) {
      std::vector< std::vector<size_t> > out(Nk*Nk, std::vector<size_t>(4));
      size_t it(0);
      for (size_t m(0); m < Nk; m++)
      for (size_t n(0); n < Nk; n++)
        out[n + m * Nk] = {idxMinusIdx(n,m), n, m, it++};
      return out;
    }
    if ( order == 4) {
      std::vector< std::vector<size_t> > out(Nk*Nk*Nk, std::vector<size_t>(5));
      size_t it(0);
      for (size_t m(0); m < Nk; m++)
      for (size_t n(0); n < Nk; n++)
      for (size_t o(0); o < Nk; o++)
        out[o + n * Nk + m * Nk * Nk] = {m, n, o, getForthIdx(o,n,m), it++};
      return out;
    }
    assert(0);
    return {};
  }

};


} // namespace ibst

#endif /* IBST_KMESH_H */
