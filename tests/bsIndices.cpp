#include "../src/bsIndices.h"
#include <iostream>
#include "../src/kmesh.h"

using ibst::BsInd;
using ibst::getIndices;

int main(){


  auto kmesh(ibst::Kmesh({3,2,1}));
  auto f(kmesh.getDefaultNzc(4));
  auto t(kmesh.getDefaultNzc(3));
  auto o(kmesh.getDefaultNzc(1));

  BsInd cind;
  std::vector< std::vector<size_t> > ref;
  std::vector<size_t> vref;
  ref = { {0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0}};
  assert(ref == kmesh.kPoints);

  for (auto tt: t){
    for (auto e: tt) printf("%ld ", e);
    printf("\n");
  }

  auto flip(t);
  std::sort(flip.begin(), flip.end(), ibst::compare({1,2,0}));
  printf("----\n");


  for (auto ff: flip){
    for (auto e: ff) printf("%ld ", e);
    printf("\n");
  }

  std::vector<size_t> remap;
  for (auto e: flip) remap.push_back(e.back());
  cind = getIndices(flip, t);
  cind = getIndices(flip, t, "bla", true);
  cind = getIndices(f, "abij", f, "baji");
  cind = getIndices(f, "abij", o, "b");
  cind = getIndices(t, "Gai", t, "Gia", remap);
  cind = getIndices(f, "abij", f, "akic", f, "cbkj");
  //cind = ibst::
  return 0;

}
