#include "../src/util.h"
#include <iostream>
#include "../src/kmesh.h"


int main(){



  ibst::checkDuplicate("ij");
  ibst::checkDuplicate("3142");
  ibst::checkDuplicate("faih");


  auto kmesh(ibst::Kmesh({3,2,1}));
  auto a(kmesh.getDefaultNzc(4));


  std::vector< std::vector<size_t> > ref;
  std::vector<size_t> vref;
  ref = { {0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0}};
  assert(ref == kmesh.kPoints);

  std::vector< std::pair<size_t, size_t> > ca;
  std::vector<size_t> c;

  // Example: C["abcd"] = A["abef"] * ...
  vref = {1, 0, 5, 5, 41};
  assert(a[41] == vref);

  ca = {{0, 0}, {1, 1}};
  c = {0, 1, 2, 3};

  auto beginA = std::find_if(a.begin(), a.end(), ibst::find(c, ca));
  auto endA = std::find_if_not(beginA, a.end(), ibst::find(c, ca));

  assert(std::distance(a.begin(), beginA) == 6);
  assert(std::distance(a.begin(), endA) == 12);


  // Example: C["abcd"] = A["aebf"] * ...
  std::vector<size_t> idx = {0,2,1,3};
  std::sort(a.begin(), a.end(), ibst::compare(idx));

  vref = {0, 2, 1, 1, 13};
  assert(a[8] == vref);

  c = {0, 3, 2, 1};
  ca = {{0, 0}, {1, 2}};

  beginA = std::find_if(a.begin(), a.end(), ibst::find(c, ca));
  endA = std::find_if_not(beginA, a.end(), ibst::find(c, ca));

  assert(std::distance(a.begin(), beginA) == 18);
  assert(std::distance(a.begin(), endA) == 24);


  // Example: C["abcd"] = B["acef"] * A["bdfe"]
  idx = {1,3,2,0};
  std::sort(a.begin(), a.end(), ibst::compare(idx));

  for (auto aa: a){
    for (auto e: aa) printf("%ld ", e);
    printf("\n");
  }

  vref = {4, 0, 1, 3, 145};
  assert(a[19] == vref);

  c = {0, 0, 1, 2};
  ca = {{1, 0}, {3, 1}};

  beginA = std::find_if(a.begin(), a.end(), ibst::find(c, ca));
  endA = std::find_if_not(beginA, a.end(), ibst::find(c, ca));

  assert(std::distance(a.begin(), beginA) == 12);
  assert(std::distance(a.begin(), endA) == 18);


  return 0;

}
