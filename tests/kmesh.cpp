#include "../src/kmesh.h"
#include <cassert>
#include <iostream>


int main(){

  auto mesh = ibst::Kmesh({3,2,1});
  std::vector< std::vector<size_t> > ref;
  std::vector<size_t> vref;
  ref = { {0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0}};

  assert(ref == mesh.kPoints);

  auto zero  = mesh.getDefaultNzc(0);
  auto one   = mesh.getDefaultNzc(1);
  auto two   = mesh.getDefaultNzc(2);
  auto three = mesh.getDefaultNzc(3);
  auto four  = mesh.getDefaultNzc(4);

  ref = {{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}, {5, 5}};
  assert(ref == one);

  ref = {{0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {3, 3, 3}, {4, 4, 4}, {5, 5, 5}};
  assert(ref == two);

  assert(three.size() == 36);
  vref = {0, 1, 1, 7};
  assert(three[7] == vref);
  vref = {3, 0, 3, 18};
  assert(three[18] == vref);
  vref = {5, 2, 4, 26};
  assert(three[26] == vref);


  assert(four.size() == 216);
  vref = {5, 3, 5, 3, 203};
  assert(four[203] == vref);
  vref = {4, 4, 2, 0, 170};
  assert(four[170] == vref);
  vref = {3, 4, 0, 1, 132};
  assert(four[132] == vref);
  vref = {1, 2, 2, 1, 50};
  assert(four[50] == vref);
  vref = {0, 2, 2, 0, 14};
  assert(four[14] == vref);

/*
  for (auto o: one){
    for (auto e: o) printf("%ld ", e);
    printf(": ");
    for (auto e: o){
      auto v(mesh.idxToK(e));
      printf("( ");
      for (auto vv: v) printf("%ld ", vv);
      printf(")");
    }
    printf("\n");
  }
  for (auto t: two){
    for (auto e: t) printf("%ld ", e);
    printf(": ");
    for (auto e: t){
      auto v(mesh.idxToK(e));
      printf("( ");
      for (auto vv: v) printf("%ld ", vv);
      printf(")");
    }
    printf("\n");
  }
  for (auto t: three){
    for (auto e: t) printf("%ld ", e);
    printf(": ");
    for (auto e: t){
      auto v(mesh.idxToK(e));
      printf("( ");
      for (auto vv: v) printf("%ld ", vv);
      printf(")");
    }
    printf("\n");
  }

  for (auto f: four){
    for (auto e: f) printf("%ld ", e);
    printf(": ");
    for (auto e: f){
      auto v(mesh.idxToK(e));
      printf("( ");
      for (auto vv: v) printf("%ld ", vv);
      printf(")");
    }
    printf("\n");
  }
  */
  return 0;
}
