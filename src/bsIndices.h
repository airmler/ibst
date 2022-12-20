#ifndef IBST_BSINDICES_H
#define IBST_BSINDICES_H

#include "util.h"
#include <numeric>
#include <set>

namespace ibst {

class BsInd {
 public:
  std::vector< std::vector<size_t> > tasks;
  std::vector< std::vector<size_t> > nzcC, nzcA, nzcB;
  std::string cs;

  BsInd() {};
  BsInd(std::string _cs) { cs = _cs;}

  std::string nzcToString(std::vector<size_t> nzc) {
    std::string out;
    for (size_t i(0); i < nzc.size() - 1; i++) out += _FORMAT("%ld ", nzc[i]);
    out += _FORMAT("(%ld)", nzc.back());
    return out;
  }

  std::string toStringContract() {
    std::string out;

    auto sizeC(std::set<std::vector<size_t>>(nzcC.begin(), nzcC.end()).size());
    auto sizeB(std::set<std::vector<size_t>>(nzcB.begin(), nzcB.end()).size());
    auto sizeA(std::set<std::vector<size_t>>(nzcA.begin(), nzcA.end()).size());

    out += _FORMAT( "%ld : %ld %ld %ld\n"
                  , tasks.size(), sizeC, sizeA, sizeB
                  );

    for (size_t i(0); i < nzcC.size(); i++){
      out += std::string(nzcToString(nzcC[i])) + " = ";
      out += std::string(nzcToString(nzcA[i])) + " x ";
      out += std::string(nzcToString(nzcB[i])) + '\n';
//      if (i+1 == nzcC.size()); // do nothing for the last element
//			else if (nzcC[i] == nzcC[i+1]) out += " + ";
//      else out += '\n' + std::string(nzcToString(nzcC[i])) + " = ";
    }
    out.pop_back();
    return out;
  }

  std::string toStringSum() {
    std::string out;

    auto sizeB(std::set<std::vector<size_t>>(nzcB.begin(), nzcB.end()).size());
    auto sizeA(std::set<std::vector<size_t>>(nzcA.begin(), nzcA.end()).size());

    out += _FORMAT("%ld : %ld %ld\n", tasks.size(), sizeB, sizeA);

    for (size_t i(0); i < nzcB.size(); i++){
      out += std::string(nzcToString(nzcB[i])) + " = ";
      out += std::string(nzcToString(nzcA[i])) + '\n';
    }
    out.pop_back();
    return out;
  }

  void print() {
    ILOG(cs);
    ILOG("Non-zero conditions:");
    if (nzcC.size()) {
      ILOG(toStringContract());
    } else {
      ILOG(toStringSum());
    }
    ILOG("+++++++");
  }
};

inline BsInd getIndices( std::vector< std::vector<size_t>> _nzcB
                       , std::vector< std::vector<size_t>> _nzcA
                       , std::string _cs = ""
                       , bool enforce = false
                       ) {
  BsInd cind(_cs);
  assert(_nzcA.size() == _nzcB.size());
  if (!enforce) {
    std::sort(_nzcA.begin(), _nzcA.end());
    std::sort(_nzcB.begin(), _nzcB.end());
  }
  cind.tasks.resize(_nzcA.size(), std::vector<size_t>(2));
  for (size_t i(0); i < _nzcB.size(); i++) {
    cind.tasks[i] = {_nzcB[i].back(), _nzcA[i].back()};
    assert(_nzcB[i].size() == _nzcA[i].size());
    if (!enforce) assert(_nzcB[i] == _nzcA[i]);
    cind.nzcB.push_back(_nzcB[i]);
    cind.nzcA.push_back(_nzcA[i]);
  }
  return cind;
}


inline BsInd getIndices( std::vector< std::vector<size_t> > _nzcB
                       , const std::string &_cidxB
                       , std::vector< std::vector<size_t> > _nzcA
                       , const std::string &_cidxA
                       , std::string _cs = ""
                       ) {
  BsInd cind(_cs);
  checkDuplicate(_cidxA);
  checkDuplicate(_cidxB);
  std::sort(_nzcB.begin(), _nzcB.end());
  // analyze the cidx of matrix A and matrix B;
  // if the rhs apears on the lhs store the position in the vector idx
  std::vector<size_t> idx;
  for (auto c: _cidxA) {
    auto it = std::find(_cidxB.begin(), _cidxB.end(), c);
    if (it != std::end(_cidxB)) {
      idx.push_back(std::distance( _cidxB.begin(), it ));
    }
  }

  // case 1: B["abc"] = A["cab"];
  // it is assumed that lhs and rhs have both "sane" non zero conditions.
  // if not, a correct behaviour is not guaranteed
  if (_cidxB.size() == _cidxA.size()) {
    // as lhs is a permutation of rhs
    // they should have the same number of non zero conditions
    assert(_nzcA.size() == _nzcB.size());
    //make sure that cidxB is any permutation of cidxA!!
    //assert( idx.size() == B.order );
    std::sort(_nzcA.begin(), _nzcA.end(), ibst::compare(idx));
    for (size_t i(0); i < _nzcA.size(); i++) {
      cind.tasks.push_back({_nzcB[i].back(), _nzcA[i].back()});
      cind.nzcB.push_back(_nzcB[i]);
      cind.nzcA.push_back(_nzcA[i]);
    }
  } else if (_cidxB.size() > _cidxA.size()) {
    // case 2: B["abij"] = A["i"];
    // idea is to take the nzc-map from B for matrix A.
    // then we remove all elements from this map which do not appear in A
    auto nzc(_nzcB);
    std::vector<int> iidx(_cidxB.size());
    std::iota(iidx.begin(), iidx.end(), 0);
    std::vector<int> toRemove;
    // get all elements in B which do not appear in A
    // In the upper example: this would be ["ab j"] ie: 0, 1, 3
    std::set_difference( iidx.begin(), iidx.end()
                       , idx.begin(), idx.end()
                       , std::back_inserter(toRemove)
                       );
    size_t id(0);
    for (auto &n: nzc) {
      // mark all elements to be removed with -1 and delete these elements
      for (auto t: toRemove) n[t] = -1;
      n.erase( std::remove(n.begin(), n.end(), -1), n.end());

      // now we have the correct nzc maps on the lhs but with the wrong index
      // match the nzc with the A.nzc and take over the last element (index)
      bool replaced(false);
      for (auto &a: _nzcA)
      if (std::equal( a.begin(), a.end() - 1, n.begin())){
        n.back() = a.back();
        replaced = true;
      }
      assert(replaced);
      cind.tasks.push_back({_nzcB[id].back(), nzc[id].back()});
      cind.nzcB.push_back(_nzcB[id]);
      cind.nzcA.push_back(nzc[id]);
      id++;
    }
  } else {  // not implemented yet: B["a"] = A["aj"] which is valid einsum code
    assert(0);
  }
  return cind;
}

inline BsInd getIndices( std::vector< std::vector<size_t> > _nzcB
                       , const std::string &_cidxB
                       , std::vector< std::vector<size_t> > _nzcA
                       , const std::string &_cidxA
                       , std::vector<size_t> remap
                       , std::string _cs = ""
                       ) {
  BsInd cind(_cs);
  checkDuplicate(_cidxB);
  checkDuplicate(_cidxA);
  assert(_cidxB.size() == _cidxA.size());
  assert(_nzcB.size()  == _nzcA.size());
  assert(_nzcA.size() == remap.size());
  cind.tasks.resize(_nzcB.size(), std::vector<size_t>(2));
  //write tensor B as it is
  for (size_t id(0); id < _nzcB.size(); id++) {
    auto p = std::find_if( _nzcB.begin()
                         , _nzcB.end()
                         , [id](std::vector<size_t> i){return i.back() == id;}
                         );
    cind.nzcB.push_back(*p);
    cind.tasks[id][0] = std::distance(_nzcB.begin(), p);
  }
  // Find nzcA in Tensor A
  for (size_t id(0); id < _nzcA.size(); id++) {
    auto t(remap[id]);
    auto p = std::find_if( _nzcA.begin()
                         , _nzcA.end()
                         , [t](std::vector<size_t> i) {return i.back() == t;}
                         );
	  assert(p != _nzcA.end());
    cind.nzcA.push_back(*p);
    cind.tasks[id][1] = std::distance(_nzcA.begin(), p);
  }

  return cind;
}

inline BsInd getIndices( std::vector< std::vector<size_t> > _nzcC
                       , std::string _cidxC
                       , std::vector< std::vector<size_t> > _nzcA
                       , std::string _cidxA
                       , std::vector< std::vector<size_t> > _nzcB
                       , std::string _cidxB
                       , std::string _cs = ""
                       ) {
  BsInd cind(_cs);

  //TODO
  // still some sane tensor contractions are not supported:
  // C["ij"] = A["ik"] * B["kkj"]; C["ijl"] = A["ikl"] * B["kjl"]
  checkDuplicate(_cidxC);
  checkDuplicate(_cidxA);
  checkDuplicate(_cidxB);
  // All nzc-maps are sorted. A && B are sorted such that the indices
  // which appear on the lhs are the slowest
  // idxA && idxB specifies how to sort
  std::vector<size_t> idxA, idxB;
  std::vector< std::pair<size_t, size_t> > ca, cb, ab;
  // Find the indices of A && B which appear on C and write it to idxA/idxB
  // ca, cb, ab is list of matching indices between two tensors
  for (size_t i(0); i < _cidxC.size(); i++) {
    for (size_t j(0); j < _cidxA.size(); j++)
    if ( _cidxA[j] == _cidxC[i]) {
      idxA.push_back(j);
      ca.push_back({i,j});
    }
    for (size_t j(0); j < _cidxB.size(); j++)
    if ( _cidxB[j] == _cidxC[i]) {
      idxB.push_back(j);
      cb.push_back({i,j});
    }
  }
  // write the indices which are contracted into idxA/idxB
  for (size_t i(0); i < _cidxA.size(); i++)
  for (size_t j(0); j < _cidxB.size(); j++)
  if ( _cidxA[i] == _cidxB[j]) {
    ab.push_back({i,j});
    idxA.push_back(i);
    idxB.push_back(j);
  }
  // create temporary copies of the nzc-maps and sort according to idxA/idxB
  std::sort(_nzcA.begin(), _nzcA.end(), ibst::compare(idxA));
  std::sort(_nzcB.begin(), _nzcB.end(), ibst::compare(idxB));
  std::sort(_nzcC.begin(), _nzcC.end());
  //loop over all nzc of C
  //for a given nzcC look for all elements of A && B
  // which have the correct indices on the lhs
  for (auto &c: _nzcC) {
    auto beginA = std::find_if(_nzcA.begin(), _nzcA.end(), ibst::find(c, ca));
    auto endA   = std::find_if_not(beginA, _nzcA.end(),    ibst::find(c, ca));
    auto beginB = std::find_if(_nzcB.begin(), _nzcB.end(), ibst::find(c, cb));
    auto endB   = std::find_if_not(beginB, _nzcB.end(),    ibst::find(c, cb));
    for (auto itA(beginA); itA != endA; ++itA)
    for (auto itB(beginB); itB != endB; ++itB) {
      // check if ALL contracting indices agree.
      bool m(true);
      for (auto &x: ab) {
        if ((*itA)[x.first] != (*itB)[x.second]) m = false;
      }
      if (!m) continue;
      cind.tasks.push_back({c.back(), (*itA).back(), (*itB).back()});
      cind.nzcC.push_back(c);
      cind.nzcA.push_back(*itA);
      cind.nzcB.push_back(*itB);
    }
  }
  return cind;
}

}  // namespace ibst

#endif /* IBST_BSINDICES_H */
