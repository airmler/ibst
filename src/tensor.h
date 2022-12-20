#ifndef IBST_TENSOR_H
#define IBST_TENSOR_H

#include <ctf.hpp>

#include "util.h"
#include "bsIndices.h"


//TODO make tensor to a template argument and remove all CTF related stuff
namespace ibst {


template <typename F=double>
class Tensor {
 public:

  Tensor( size_t const                               _order
        , std::vector<size_t> const                  _lens
        , std::vector< std::vector<size_t> > const   _nzc
        , CTF::World                                *_world
        , std::string                        const   _name = "dummy"
        ) {
    order = _order;
    lens = _lens;
    nzc = _nzc;
    nBlocks = nzc.size();
    name = _name;
    world = _world;
    elements = std::accumulate( this->lens.begin()
                              , this->lens.end()
                              , 1UL
                              , std::multiplies<size_t>()
                              );
    //, world(_world) {
    // Do some santiy checks of the provided tensor data
    assert(lens.size() == order);
    assert(nBlocks > 0);
    std::vector<size_t> idx, ref(nBlocks);
    idx.reserve(nBlocks);
    std::iota(ref.begin(), ref.end(), 0);
    std::vector<int> sym(order);
    for (auto n: nzc) {
      assert(n.size() == order + 1);
      idx.push_back(n.back());
      tensors.push_back(
        new CTF::Tensor<F>( (int) order, (int64_t const *) lens.data()
                          , sym.data(), *_world, name.c_str()
                          )
      );
    }
    std::sort(idx.begin(), idx.end());
    assert(idx.size() == ref.size());
    assert(idx == ref);
  };

  ~Tensor() {
    for (auto &t: tensors) delete t;
  }

  // TODO read/write without an (blockSparse-Block) offset in the file.
  void read_dense_from_file(MPI_File & file, int64_t block = -1) {

    // Write all blocks (-1) or just one specific nzc-Block of the whole Tensor
    if (block < 0) {
      int64_t offset(0L);
      for (auto &t: this->tensors) {
        t->read_dense_from_file(file, offset*sizeof(F));
        offset += elements;
      }
    }
    else {
      assert(block < this->nBlocks);
      this->tensors[block]->read_dense_from_file(
        file, block*elements*sizeof(F)
      );
    }
  }

  // TODO read/write without an (blockSparse-Block) offset in the file.
  void write_dense_to_file(MPI_File & file, int64_t block = -1) {
    // Write all blocks (-1) or just one specific nzc-Block of the whole Tensor
    if (block < 0) {
      int64_t offset(0L);
      for (auto &t: this->tensors) {
        t->write_dense_to_file(file, offset*sizeof(F));
        offset += elements;
      }
    }
    else {
      assert(block < this->nBlocks);
      this->tensors[block]->write_dense_to_file(file, block*elements*sizeof(F));
    }
  }

  void write(size_t                    npair,
             std::vector<size_t> const global_idx,
             std::vector<F>      const data,
             int64_t const             block = -1)
  {
    if (block < 0) {
      // global_idx is identical in all blocks!!!
      size_t i(0);
      for (auto &t: this->tensors) {
        t->write( (int64_t ) npair
                , (const int64_t *) global_idx.data()
                , (const F *) &data[(i++)*npair]
                );
      }
    } else {
      assert(block < this->nBlocks);
      this->tensors[block]->write( (int64_t) npair
                                 , (const int64_t *) global_idx.data()
                                 , (const F*) data.data()
                                 );
    }
  }


  void read(size_t                     npair,
            std::vector<size_t> const  global_idx,
            std::vector<F>             data,
            int64_t const              block = -1)
  {
    //TODO check if the data pointer is large enough!
    if (block < 0) {
      // global_idx is identical in all blocks!!!
      size_t i(0);
      for (auto &t: this->tensors)
        t->read( (int64_t) npair
               , (const int64_t *) global_idx.data()
               , (F *) &data[(i++)*npair]
               );
    } else {
      this->tensors[block]->read( (int64_t) npair
                                , (const int64_t *) global_idx.data()
                                , (F *) data.data()
                                );
    }
  }

  void read(F &data) {
    assert(nBlocks == 1);
    assert(elements == 1);
    std::vector<int64_t> global_idx(1);
    std::vector<F> vdata(1);
    int64_t npair(1);
    this->tensors[0]->read_all(vdata.data());
    data = vdata[0];
  }



  void sum( F                          alpha
          , const Tensor              &A
          , const std::string          cidxA
          , F                          beta
          , const std::string          cidxB
          , const std::function<F(const F)> &fseq
          , bool                       verbose = false
          ) {
    std::string cs(name + "[" + cidxB + "] = f(" + A.name + "[" + cidxA + "])");
    auto cind = getIndices(this->nzc, cidxB, A.nzc, cidxA, cs);
    if (verbose) cind.print();
    for (auto t: cind.tasks) {
      this->tensors[t[0]]->sum( alpha
                             , *(A.tensors[t[1]])
                             , cidxA.c_str()
                             , beta
                             , cidxB.c_str()
                             , CTF::Univar_Function<F>(fseq)
                             );
    }

  }

  void sum( F                 alpha
          , const Tensor     &A
          , const std::string cidxA
          , F                 beta
          , const std::string cidxB
          , bool              verbose = false
          ) {
    std::string cs(name + "[" + cidxB + "] = " + A.name + "[" + cidxA + "]");
    auto cind = getIndices(this->nzc, cidxB, A.nzc, cidxA, cs);
    if (verbose) cind.print();
    for (auto t: cind.tasks) {
      this->tensors[t[0]]->sum( alpha
                             , *(A.tensors[t[1]])
                             , cidxA.c_str()
                             , beta
                             , cidxB.c_str()
                             );
    }
  }

  void sum( F                          alpha
          , const Tensor              &A
          , const std::string          cidxA
          , F                          beta
          , const std::string          cidxB
          , std::vector<size_t>        remap
          , const std::function<F(const F)> &fseq
          , bool                       verbose = false
          ) {
    std::string cs(name + "[" + cidxB + "] = f(" + A.name + "[" + cidxA + "])");
    auto cind = getIndices(this->nzc, cidxB, A.nzc, cidxA, remap, cs);
    if (verbose) cind.print();

    for (auto t: cind.tasks) {
      this->tensors[t[0]]->sum( alpha
                             , *(A.tensors[t[1]])
                             , cidxA.c_str()
                             , beta
                             , cidxB.c_str()
                             , CTF::Univar_Function<F>(fseq)
                             );
    }

  }

  void sum( F                   alpha
          , const Tensor       &A
          , const std::string   cidxA
          , F                   beta
          , const std::string   cidxB
          , std::vector<size_t> remap
          , bool                verbose = false
          ) {
    std::string cs(name + "[" + cidxB + "] = " + A.name + "[" + cidxA + "]");
    auto cind = getIndices(this->nzc, cidxB, A.nzc, cidxA, remap, cs);
    if (verbose) cind.print();
    for (auto t: cind.tasks) {
      this->tensors[t[0]]->sum( alpha
                             , *(A.tensors[t[1]])
                             , cidxA.c_str()
                             , beta
                             , cidxB.c_str()
                             );

    }
  }

  void contract( F alpha
               , const Tensor &A
               , const std::string cidxA
               , const Tensor &B
               , const std::string cidxB
               , F beta
               , const std::string cidxC
               , bool verbose = false
               ) {
    std::string cs( name + "[" + cidxC + "] = " + A.name + "[" + cidxA + "] * "
                  + B.name + "[" + cidxB + "]");
    auto cind = getIndices(this->nzc, cidxC, A.nzc, cidxA, B.nzc, cidxB, cs);
    if (verbose) cind.print();
    // If beta is not equal one. We scale the target Tensor once
    // and set beta to one
    if (beta != (F) 1) {
      for (auto t: this->tensors) t->set_zero();
      beta = (F) 1;
    }
    for (auto t: cind.tasks) {
      tensors[t[0]]->contract( alpha
                                  , *(A.tensors[t[1]])
                                  , cidxA.c_str()
                                  , *(B.tensors[t[2]])
                                  , cidxB.c_str()
                                  , beta
                                  , cidxC.c_str()
                                  );
    }
  }

  void slice(
    const std::vector<size_t> begins,
    const std::vector<size_t> ends,
    F                         beta,
    Tensor                   &A,
    const std::vector<size_t> aBegins,
    const std::vector<size_t> aEnds,
    F                         alpha,
    bool                      verbose = false
  ) {
    assert(nBlocks == A.nBlocks);
    for (int i(0); i < this->nBlocks; i++) assert(this->nzc[i] == A.nzc[i]);

    std::string cs(name + " = " + A.name);
    auto cind = getIndices(this->nzc, A.nzc, cs);
    if (verbose) cind.print();
    for (auto t: cind.tasks) {
      this->tensors[t[0]]->slice( (const int64_t *) begins.data()
                                , (const int64_t *) ends.data()
                                , beta
                                , *(A.tensors[t[1]])
                                , (const int64_t *) aBegins.data()
                                , (const int64_t *) aEnds.data()
                                , alpha
                                );
    }
  }


  void relabelBlocks( std::vector< std::vector<size_t> > nzcOut
                    , bool              verbose = false
                    ) {
    std::string cs(name);
    auto cind = getIndices(this->nzc, nzcOut, cs, true);
    if (verbose) cind.print();
    for (auto t: cind.tasks) this->nzc[t[0]] = nzcOut[t[1]];
  }

  //TODO make tensor to a template argument
  std::vector<CTF::Tensor<F> * > tensors;
  CTF::World *world;
  std::vector< std::vector<size_t> > nzc;
  std::vector<size_t> lens;
  size_t nBlocks;
  size_t order;
  size_t elements;
  std::string name;
};
} // namespace ibst

#endif /* IBST_TENSOR_H */

