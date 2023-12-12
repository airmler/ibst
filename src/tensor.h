#ifndef IBST_TENSOR_H
#define IBST_TENSOR_H

#include "bs_indices.h"
#include "tracker.h"
#include "util.h"
#include "world.h"

#include "ctf_tensor.h"

namespace ibst {

template <typename F = double>
class Tensor {
public:
  Tensor(int64_t const _order,
         std::vector<int64_t> const _lens,
         std::vector<std::vector<int64_t>> const _nzc,
         World *_world,
         std::string const _name = "dummy")
      : order(_order)
      , lens(_lens)
      , nzc(_nzc)
      , name(_name) {
    nBlocks = nzc.size();
    world = _world;
    elements = std::accumulate(this->lens.begin(),
                               this->lens.end(),
                               1L,
                               std::multiplies<int64_t>());
    int64_t allocMem = elements * nBlocks * (traits::is_complex<F>() ? 16 : 8);

    std::string line(tensor_to_line(name, lens, nzc.size()));
    auto event =
        Event().operation("Initialize").memoryUsage(allocMem).line(line);
    world->tracker->add(event);
    // Do some santiy checks of the provided tensor data
    assert(lens.size() == order);
    assert(nBlocks > 0);
    std::vector<int64_t> idx, ref(nBlocks);
    idx.reserve(nBlocks);
    std::iota(ref.begin(), ref.end(), 0);
    std::vector<int> sym(order);
    for (auto n : nzc) {
      assert(n.size() == order + 1);
      idx.push_back(n.back());
      tensors.push_back(new Blocked_tensor<F>((int)order,
                                              lens.data(),
                                              sym.data(),
                                              *(_world->wrld),
                                              name.c_str()));
    }
    std::sort(idx.begin(), idx.end());
    assert(idx.size() == ref.size());
    assert(idx == ref);
  }

  ~Tensor() {
    int64_t freeMem = -elements * nBlocks * (traits::is_complex<F>() ? 16 : 8);
    std::string line(tensor_to_line(name, lens, nzc.size()));
    auto event = Event().operation("Destroy").memoryUsage(freeMem).line(line);
    world->tracker->add(event);
    for (auto &t : tensors) delete t;
  }

  // TODO read/write without an (blockSparse-Block) offset in the file.
  void read_dense_from_file(MPI_File &file, int64_t block = -1) {
    auto fc = (block < 0) ? elements * nBlocks : elements;

    auto start = chrono::high_resolution_clock::now();

    // Write all blocks (-1) or just one specific nzc-Block of the whole Tensor
    if (block < 0) {
      int64_t offset(0L);
      for (auto &t : this->tensors) {
        t->read_dense_from_file(file, offset * sizeof(F));
        offset += elements;
      }
    } else {
      assert(block < this->nBlocks);
      this->tensors[block]->read_dense_from_file(file,
                                                 block * elements * sizeof(F));
    }
    auto end = chrono::high_resolution_clock::now();
    auto time =
        chrono::duration_cast<chrono::microseconds>(end - start).count() * 1e-6;
    std::string line(tensor_to_line(name, lens, nzc.size()));
    auto event =
        Event().operation("ReadFile").flopCount(fc).time(time).line(line);
    world->tracker->add(event);
  }

  // TODO read/write without an (blockSparse-Block) offset in the file.
  void write_dense_to_file(MPI_File &file, int64_t block = -1) {
    auto fc = (block < 0) ? elements * nBlocks : elements;
    ///    Event event{.operation = "WriteToFile", .flopCount = fc};
    auto start = chrono::high_resolution_clock::now();

    // Write all blocks (-1) or just one specific nzc-Block of the whole Tensor
    if (block < 0) {
      int64_t offset(0L);
      for (auto &t : this->tensors) {
        t->write_dense_to_file(file, offset * sizeof(F));
        offset += elements;
      }
    } else {
      assert(block < this->nBlocks);
      this->tensors[block]->write_dense_to_file(file,
                                                block * elements * sizeof(F));
    }
    auto end = chrono::high_resolution_clock::now();
    auto time =
        chrono::duration_cast<chrono::microseconds>(end - start).count() * 1e-6;
    std::string line(tensor_to_line(name, lens, nzc.size()));
    auto event = Event().operation("WriteFile").flopCount(fc).time(time).line(line);
    world->tracker->add(event);
  }

  void write(int64_t npair,
             std::vector<int64_t> const global_idx,
             std::vector<F> const data,
             int64_t const block = -1) {
    auto fc = (block < 0) ? npair * nBlocks : npair;
    auto start = chrono::high_resolution_clock::now();

    if (block < 0) {
      // global_idx is identical in all blocks!!!
      int64_t i(0);
      for (auto &t : this->tensors) {
        t->write(npair, global_idx.data(), (const F *)&data[(i++) * npair]);
      }
    } else {
      assert(block < this->nBlocks);
      this->tensors[block]->write(npair,
                                  global_idx.data(),
                                  (const F *)data.data());
    }
    auto end = chrono::high_resolution_clock::now();
    auto time =
        chrono::duration_cast<chrono::microseconds>(end - start).count() * 1e-6;
    std::string line(tensor_to_line(name, lens, nzc.size()));
    auto event = Event().operation("Write").flopCount(fc).time(time).line(line);
    world->tracker->add(event);
  }

  void read(int64_t npair,
            std::vector<int64_t> const global_idx,
            std::vector<F> data,
            int64_t const block = -1) {
    auto fc = (block < 0) ? npair * nBlocks : npair;

    auto start = chrono::high_resolution_clock::now();
    // TODO check if the data pointer is large enough!
    if (block < 0) {
      // global_idx is identical in all blocks!!!
      int64_t i(0);
      for (auto &t : this->tensors)
        t->read(npair, global_idx.data(), (F *)&data[(i++) * npair]);
    } else {
      this->tensors[block]->read(npair, global_idx.data(), (F *)data.data());
    }
    auto end = chrono::high_resolution_clock::now();
    auto time =
        chrono::duration_cast<chrono::microseconds>(end - start).count() * 1e-6;
    std::string line(tensor_to_line(name, lens, nzc.size()));
    auto event = Event().operation("Read").flopCount(fc).time(time).line(line);
    world->tracker->add(event);
  }

  void read(F &data) {
    assert(nBlocks == 1);
    assert(elements == 1);
    std::vector<int64_t> global_idx(1);
    std::vector<F> vdata(1);
    int64_t npair(1);
    ///    Event event{.operation = "Read", .flopCount = 1};
    this->tensors[0]->read_all(vdata.data());
    data = vdata[0];
  }

  void sum(F alpha,
           const Tensor &A,
           const std::string cidxA,
           F beta,
           const std::string cidxB,
           const std::function<F(const F)> &fseq,
           bool verbose = false) {
    std::string cs(name + "[" + cidxB + "] = f(" + A.name + "[" + cidxA + "])");
    auto cind = getIndices(this->nzc, cidxB, A.nzc, cidxA, cs);
    if (verbose)
      cind.print();
    // if we dont allow dublicates all elements of the tensor are involved
    double ops(traits::is_complex<F>() ? 2.0 : 1.0);
    ops *= elements;
    ops *= cind.tasks.size();

    auto start = chrono::high_resolution_clock::now();
    for (auto t : cind.tasks) {
      this->tensors[t[0]]->sum(alpha,
                               *(A.tensors[t[1]]),
                               cidxA.c_str(),
                               beta,
                               cidxB.c_str(),
                               (const Univar_function<F>)fseq);
    }
    auto end = chrono::high_resolution_clock::now();
    auto time =
        chrono::duration_cast<chrono::microseconds>(end - start).count() * 1e-6;
    auto event = Event().operation("Sum").flopCount(ops).time(time).line(cs);
    world->tracker->add(event);
  }

  void sum(F alpha,
           const Tensor &A,
           const std::string cidxA,
           F beta,
           const std::string cidxB,
           bool verbose = false) {
    std::string cs(name + "[" + cidxB + "] = " + A.name + "[" + cidxA + "]");
    auto cind = getIndices(this->nzc, cidxB, A.nzc, cidxA, cs);
    if (verbose)
      cind.print();
    // if we dont allow dublicates all elements of the tensor are involved
    double ops(traits::is_complex<F>() ? 2.0 : 1.0);
    ops *= elements;
    ops *= cind.tasks.size();

    auto start = chrono::high_resolution_clock::now();
    for (auto t : cind.tasks) {
      this->tensors[t[0]]->sum(alpha,
                               *(A.tensors[t[1]]),
                               cidxA.c_str(),
                               beta,
                               cidxB.c_str());
    }
    auto end = chrono::high_resolution_clock::now();
    auto time =
        chrono::duration_cast<chrono::microseconds>(end - start).count() * 1e-6;
    auto event = Event().operation("Sum").flopCount(ops).time(time).line(cs);
    world->tracker->add(event);
  }

  void sum(F alpha,
           const Tensor &A,
           const std::string cidxA,
           F beta,
           const std::string cidxB,
           std::vector<int64_t> remap,
           const std::function<F(const F)> &fseq,
           bool verbose = false) {
    std::string cs(name + "[" + cidxB + "] = f(" + A.name + "[" + cidxA + "])");
    auto cind = getIndices(this->nzc, cidxB, A.nzc, cidxA, remap, cs);
    if (verbose)
      cind.print();
    // if we dont allow dublicates all elements of the tensor are involved
    double ops(traits::is_complex<F>() ? 2.0 : 1.0);
    ops *= elements;
    ops *= cind.tasks.size();

    auto start = chrono::high_resolution_clock::now();
    for (auto t : cind.tasks) {
      this->tensors[t[0]]->sum(alpha,
                               *(A.tensors[t[1]]),
                               cidxA.c_str(),
                               beta,
                               cidxB.c_str(),
                               (const Univar_function<F>)fseq);
    }
    auto end = chrono::high_resolution_clock::now();
    auto time =
        chrono::duration_cast<chrono::microseconds>(end - start).count() * 1e-6;
    auto event = Event().operation("Sum").flopCount(ops).time(time).line(cs);
    world->tracker->add(event);
  }

  void sum(F alpha,
           const Tensor &A,
           const std::string cidxA,
           F beta,
           const std::string cidxB,
           std::vector<int64_t> remap,
           bool verbose = false) {
    std::string cs(name + "[" + cidxB + "] = " + A.name + "[" + cidxA + "]");
    auto cind = getIndices(this->nzc, cidxB, A.nzc, cidxA, remap, cs);
    if (verbose)
      cind.print();
    // if we dont allow dublicates all elements of the tensor are involved
    double ops(traits::is_complex<F>() ? 2.0 : 1.0);
    ops *= elements;
    ops *= cind.tasks.size();

    auto start = chrono::high_resolution_clock::now();
    for (auto t : cind.tasks) {
      this->tensors[t[0]]->sum(alpha,
                               *(A.tensors[t[1]]),
                               cidxA.c_str(),
                               beta,
                               cidxB.c_str());
    }
    auto end = chrono::high_resolution_clock::now();
    auto time =
        chrono::duration_cast<chrono::microseconds>(end - start).count() * 1e-6;
    auto event = Event().operation("Sum").flopCount(ops).time(time).line(cs);
    world->tracker->add(event);
  }

  void contract(F alpha,
                const Tensor &A,
                const std::string cidxA,
                const Tensor &B,
                const std::string cidxB,
                F beta,
                const std::string cidxC,
                bool verbose = false) {
    assert(cidxC.size() == this->lens.size() && cidxA.size() == A.lens.size() && cidxB.size() == B.lens.size());
    std::string cs(name + "[" + cidxC + "] = " + A.name + "[" + cidxA + "] * "
                   + B.name + "[" + cidxB + "]");
    auto cind = getIndices(this->nzc, cidxC, A.nzc, cidxA, B.nzc, cidxB, cs);
    if (verbose)
      cind.print();
    // currently we dont allow dublicated indices on one tensor
    double ops(traits::is_complex<F>() ? 8 : 2);
    //ops *= traits::is_complex<F>() ? 8 : 2;
    std::string idx(cidxA+cidxB+cidxC);
    std::sort(idx.begin(), idx.end());
    idx.erase( std::unique( idx.begin(), idx.end() ), idx.end() );
    std::map<char, int64_t> dimMap;
    for (size_t c(0); c < lens.size(); c++)   dimMap[cidxC[c]] = lens[c];
    for (size_t a(0); a < A.lens.size(); a++) dimMap[cidxA[a]] = A.lens[a];
    for (size_t b(0); b < B.lens.size(); b++) dimMap[cidxB[b]] = B.lens[b];

    for (auto c: idx) ops *= dimMap[c];
    ops *= cind.tasks.size();

    auto start = chrono::high_resolution_clock::now();
    // If beta is not equal one. We scale the target Tensor once
    // and set beta to one. TODO: not correct for all beta
    if (beta != (F)1) {
      for (auto t : this->tensors) t->set_zero();
      beta = (F)1;
    }
    for (auto t : cind.tasks) {
      tensors[t[0]]->contract(alpha,
                              *(A.tensors[t[1]]),
                              cidxA.c_str(),
                              *(B.tensors[t[2]]),
                              cidxB.c_str(),
                              beta,
                              cidxC.c_str());
    }
    auto end = chrono::high_resolution_clock::now();
    auto time =
        chrono::duration_cast<chrono::microseconds>(end - start).count() * 1e-6;
    auto event = Event().operation("Contract").flopCount(ops).time(time).line(cs);
    world->tracker->add(event);
  }

  void slice(const std::vector<int64_t> begins,
             const std::vector<int64_t> ends,
             F beta,
             Tensor &A,
             const std::vector<int64_t> aBegins,
             const std::vector<int64_t> aEnds,
             F alpha,
             bool verbose = false) {
    assert(nBlocks == A.nBlocks);
    for (int i(0); i < this->nBlocks; i++) assert(this->nzc[i] == A.nzc[i]);

    std::string cs(name + " = " + A.name);
    auto cind = getIndices(this->nzc, A.nzc, cs);
    if (verbose)
      cind.print();
    double ops(traits::is_complex<F>() ? 2.0 : 1.0);
    for (size_t i(0); i < begins.size(); i++) ops *= ends[i] - begins[i];
    ops *= cind.tasks.size();

    auto start = chrono::high_resolution_clock::now();
    for (auto t : cind.tasks) {
      this->tensors[t[0]]->slice(begins.data(),
                                 ends.data(),
                                 beta,
                                 *(A.tensors[t[1]]),
                                 aBegins.data(),
                                 aEnds.data(),
                                 alpha);
    }
    auto end = chrono::high_resolution_clock::now();
    auto time =
        chrono::duration_cast<chrono::microseconds>(end - start).count() * 1e-6;
    auto event = Event().operation("Slice").flopCount(ops).time(time).line(cs);
    world->tracker->add(event);
  }

  void relabelBlocks(std::vector<std::vector<int64_t>> nzcOut,
                     bool verbose = false) {
    std::string cs(name);
    auto cind = getIndices(this->nzc, nzcOut, cs, true);
    if (verbose)
      cind.print();
    for (auto t : cind.tasks) this->nzc[t[0]] = nzcOut[t[1]];
  }

  std::vector<Blocked_tensor<F> *> tensors;
  std::vector<std::vector<int64_t>> nzc;
  std::vector<int64_t> lens;
  int64_t nBlocks;
  int64_t order;
  int64_t elements;
  std::string name;
  World *world;
};
} // namespace ibst

#endif /* IBST_TENSOR_H */
