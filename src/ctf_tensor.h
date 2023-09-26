#pragma once

#include <ctf.hpp>
using Machine_world = CTF::World;

template <typename F>
using Blocked_tensor = typename CTF::Tensor<F>;

template <typename F>
using Univar_function = typename CTF::Univar_Function<F>;

/*
template <typename F>
class Blocked_tensor {
  public:
    CTF::Tensor<F> *machine_tensor;

    Blocked_tensor(int64_t                                   order,
                   std::vector<int64_t> const                len,
                   CTF::World                                wrld,
                   std::string const                         name) {
      std::vector<int> sym(order);
      machine_tensor = new CTF::Tensor<F>(
        (int) order, len.data(), sym.data(), wrld, name
      );
    }

    ~Blocked_tensor() { delete machine_tensor; }

    void
    contract(F          alpha,
             Blocked_tensor     &A,
             std::string const idx_A,
             Blocked_tensor     &B,
             std::string const idx_B,
             F          beta,
             std::string const idx_C
            ) {

      machine_tensor->contract(alpha,
                               *A.machine_tensor,
                               idx_A.c_str(),
                               *B.machine_tensor,
                               idx_B.c_str(),
                               beta,
                               idx_C.c_str()
                              );

    }

    void sum(F                                   alpha,
             Blocked_tensor                    & A,
             std::string const                   idx_A,
             F                                   beta,
             std::string const                   idx_B,
             std::function<F(const F)>         & fseq
            ) {

      machine_tensor->sum(alpha,
                          *A.machine_tensor,
                          idx_A.c_str(),
                          beta,
                          idx_B.c_str(),
                          CTF::Univar_Function<F>(fseq)
                         );

    }

    void sum(F                                   alpha,
             Blocked_tensor                    & A,
             std::string const                   idx_A,
             F                                   beta,
             std::string const                   idx_B
            ) {

      machine_tensor->sum(alpha,
                          *A.machine_tensor,
                          idx_A.c_str(),
                          beta,
                          idx_B.c_str()
                         );

    }

    void read_dense_from_file(MPI_File &file) {
      machine_tensor->read_dense_from_file(file);
    }

    void write(int64_t              const npair,
               std::vector<int64_t> const global_idx,
               std::vector<F>       const data) {
      machine_tensor->write(npair, global_idx.data(), data.data());
    }

    void read_all(F *data) {
      machine_tensor->read_all(data);
    }

    void slice(std::vector<int64_t> const offsets,
               std::vector<int64_t> const ends,
               F                          beta,
               Blocked_tensor       const &A,
               std::vector<int64_t> const offsets_A,
               std::vector<int64_t> const ends_A,
               F                          alpha ) {

      machine_tensor->slice( offsets.data()
                           , ends.data()
                           , beta
                           , *A.machine_tensor
                           , offsets_A.data()
                           , ends_A.data()
                           , alpha);
    }

};
*/
