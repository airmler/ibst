#include "../src/tensor.h"
#include "../src/kmesh.h"

using ibst::Tensor;

double inv(double x){ return 1.0/x; }

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  auto kmesh(ibst::Kmesh({3,2,1}));

  auto dw = new CTF::World();

  Tensor<double> G(3, {300, 100, 10}, kmesh.getDefaultNzc(3), dw, "G");
  Tensor<double> H(3, {300, 100, 10}, kmesh.getDefaultNzc(3), dw, "H");
  Tensor<double> C(3, {300, 100, 10}, kmesh.getDefaultNzc(3), dw, "C");
  Tensor<double> K(3, {300, 10, 100}, kmesh.getDefaultNzc(3), dw, "K");
  Tensor<double> T(4, {100, 100, 10, 10}, kmesh.getDefaultNzc(4), dw, "T");
  Tensor<double> D(4, {100, 100, 10, 10}, kmesh.getDefaultNzc(4), dw, "D");
  Tensor<double> E(1, {110}, kmesh.getDefaultNzc(1), dw, "E");
  Tensor<double> O(1, {10}, kmesh.getDefaultNzc(1), dw, "O");
  Tensor<double> U(1, {100}, kmesh.getDefaultNzc(1), dw, "U");

  for (auto &g: G.tensors) g->fill_random(0., 1.);

  if (!dw->rank) std::ofstream("testG");
  MPI_File file;
  MPI_File_open(dw->comm, "testG", MPI_MODE_RDWR, MPI_INFO_NULL, &file);
  G.write_dense_to_file(file);
  H.read_dense_from_file(file);
  MPI_File_close(&file);
  if (!dw->rank) {
    int ioerr = remove("testG");
    assert(!ioerr);
  }

  std::vector<double> v(100);
  for (auto &vv: v) vv = 3.41;
  std::vector<size_t> idx(100);
  std::iota(idx.begin(), idx.end(), 0);
  if (dw->rank) idx.resize(0);

  G.write(idx.size(), idx, v);
  G.read(idx.size(), idx, v);

  G.write(idx.size(), idx, v, 4);
  G.read(idx.size(), idx, v, 4);

  O.slice({0}, {10}, 0.0, E, {0}, {10}, 1.0);
  U.slice({0}, {100}, 0.0, E, {10}, {110}, 1.0);


  std::function<double(const double)> fInv(&inv);

  std::vector<size_t> remap(O.nzc.size());
  for (size_t i(0); i < remap.size(); i++) remap[i] = (i*17)%remap.size();

  D.sum(1.0, O, "i", 0.0, "abij");
  D.sum(1.0, O, "j", 1.0, "abij");
  D.sum(1.0, U, "a", 1.0, "abij");
  D.sum(1.0, U, "b", 1.0, "abij");
  D.sum(1.0, D, "abij", 0.0, "abij", fInv);
  T.sum(1.0, T, "baji", 1.0, "abij");

  O.sum(1.0, O, "a", 0.0, "a", remap);
  remap.resize(G.nzc.size());
  auto nzc(kmesh.getDefaultNzc(3));
  std::sort(nzc.begin(), nzc.end(), ibst::compare({1,2,0}));
  for (size_t i(0); i < nzc.size(); i++) remap[i] = nzc[i].back();
  C.sum(1.0, K, "Gia", 0.0, "Gai", remap, fInv);

  nzc = C.nzc;
  for (auto &n: nzc) n[0] = kmesh.getMinusIdx(n[0]);
  C.relabelBlocks(nzc);
  T.contract(1.0, C, "Gai", G, "Gbj", 0.0, "abij");
  T.contract(1.0, T, "acik", T, "cbkj", 0.0, "abij");
  MPI_Finalize();
  return 0;
}
