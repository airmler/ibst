#ifndef IBST_UTIL_H
#define IBST_UTIL_H

#include <algorithm>
#include <cassert>
#include <chrono>
#include <functional>
#include <mpi.h>
#include <string>
#include <vector>
#include <complex>
#include <iostream>

namespace ibst {

// TODO MPI command in a library is maybe a problem
//#define ASSERT(...)                \
//  do { if (!(__VA_ARGS__)){ int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank); if
//  (rank == 0){ printf("CTF ERROR: %s:%d, ASSERT(%s)
//  failed\n",__FILE__,__LINE__,#__VA_ARGS__); } CTF_int::handler();
//  assert(__VA_ARGS__); } } while (0)

#define ILOG(name)                                                             \
  do {                                                                         \
    int rank;                                                                  \
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);                                      \
    if (!rank)                                                                 \
      std::cout << name << std::endl;                                          \
  } while (0)

#define _FORMAT(_fmt, ...)                                                     \
  ([&](void) -> std::string {                                                  \
    int _sz = std::snprintf(nullptr, 0, _fmt, __VA_ARGS__);                    \
    std::vector<char> _out(_sz + 1);                                           \
    std::snprintf(&_out[0], _out.size(), _fmt, __VA_ARGS__);                   \
    return std::string(_out.data());                                           \
  })()

inline void checkDuplicate(std::string t) {
  std::sort(t.begin(), t.end());
  auto it = std::unique(t.begin(), t.end());
  if (it != t.end())
    assert(0);
}

// p provides a the order for sorting the input array
// p[0] indicates the slowest index....p.back() is the fastest index
std::function<bool(const std::vector<int64_t> &,
                   const std::vector<int64_t>
                       &)> inline compare(const std::vector<int64_t> p) {
  return [p](const std::vector<int64_t> &a,
             const std::vector<int64_t> &b) -> bool {
    size_t n(p.size());
    std::vector<int64_t> c(n), d(n);
    for (size_t i(0); i < n; i++) {
      c[i] = a[p[i]];
      d[i] = b[p[i]];
    }
    return c < d;
  };
}

std::function<bool(const std::vector<int64_t> &)> inline find(
    const std::vector<int64_t> c,
    const std::vector<std::pair<int64_t, int64_t>> ca) {
  return [c, ca](const std::vector<int64_t> &a) -> bool {
    for (size_t i(0); i < ca.size(); i++)
      for (auto &x : ca)
        if (a[x.second] != c[x.first])
          return false;
    return true;
  };
}

inline std::string hash_me(std::string s) {
  std::string o;
  // transform string to a size_t hash and reduce it to int-size
  int h(std::hash<std::string>{}(s) % 2147483648);
  // a integer number can be written as a 6 digit alphanumeric object
  std::string an = "abcdefghijklmnopqrstuvwxyz0123456789";
  for (int i(0); i < 6; i++) {
    o += an[h % 36];
    h /= 36;
  }
  return o;
}

inline std::string tensor_to_line(std::string _name, std::vector<int64_t> _lens, int64_t _nzc_dim) {
  std::string line = _name + "[";
  for (auto l: _lens) line += std::to_string(l) + ",";
  if (_lens.size()) line.pop_back();
  line += "](";
  line += std::to_string(_nzc_dim);
  line += ")";

  return line;
}

namespace traits {

template <typename F>
inline bool is_complex() {
  return false;
}

template <>
inline bool is_complex<double>() {
  return false;
}

template <>
inline bool is_complex< std::complex<double> >() {
  return true;
}

}



} // namespace ibst
#endif /* IBST_UTIL_H */
