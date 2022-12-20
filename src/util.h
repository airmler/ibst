#ifndef IBST_UTIL_H
#define IBST_UTIL_H

#include <mpi.h>
#include <cassert>
#include <string>
#include <functional>
#include <vector>
#include <algorithm>

namespace ibst {

//TODO MPI command in a library is maybe a problem
//#define ASSERT(...)                \
//  do { if (!(__VA_ARGS__)){ int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank); if (rank == 0){ printf("CTF ERROR: %s:%d, ASSERT(%s) failed\n",__FILE__,__LINE__,#__VA_ARGS__); } CTF_int::handler(); assert(__VA_ARGS__); } } while (0)

#define ILOG(name)                                \
  do { int rank;                                  \
       MPI_Comm_rank(MPI_COMM_WORLD, &rank);      \
       if (!rank) std::cout << name << std::endl; \
     } while(0)


#define _FORMAT(_fmt, ...)                                    \
  ([&] (void) -> std::string {                                \
     int _sz = std::snprintf(nullptr, 0, _fmt, __VA_ARGS__);  \
     std::vector<char>  _out(_sz  +  1);                      \
     std::snprintf(&_out[0], _out.size(), _fmt, __VA_ARGS__); \
     return std::string(_out.data());                         \
   })()



inline void checkDuplicate(std::string t) {
  std::sort(t.begin(), t.end());
  auto it = std::unique(t.begin(), t.end());
  if (it != t.end()) assert(0);
}

// p provides a the order for sorting the input array
// p[0] indicates the slowest index....p.back() is the fastest index
std::function<bool(const std::vector<size_t> &, const std::vector<size_t> &)>
inline compare(const std::vector<size_t> p)
{
  return [p] (const std::vector<size_t> &a, const std::vector<size_t> &b) -> bool {
    size_t n(p.size());
    std::vector<size_t> c(n), d(n);
    for (size_t i(0); i < n; i++) {
      c[i] = a[p[i]]; d[i] = b[p[i]];
    }
    return c < d;
  };
}

std::function<bool(const std::vector<size_t> &)>
inline find(const std::vector<size_t> c, const std::vector< std::pair<size_t,size_t> > ca)
{
  return [c, ca] (const std::vector<size_t> &a) -> bool {
    for (size_t i(0); i < ca.size(); i++)
    for (auto &x: ca)  if (a[x.second] != c[x.first]) return false;
    return true;
  };
}

}  //namespace
#endif /* IBST_UTIL_H */
