#ifndef IBST_WORLD_H
#define IBST_WORLD_H

#include "ctf_tensor.h"
#include "tracker.h"
#include <mpi.h>

namespace ibst {

class World {
public:
  World() {
    dryRun = false;
    if (!dryRun) {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &np);
      comm = MPI_COMM_WORLD;
      wrld = new Machine_world();
      tracker = new Tracker(rank, np);
    } else {
      rank = 0;
      np = 1;
    }
  }

  ~World() {
    delete wrld;
    delete tracker;
  }

  int rank;
  int np;
  bool dryRun;
  MPI_Comm comm;
  Machine_world *wrld;
  Tracker *tracker;
};

} // namespace ibst

#endif /* IBST_WORLD_H */
