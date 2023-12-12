#ifndef IBST_TRACKER_H
#define IBST_TRACKER_H

#include <fstream>
#include <string>

namespace ibst {

class Event {
public:
  /* possible Operations:
   * Initialize, Destroy, Contract, Sum, Slice,
   * Read, Write, ReadFromFile, WriteToFile
   */
  std::string _operation;

  // In case of Initialize/Destroy: add/substract memory footprint (# words)
  int64_t _memoryUsage = 0;

  // In case of all other operations: provide number of operations
  double _flopCount = 0.0;

  // In case of dryRun/calculation: provide estimate/elapsed time
  double _time = 0.0;

  // Provide operation information
  std::string _line = "";

  // A 6 digit alphanumerical hash of the operation information
  std::string _hash = "";

  Event &operation(std::string v) {
    _operation = v;
    return *this;
  }
  Event &memoryUsage(int64_t v) {
    _memoryUsage = v;
    return *this;
  }
  Event &flopCount(double v) {
    _flopCount = v;
    return *this;
  }
  Event &time(double v) {
    _time = v;
    return *this;
  }

  Event &line(std::string v) {
    _line = v;
    _hash = hash_me(v);
    return *this;
  }

  std::string print(int np, int64_t hwm) {
    assert(np);
    std::string out;
    out += _FORMAT("%-10s: ", _operation.c_str());
    if (_line.size()) {
      out += _hash + " ";
      out += _FORMAT("%-48s: ", _line.c_str());
    }
    if (_memoryUsage != 0) {
      out += _FORMAT("mem %6.4f (GB/rank), ", 9.3132e-10 *  _memoryUsage / np );
      out += _FORMAT("hwm %6.4f (GB/rank)" , 9.3132e-10 * hwm / np );
    }
    if (_flopCount > 0.0) {
      out += _FORMAT("%4.2e GF/rank, ", _flopCount / np / 1e9);
    }
    if (_time > 0.0) {
      out += _FORMAT("%7.4lf s, ", _time);
      out += _FORMAT("%4.2e GF/rank/s", _flopCount / np / 1e9 / _time);
    }
    return out;
  }

  /*

    Event event
          .operations("shit in shit")
          .memoryUsage(546.654);

  */
};

class Tracker {
public:
  std::vector<Event> events;
  int rank;
  int np;
  size_t bufSize;
  size_t bufCount;
  size_t hwm;
  std::string fileName;
  std::ofstream logFile;

  void add(Event &_event) {
    if (rank)
      return;
    if (_event._memoryUsage != 0)
      hwm += _event._memoryUsage;
    events.push_back(_event);
    bufCount++;
    if (bufCount == bufSize) {
      bufCount = 0;
      // flush
      logFile << _event.print(np, hwm) << std::endl;
    }
  }
  Tracker(int _rank, int _np)
      : rank(_rank)
      , np(_np) {

    fileName = "ibst.log";
    if (!rank) {
      logFile = std::ofstream(fileName.c_str(),
                              std::ofstream::out | std::ofstream::trunc);
    }
    bufSize = 1;
    bufCount = 0;
    hwm = 0;
    if (!rank) {
      logFile << "Operation :  hash ";
      logFile << " __  __  __   __  __  __  __  __  __  __  __  __ :";
      logFile << " Δmem | hwm .or. ops | time | perf ";
      logFile << std::endl;
      logFile << "-----\n";
    }
  }
};

} // namespace ibst

#endif /* IBST_TRACKER_H */
