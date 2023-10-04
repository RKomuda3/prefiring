#ifndef UserCode_OmtfAnalysis_SynchroCounts_H
#define UserCode_OmtfAnalysis_SynchroCounts_H

#include <vector>
#include <string>

class SynchroCounts {
public:
  SynchroCounts() : theCounts(std::vector<unsigned int>(8,0)) {}
  SynchroCounts(const unsigned int *hits) : theCounts(std::vector<unsigned int>(hits,hits+8)) {}

  void increment(unsigned int bxDiff);
  void set(unsigned int bxDiff);
  unsigned int firstHit() const;

  double rms() const;
  double mean() const;
  unsigned int sum() const { return mom0(); }

  std::string print() const;
  const std::vector<unsigned int> & counts() const { return theCounts; }
  bool operator==(const SynchroCounts &) const;
  SynchroCounts & operator+=(const SynchroCounts &rhs);
private:
  unsigned int  mom0() const;
  double mom1() const;
  std::vector<unsigned int> theCounts;
};
#endif

