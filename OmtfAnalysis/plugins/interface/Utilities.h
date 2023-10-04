#ifndef OmtfAnalysis_Utilities_H
#define OmtfAnalysis_Utilities_H

#include <vector>
#include <map>
#include <iomanip>
#include <iostream>

#include "L1Trigger/L1TMuonOverlap/interface/OmtfName.h"
class L1Obj;

namespace omtfUtilities {
  OmtfName makeName(const L1Obj & obj);

  const int nEtaBins = 14;
  extern const int  etaBinVal[nEtaBins]; // = {73, 75, 78, 79, 85, 90, 92, 94, 95, 99, 103, 110, 115, 121};

  double code2pt(int pt);

  template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); } 

  int code2HistoBin (int code);

  extern std::vector<std::string> layerNames;


} // end omtfUtilities


// variable-size PT bins corresponding to L1 trigger scale. 
// Lower bound of bin ipt corresponds to L1Pt(ipt),
// with exception of ipt=0 and 1 (which should be NAN and 0. instead of 0. and 0.1 respectively)
class L1PtScale {
public:
  static unsigned int ptCode(float ptValue);
  static float ptValue(unsigned int ptCode);
  static const unsigned int nPtBins = 38;
  static double ptBins[nPtBins + 1];
};

// variable-size ETA bins corresponding to PAC trigger granularity
class L1RpcEtaScale {
public:
  static int etaCode(float etaValue); // tower for given eta, +-17 for over/undeflows
  static float etaValue(int etaCode); // central eta value of a tower, +-2.4 for over/underrflows 
  static  const unsigned int nEtaBins=33;
  static double etaBins[nEtaBins+1];
};

class RunEffMap {
public:
  struct EffAndErr : public std::pair<double,double> {
    EffAndErr(double aEff, double aErr) :  std::pair<double,double>(std::make_pair(aEff,aErr)) {}
    double    eff() const { return first; }
    double effErr() const { return second; }
    friend std::ostream & operator<< (std::ostream &out, const EffAndErr &o) {
       return (out << std::setprecision(4) <<o.eff() <<" +/- "<< o.effErr() );
    }
  };
  std::vector<unsigned int> runs() const;
  bool hasRun(unsigned int run) const;
  void addEvent(unsigned int run, bool fired);
  EffAndErr effAndErr(unsigned int run) const;
  EffAndErr effAndErr() const;

  std::pair<unsigned int, unsigned int> stat(unsigned int run) const;

private:
   typedef std::map< unsigned int, std::pair<unsigned int, unsigned int> > Map;
   EffAndErr effAndErr(unsigned int nFired, unsigned int nEvents) const;
   Map theRunMap;
};

#endif
