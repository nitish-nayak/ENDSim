#include <End.hh>
#include <RAT/DB.hh>
#include "MyDAQProc.hh"

namespace END {
End::End(RAT::AnyParse* parser, int argc, char** argv) : Rat(parser, argc, argv) {
  // Append an additional data directory (for ratdb and geo)
  char* enddata = getenv("ENDDATA");
  if (enddata != NULL) {
    ratdb_directories.insert(static_cast<std::string>(enddata) + "/ratdb");
    model_directories.insert(static_cast<std::string>(enddata) + "/models");
  }
  // AppendProcessor constructs a throwaway instance to read each processor's
  // name, and OutNtupleProc's constructor reads the IO table. Configure() does
  // not load the database until Begin(), so load defaults now (idempotent: the
  // later Configure() call overwrites the same table entries).
  char* ratshare = getenv("RATSHARE");
  if (ratshare != NULL) {
    ratdb_directories.insert(static_cast<std::string>(ratshare) + "/ratdb");
    model_directories.insert(static_cast<std::string>(ratshare) + "/models");
  }
  RAT::DB::Get()->LoadDefaults();
  // Initialize a geometry factory
  new GeoEndFactory();
  new RAT::GeoGdmlFactory();
  //new DichroiconArrayFactory();
#if TENSORFLOW_Enabled && NLOPT_Enabled
  RAT::ProcBlockManager::AppendProcessor<HitmanProc>();
#endif
  RAT::ProcBlockManager::AppendProcessor<NtupleProc>();
  RAT::ProcBlockManager::AppendProcessor<END::MyDAQProc>();
  // Include a new type of processor
  // Add a unique component to the datastructure
  // Register generators
  RAT::GlobalFactory<GLG4Gen>::Register("laserball", new RAT::Alloc<GLG4Gen, LaserballGenerator>);
}
}  // namespace END
