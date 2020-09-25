#include "gen.h"


thread_local Gen::Ctxt Gen::ctxt;


using Bool = Gen::Bool;


Bool Bool::operator&(const Bool& o) const {
  auto& ctxt = Gen::ctxt;
};
