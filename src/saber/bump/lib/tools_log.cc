/*
 * (C) Copyright 2022 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <cstdint>
#include <cstring>

#include "eckit/log/Channel.h"
#include "eckit/log/Log.h"
#include "eckit/runtime/Main.h"

using int32 = std::int32_t;

namespace bump_lib {

// -----------------------------------------------------------------------------

extern "C" {
  void log__write_log(eckit::Channel * channel, char* msg, int32 newl, int32 flush) {
    if (eckit::Main::ready()) {
      if (::strlen( msg ) )
        *channel << msg;
      else
        *channel << " ";
      if ( newl )
        *channel << eckit::newl;
      if ( flush )
        *channel << std::flush;
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace bump_lib
