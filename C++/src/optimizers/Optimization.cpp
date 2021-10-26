#include "Optimization.h"

std::string
getNLOPTResultDescription(nlopt::result resultcode)
{
  switch (resultcode)
  {
    case nlopt::result::FAILURE:
      return "Failed - Generic result code";
      break;

    case nlopt::result::INVALID_ARGS:
      return "Failed - Invalid arguments";
      break;

    case nlopt::result::OUT_OF_MEMORY:
      return "Failed - Out of memory";
      break;

    case nlopt::result::ROUNDOFF_LIMITED:
      return "Failed - Round-off limited";
      break;

    case nlopt::result::FORCED_STOP:
      return "Failed - Forcefully stopped";
      break;

    case nlopt::result::SUCCESS:
      return "Success - Generic result code";
      break;

    case nlopt::result::STOPVAL_REACHED:
      return "Success - Stop value reached";
      break;

    case nlopt::result::FTOL_REACHED:
      return "Success - Relative f-tolerance reached";
      break;

    case nlopt::result::XTOL_REACHED:
      return "Success - Relative x-tolerance reached";
      break;

    case nlopt::result::MAXEVAL_REACHED:
      return "Success - Maximum evaluation count reached";
      break;

    case nlopt::result::MAXTIME_REACHED:
      return "Success - Maximum time reached";
      break;

    default:
      break;
  }

  return "Unknown result code";
}
