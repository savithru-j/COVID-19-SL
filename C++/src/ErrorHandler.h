#ifndef ERRORHANDLER_H
#define ERRORHANDLER_H

#include <string>

void throwError(const std::string& msg)
{
  std::cout << std::endl << "Error: " << msg << std::endl;
  exit(EXIT_FAILURE);
}

#endif //ERRORHANDLER_H
