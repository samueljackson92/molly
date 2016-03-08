
/* Catch requires CATCH_CONFIG_MAIN to be defined *exactly* once. This file
 includes the define statement and imports the header file. This prevents the
 recompilation catch for every test case.
*/
#ifndef CATCH_CONFIG_MAIN
#define CATCH_CONFIG_MAIN
#endif

#include "catch.hpp"
