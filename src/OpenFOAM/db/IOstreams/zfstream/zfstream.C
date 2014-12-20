#if defined(__GNUC__) && __GNUC__ < 3
#   include "zfstream2.C"
#else
#   include "std_zfstream.C"
#endif
