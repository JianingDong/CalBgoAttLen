#include <boost/python.hpp>
#include "DmpAlgBgoAttLen.h"

BOOST_PYTHON_MODULE(libDmpBgoAttLen){
  using namespace boost::python;

  class_<DmpAlgBgoAttLen,boost::noncopyable,bases<DmpVAlg> >("DmpAlgBgoAttLen",init<>());
}
