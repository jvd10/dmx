cmake_minimum_required (VERSION 2.6)
project (sandbox_jvd_apps_dmx)

SET(CMAKE_BUILD_TYPE Release CACHE STRING "Release" FORCE)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(dmx dmx.cpp dmxCore.cpp dmxIO.cpp dmxRead.cpp dmxBarcode.cpp)


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")


include_directories(/usr/include /home/ghedin/common/sl/bld/tbb/tbb40_297oss/include /home/ghedin/common/sl/include/ltilib)
link_directories(/usr/lib64/ /home/ghedin/common/sl/bld/tbb/tbb40_297oss/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21 /home/ghedin/common/sl/lib/ltilib)
target_link_libraries(dmx z boost_iostreams /home/ghedin/common/sl/bld/tbb/tbb40_297oss/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21/libtbb.so /home/ghedin/common/sl/lib/ltilib/libltid.a /home/ghedin/common/sl/lib/ltilib/libltinvd.a /home/ghedin/common/sl/lib/ltilib/libltinvr.a /home/ghedin/common/sl/lib/ltilib/libltir.a)
