cmake_minimum_required(VERSION 2.8...3.23)
project(qpgen)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread -std=c++14")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")

find_path(HTS_INCLUDE_DIRS htslib/bgzf.h ../htslib/bgzf.h HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ submodules/htslib ../htslib)
if(NOT HTS_INCLUDE_DIRS )
    message(FATAL_ERROR "libhts HTS_INCLUDE_DIRS not found")
endif()

message(STATUS ${HTS_INCLUDE_DIRS})

find_library(HTS_LIBRARIES libhts.a HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ submodules/htslib ../htslib/)
if(NOT HTS_LIBRARIES)
    message(FATAL_ERROR "libhts HTS_LIBRARIES not found")
endif()

message(STATUS ${HTS_LIBRARIES})

find_path(QGEN_INCLUDE_DIRS qgenlib/qgen_error.h  HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ submodules/qgenlib ../qgenlib)
if(NOT QGEN_INCLUDE_DIRS )
    message(FATAL_ERROR "libqgen QGEN_INCLUDE_DIRS not found")
endif()

message(STATUS ${QGEN_INCLUDE_DIRS})

find_library(QGEN_LIBRARIES libqgen.a  HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ submodules/qgenlib/lib/ ../qgenlib/lib/)
if(NOT QGEN_LIBRARIES)
    message(FATAL_ERROR "libqgen QGEN_LIBRARIES not found")
endif()

message(STATUS ${QGEN_LIBRARIES})

# Eigen
set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/submodules/eigen")
if (NOT EXISTS "${EIGEN3_INCLUDE_DIR}/Eigen/Dense")
    message(FATAL_ERROR "Eigen submodule missing. Did you run 'git submodule update --init'?")
endif()
add_library(Eigen3::Eigen INTERFACE IMPORTED)
set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")


include_directories(. simde/ ${HTS_INCLUDE_DIRS} ${QGEN_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})

if(APPLE)
    message(STATUS "Configuring for macOS")
    set(CMAKE_C_COMPILER "/opt/homebrew/opt/llvm/bin/clang")
    set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
endif()


set(APP_FILES
    qpgentools.cpp
    qpgentools.h
    cmd_test_qpgen.cpp
    cmd_match_plp_geno.cpp
    cmd_join_plp_geno.cpp
    cmd_geno2tsv.cpp
    cmd_pair_assoc.cpp
)

set(SOURCE_FILES
    qpgen.cpp
    qpgen.h
    assoc_utils.h
    assoc_utils.cpp
    include/pgenlib_misc.cc
    include/pgenlib_misc.h
    include/pgenlib_read.cc
    include/pgenlib_read.h
    include/plink2_base.cc
    include/plink2_base.h
    include/plink2_bits.cc
    include/plink2_bits.h
    pgenlib_ffi_support.cpp
    pgenlib_ffi_support.h
    pgenlibr.cpp
    pgenlibr.h
    pvar_ffi_support.cc
    pvar_ffi_support.h
    simde/arm/neon.h
    simde/arm/neon/aba.h
    simde/arm/neon/abd.h
    simde/arm/neon/abdl.h
    simde/arm/neon/abs.h
    simde/arm/neon/add.h
    simde/arm/neon/addhn.h
    simde/arm/neon/addl.h
    simde/arm/neon/addl_high.h
    simde/arm/neon/addlv.h
    simde/arm/neon/addv.h
    simde/arm/neon/addw.h
    simde/arm/neon/addw_high.h
    simde/arm/neon/and.h
    simde/arm/neon/bcax.h
    simde/arm/neon/bic.h
    simde/arm/neon/bsl.h
    simde/arm/neon/cage.h
    simde/arm/neon/cagt.h
    simde/arm/neon/ceq.h
    simde/arm/neon/ceqz.h
    simde/arm/neon/cge.h
    simde/arm/neon/cgez.h
    simde/arm/neon/cgt.h
    simde/arm/neon/cgtz.h
    simde/arm/neon/cle.h
    simde/arm/neon/clez.h
    simde/arm/neon/cls.h
    simde/arm/neon/clt.h
    simde/arm/neon/cltz.h
    simde/arm/neon/clz.h
    simde/arm/neon/cmla.h
    simde/arm/neon/cmla_rot180.h
    simde/arm/neon/cmla_rot270.h
    simde/arm/neon/cmla_rot90.h
    simde/arm/neon/cnt.h
    simde/arm/neon/combine.h
    simde/arm/neon/create.h
    simde/arm/neon/cvt.h
    simde/arm/neon/dot.h
    simde/arm/neon/dot_lane.h
    simde/arm/neon/dup_lane.h
    simde/arm/neon/dup_n.h
    simde/arm/neon/eor.h
    simde/arm/neon/ext.h
    simde/arm/neon/fma.h
    simde/arm/neon/fma_lane.h
    simde/arm/neon/fma_n.h
    simde/arm/neon/get_high.h
    simde/arm/neon/get_lane.h
    simde/arm/neon/get_low.h
    simde/arm/neon/hadd.h
    simde/arm/neon/hsub.h
    simde/arm/neon/ld1.h
    simde/arm/neon/ld1_dup.h
    simde/arm/neon/ld1_lane.h
    simde/arm/neon/ld2.h
    simde/arm/neon/ld3.h
    simde/arm/neon/ld4.h
    simde/arm/neon/ld4_lane.h
    simde/arm/neon/max.h
    simde/arm/neon/maxnm.h
    simde/arm/neon/maxv.h
    simde/arm/neon/min.h
    simde/arm/neon/minnm.h
    simde/arm/neon/minv.h
    simde/arm/neon/mla.h
    simde/arm/neon/mla_n.h
    simde/arm/neon/mlal.h
    simde/arm/neon/mlal_high.h
    simde/arm/neon/mlal_high_n.h
    simde/arm/neon/mlal_lane.h
    simde/arm/neon/mlal_n.h
    simde/arm/neon/mls.h
    simde/arm/neon/mls_n.h
    simde/arm/neon/mlsl.h
    simde/arm/neon/mlsl_high.h
    simde/arm/neon/mlsl_high_n.h
    simde/arm/neon/mlsl_lane.h
    simde/arm/neon/mlsl_n.h
    simde/arm/neon/movl.h
    simde/arm/neon/movl_high.h
    simde/arm/neon/movn.h
    simde/arm/neon/movn_high.h
    simde/arm/neon/mul.h
    simde/arm/neon/mul_lane.h
    simde/arm/neon/mul_n.h
    simde/arm/neon/mull.h
    simde/arm/neon/mull_high.h
    simde/arm/neon/mull_lane.h
    simde/arm/neon/mull_n.h
    simde/arm/neon/mvn.h
    simde/arm/neon/neg.h
    simde/arm/neon/orn.h
    simde/arm/neon/orr.h
    simde/arm/neon/padal.h
    simde/arm/neon/padd.h
    simde/arm/neon/paddl.h
    simde/arm/neon/pmax.h
    simde/arm/neon/pmin.h
    simde/arm/neon/qabs.h
    simde/arm/neon/qadd.h
    simde/arm/neon/qdmulh.h
    simde/arm/neon/qdmulh_lane.h
    simde/arm/neon/qdmulh_n.h
    simde/arm/neon/qdmull.h
    simde/arm/neon/qmovn.h
    simde/arm/neon/qmovn_high.h
    simde/arm/neon/qmovun.h
    simde/arm/neon/qneg.h
    simde/arm/neon/qrdmulh.h
    simde/arm/neon/qrdmulh_lane.h
    simde/arm/neon/qrdmulh_n.h
    simde/arm/neon/qrshrn_n.h
    simde/arm/neon/qrshrun_n.h
    simde/arm/neon/qshl.h
    simde/arm/neon/qshlu_n.h
    simde/arm/neon/qshrn_n.h
    simde/arm/neon/qshrun_n.h
    simde/arm/neon/qsub.h
    simde/arm/neon/qtbl.h
    simde/arm/neon/qtbx.h
    simde/arm/neon/rbit.h
    simde/arm/neon/recpe.h
    simde/arm/neon/recps.h
    simde/arm/neon/reinterpret.h
    simde/arm/neon/rev16.h
    simde/arm/neon/rev32.h
    simde/arm/neon/rev64.h
    simde/arm/neon/rhadd.h
    simde/arm/neon/rnd.h
    simde/arm/neon/rndi.h
    simde/arm/neon/rndm.h
    simde/arm/neon/rndn.h
    simde/arm/neon/rndp.h
    simde/arm/neon/rshl.h
    simde/arm/neon/rshr_n.h
    simde/arm/neon/rshrn_n.h
    simde/arm/neon/rsqrte.h
    simde/arm/neon/rsqrts.h
    simde/arm/neon/rsra_n.h
    simde/arm/neon/set_lane.h
    simde/arm/neon/shl.h
    simde/arm/neon/shl_n.h
    simde/arm/neon/shll_n.h
    simde/arm/neon/shr_n.h
    simde/arm/neon/shrn_n.h
    simde/arm/neon/sqadd.h
    simde/arm/neon/sra_n.h
    simde/arm/neon/sri_n.h
    simde/arm/neon/st1.h
    simde/arm/neon/st1_lane.h
    simde/arm/neon/st2.h
    simde/arm/neon/st2_lane.h
    simde/arm/neon/st3.h
    simde/arm/neon/st3_lane.h
    simde/arm/neon/st4.h
    simde/arm/neon/st4_lane.h
    simde/arm/neon/sub.h
    simde/arm/neon/subhn.h
    simde/arm/neon/subl.h
    simde/arm/neon/subl_high.h
    simde/arm/neon/subw.h
    simde/arm/neon/subw_high.h
    simde/arm/neon/tbl.h
    simde/arm/neon/tbx.h
    simde/arm/neon/trn.h
    simde/arm/neon/trn1.h
    simde/arm/neon/trn2.h
    simde/arm/neon/tst.h
    simde/arm/neon/types.h
    simde/arm/neon/uqadd.h
    simde/arm/neon/uzp.h
    simde/arm/neon/uzp1.h
    simde/arm/neon/uzp2.h
    simde/arm/neon/xar.h
    simde/arm/neon/zip.h
    simde/arm/neon/zip1.h
    simde/arm/neon/zip2.h
    simde/arm/sve.h
    simde/arm/sve/add.h
    simde/arm/sve/and.h
    simde/arm/sve/cmplt.h
    simde/arm/sve/cnt.h
    simde/arm/sve/dup.h
    simde/arm/sve/ld1.h
    simde/arm/sve/ptest.h
    simde/arm/sve/ptrue.h
    simde/arm/sve/qadd.h
    simde/arm/sve/reinterpret.h
    simde/arm/sve/sel.h
    simde/arm/sve/st1.h
    simde/arm/sve/sub.h
    simde/arm/sve/types.h
    simde/arm/sve/whilelt.h
    simde/check.h
    simde/debug-trap.h
    simde/hedley.h
    simde/simde-arch.h
    simde/simde-common.h
    simde/simde-constify.h
    simde/simde-detect-clang.h
    simde/simde-diagnostic.h
    simde/simde-features.h
    simde/simde-math.h
    simde/x86/avx.h
    simde/x86/avx2.h
    simde/x86/avx512.h
    simde/x86/avx512/abs.h
    simde/x86/avx512/add.h
    simde/x86/avx512/adds.h
    simde/x86/avx512/and.h
    simde/x86/avx512/andnot.h
    simde/x86/avx512/avg.h
    simde/x86/avx512/blend.h
    simde/x86/avx512/broadcast.h
    simde/x86/avx512/cast.h
    simde/x86/avx512/cmp.h
    simde/x86/avx512/cmpeq.h
    simde/x86/avx512/cmpge.h
    simde/x86/avx512/cmpgt.h
    simde/x86/avx512/cmple.h
    simde/x86/avx512/cmplt.h
    simde/x86/avx512/copysign.h
    simde/x86/avx512/cvt.h
    simde/x86/avx512/cvts.h
    simde/x86/avx512/div.h
    simde/x86/avx512/extract.h
    simde/x86/avx512/fmadd.h
    simde/x86/avx512/fmsub.h
    simde/x86/avx512/fnmadd.h
    simde/x86/avx512/fnmsub.h
    simde/x86/avx512/insert.h
    simde/x86/avx512/load.h
    simde/x86/avx512/loadu.h
    simde/x86/avx512/lzcnt.h
    simde/x86/avx512/max.h
    simde/x86/avx512/min.h
    simde/x86/avx512/mov.h
    simde/x86/avx512/mov_mask.h
    simde/x86/avx512/movm.h
    simde/x86/avx512/mul.h
    simde/x86/avx512/mulhi.h
    simde/x86/avx512/mulhrs.h
    simde/x86/avx512/mullo.h
    simde/x86/avx512/negate.h
    simde/x86/avx512/or.h
    simde/x86/avx512/packs.h
    simde/x86/avx512/packus.h
    simde/x86/avx512/permutex2var.h
    simde/x86/avx512/permutexvar.h
    simde/x86/avx512/sad.h
    simde/x86/avx512/set.h
    simde/x86/avx512/set1.h
    simde/x86/avx512/set4.h
    simde/x86/avx512/setone.h
    simde/x86/avx512/setr.h
    simde/x86/avx512/setr4.h
    simde/x86/avx512/setzero.h
    simde/x86/avx512/shuffle.h
    simde/x86/avx512/sll.h
    simde/x86/avx512/slli.h
    simde/x86/avx512/sllv.h
    simde/x86/avx512/sqrt.h
    simde/x86/avx512/sra.h
    simde/x86/avx512/srai.h
    simde/x86/avx512/srav.h
    simde/x86/avx512/srl.h
    simde/x86/avx512/srli.h
    simde/x86/avx512/srlv.h
    simde/x86/avx512/store.h
    simde/x86/avx512/storeu.h
    simde/x86/avx512/sub.h
    simde/x86/avx512/subs.h
    simde/x86/avx512/test.h
    simde/x86/avx512/types.h
    simde/x86/avx512/unpackhi.h
    simde/x86/avx512/unpacklo.h
    simde/x86/avx512/xor.h
    simde/x86/avx512/xorsign.h
    simde/x86/fma.h
    simde/x86/gfni.h
    simde/x86/mmx.h
    simde/x86/sse.h
    simde/x86/sse2.h
    simde/x86/sse3.h
    simde/x86/sse4.1.h
    simde/x86/sse4.2.h
    simde/x86/ssse3.h
    simde/x86/svml.h
)

ADD_LIBRARY(qpgen STATIC ${SOURCE_FILES})

SET(APP_EXE qpgentools)

#find_package(OpenMP REQUIRED)
ADD_EXECUTABLE(${APP_EXE} ${SOURCE_FILES} ${APP_FILES})

find_library(ZLIB z HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT ZLIB)
    message(FATAL_ERROR "libz library not found")
endif()

find_library(CURLLIB curl HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CURLLIB)
    message(FATAL_ERROR "libcurl library not found")
endif()

find_package (BZip2)
if (NOT BZIP2_FOUND)
    message(FATAL_ERROR "libbz2 library not found")
else()
    include_directories(${BZIP2_INCLUDE_DIRS})
    set(BZIP2 ${BZIP2_LIBRARIES})
endif()

find_library(LZMA lzma HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT LZMA)
    message(WARNING "liblzma library not found, if you specified --disable-lzma when compiling libhts, please ignore this warning")
    set(LZMA "")
else()
    set(LZMA ${LZMA_LIBRARIES})
endif()

find_library(DEFLATELIB deflate HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT DEFLATELIB)
    message(FATAL_ERROR "libdeflate library not found")
endif()

find_library(CRYPTOLIB crypto HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CRYPTOLIB)
    message(FATAL_ERROR "libcrypto library not found")
endif()

#TARGET_LINK_LIBRARIES(${APP_EXE} PRIVATE ${QGEN_LIBRARIES} ${HTS_LIBRARIES} ${ZLIB} ${LZMA} ${BZIP2} ${CURLLIB} ${DEFLATELIB} ${CRYPTOLIB} "${OpenMP_CXX_FLAGS}")
TARGET_LINK_LIBRARIES(${APP_EXE} PRIVATE Eigen3::Eigen ${QGEN_LIBRARIES} ${HTS_LIBRARIES} ${ZLIB} ${LZMA} ${BZIP2} ${CURLLIB} ${DEFLATELIB} ${CRYPTOLIB} "${OpenMP_CXX_FLAGS}")
#target_compile_options(${APP_EXE} PRIVATE "${OpenMP_CXX_FLAGS}")


install(TARGETS ${APP_EXE} RUNTIME DESTINATION bin)
