#############################################################################
# Initialization
AC_PREREQ([2.69])
AC_INIT([graph_libs_eval], [0.1])
AC_CONFIG_AUX_DIR([build-aux])
AC_PRESERVE_HELP_ORDER
AC_LANG([C++])

#############################################################################
# Main makefile
AC_CONFIG_FILES([Makefile])

#############################################################################
# Set the compiler $CC and $CXX. Prefer clang over gcc
old_CFLAGS="${CFLAGS}"; old_CXXFLAGS="${CXXFLAGS}" dnl Do not automatically set -g -O2
AC_PROG_CC([clang cc gcc icc])
AC_PROG_CXX([clang++ c++ cxx g++ icpc])
CFLAGS="${old_CFLAGS}"; unset old_CFLAGS; CXXFLAGS="${old_CXXFLAGS}"; unset old_CXXFLAGS;
AX_CXX_COMPILE_STDCXX_17 dnl Ask for C++17, thanks!

#############################################################################
# Check whether the user has explicitly set CPPFLAGS, CFLAGS and CXXFLAGS. If 
# so we try to avoid polluting these flags and respect the user setting
m4_divert_push([INIT_PREPARE]) 
if test "x${CPPFLAGS}" != "x"; then ac_user_cppflags="yes"; fi
if test "x${CFLAGS}" != "x"; then ac_user_cflags="yes"; fi
if test "x${CXXFLAGS}" != "x"; then ac_user_cxxflags="yes"; fi
m4_divert_pop([INIT_PREPARE])

#############################################################################
# cmake, dependency for libcommon
AC_CHECK_PROGS(cmake, [cmake], [no])
if( test x"${cmake}" = x"no" ); then
   AC_MSG_ERROR([missing prerequisite: this program requires cmake (dependency for libcommon)])
fi

#############################################################################
# git, dependency for libcommon
AC_CHECK_PROGS(git, [git], [no])
if( test x"${git}" = x"no" ); then
   AC_MSG_ERROR([missing prerequisite: this program requires git (dependency for libcommon)])
fi

#############################################################################
# OpenMP
have_openmp="yes";
m4_foreach_w([lang], [C, C++], [AC_LANG_PUSH(lang) AX_OPENMP([], [have_openmp="no"]) AC_LANG_POP(lang)])
CFLAGS="${CFLAGS} ${OPENMP_CFLAGS}"
CXXFLAGS="${CXXFLAGS} ${OPENMP_CXXFLAGS}"
LIBS="${LIBS} ${OPENMP_CXXFLAGS}"
if test x"${have_openmp}" == x"yes"; then
	CPPFLAGS="${CPPFLAGS} -DHAVE_OPENMP"
fi

#############################################################################
# Intel TBB
have_tbb="unknown"; # Perform the check only once when invoking MY_INTEL_TBB()
AC_DEFUN([MY_INTEL_TBB], [
	if test x"${have_tbb}" == x"unknown"; then
		have_tbb="yes" # assume so
		AC_CHECK_HEADER([tbb/tbb.h], [
			AC_SEARCH_LIBS([TBB_runtime_interface_version], [tbb], 
			[CPPFLAGS="${CPPFLAGS} -DHAVE_TBB"], 
			[have_tbb="no"])
		], [have_tbb="no"], [ [/* avoid default includes */] ])
	fi
	
	if test x"${have_tbb}" == x"no"; then
		AC_MSG_ERROR(["The dependency library Intel TBB is not available in this machine"]);
	fi
])

#############################################################################
# libevent >2.1.1
have_libevent="unknown"; # Perform the check only once when invoking MY_LIBEVENT()
AC_DEFUN([MY_LIBEVENT], [
    if test x"${have_libevent}" == x"unknown"; then
        for header in "event2/event.h" "event2/thread.h"; do
            AC_CHECK_HEADER([${header}], [],
                [AC_MSG_ERROR([missing prerequisite: this program requires the headers for libevent v2.1.x with support for pthreads]) ],
                [ [/* avoid default includes */] ])
        done
        AC_SEARCH_LIBS([libevent_global_shutdown], [event_core event], [],
            [ AC_MSG_ERROR([missing prerequisite: this program requires libevent]) ])
        AC_SEARCH_LIBS([evthread_use_pthreads], [event_pthreads event], [],
            [ AC_MSG_ERROR([missing prerequisite: this program requires libevent with support for pthreads]) ])
        have_libevent="yes";
    fi
])
    
#############################################################################
# pthreads
AC_SEARCH_LIBS([pthread_create], [pthread], [],
    [ AC_MSG_ERROR([missing prerequisite: this program requires pthreads to work (dependency for sqlite3)]) ])
    
#############################################################################
# Link with tcmalloc ?
AC_ARG_WITH([tcmalloc], AS_HELP_STRING([--with-tcmalloc], [Whether to link the program with tcmalloc, a replacement of the glibc malloc implementation.]))
AS_IF(
    [test "x${with_tcmalloc}" = "xyes"], [tcmalloc_system=yes; tcmalloc_bundle=yes;],
    [test "x${with_tcmalloc}" = "xsystem"], [tcmalloc_system=yes;],
    [test "x${with_tcmalloc}" = "xbundled" -o "x${with_tcmalloc}" = "xbundle"], [tcmalloc_bundle=yes;],
    [test "x${with_tcmalloc}" = "x"], [], dnl nop
    [AC_MSG_FAILURE([Invalid argument for the argument --with-tcmalloc="${with_tcmalloc}". Recognised arguments are: yes, system, bundled.])]
)

# Check whether tcmalloc is already available in the system
if( test "x${tcmalloc_system}" = "xyes" ); then
    AC_SEARCH_LIBS([tc_malloc], [tcmalloc_minimal tcmalloc], [tcmalloc_avail=yes])
fi
if( test "x${tcmalloc_bundle}" = "xyes" -a "x${tcmalloc_built}" != "xyes" ); then
    abs_top_builddir=$(realpath .)
    AX_SUBDIRS_CONFIGURE([third-party/gperftools-2.7], [[--enable-minimal], [--prefix=${abs_top_builddir}/third-party/gperftools-2.7/install]], [[--disable-shared], [--enable-static]])
    EXTRA_LDFLAGS="${EXTRA_LDFLAGS} -L${abs_top_builddir}/third-party/gperftools-2.7/install/lib -ltcmalloc_minimal"
    tcmalloc_avail=yes;
fi

if( test "x${tcmalloc_avail}" != "xyes" ); then
    if( test "x${tcmalloc_system}" = "xyes" -o "x${tcmalloc_bundle}" = "xyes" ); then
        AC_MSG_FAILURE([Unable to link the library tcmalloc])
    fi
fi

# Only to be seen by autoreconf -iv
if false; then
    AC_CONFIG_SUBDIRS([third-party/gperftools-2.7])
fi
    
#############################################################################
# libnuma
have_libnuma="yes"
AC_CHECK_HEADERS([numaif.h numa.h], [], [have_libnuma="no"; break;], [ [/* avoid default includes */] ])
AS_IF([test x"${have_libnuma}" == x"yes"], [AC_SEARCH_LIBS([numa_available], [numa], [], [have_libnuma="no"])])
if test x"${have_libnuma}" == x"yes"; then
    AC_MSG_NOTICE([libnuma support enabled...])
    CPPFLAGS="${CPPFLAGS} -DHAVE_LIBNUMA";
else
    AC_MSG_WARN([libnuma support disabled...])
fi

#############################################################################
# libpapi (profiler support)
have_libpapi="yes"
AC_CHECK_HEADERS([papi.h], [], [have_libpapi="no"; break;], [ [/* avoid default includes */] ])
AS_IF([test x"${have_libpapi}" == x"yes"], [AC_SEARCH_LIBS([PAPI_library_init], [papi], [], [have_libpapi="no"])])
if test x"${have_libpapi}" == x"yes"; then
    AC_MSG_NOTICE([libpapi support enabled...])
    CPPFLAGS="${CPPFLAGS} -DHAVE_LIBPAPI";
else
    AC_MSG_WARN([libpapi support disabled...])
fi

#############################################################################
# zlib
AC_CHECK_HEADERS([zlib.h], [], [ AC_MSG_ERROR([zlib.h not found. The library zlib is a required dependency.]) ], [ [/* avoid default includes */] ])
AC_SEARCH_LIBS([inflate], [z], [ ], [ AC_MSG_ERROR([Library zlib not found. This library is a required dependency.]) ], [ ]):

#############################################################################
# sqlite3, run-time support for libcommon. If not present, libcommon uses its own version bundled.
AX_LIB_SQLITE3()
LIBS="${SQLITE3_LDFLAGS} ${LIBS}"

#############################################################################
# Debug flags (-g)
MY_ARG_ENABLE([debug], 
    [Whether to enable the debug flags], 
    [yes no], [yes]) 
dnl first argument is the variable with the flags, the second argument is the language
m4_defun([_my_set_debug_flags], [
    m4_pushdef([_FLAGS], [m4_translit([$1], [+], [X])FLAGS]) dnl => CFLAGS, CXXFLAGS
    [if test -n "${ac_user_]m4_tolower(_FLAGS)[}"; then]
         AC_MSG_WARN([Action --enable-debug ignored as _FLAGS has been explicitly set through command line])
    else
        # Append either -g or -g3
        AX_CHECK_COMPILE_FLAG([[-g3]], [AS_VAR_APPEND([_FLAGS], " -g3")], [AS_VAR_APPEND([_FLAGS], " -g")] )

        # Force clang to emit the whole debug information
        AC_LANG_PUSH([$1])
        MY_SET_CC_FLAG([_FLAGS], [-fno-limit-debug-info])
        MY_SET_CC_FLAG([_FLAGS], [-fno-omit-frame-pointer])
        AC_LANG_POP([$1])
        
    fi
    m4_popdef([_FLAGS])
])
if( test x"${enable_debug}" = x"yes" ); then
    _my_set_debug_flags([C])
    _my_set_debug_flags([C++])
fi
m4_undefine([_my_set_debug_flags])

#############################################################################
# Warning flags (-Wall)
MY_ARG_ENABLE([warnings], 
    [Whether to enable all warnings (-Wall)], 
    [yes no], [yes])
m4_defun([_my_set_warnings], [
    m4_pushdef([_FLAGS], [m4_translit([$1], [+], [X])FLAGS]) dnl => CFLAGS, CXXFLAGS
    [if test -n "${ac_user_]m4_tolower(_FLAGS)[}"; then]
        AC_MSG_WARN([Action --enable-warnings ignored as _FLAGS has been explicitly set through command line])
    [else]
        AS_VAR_APPEND([_FLAGS], [" -Wall"])
    [fi]
    m4_popdef([_FLAGS])
])
if( test x"${enable_warnings}" = x"yes" ); then
    _my_set_warnings([C])
    _my_set_warnings([CXX])
fi
m4_undefine([_my_set_warnings])

#############################################################################
# Optimization flags (-O3)
MY_ARG_ENABLE([optimize], [Whether to enable the optimization flags], [yes no], [no])

m4_defun([_my_set_optimization_flags], [
    m4_pushdef([_FLAGS], [m4_translit([$1], [+], [X])FLAGS]) dnl => CFLAGS, CXXFLAGS
    [if test -n "${ac_user_]m4_tolower(_FLAGS)[}"; then]
        AC_MSG_WARN([Action --enable-optimize ignored as _FLAGS has been explicitly set through command line])
    [else]
        if( test x"${enable_optimize}" = x"yes" ); then
            AS_VAR_APPEND([_FLAGS], [[" -O3"]])
            AC_LANG_PUSH([$1])
            MY_SET_CC_FLAG([_FLAGS], [-march=native])
            MY_SET_CC_FLAG([_FLAGS], [-mtune=native])
            MY_SET_CC_FLAG([_FLAGS], [-fno-stack-protector])
            AC_LANG_POP([$1])
        else
            AS_VAR_APPEND([_FLAGS], [[" -O0"]])
        fi
    [fi]
    m4_popdef([_FLAGS])
])
_my_set_optimization_flags([C])
_my_set_optimization_flags([C++])
m4_undefine([_my_set_optimization_flags])

#############################################################################
# Assertions. Possible values:
# yes => nop
# no => CPPFLAGS += -DNDEBUG
# auto => yes if the optimize mode is not enabled, no otherwise
MY_ARG_ENABLE([assert], 
    [Whether to enable assertions. The option 'auto' implies the assertions are enabled when --enable-optimize is not specified], 
    [yes no auto], [auto])
    
if (test x"${enable_assert}" = x"auto"); then
    if (test x"${enable_optimize}" != x"yes"); then
        enable_assert=yes
    else
        enable_assert=no
    fi
fi
if (test x"${enable_assert}" = x"yes"); then
    : ; # nop
elif (test x"${enable_assert}" = x"no"); then
    CPPFLAGS="${CPPFLAGS} -DNDEBUG"
else
    AC_MSG_ERROR([Invalid value for --enable-assert: ${enable_assert}])
fi

#############################################################################
# CMake build type (for libcommon)
CMAKE_BUILD_TYPE="Debug" #default
if( test x"${enable_optimize}" = x"yes" ); then
    if( test x"${enable_debug}" = x"yes" ); then
        CMAKE_BUILD_TYPE="RelWithDebInfo"
    else # debug flags not set
        CMAKE_BUILD_TYPE="Release"
    fi
fi
AC_SUBST([CMAKE_BUILD_TYPE])

#############################################################################
# Switch to LLVM libc++ (-stdlib=libc++)
#MY_CHECK_STDLIB_LIBCXX([CXX="$CXX -stdlib=libc++"])

#############################################################################
# Support for the Teseo
# Usage:
# --with-teseo=/path/to/teseo/build
#

AC_ARG_WITH([teseo], [AS_HELP_STRING([--with-teseo@<:@=ARG@:>@], [
    Enable the support to Teseo. The argument needs to be the path where the library has been built.
])], [])

if ( test x"${with_teseo}" != x"" && test x"${with_teseo}" != "no" ); then
    if( test x"${with_teseo}" == x"yes" ); then
        AC_MSG_ERROR(["the full path to the Teseo build needs to be specified, e.g. --with-teseo=/path/to/teseo/build"]);
    fi
    
    MY_LIBEVENT()
    
    path_teseo=$(realpath "${with_teseo}") 
    path_teseo_library="${path_teseo}/libteseo.a"
    AC_CHECK_FILE([${path_teseo_library}], [  ],  [ AC_MSG_ERROR([file not found: ${path_teseo_library}]) ])
    
    path_teseo_makefile="${path_teseo}/Makefile"
    AC_CHECK_FILE([${path_teseo_makefile}], [  ],  [ AC_MSG_ERROR(["file not found: ${path_teseo_makefile}"]) ])
    
    path_teseo_includedir="$(realpath "${with_teseo}/$(make -C $(dirname ${path_teseo_makefile}) -n -p | awk '/^srcdir/ {print $3}')")/include" 
    AC_MSG_NOTICE([Teseo include path: ${path_teseo_includedir}])
    
    AS_VAR_APPEND([CPPFLAGS], [" -I${path_teseo_includedir} -I${path_teseo}/include"])
    AS_VAR_SET([LIBS], ["-L${path_teseo} -lteseo ${LIBS}"])
    
    # Try to compile a program with Teseo
    AC_LINK_IFELSE([
        AC_LANG_PROGRAM([[
            @%:@include "teseo.hpp"
        ]], [[
            teseo::Teseo instance;
        ]])
    ], [], [
        AC_MSG_FAILURE([unable to compile and link an executable with teseo])
    ])
    
    AC_MSG_NOTICE([Teseo support enabled])
    AS_VAR_APPEND([CPPFLAGS], [" -DHAVE_TESEO"])
fi

#############################################################################
# Support for the library Stinger
# Usage:
# default (auto) => attempts to link with stinger, but it does not abort in case it cannot 
# --without-stinger => disable support for stinger
# --with-stinger => attempt to link with stinger and aborts in case of error
# --with-stinger=/path/to/stinger/build => add the given path into the list of libs and include paths for the $CC, it fails in case it cannot link with stinger 
#
AC_ARG_WITH([stinger], [AS_HELP_STRING([--with-stinger@<:@=ARG@:>@], [
    Link with the (externally built) Stinger library. The argument needs to be the path where stinger is located. @<:@Default: auto@:>@
    ])], [])    
m4_divert_push([DEFAULTS])
    AS_VAR_SET([with_stinger], ["auto"])
m4_divert_pop([DEFAULTS])

if ( test -z "${with_stinger}" ); then with_stinger="no"; dnl transform the empty string '' into 'no'
elif ( test x"${with_stinger}" != x"yes" && test x"${with_stinger}" != x"auto" ); then dnl --with-stinger=/path/to/stinger/build
    AS_VAR_APPEND([CPPFLAGS], [" -I${with_stinger}/include"])
    AS_VAR_APPEND([LIBS], [" -L${with_stinger}/lib -Wl,-rpath=${with_stinger}/lib"])
fi
    
dnl first argument is the function, second argument is the library
m4_defun([_check_stinger_lib], [
    if ( test x"${with_stinger}" != x"no" ); then
        AC_SEARCH_LIBS([$1], [$2], [], [
            if ( test x"${with_stinger}" = x"auto" ); then
                AC_MSG_NOTICE([stinger support disabled, missing -l$2])
                with_stinger="no";
            else
                AC_MSG_ERROR([cannot link with the library $2])
            fi
        ])
    fi
])
_check_stinger_lib([stinger_new], [stinger_core])
_check_stinger_lib([load_metisish_graph], [stinger_utils])
_check_stinger_lib([page_rank], [stinger_alg])
m4_undefine([_check_stinger_lib])
    
dnl check whether we can compile & link a simple program invoking functions from the stinger library
if ( test x"${with_stinger}" != x"no" ); then
    AC_LINK_IFELSE([
        AC_LANG_PROGRAM([[
            @%:@include "stinger_alg/shortest_paths.h"
            @%:@include "stinger_core/stinger.h"
            @%:@include "stinger_utils/metisish_support.h"
        ]], [[
            auto graph = stinger_new();
            load_metisish_graph(graph, "my_test.metis"); // the file my_test.metis does not exist, just to check whether it compiles
            a_star(graph, stinger_max_active_vertex(graph) +1, 0, 1, true); // bfs, to check stinger_alg
        ]])
    ], [], [
        if ( test x"${with_stinger}" = x"auto" ); then
            AC_MSG_NOTICE([stinger support disabled, cannot compile and link the test program])
            with_stinger="no";
        else
            AC_MSG_FAILURE([unable to compile and link an executable with stinger])
        fi
    ])
fi

dnl add the definition HAVE_STINGER into the list of preprocessor variables $CPPFLAGS
if( test x"${with_stinger}" != x"no" ); then
    AC_MSG_NOTICE([stinger support enabled])
    AS_VAR_APPEND([CPPFLAGS], [" -DHAVE_STINGER"])
fi

#############################################################################
# Support for the LLAMA framework
# Usage:
# default (auto) => it attempts to support llama, but it does not abort in case it cannot 
# --without-llama => disable support for llama
# --with-llama => attempt to support llama and aborts in case of error
# --with-llama=/path/to/llama_repo => add the given path into the list of include paths for the $CC, it fails in case it cannot enable llama 
#
AC_ARG_WITH([llama], [AS_HELP_STRING([--with-llama@<:@=ARG@:>@], [
    Add support for the LLAMA library. The argument needs to be the path where the llama library is located. @<:@Default: auto@:>@
    ])], [])    
m4_divert_push([DEFAULTS])
    AS_VAR_SET([with_llama], ["auto"])
m4_divert_pop([DEFAULTS])


if ( test -z "${with_llama}" ); then with_llama="no"; dnl transform the empty string '' into 'no'
elif ( test x"${with_llama}" != x"yes" && test x"${with_llama}" != x"auto" ); then dnl --with-llama=/path/to/llama_repo/llama/include
    with_llama=${with_llama%/} # remove the trailing slash if present
    llama_includedir="";

    AC_CHECK_FILE([${with_llama}/llama.h], [llama_includedir="${with_llama}"])
    
    if( test x"${llama_includedir}" == x"" ); then
        AC_CHECK_FILE([${with_llama}/include/llama.h], [llama_includedir="${with_llama}/include"])
    fi
    if( test x"${llama_includedir}" == x"" ); then
        AC_CHECK_FILE([${with_llama}/llama/include/llama.h], [llama_includedir="${with_llama}/llama/include"])
    fi

    if( test x"${llama_includedir}" == x"" ); then
        AC_MSG_FAILURE([Cannot find the header llama.h inside ${with_llama}])
    fi
    
    AS_VAR_APPEND([CPPFLAGS], [" -I${llama_includedir}"])
fi
    
dnl check whether we can compile & link a program that makes usage of omp.h
if ( test "${have_openmp}" == "no" ); then
    if ( test x"${with_llama}" == x"auto" ); then
        AC_MSG_NOTICE([llama support disabled: the OpenMP header omp.h is missing.])
        with_llama="no";
    else
        AC_MSG_FAILURE([llama: missing prerequisite header omp.h (OpenMP runtime library)])
    fi
fi

dnl check whether we have Intel TBB
if ( test x"${with_llama}" != x"no" &&  test x"${with_llama}" != x"auto" ); then
	AC_MSG_NOTICE([llama: checking whether Intel TBB is present in the machine...])
	MY_INTEL_TBB() dnl it raises an error if not satisfied
	AS_VAR_APPEND([CPPFLAGS], [" -DLLAMA_HASHMAP_WITH_TBB"])
fi

dnl check whether we can compile & link a simple program invoking functions from the llama library
if ( test x"${with_llama}" != x"no" ); then
    AC_MSG_CHECKING([whether we can compile and link against the llama library])
    AC_LINK_IFELSE([
        AC_LANG_PROGRAM([[
            #define LL_MEMORY_ONLY 
            #define LL_DELETIONS
            #define register /* disable the usage of register */
            @%:@include "llama.h"
            
            ll_memory_pool __w_pool; /* required global */
        ]], [[
            ll_database llama{"/tmp/llama"};
        ]])
    ], [
        AC_MSG_RESULT([yes])
    ], [
        AC_MSG_RESULT([no])
        if ( test x"${with_llama}" = x"auto" ); then
            AC_MSG_NOTICE([llama support disabled])
            with_llama="no";
        else
            AC_MSG_FAILURE([unable to compile and link an executable with llama])
        fi
    ])
fi

dnl check whether we have Intel TBB (auto)
if (test x"${with_llama}" == x"auto" ); then
    AC_MSG_NOTICE([llama: checking whether Intel TBB is present in the machine...])
    MY_INTEL_TBB() dnl it raises an error if not satisfied
    AS_VAR_APPEND([CPPFLAGS], [" -DLLAMA_HASHMAP_WITH_TBB"])
fi

dnl add the definition HAVE_LLAMA into the list of preprocessor variables $CPPFLAGS
if( test x"${with_llama}" != x"no" ); then
    AC_MSG_NOTICE([llama support enabled])
    AS_VAR_APPEND([CPPFLAGS], [" -DHAVE_LLAMA"])
    AS_VAR_APPEND([CPPFLAGS], [" -DLLAMA_FAIR_SHARED_MUTEX"])
fi

#############################################################################
# Support for the GraphOne framework
# Usage:
# --with-graphone=/path/to/graphone/build
#
AC_ARG_WITH([graphone], [AS_HELP_STRING([--with-graphone@<:@=ARG@:>@], [
    Link to the GraphOne library. The argument needs to be the path where the library has been built.
])], [])

if ( test x"${with_graphone}" != x"" && test x"${with_graphone}" != "no" ); then
	dnl check whether we can compile & link a program that makes usage of omp.h
	if ( test "${have_openmp}" == "no" ); then
        AC_MSG_FAILURE([GraphOne: missing prerequisite header omp.h (OpenMP runtime library)])
	fi
	
	AC_MSG_NOTICE([GraphOne: checking whether Intel TBB is present in the machine...])
	MY_INTEL_TBB()
	
	if( test x"${with_graphone}" == x"yes" ); then
		AC_MSG_ERROR(["the full path to the GraphOne build needs to be specified, e.g. --with-graphone=/path/to/graphone/build"]);
	fi
	
	path_onedata=$(realpath "${with_graphone}/onedata/libonedata64.a")
		AC_CHECK_FILE([${path_onedata}], 
		[ LIBS="${path_onedata} ${LIBS}" ], 
		[ AC_MSG_ERROR(["file not found: ${path_onedata}"]) ]
	)
		
	path_graphone_library=$(realpath "${with_graphone}/src/libsrc64.a")
	AC_CHECK_FILE([${path_graphone_library}], 
		[ LIBS="${path_graphone_library} ${LIBS}" ], 
		[ AC_MSG_ERROR(["file not found: ${path_graphone_library}"]) ]
	)

	path_graphone_cmakecache=$(realpath "${with_graphone}/CMakeCache.txt")
	AC_CHECK_FILE([${path_graphone_cmakecache}], [], [AC_MSG_ERROR(["file not found: ${path_graphone_cmakecache}. This file is used to infer the path to the source directory of GraphOne"])])
	path_graphone_srcdir=$(awk -F= '$1 ~ /CMAKE_HOME_DIRECTORY/ {print $2}' ${path_graphone_cmakecache})
	AC_MSG_NOTICE([GraphOne source directory: ${path_graphone_srcdir}])
    AS_VAR_APPEND([CPPFLAGS], [" -I${path_graphone_srcdir}/src -I${path_graphone_srcdir}/include -I${path_graphone_srcdir}/gview -I${path_graphone_srcdir}/analytics -I${path_graphone_srcdir}/onedata"])
    
	dnl check whether we can compile & link a simple program invoking functions from the graphone library
    AC_MSG_CHECKING([whether we can compile and link against the GraphOne library])
    AC_LINK_IFELSE([
        AC_LANG_PROGRAM([[
            #define OVER_COMMIT 
            #define TBB
           	#define PLAIN_GRAPH
           	#define B64
           	#define DEL
        	@%:@include <omp.h>
            @%:@include "graph.h"
			
            /* required globals */
            graph* g = new graph; 
			int THD_COUNT = omp_get_max_threads() -1;
        ]], [[
            /* empty main, it uses the global variable g as a global */
        ]])
    ], [
        AC_MSG_RESULT([yes])
    ], [
        AC_MSG_RESULT([no])
        AC_MSG_FAILURE([unable to compile and link an executable with GraphOne])
    ])
    
    AC_MSG_NOTICE([GraphOne support enabled])
    AS_VAR_APPEND([CPPFLAGS], [" -DHAVE_GRAPHONE"])
fi

#############################################################################
# Remove extra blanks from our variables
EXTRA_CPPFLAGS=$(echo ${EXTRA_CPPFLAGS} | xargs)
CPPFLAGS=$(echo ${CPPFLAGS} | xargs);
CFLAGS=$(echo ${CFLAGS} | xargs);
EXTRA_CFLAGS=$(echo ${EXTRA_CFLAGS} | xargs);
CXXFLAGS=$(echo ${CXXFLAGS} | xargs);
EXTRA_CXXFLAGS=$(echo ${EXTRA_CXXFLAGS} | xargs);
EXTRA_LDFLAGS=$(echo ${EXTRA_LDFLAGS} | xargs);
# these two variables are only for presentation, overriding won't achieve much
ALL_CFLAGS=$(echo ${EXTRA_CPPFLAGS} ${CPPFLAGS} ${EXTRA_CFLAGS} ${CFLAGS} | xargs)
ALL_CXXFLAGS=$(echo ${EXTRA_CPPFLAGS} ${CPPFLAGS} ${EXTRA_CXXFLAGS} ${CXXFLAGS} | xargs)
LDFLAGS="${LDFLAGS} ${LIBS} ${EXTRA_LDFLAGS}"

#############################################################################
# CC, CXX and linker additional output variables
AC_SUBST([EXTRA_CPPFLAGS])
AC_SUBST([EXTRA_CFLAGS])
AC_SUBST([EXTRA_CXXFLAGS])

#############################################################################
# Create the configure script
AC_OUTPUT

#############################################################################
# Final summary
echo \
"-------------------------------------------------
${PACKAGE_NAME} version ${PACKAGE_VERSION}
Compiler C..........: ${CC} ${ALL_CFLAGS}
Compiler C++........: ${CXX} ${ALL_CXXFLAGS}
Linker..............: ${CXX} ${LDFLAGS}
Enable assertions...: ${enable_assert}
Enable debug........: ${enable_debug}
Enable optimize.....: ${enable_optimize}

Now type 'make -j'
--------------------------------------------------"
