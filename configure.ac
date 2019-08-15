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
AC_PROG_CC([clang gcc icc cc])
AC_PROG_CXX([clang++ g++ icpc c++ cxx])
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
# git, dependency for libcommon
AC_CHECK_PROGS(git, [git], [no])
if( test x"${git}" = x"no" ); then
   AC_MSG_ERROR([missing prerequisite: this program requires git (dependency for libcommon)])
fi

#############################################################################
# OpenMP
m4_foreach_w([lang], [C, C++], [AC_LANG_PUSH(lang) AX_OPENMP AC_LANG_POP(lang)])
CFLAGS="${CFLAGS} ${OPENMP_CFLAGS}"
CXXFLAGS="${CXXFLAGS} ${OPENMP_CXXFLAGS}"
LIBS="${LIBS} ${OPENMP_CXXFLAGS}"

#############################################################################
# pthreads
AC_SEARCH_LIBS([pthread_create], [pthread], [],
    [ AC_MSG_ERROR([missing prerequisite: this program requires pthreads to work (dependency for sqlite3)]) ])

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
# GMP GNU Multiple Precision Arithmetic Library
# check if GMP is already available in the system
gmp_system_avail="yes"
AC_CHECK_HEADERS([gmpxx.h], [], [gmp_system_avail="no"; break;], [ [/* avoid default includes */] ])
AS_IF([test x"${gmp_system_avail}" == x"yes"], [ AC_CHECK_LIB([gmp], [__gmpz_init], [ nop="nop" ], [gmp_system_avail="no";]) ])
AS_IF([test x"${gmp_system_avail}" == x"yes"], [ AC_CHECK_LIB([gmpxx], [_ZlsRSoPK12__mpz_struct], [ nop="nop" ], [gmp_system_avail="no";]) ])

dnl Use the bundled package of GMP inside third-party/gmp-*
AS_IF([test x"${gmp_system_avail}" == x"no"], [
    m4_define([GMP_BUNDLED_VERSION], [6.1.2])
    AC_MSG_NOTICE([Configuring the bundled package GNU GMP v]GMP_BUNDLED_VERSION)
    abs_top_builddir=$(realpath .)
    AX_SUBDIRS_CONFIGURE([third-party/gmp-]GMP_BUNDLED_VERSION, [[--enable-cxx], [--prefix=${abs_top_builddir}/third-party/gmp-]GMP_BUNDLED_VERSION[/build], [CC=${CC}], [CXX=${CXX}]], [], [], [])
    EXTRA_CPPFLAGS="${EXTRA_CPPFLAGS} -isystem ${abs_top_builddir}/third-party/gmp-GMP_BUNDLED_VERSION/include"
    EXTRA_LDFLAGS="${EXTRA_LDFLAGS} -L${abs_top_builddir}/third-party/gmp-GMP_BUNDLED_VERSION/lib -lgmpxx -lgmp"
], dnl else 
[ LIBS="-lgmpxx -lgmp ${LIBS}" ]
) dnl end if


#AC_ARG_WITH([gmp], AS_HELP_STRING([--with-gmp], [Link with]))
#AS_IF(
#    [test "x${with_tcmalloc}" = "xyes"], [tcmalloc_system=yes; tcmalloc_bundle=yes;],
#    [test "x${with_tcmalloc}" = "xsystem"], [tcmalloc_system=yes;],
#    [test "x${with_tcmalloc}" = "xbundled" -o "x${with_tcmalloc}" = "xbundle"], [tcmalloc_bundle=yes;],
#    [test "x${with_tcmalloc}" = "x"], [], dnl nop
#    [AC_MSG_FAILURE([Invalid argument for the argument --with-tcmalloc="${with_tcmalloc}". Recognised arguments are: yes, system, bundled.])]
#)

# Check whether tcmalloc is already available in the system
#if( test "x${tcmalloc_system}" = "xyes" ); then
#    AC_SEARCH_LIBS([tc_malloc], [tcmalloc_minimal tcmalloc], [tcmalloc_avail=yes])
#fi
#if( test "x${tcmalloc_bundle}" = "xyes" -a "x${tcmalloc_built}" != "xyes" ); then
#    abs_top_builddir=$(realpath .)
#    AX_SUBDIRS_CONFIGURE([third-party/gperftools-2.7], [[--enable-minimal], [--prefix=${abs_top_builddir}/third-party/gperftools-2.7]], [], [], [])
#    EXTRA_LDFLAGS="${EXTRA_LDFLAGS} -L${abs_top_builddir}/third-party/gperftools-2.7/lib -Wl,-rpath=${abs_top_builddir}/third-party/gperftools-2.7/lib  -ltcmalloc_minimal"
#    tcmalloc_avail=yes;
#fi
#
#if( test "x${tcmalloc_avail}" != "xyes" ); then
#    if( test "x${tcmalloc_system}" = "xyes" -o "x${tcmalloc_bundle}" = "xyes" ); then
#        AC_MSG_FAILURE([Unable to link the library tcmalloc])
#    fi
#fi


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
        CMAKE_BUILD_TYPE="RELWITHDEBINFO"
    else # debug flags not set
        CMAKE_BUILD_TYPE="Release"
    fi
fi
AC_SUBST([CMAKE_BUILD_TYPE])

#############################################################################
# Switch to LLVM libc++ (-stdlib=libc++)
MY_CHECK_STDLIB_LIBCXX([CXX="$CXX -stdlib=libc++"])

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

LIBS="${LIBS} ${EXTRA_LDFLAGS}"

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
Linker..............: ${LIBS}
Enable assertions...: ${enable_assert}
Enable debug........: ${enable_debug}
Enable optimize.....: ${enable_optimize}

Now type 'make -j'
--------------------------------------------------"
