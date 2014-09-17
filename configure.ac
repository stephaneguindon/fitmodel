#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

dnl AC_PREREQ(2.61)
AC_INIT([Fitmodel],esyscmd([sh -c "date \"+%Y%m%d\" | tr -d '\n'"]),[s.guindon@auckland.ac.nz])
dnl AM_SILENT_RULES([yes])
AM_INIT_AUTOMAKE([foreign])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES])

dnl AC_INIT([PhyML],[beta],[guindon@stat.auckland.ac.nz])
dnl AC_CONFIG_SRCDIR([simu.c])
AC_CONFIG_SRCDIR([utilities.c],[example/p1.nxs])
AC_CONFIG_HEADERS([config.h])

AC_DEFINE([UNIX],[1],[Unix tag on])
AC_DEFINE([DEBUG],[1],[Debug tag on])

AM_INIT_AUTOMAKE
AC_CANONICAL_HOST

AC_PROG_CC([gcc])

# Checks for libraries.
AC_CHECK_LIB([m], [log])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_SETVBUF_REVERSED
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([floor pow rint sqrt strchr strstr])

# Thanks to Shane Sturrock from Biomatters
#case "${host}" in
#*darwin*)
#AM_INIT_AUTOMAKE([no-dependencies])
#ARCH_flag="-arch i386 -arch ppc -arch x86_64 -mmacosx-version-min=10.4"
#;;
#esac


dnl CFLAGS="-O3 -fomit-frame-pointer -funroll-loops -Wall ${ARCH_flag}"
dnl CFLAGS="-O3 -Wfloat-equal -fomit-frame-pointer -funroll-loops"
dnl CFLAGS="-O3 -Wfloat-equal"
dnl CFLAGS="-g"

LDFLAGS="${ARGC_flag}"

dnl LT_INIT
dnl AC_PROG_LIBTOOL
dnl Add option for compiling with debugging informations

AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [Remove optimization options and add debug informations.])])
AS_IF([test "x$enable_debug" = "xyes"],
      [CFLAGS="-ansi -pedantic -Wall -std=c99 -O0 -g"],
      [CFLAGS="-Wall -O2 -msse -fomit-frame-pointer -funroll-loops ${ARCH_flag}"]) 

AC_ARG_ENABLE([gprof],
              [AS_HELP_STRING([--enable-gprof],
                              [Remove optimization options and add profiling informations.])])
AS_IF([test "x$enable_gprof" = "xyes"],
      [CFLAGS="-ansi -pedantic -Wall -std=c99 -O0 -g -pg"],
      [AS_IF([test "x$enable_debug" = "xno"],[CFLAGS="-Wall -O2 -msse -fomit-frame-pointer -funroll-loops ${ARCH_flag}"])]) 

AC_ARG_ENABLE([win],
              [AS_HELP_STRING([--enable-win],
                              [Compile with mingw instead of gcc.])])
AS_IF([test "x$enable_win" = "xyes"],[CC="i586-mingw32msvc-cc"])


AC_ARG_ENABLE([evolve],[AS_HELP_STRING([--enable-evolve],[Compile evolve])],[evolve=yes],[evolve=no])
AM_CONDITIONAL([WANT_EVOLVE], [test "$evolve" = yes]) 
if test "$evolve" = yes; then
   AC_DEFINE([EVOLVE],[1],[EVOLVE tag on])
fi

AC_ARG_ENABLE([coltree],[AS_HELP_STRING([--enable-coltree],[Compile coltree])],[coltree=yes],[coltree=no])
AM_CONDITIONAL([WANT_COLTREE], [test "$coltree" = yes]) 
if test "$coltree" = yes; then
   AC_DEFINE([COLTREE],[1],[COLTREE tag on])
fi


if test "$coltree" = no; then
if test "$evolve" = no; then
   AC_DEFINE([FITMODEL],[1],[FITMODEL tag on])
fi
fi
AC_CONFIG_FILES([Makefile])
dnl AC_CONFIG_FILES([Makefile])
AC_OUTPUT
