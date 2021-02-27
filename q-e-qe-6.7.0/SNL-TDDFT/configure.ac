# Autoconf template for snl-tddft
AC_INIT(snl-tddft, 6.4.1, adbacze@sandia.gov, snl-tddft)

# Option: Quantum-Espresso source directory (QE_SOURCE)
AC_ARG_WITH(qe-source,
   [AC_HELP_STRING([--with-qe-source=DIR],[Quantum-Espresso source directory (QE_SOURCE)])],
   [qe_source=$withval],
   [qe_source=$QE_SOURCE])

# make path absolute
abspath=`readlink -f ${qe_source} 2>/dev/null` && qe_source=$abspath

# Look for make.sys
AC_CHECK_FILE(${qe_source}/make.inc, [], [AC_MSG_ERROR([Cannot find make.inc in ${qe_source}])])

# Look for Modules/version.f90
AC_CHECK_FILE(${qe_source}/Modules/version.f90, [], [AC_MSG_ERROR([Cannot find version.f90 in ${qe_source}/Modules])])

# Extract version
AC_MSG_CHECKING([Quantum-Espresso version])
version=`grep version_number ${qe_source}/Modules/version.f90 | cut -d\' -f2`
AC_MSG_RESULT([$version])
case $version in
  6.4.1   ) ;;
  *) AC_MSG_ERROR([Cannot compile against this version of Quantum-Espresso]) ;;
esac

# export settings to generated files
AC_SUBST(qe_source)
AC_SUBST(dflags_tddft)
AC_OUTPUT([makedeps.sh src/Makefile], [chmod +x makedeps.sh])

# default commands
./makedeps.sh
