AC_INIT(origami,0.1)

AC_LANG_PUSH([C++])

AC_MSG_NOTICE([Checking for software tools needed by origami])

AC_CHECK_TOOLS([CUTADAPTOR],[cutadapt],[no])
AC_CHECK_TOOLS([MACS2],[macs2],[no])
AC_CHECK_TOOLS([BOWTIE],[bowtie],[no])
AC_CHECK_TOOLS([BEDTOOLS],[bedtools],[no])
AC_CHECK_TOOLS([REXE],[R],[no])
AC_PROG_CXX([g++ gcc])

if test "$ac_cv_prog_CUTADAPTOR" == no
then
	AC_MSG_ERROR([Cannot find cutadapt in PATH, please either add cutadapt to PATH or set CUTADAPTOR=/path/to/cutadapt manually])
fi

if test "$ac_cv_prog_MACS2" == no
then
  AC_MSG_ERROR([Cannot find macs2 in PATH, please either add macs2 to PATH or set MACS2=/path/to/macs2 manually])
fi

if test "$ac_cv_prog_BOWTIE" == no
then
  AC_MSG_ERROR([Cannot find bowtie in PATH, please either add bowtie to PATH or set BOWTIE=/path/to/bowtie manually])
fi

if test "$ac_cv_prog_BEDTOOLS" == no
then
  AC_MSG_ERROR([Cannot find bedtools in PATH, please either add bedtools to PATH or set BEDTOOLS=/path/to/bedtools manually])
fi

if test "$ac_cv_prog_REXE" == no
then
  AC_MSG_ERROR([Cannot find R in PATH, please either add R to PATH or set REXE=/path/to/R manually])
fi

AC_MSG_NOTICE([Checking for C++ libraries needed by origami])
AC_CHECK_HEADER([api/BamReader.h])
AC_CHECK_HEADER([api/BamWriter.h])

AC_MSG_NOTICE([Checking for R libraries needed by origami])
### GenomicRanges,

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([bin/origami-alignment],[chmod 775 bin/origami-alignment])

AC_OUTPUT


