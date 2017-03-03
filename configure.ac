AC_INIT(origami,0.1)

####  

## params
## $1 -- library to check
AC_DEFUN([R_LIBRARY_TEST],[
  AC_MSG_CHECKING([if R library $1 is installed])
  VAL=`$ac_cv_prog_ac_ct_REXE --no-save &> /dev/null <<END 
library($1)
END
`

  
  if test $? -eq 0
  then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_RESULT([no])
    AC_MSG_WARN([R library $1 is not installed (or at least accessible by the identified R executable in REXE), please install before using the package])
  fi
])

####

AC_LANG_PUSH([C++])

AC_MSG_NOTICE([Checking for software tools needed by origami])

AC_ARG_VAR([BEDTOOLS],[Path to the bedtools executable])

AC_CHECK_TOOLS([CUTADAPTOR],[cutadapt],[no])
AC_CHECK_TOOLS([MACS2],[macs2],[no])
AC_CHECK_TOOLS([BOWTIE],[bowtie],[no])
AC_CHECK_TOOLS([BEDTOOLS],[bedtools],[no])
AC_CHECK_TOOLS([REXE],[R],[no])
AC_CHECK_TOOLS([RSCRIPTEXE],[Rscript],[no])
AC_CHECK_TOOLS([SAMTOOLS],[samtools],[no])
AC_CHECK_TOOLS([PERL],[perl],[no])
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

if test "$ac_cv_prog_RSCRIPTEXE" == no
then
  AC_MSG_ERROR([Cannot find Rscript in PATH, please either add Rscript to PATH or set RSCRIPTEXE=/path/to/Rscript manually])
fi

if test "$ac_cv_prog_SAMTOOLS" == no
then
  AC_MSG_ERROR([Cannot find samtools in PATH, please either add samtools to PATH or set SAMTOOLS=/path/to/samtools manually])
fi

if test "$ac_cv_prog_PERL" == no
then
  AC_MSG_ERROR([Cannot find perl in PATH, please either add perl to PATH or set PERL=/path/to/perl manually])
fi

AC_MSG_NOTICE([Checking for bedtools version])

if test `bedtools --version | cut -d ' ' -f 2 | cut -d '.' -f 2` -lt 24
then
  BEDTOOLSCOVERAGE="$BEDTOOLS coverage -counts -abam \$READSFILE -b \$TMPFILE > \$PREFIX-raw-counts.out"
else
  BEDTOOLSCOVERAGE="$BEDTOOLS coverage -sorted -counts -a \$TMPFILE -b \$SORTEDREADSFILE  > \$PREFIX-raw-counts.out"
fi

AC_SUBST([BEDTOOLSCOVERAGE])

#AC_MSG_NOTICE([Checking for C++ libraries needed by origami])
#AC_CHECK_HEADER([api/BamReader.h])
#AC_CHECK_HEADER([api/BamWriter.h])

AC_MSG_NOTICE([Checking for R libraries needed by origami])
R_LIBRARY_TEST([GenomicRanges])
R_LIBRARY_TEST([matrixStats])

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([bin/origami-alignment],[chmod 775 bin/origami-alignment])
AC_CONFIG_FILES([bin/origami-analysis],[chmod 775 bin/origami-analysis])


AC_OUTPUT


