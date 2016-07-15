AC_INIT(origami,0.1)

AC_CHECK_TOOLS([CUTADAPTOR],[cutadapt],[no])

if test "$ac_cv_prog_CUTADAPTOR" == no
then
	AC_MSG_ERROR([Cannot find cutadapt in PATH, please either add cutadapt to PATH or set CUTADAPTOR=/path/to/cutadapt manually],[1])
fi

AC_CHECK_TOOLS([MACS2],[macs2],[:])
AC_CHECK_TOOLS([BOWTIE],[bowtie],[:])
AC_CHECK_TOOLS([BEDTOOLS],[bedtools],[:])


AC_CHECK_HEADERS([api/BamReader.h api/BamWriter.h])


AC_CONFIG_FILES([Makefile])

AC_OUTPUT
