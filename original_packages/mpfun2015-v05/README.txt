MPFUN2015: A thread-safe Fortran arbitrary precision computation package
Version date:  19 May 2015

AUTHOR:
   David H. Bailey
   Lawrence Berkeley National Lab (retired) and University of California, Davis
   Email: dhbailey@lbl.gov
   
COPYRIGHT AND DISCLAIMER:
  All software in this package (c) 2015 David H. Bailey.
  By downloading or using this software you agree to the copyright, disclaimer
  and license agreement in the accompanying file DISCLAIMER.txt.

1. PURPOSE OF PACKAGE:
  This package permits one to perform floating-point computations (real and
  complex) to arbitrarily high numeric precision, by making only relatively
  minor changes to existing Fortran-90 programs (mostly changes to type
  statements).  All basic arithmetic operations and transcendental functions
  are supported.  Advanced techniques, including FFT-based multiplication and
  quadratically convergent transcendental algorithms, are employed.
  
  In addition to fast execution times, one key feature of this package is a
  100% THREAD-SAFE design, which means that user-level applications can be
  easily converted for parallel execution, say using a threaded parallel
  environment such as OpenMP.  There are NO global shared variables (except
  static compile-time data), and NO initialization is necessary unless
  extremely high precision (> 19,500 digits) is required.

2. DOCUMENTATION:
   A detailed description of this package, and instructions for compiling
   and testing this program on various specific systems are included in the
   README file accompanying this package, and, in more detail, in the
   following technical paper:
   
   David H. Bailey, "MPFUN2015: A thread-safe arbitrary precision package," 
   available at http://www.davidhbailey.com/dhbpapers/mpfun2015.pdf.
   
3. COMPILATION AND LINKING:

  Compilation and linking is simple, provided that you have a command-line
  Unix-based system, such as Linux or Mac OSX.  
  
  For Mac OSX systems, you first must install the "Command Line Tools"
  package, which is available (for free) from the Apple Developer website.
  Login (or register, if this your first access) at the URL:
    https://developer.apple.com/devcenter/mac/index.action
  Then click on "View all downloads" and select "Command Line Tools" for
  your particular version of the MAC OSX operating system.

  The gfortran compiler, which is highly recommended for this package, is
  available for a variety of systems at the website
    https://gcc.gnu.org/wiki/GFortranBinaries
  The package also works with IBM's xlf_r compiler, Intel's ifort compiler
  and Portland Group's pgf90 compiler.

  There are actually four versions of this package, all of which are included
  in the distribution file:
  
  Version 1: This is recommended for most applications, particularly those
    that do not dynamically change the precision level.
  Version 2: This is recommended for more sophisticated applications
    that dynamically change the precision level (see below).
  Version Q1: This is the same as version 1, except that it includes
    interfaces for the real*16 (quad precision) datatype, which is now
    supported by the gfortran compiler on some platforms.
  Version Q2: This is the same as version 2, except that it includes
    interfaces for the real*16 (quad precision) datatype, which is now
    supported by the gfortran compiler on some platforms.

  The different four versions of the packages are merely four different versions
  of module MPFUNF, the high-level language interface module, and two different
  versions of module MPFUNB.  Compile/link scripts are available in the f90
  directory for each of the four compilers mentioned above (GNU's gfortran,
  IBM's xlf_r, Intel's ifort and Portland Group's pgf90), with the "Q" versions
  invoked for gfortran and ifort.

  For example, to compile version Q1 of the library with the GNU gfortran
  compiler, type
    ./gnu-complib1.scr
  and to compile and link the application program prog.f90 for version Q1,
  producing the executable file prog, type
    ./gnu-complink1.scr prog

  These scripts can be easily changed for other compilers or systems, since
  the MPFUN2015 programs do not invoke any special options or system-dependent
  features.

4. CODING INSTRUCTIONS AND USAGE:

  Here is a brief summary of Fortran coding instructions.  For full details,
  see the documentation paper mentioned above.
  
a. General instructions:

  First set the parameter mpipl, the "default" precision level in digits,
  which is the maximum precision level to be used for subsequent computation,
  and is used to set the amount of storage required for multiprecision data.  
  mpipl is set in a parameter statement at the start of module MPFUNF.  This
  module is in file mpfunf1.f90, mpfunq1.f90, mpfunf2.f90, or mpfunq2.f90,
  depending on which of the four versions is being used -- see the previous
  section.  By default, mpipl is set to 1000 digits, and all computations are
  performed to this precision unless otherwise specified by the user.  mpipl
  is automatically converted to mantissa words by the formula mpwds = 
  int (mpipl/log(2^{48}) + 2).  The parameter mpwds is the the internal default
  precision level for the MPFUN2015 software.

  Next, to invoke the package, place this line in every subprogram that contains
  a multiprecision variable or array, in the declaration section before
  any implicit or type statements: 
    use mpmodule

  To designate a variable or array as multiprecision real (MPR) in your
  application code, use the Fortran-90 type statement with "mp_real", e.g.,
    type (mp_real) a, b(m), c(m,n)
  Similarly, to designate a variable or array as multiprecision complex
  (MPC), use the type statement with "mp_complex".

  Thereafter when one of these variables or arrays appears in code, e.g.,
     d = a + b(i) * sqrt(3.d0 - c(i,j))
  the proper multiprecision routines are automatically called by the
  Fortran compiler.
  
  All usual mixed-mode combinations (operations, comparisons and assignments)
  involving MPR and double precision (DP) entities are supported, but not
  mixed-mode combinations of MPR and single-precision or integer data.
  Similarly, all usual mixed mode combinations of MPC and double complex (DC)
  entities are supported, but not mixed-mode combinations of MPC and 
  single-precision, double-precision or integer entities.  See the note
  below about DP and DC constants and expressions.

  Input/output of MP variables or array elements is done using the
  subroutines mpread and mpwrite.  See documentation for details.
  
  The above instructions apply if the precision level, namely mpipl, is 19,500
  digits or less.  For higher precision, in addition to changing mpipl to this
  higher level, one must call mpinit at the start of execution, before any
  multiprecision computation is done.  If this is a multithreaded application,
  this initialization must be done in single-threaded mode.  With version 1 or
  Q1, subroutine mpinit has an optional argument, which is the maximum precision
  level, in words; if not present, the default precision, namely mpwds words
  (which corresponds to mpipl digits), is assumed.  In version 2 or Q2, this
  argument is required.  See documentation for details.

b. Functions and subroutines:

  The following Fortran intrinsics are supported with multiprecision real (MPR)
  arguments, and they operate similarly to the standard double precision (DP)
  equivalents:
    abs, acos, aint, anint, asin, atan, atan2, cos, cosh, dble, exp, log,
    max, min, sign, sin, sqrt, tan and tanh.
  The dble intrinsic returns a DP approximation for a MPR argument.
    
  The following Fortran intrinsics are supported with multiprecision complex
  (MPC) arguments:
    abs, aimag, conjg, dcmplx, exp, log and sqrt.

  The mpreal function may be used to obtain the real part of a MPC argument
  argument (see table below), and the aimag intrinsic may be used to obtain
  the imaginary part.  The dcmplx intrinsic returns a DC approximation for a
  MPC argument.

  Here is a table of special functions and subroutines provided in this package.
  In this table, "F" denotes function, "S" denotes subroutine, "MPR" denotes
  multiprecision real, "MPC" denotes multiprecision complex, "DP" denotes
  double precision, "DC" denotes double complex, "Int" denotes integer and
  "Q" denotes quad precision or real*16.  The variable names r1, r2 and r3 are 
  multiprecision real, d1 and d2 are double precision, dc1 is double complex,
  s1 is type character*1, str1 is type character*(*), and iu, n1 and n2
  denotes integers.  The two items at the end are provided with versions Q1
  and Q2.

  F(MPR) : bessj (d1, r1) :  BesselJ function (for integer or half-integer d1).
  F(MPR) : gamma (r1) : Gamma function.
  F(MPC) : mpcmplx (r1,r2) :  Converts (r1,r2) to MPC. [1]   
  F(MPC) : mpcmplx (dc1) :  Converts dc1 to MPC. [1]
  F(MPC) : mpcmplx (z1) :  Converts z1 to MPC. [1]
  F(MPC) : mpcmplxdc (dc1) : Converts dc1 to MPC. [1, 2]
  S      : mpcssh (r1,r2,r3) : Returns both cosh and sinh of r1, in the same
         :   time as calling just cosh or just sinh.
  S      : mpcssn (r1,r2,r3) :  Returns both cos and sin of r1, in the same
         :   time as calling just cos or just sin.
  S      : mpeform (r1,n1,n2,s1) :  Returns r1 in decimal array s1 of length n1
         :   in En1.n2 format, suitable for output. [3]
  S      : mpfform (r1,n1,n2,s1) :  Returns r1 in decimal array s1 of length n1
         :   in Fn1.n2 format, suitable for output. [3]
  S      : mpinit :  Initializes for higher levels of precision.
  F(MPR) : mplog2 () :  Returns log (2). [1]
  F(MPR) : mpnrt (r1,n1) :  Returns the n1-th root of r1.
  F(MPR) : mppi () :  Returns pi. [1]
  F(MPR) : mpprodd (r1,d1) :  Returns the product of r1 and d1. [2]
  F(MPR) : mpquotd (r1, d1) :  Returns the quotient of r1 and d1. [2]
  S      : mpread (iu,r1) : Inputs r1 from Fortran unit iu.  Up to five
         :   MPR arguments may be listed. [1,3]
  S      : mpread (iu,z1) : Inputs z1 from Fortran unit iu.  Up to five
         :   MPC arguments may be listed. [1,3]
  F(MPR) : mpreal (r1) :  Converts r1 to MPR. [1]
  F(MPR) : mpreal (z1) :  Converts z1 to MPR. [1]
  F(MPR) : mpreal (d1) :  Converts d1 to MPR. [1, 2]
  F(MPR) : mpreal (str1) :  Converts str1 to MPR. [1, 2]
  F(MPR) : mpreald (d1) :  Converts d1 to MPR. [1, 2]
  S      : mpstrcon (s1,n1,r1)} & Converts s1, of length n1, to MPR [3].
  F(Int) : mpwprec (r1) :  Working precision level assigned to r1.
  F(Int) : mpwprec (z1) :  Working precision level assigned to z1.
  S      : mpwrite (iu,n1,n2,r1) :  Outputs r1 in En1.n2 format to unit iu.
         :   Up to five MPR arguments may be listed. [3]
  S      : mpwrite (iu,n1,n2,z1) :  Outputs z1 in En1.n2 format to unit iu.
         :   Up to five MPC arguments may be listed. [3]

  F(MPR) : mpreal (q1) :  Converts the real*16 argument to MPR. [1, 2]
  F(Q)   : qreal (r1) :  Converts MPR to real*16.

  Notes to above list:

  [1]: An additional integer argument may be added here as the final
  argument, which is the working precision level (in words) to be assigned
  to the result.  This argument is optional for versions 1 and Q1, but is
  mandatory for versions 2 and Q2.  See the note below on dynamically changing
  precision for additional details.

  [2]: These routines do NOT check DP or DC values.  See the note
  immediately below on double precision constants and expressions.
    
  [3]: These routines are for input/output of multiprecision data.  See
  documentation for details on their usage.
    
c. Note about double precision constants and expressions:

  While mixed-mode operations involving DP and MPR entities, or between
  DC and MPC entities, are permitted, there are some hazards:

  For example, the code r1 = 3.14159d0, where r1 is MPR, does NOT produce
  the true multiprecision equivalent of 3.14159, unless the numerical
  value is a modest-sized whole number or exact binary fraction.  To avoid
  problems, write instead r1 = '3.14159' or r1 = mpreal ('3.14159').  The
  usage of a literal here is a cue to the compiler to convert these
  constants to full multiprecision.  See also the note about dynamically
  varying precision below.  

  Similarly, the code r2 = r1 + 3.d0 * sqrt (2.d0), where r1 and r2 are
  MPR, does NOT produce the true multiprecision value, since the expression
  3.d0 * sqrt (2.d0) will be performed in double precision (according to
  standard Fortran-90 precedence rules).  To avoid problems, write instead
  something like r2 = r1 + 3.d0 * sqrt (mpreal (2.d0)), which forces all
  operations to be done using the multiprecision routines.

  To help avoid such problems, the MPFUN2015 low-level software checks
  *every* double precision value (constants, variables and expression values)
  in a multiprecision statement at *execution time* to see if it has more
  than 40 significant bits (i.e., if it has nonzero bits in the last 13 bits
  of the 53-bit mantissa).  If so, it is flagged as an error, since very
  likely such usage represented an unintended inaccurate value in the
  application program.  This feature catches 99.99% of accuracy loss problems
  due to double precision usage.

  On the other hand, some applications contain legitimate double precision
  constants that are trapped by this test.  In the author's two-level
  multipair PSLQ code, for example, double precision multipliers can arise
  that are greater than 40 bits in size.  In order to permit such usage, four
  special functions have been provided: mpprodd, mpquotd, mpreald, mpcmplxdc.
  The first and second return the product and quotient, respectively, of a MPR
  argument and a DP argument; the third converts a DP value to MPR (with an
  optional precision level parameter -- see next subsection); and the fourth
  converts a DC value to MPC (with an optional precision level parameter --
  see next subsection).  These routines do *not* check the double precision
  argument to see if it has more than 40 significant bits.

d. Note for applications that dynamically change precision level:

  Different applications have different requirements for language support.
  One distinction that immediately arises is between applications that do not
  need to change the working precision from the initially-defined default
  level (or change it only rarely) and those which, often for performance
  reasons, require that the working precision be changed frequently.

  Accordingly, there are four versions of the language interface module
  MPFUNF:
  
  Version 1 (in file mpfunf1.f90) and version Q1 (in file mpfunfq1.f90):
    Precision specifications are optional.  This version is recommended
    for applications that do not change the precision level, or do so
    only sparingly.
  Version 2 (in file mpfunf2.f90) and version Q2 (in file mpfunfq2.f90):
    Precision specifications are required.  This version is recommended
    for applications that frequently change the working precision level.

  See documentation for full details on the differences between these
  two versions.

  Note that the mpreal function, with the precision level (in words) as the
  second argument, can be used to convert an MPR argument to an MPR value
  with a different working precision level.  The same is true of mpcmplx.
  The working precision currently assigned to any MP variable or array
  element may be obtained by using the function mpwprec -- see Table above.
  
e. Sample application programs and output files:

  Six application programs (tpslq1.f90, tpslqm1.f90, tpslqm2.f90,
  tpslqm3.f90, tquadts.f90 and tquadtsp.f90 are included, together with
  corresponding output files for comparison with user results.  If, after
  compiling the library and running each of these programs, the results in
  these reference output files can be reproduced (except for timings,
  iteration counts, etc.), then one can be fairly confident that the
  software is working properly.  Full descriptions of these six programs
  are included in the documentation paper.


