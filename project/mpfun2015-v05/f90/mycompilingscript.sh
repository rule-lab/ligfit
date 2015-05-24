gfortran -O2 -c mpmodule.f90 mpfuna.f90 mpfunbq.f90 mpfunc.f90 \
  mpfund.f90 mpfune.f90 mpfunfq1.f90

echo
echo " ---- First compile command complete, starting second compile command ---- "
echo

gfortran -O2 -c mpmodule.f90 mpfuna.f90 mpfunbq.f90 mpfunc.f90 \
  mpfund.f90 mpfune.f90 mpfunfq1.f90
