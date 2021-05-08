#!/bin/sh
# compute dependencies for the PWscf directory tree

# make sure there is no locale setting creating unneeded differences.
export LC_ALL=C

QE_SOURCE=/pscratch/adbacze/Repos/qe-tddft/q-e-qe-6.7.0
DEPENDS="${QE_SOURCE}/include ${QE_SOURCE}/upflib 
         ${QE_SOURCE}/Modules ${QE_SOURCE}/PW/src
         ${QE_SOURCE}/FFTXlib ${QE_SOURCE}/UtilXlib"

cd src

${QE_SOURCE}/install/moduledep.sh $DEPENDS > make.depend
${QE_SOURCE}/install/includedep.sh $DEPENDS >> make.depend

# handle special cases
sed -i '/@\/cineca\/prod\/hpm\/include\/f_hpm.h@/d' make.depend
sed -i '/@iso_c_binding@/d;/@ifcore@/d' make.depend
    
# check for missing dependencies 
if grep @ make.depend
then
  echo WARNING: dependencies not found in directory 'src'
  exit 1
fi
