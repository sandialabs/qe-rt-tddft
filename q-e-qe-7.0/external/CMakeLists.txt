###########################################################
# FoX
###########################################################
include(fox.cmake)

###########################################################
# WANNIER90
###########################################################
include(wannier90.cmake)

###########################################################
# MBD
###########################################################
include(mbd.cmake)

###########################################################
# Eigensolver_GPU
###########################################################
if(QE_ENABLE_CUDA AND NOT QE_ENABLE_LAXLIB_CUSOLVER)
    include(eigensolver_gpu.cmake)
endif()

###########################################################
# DeviceXlib
###########################################################
include(devxlib.cmake)
