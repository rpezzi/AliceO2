# Copyright CERN and copyright holders of ALICE O2. This software is distributed
# under the terms of the GNU General Public License v3 (GPL Version 3), copied
# verbatim in the file "COPYING".
#
# See http://alice-o2.web.cern.ch/license for full licensing information.
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.

message(STATUS "Building ITSMFT Common CUDA Library")
o2_add_library(ITSMFTCUDA 
                   SOURCES src/Context.cu
                   SOURCES src/Stream.cu
                   SOURCES src/Utils.cu
                   PUBLIC_LINK_LIBRARIES O2::GPUCommon
                   TARGETVARNAME targetName)


target_compile_definitions(
  ${targetName} PRIVATE $<TARGET_PROPERTY:O2::ITStracking,COMPILE_DEFINITIONS>
  ${targetName} PRIVATE $<TARGET_PROPERTY:O2::ITStrackingCUDA,COMPILE_DEFINITIONS>)
