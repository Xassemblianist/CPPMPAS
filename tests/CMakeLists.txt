# SPDX-License-Identifier: GPL-3.0-or-later
# Tests on the NCAR test-case meshes. Data locations are cache variables; a test whose data is missing fails
# (a check that did not run is a failure).
set(CPPMPAS_CASES_DIR "" CACHE PATH "extracted NCAR test-case tarballs (jw_baroclinic_wave, supercell, mountain_wave)")
set(CPPMPAS_TEST_WORK "${CMAKE_BINARY_DIR}/test_work" CACHE PATH "scratch directory for test files")
set(CPPMPAS_MPAS_COMPARE "${CMAKE_SOURCE_DIR}/tools/mpas_compare.py" CACHE FILEPATH "tools/mpas_compare.py")
set(CPPMPAS_MPAS_REGISTRY "" CACHE FILEPATH "MPAS v8.4.2 src/core_atmosphere/Registry.xml")
find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_executable(test_decomp test_decomp.cpp)
target_link_libraries(test_decomp PRIVATE cppmpas)
add_test(NAME unit_decomp COMMAND test_decomp)

set(PY ${Python3_EXECUTABLE} -I)
set(T ${CMAKE_CURRENT_SOURCE_DIR})
set(EXE $<TARGET_FILE:cppmpas_cli>)
set(W ${CPPMPAS_TEST_WORK})
set(C ${CPPMPAS_CASES_DIR})
set(MW ${C}/mountain_wave/mtn_wave_grid.nc)
set(SC ${C}/supercell/supercell_grid.nc)
set(JW ${C}/jw_baroclinic_wave/x1.40962.grid.nc)

add_test(NAME make_inputs COMMAND ${PY} ${T}/make_inputs.py ${C} ${W})
set_tests_properties(make_inputs PROPERTIES FIXTURES_SETUP inputs)

# local edge / vertex order and in-memory arrays against an independent re-implementation
foreach(case MW SC JW)
  add_test(NAME order_${case} COMMAND ${CMAKE_COMMAND} -E env sh -c
    "'${EXE}' dump '${${case}}' '${W}/dump_${case}.nc' && ${Python3_EXECUTABLE} -I '${T}/check_order.py' '${${case}}' '${W}/dump_${case}.nc' ${case}")
endforeach()
add_test(NAME order_init_MW COMMAND ${CMAKE_COMMAND} -E env sh -c
  "'${EXE}' dump '${W}/mw_init.nc' '${W}/dump_mw_init.nc' --init && ${Python3_EXECUTABLE} -I '${T}/check_order.py' '${W}/mw_init.nc' '${W}/dump_mw_init.nc' MW")
set_tests_properties(order_init_MW PROPERTIES FIXTURES_REQUIRED inputs)

# read -> write -> tools/mpas_compare.py identical
foreach(case MW SC JW)
  add_test(NAME roundtrip_grid_${case} COMMAND ${PY} ${T}/roundtrip.py ${EXE} ${CPPMPAS_MPAS_COMPARE} ${${case}} ${W} grid)
endforeach()
add_test(NAME roundtrip_init_MW COMMAND ${PY} ${T}/roundtrip.py ${EXE} ${CPPMPAS_MPAS_COMPARE} ${W}/mw_init.nc ${W} init qv)
add_test(NAME roundtrip_init_SC COMMAND ${PY} ${T}/roundtrip.py ${EXE} ${CPPMPAS_MPAS_COMPARE} ${W}/sc_init.nc ${W} init qv,qc,qr)
set_tests_properties(roundtrip_init_MW roundtrip_init_SC PROPERTIES FIXTURES_REQUIRED inputs)

# history writer layout and values
set(HIST_FIELDS "latCell,lonCell,xCell,yCell,zCell,indexToCellID,latEdge,lonEdge,xEdge,yEdge,zEdge,indexToEdgeID,latVertex,lonVertex,xVertex,yVertex,zVertex,indexToVertexID,cellsOnEdge,nEdgesOnCell,nEdgesOnEdge,edgesOnCell,edgesOnEdge,weightsOnEdge,dvEdge,dcEdge,angleEdge,areaCell,areaTriangle,cellsOnCell,verticesOnCell,verticesOnEdge,edgesOnVertex,cellsOnVertex,kiteAreasOnVertex,meshDensity,zgrid,fzm,fzp,rdzw,cf1,cf2,cf3,zz,xtime,u,w,rho,theta,surface_pressure")
add_test(NAME history_MW COMMAND ${PY} ${T}/check_history.py ${EXE} ${W}/mw_init.nc ${CPPMPAS_MPAS_REGISTRY} ${W} "${HIST_FIELDS},qv" qv)
add_test(NAME history_SC COMMAND ${PY} ${T}/check_history.py ${EXE} ${W}/sc_init.nc ${CPPMPAS_MPAS_REGISTRY} ${W} "${HIST_FIELDS},qv,qc,qr" qv,qc,qr)
set_tests_properties(history_MW history_SC PROPERTIES FIXTURES_REQUIRED inputs)

# unhandled input is an error with a clear message
function(cppmpas_error_test name text)
  add_test(NAME err_${name} COMMAND ${PY} ${T}/expect_error.py "${text}" ${EXE} ${ARGN})
  set_tests_properties(err_${name} PROPERTIES FIXTURES_REQUIRED inputs)
endfunction()
cppmpas_error_test(periodic "is_periodic = YES" info ${W}/err_periodic.nc)
cppmpas_error_test(missing_var "required mesh variable dvEdge is missing" info ${W}/err_missing_dvEdge.nc)
cppmpas_error_test(float_var "variable dvEdge has NetCDF type 5" info ${W}/err_float_dvEdge.nc)
cppmpas_error_test(netcdf4 "NetCDF format netCDF-4/HDF5 is not handled" info ${W}/err_netcdf4.nc)
cppmpas_error_test(index_ids "indexToEdgeID(1) = 2" info ${W}/err_index_ids.nc)
cppmpas_error_test(boundary "lies on the mesh boundary (cellsOnCell = 0)" info ${W}/err_boundary.nc)
cppmpas_error_test(bdymask "limited-area meshes are not handled" info --init ${W}/err_bdymask.nc)
cppmpas_error_test(two_records "has 2 time records" info --init ${W}/err_two_records.nc)
cppmpas_error_test(missing_init_var "required variable zz of the input stream is missing" info --init ${W}/err_missing_zz.nc)
cppmpas_error_test(missing_scalar "active scalar qr is missing" info --init --scalars qv,qc,qr ${W}/err_missing_qr.nc)
cppmpas_error_test(bad_dims "variable u has dimensions (Time, nVertLevels, nEdges)" info --init ${W}/err_bad_dims_u.nc)
cppmpas_error_test(unknown_scalar "scalar 'qz' is not a constituent" info --init --scalars qv,qz ${W}/mw_init.nc)
cppmpas_error_test(unported_field "variable 'pressure' has no CPPMPAS port" history ${W}/mw_init.nc ${W}/x.nc --fields u,pressure)
cppmpas_error_test(unknown_field "variable 'notavar' is not defined in MPAS-Atmosphere's Registry" history ${W}/mw_init.nc ${W}/x.nc --fields notavar)
