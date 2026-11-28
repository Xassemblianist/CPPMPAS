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
add_executable(test_io_core test_io_core.cpp)
target_link_libraries(test_io_core PRIVATE cppmpas)

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
add_test(NAME unit_io_core COMMAND ${CMAKE_COMMAND} -E env sh -c
  "mkdir -p '${W}' && '$<TARGET_FILE:test_io_core>' '${W}' '${MW}' '${W}/mw_init_long_parent.nc'")
set_tests_properties(unit_io_core PROPERTIES FIXTURES_REQUIRED inputs)

# tools/mpas_compare.py: its unit tests and its mutation check (tests/compare)
add_test(NAME compare_unit COMMAND ${Python3_EXECUTABLE} -m unittest discover -s ${T}/compare)
add_test(NAME compare_mutants COMMAND ${Python3_EXECUTABLE} ${T}/compare/mutate.py)
set_tests_properties(compare_unit compare_mutants PROPERTIES ENVIRONMENT "MPAS_COMPARE_TMP=${W};PYTHONDONTWRITEBYTECODE=1"
  FIXTURES_REQUIRED work_dir)
add_test(NAME make_work_dir COMMAND ${CMAKE_COMMAND} -E make_directory ${W})
set_tests_properties(make_work_dir PROPERTIES FIXTURES_SETUP work_dir)

# src/io/registry_atm_table.inc is what tools/gen_registry_table.py makes of MPAS's Registry.xml
add_test(NAME registry_table_current COMMAND ${CMAKE_COMMAND} -E env sh -c
  "${Python3_EXECUTABLE} -I '${CMAKE_SOURCE_DIR}/tools/gen_registry_table.py' '${CPPMPAS_MPAS_REGISTRY}' > '${W}/registry_atm_table.inc' && ${CMAKE_COMMAND} -E compare_files '${W}/registry_atm_table.inc' '${CMAKE_SOURCE_DIR}/src/io/registry_atm_table.inc'")
set_tests_properties(registry_table_current PROPERTIES FIXTURES_REQUIRED work_dir)

# flags that change floating-point results are refused at configure time, also in CMAKE_CXX_FLAGS_<CONFIG>
foreach(flag Ofast ffast-math)
  add_test(NAME cmake_rejects_${flag} COMMAND ${CMAKE_COMMAND} -E env sh -c
    "rm -rf '${W}/cmake_${flag}'; out=$(${CMAKE_COMMAND} -S '${CMAKE_SOURCE_DIR}' -B '${W}/cmake_${flag}' -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS_RELEASE='-O2 -${flag}' 2>&1); rc=$?; rm -rf '${W}/cmake_${flag}'; echo \"$out\" | tail -5; [ $rc -ne 0 ] && echo \"$out\" | grep -q 'changes the floating-point results'")
  set_tests_properties(cmake_rejects_${flag} PROPERTIES FIXTURES_REQUIRED work_dir)
endforeach()

# the CPU build has none of the GPU's capacity limits (work arrays, columns, RBF points)
add_executable(test_backend test_backend.cpp)
target_link_libraries(test_backend PRIVATE cppmpas)
add_test(NAME unit_backend COMMAND test_backend)

# local edge / vertex order and in-memory arrays against an independent re-implementation
foreach(case MW SC JW)
  add_test(NAME order_${case} COMMAND ${CMAKE_COMMAND} -E env sh -c
    "rm -f '${W}/dump_${case}.nc' && '${EXE}' dump '${${case}}' '${W}/dump_${case}.nc' && ${Python3_EXECUTABLE} -I '${T}/check_order.py' '${${case}}' '${W}/dump_${case}.nc' ${case}")
  set_tests_properties(order_${case} PROPERTIES FIXTURES_REQUIRED work_dir)
endforeach()
add_test(NAME order_init_MW COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -f '${W}/dump_mw_init.nc' && '${EXE}' dump '${W}/mw_init.nc' '${W}/dump_mw_init.nc' --init && ${Python3_EXECUTABLE} -I '${T}/check_order.py' '${W}/mw_init.nc' '${W}/dump_mw_init.nc' MW --init qv")
add_test(NAME order_init_SC COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -f '${W}/dump_sc_init.nc' && '${EXE}' dump '${W}/sc_init.nc' '${W}/dump_sc_init.nc' --init --scalars qv,qc,qr && ${Python3_EXECUTABLE} -I '${T}/check_order.py' '${W}/sc_init.nc' '${W}/dump_sc_init.nc' SC --init qv,qc,qr")
set_tests_properties(order_init_MW order_init_SC PROPERTIES FIXTURES_REQUIRED inputs)
# an init file read as a grid: its vertical fields (nVertLevels and nVertLevelsP1) are left out
add_test(NAME grid_read_of_init COMMAND ${EXE} info ${W}/mw_init.nc)
add_test(NAME global_atts COMMAND ${PY} ${T}/check_global_atts.py ${EXE} ${MW} ${W}/mw_init_long_parent.nc ${W})
set_tests_properties(grid_read_of_init global_atts PROPERTIES FIXTURES_REQUIRED inputs)

# read -> write -> tools/mpas_compare.py identical
foreach(case MW SC JW)
  add_test(NAME roundtrip_grid_${case} COMMAND ${PY} ${T}/roundtrip.py ${EXE} ${CPPMPAS_MPAS_COMPARE} ${${case}} ${W} grid)
endforeach()
add_test(NAME roundtrip_init_MW COMMAND ${PY} ${T}/roundtrip.py ${EXE} ${CPPMPAS_MPAS_COMPARE} ${W}/mw_init.nc ${W} init qv)
add_test(NAME roundtrip_init_SC COMMAND ${PY} ${T}/roundtrip.py ${EXE} ${CPPMPAS_MPAS_COMPARE} ${W}/sc_init.nc ${W} init qv,qc,qr)
set_tests_properties(roundtrip_init_MW roundtrip_init_SC PROPERTIES FIXTURES_REQUIRED inputs)

# history writer layout and values
# every field the reader handles (fields of several ranks and kinds; the order mixes them so the dimension order
# MPAS defines is exercised)
set(HIST_FIELDS "xtime,u,latCell,lonCell,xCell,yCell,zCell,indexToCellID,latEdge,lonEdge,xEdge,yEdge,zEdge,indexToEdgeID,latVertex,lonVertex,xVertex,yVertex,zVertex,indexToVertexID,cellsOnEdge,nEdgesOnCell,nEdgesOnEdge,edgesOnCell,edgesOnEdge,weightsOnEdge,dvEdge,dcEdge,angleEdge,areaCell,areaTriangle,cellsOnCell,verticesOnCell,verticesOnEdge,edgesOnVertex,cellsOnVertex,kiteAreasOnVertex,meshDensity,fEdge,fVertex,edgeNormalVectors,localVerticalUnitVectors,cellTangentPlane,nominalMinDc,bdyMaskCell,bdyMaskEdge,bdyMaskVertex,cf1,cf2,cf3,zgrid,rdzw,dzu,rdzu,fzm,fzp,zxu,zz,zb,zb3,dss,deriv_two,coeffs_reconstruct,initial_time,u_init,v_init,t_init,w,rho,theta,rho_base,theta_base,surface_pressure")
add_test(NAME history_MW COMMAND ${PY} ${T}/check_history.py ${EXE} ${W}/mw_init.nc ${CPPMPAS_MPAS_REGISTRY} ${W} "${HIST_FIELDS},qv" qv 0123456789)
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
# a constituent the input file lacks keeps its default 0 with MPAS's warning (the RG input has no qi, qs, qg)
add_test(NAME err_missing_scalar COMMAND ${PY} ${T}/expect_error.py --warning "Variable qr not in input file."
  ${EXE} info --init --scalars qv,qc,qr ${W}/err_missing_qr.nc)
set_tests_properties(err_missing_scalar PROPERTIES FIXTURES_REQUIRED inputs)
cppmpas_error_test(bad_dims "variable u has dimensions (Time, nVertLevels, nEdges)" info --init ${W}/err_bad_dims_u.nc)
cppmpas_error_test(unknown_scalar "scalar 'qz' is not a constituent" info --init --scalars qv,qz ${W}/mw_init.nc)
cppmpas_error_test(unported_field "variable 'relhum' has no CPPMPAS port" history ${W}/mw_init.nc ${W}/x.nc --fields u,relhum)
cppmpas_error_test(boundary_edge "edge 8 lies on the mesh boundary (cellsOnEdge = 0)" info ${W}/err_boundary_edge.nc)
cppmpas_error_test(boundary_vertex "vertex 10 lies on the mesh boundary (cellsOnVertex = 0)" info ${W}/err_boundary_vertex.nc)
cppmpas_error_test(cellsOnEdge_range "cellsOnEdge(2, 8) = 1601 outside 1..1600" info ${W}/err_cellsOnEdge_range.nc)
cppmpas_error_test(cellsOnCell_neg "cellsOnCell(3, 6) = -3 outside 1..1600" info ${W}/err_cellsOnCell_neg.nc)
cppmpas_error_test(verticesOnEdge_zero "verticesOnEdge(1, 5) = 0 outside 1..3200" info ${W}/err_verticesOnEdge_zero.nc)
cppmpas_error_test(edgesOnVertex_range "edgesOnVertex(2, 4) = 4805 outside 1..4800" info ${W}/err_edgesOnVertex_range.nc)
cppmpas_error_test(edgesOnEdge_range "edgesOnEdge(1, 7) = 4801 outside 1..4800" info ${W}/err_edgesOnEdge_range.nc)
cppmpas_error_test(nEdgesOnEdge "nEdgesOnEdge(11) = 13 outside 0..maxEdges2" info ${W}/err_nEdgesOnEdge.nc)
cppmpas_error_test(nEdgesOnCell "nEdgesOnCell(4) = 7 outside 1..maxEdges" info ${W}/err_nEdgesOnCell.nc)
cppmpas_error_test(missing_nvl "dimension nVertLevels is missing" info --init ${W}/err_missing_nvl.nc)
cppmpas_error_test(zero_dim "dimension nVertices = 0 outside 1..2147483646" info ${W}/err_zero_dim.nc)
cppmpas_error_test(two "dimension TWO is not 2" info ${W}/err_two.nc)
cppmpas_error_test(long_att "global attribute mesh_spec has 513 characters" info --init ${W}/err_long_mesh_spec.nc)
cppmpas_error_test(strlen "StrLen = 32, MPAS writes 64" info --init ${W}/err_strlen.nc)
cppmpas_error_test(nvlp1 "nVertLevelsP1 != nVertLevels + 1" info --init ${W}/err_nvlp1.nc)
cppmpas_error_test(no_qv "qv is not among the active scalars" info --init --scalars qc ${W}/mw_init.nc)
cppmpas_error_test(dup_scalar "scalar 'qv' is given twice" info --init --scalars qv,qv ${W}/mw_init.nc)
cppmpas_error_test(empty_name "empty name in the list 'u,,w'" history ${W}/mw_init.nc ${W}/x.nc --fields u,,w)
add_test(NAME err_empty_list COMMAND ${PY} ${T}/expect_error.py "empty name in the list ''" sh -c
  "'${EXE}' roundtrip '${W}/mw_init.nc' '${W}/x.nc' --init --fields ''")
set_tests_properties(err_empty_list PROPERTIES FIXTURES_REQUIRED inputs)
cppmpas_error_test(trailing_comma "empty name in the list 'u,'" history ${W}/mw_init.nc ${W}/x.nc --fields u,)
cppmpas_error_test(unknown_field "variable 'notavar' is not defined in MPAS-Atmosphere's Registry" history ${W}/mw_init.nc ${W}/x.nc --fields notavar)
# a field that is neither read nor computed is refused before the output file is created (the output directory does
# not exist: creating the file would fail with another message)
add_test(NAME err_unheld_field COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/no_such_dir'; ${Python3_EXECUTABLE} -I '${T}/expect_error.py' 'write ru: the field was neither read nor computed' '${EXE}' history '${W}/mw_init.nc' '${W}/no_such_dir/unheld.nc' --fields u,ru")
set_tests_properties(err_unheld_field PROPERTIES FIXTURES_REQUIRED inputs)

# tests/compare_arrays.py itself: value bytes, variable types and attribute types are compared
add_test(NAME compare_arrays_unit COMMAND ${PY} ${T}/test_compare_arrays.py ${W})
set_tests_properties(compare_arrays_unit PROPERTIES FIXTURES_REQUIRED work_dir)

# ---------------------------------------------------------------------------------------------------------------
# libm / libmvec entries of the reference (src/core/libm.hpp) against this machine's glibc 2.41 (FMA, AVX2)
add_executable(test_libm test_libm.cpp)
target_link_libraries(test_libm PRIVATE cppmpas mvec m)
add_test(NAME unit_libm COMMAND test_libm 2000000)

# src/core/namelist_table.inc is what tools/gen_namelist_table.py makes of the Registry and the reference namelists,
# and of the namelist of the MW restart run (tests/namelists/MW_restart: the MW forecast namelist restarted from the
# reference's restart file at 00:30:06, verified by the dycore_restart_MW_* tests)
add_test(NAME namelist_table_current COMMAND ${CMAKE_COMMAND} -E env sh -c
  "${Python3_EXECUTABLE} -I '${CMAKE_SOURCE_DIR}/tools/gen_namelist_table.py' '${CPPMPAS_MPAS_REGISTRY}' ${CPPMPAS_REFS_DIR}/mpas842-dbl-mw-001/MW/steps/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-mw-001/MW/forecast/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-sc-001/SC/steps/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-sc-001/SC/forecast/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-jw-001/JW/steps/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-jw-001/JW/forecast/namelist.atmosphere '${CMAKE_SOURCE_DIR}/tests/namelists/MW_restart/namelist.atmosphere' ${CPPMPAS_REFS_DIR}/mpas842-dbl-rg-004/RG/steps/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-rg-004/RG/forecast/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-rg24-001/RG24/steps/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-rg24-001/RG24/forecast/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-rgcp-001/RGCP/steps/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-rgcp-001/RGCP/forecast/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-rgcp24-001/RGCP24/steps/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-rgcp24-001/RGCP24/forecast/namelist.atmosphere > '${W}/namelist_table.inc' && ${CMAKE_COMMAND} -E compare_files '${W}/namelist_table.inc' '${CMAKE_SOURCE_DIR}/src/core/namelist_table.inc'")
set_tests_properties(namelist_table_current PROPERTIES FIXTURES_REQUIRED work_dir)

# ---------------------------------------------------------------------------------------------------------------
# The MPAS v8.4.2 FP64 references (Xassemblianist/cppmpas-data releases mpas842-dbl-{mw,sc,jw}-001, extracted under
# CPPMPAS_REFS_DIR): reader on the real init files, and atm_mpas_init_block against what MPAS wrote.
set(CPPMPAS_REFS_DIR "" CACHE PATH "directory holding mpas842-dbl-{mw,sc,jw}-001/ (extracted reference releases)")
set(REF_MW ${CPPMPAS_REFS_DIR}/mpas842-dbl-mw-001/MW)
set(REF_SC ${CPPMPAS_REFS_DIR}/mpas842-dbl-sc-001/SC)
set(REF_JW ${CPPMPAS_REFS_DIR}/mpas842-dbl-jw-001/JW)
set(INIT_MW ${REF_MW}/init/mtn_wave_init.nc)
set(INIT_SC ${REF_SC}/init/supercell_init.nc)
set(INIT_JW ${REF_JW}/init/x1.40962.init.nc)
set(RST_MW ${REF_MW}/forecast/restart.0000-01-01_00.30.06.nc)
set(RST_SC ${REF_SC}/forecast/restart.0000-01-01_00.30.03.nc)
set(RST_JW ${REF_JW}/forecast/restart.0000-01-02_00.07.30.nc)
set(SCAL_MW qv)
set(SCAL_SC qv,qc,qr)
set(SCAL_JW qv)
# arrays atm_core_init computes that MPAS writes: time-invariant ones in the restart stream (Registry.xml:
# 535-537, 564, 566, 680-682; nothing changes them after the init block; nominalMinDc is set from config_len_disp
# before it), and the time-level-1 state, pressure_p and xtime (= config_start_time) in the first (t = 0) record of
# the every-step history (written before step 1)
set(INIT_RESTART_FIELDS edgeNormalVectors,localVerticalUnitVectors,cellTangentPlane,coeffs_reconstruct,dss,exner_base,rtheta_base,pressure_base,nominalMinDc)
set(INIT_STEP0_FIELDS rho_zz,theta_m,pressure_p,xtime)
foreach(case MW SC JW)
  # the reader on the reference's own init files (round trip and in-memory order)
  add_test(NAME roundtrip_refinit_${case} COMMAND ${PY} ${T}/roundtrip.py ${EXE} ${CPPMPAS_MPAS_COMPARE} ${INIT_${case}}
    ${W}/ref_${case} init ${SCAL_${case}})
  add_test(NAME order_refinit_${case} COMMAND ${CMAKE_COMMAND} -E env sh -c
    "mkdir -p '${W}' && rm -f '${W}/dump_refinit_${case}.nc' && '${EXE}' dump '${INIT_${case}}' '${W}/dump_refinit_${case}.nc' --init --scalars ${SCAL_${case}} && ${Python3_EXECUTABLE} -I '${T}/check_order.py' '${INIT_${case}}' '${W}/dump_refinit_${case}.nc' ${case} --init ${SCAL_${case}}")
  # the reader on the reference restart file: its input-stream subset, values and in-memory order
  add_test(NAME readrestart_${case} COMMAND ${CMAKE_COMMAND} -E env sh -c
    "mkdir -p '${W}' && rm -f '${W}/dump_rst_${case}.nc' && '${EXE}' dump '${RST_${case}}' '${W}/dump_rst_${case}.nc' --init --scalars ${SCAL_${case}} && ${Python3_EXECUTABLE} -I '${T}/check_order.py' '${RST_${case}}' '${W}/dump_rst_${case}.nc' ${case} --init ${SCAL_${case}}")
  # atm_mpas_init_block on the reference init file with the reference namelist
  add_test(NAME init_run_${case} COMMAND ${CMAKE_COMMAND} -E env sh -c
    "mkdir -p '${W}' && rm -f '${W}/init_${case}.nc' && '${EXE}' init '${INIT_${case}}' '${W}/init_${case}.nc' --namelist '${REF_${case}}/steps/namelist.atmosphere' --scalars ${SCAL_${case}} --fields ${INIT_RESTART_FIELDS},${INIT_STEP0_FIELDS}")
  set_tests_properties(init_run_${case} PROPERTIES FIXTURES_SETUP init_${case})
  add_test(NAME init_restart_${case} COMMAND ${PY} ${T}/compare_arrays.py ${RST_${case}} ${W}/init_${case}.nc
    ${INIT_RESTART_FIELDS})
  add_test(NAME init_step0_${case} COMMAND ${PY} ${T}/compare_arrays.py ${REF_${case}}/steps/step.0000-01-01_00.00.00.nc
    ${W}/init_${case}.nc ${INIT_STEP0_FIELDS})
  set_tests_properties(init_restart_${case} init_step0_${case} PROPERTIES FIXTURES_REQUIRED init_${case})
endforeach()
# the history writer on MPAS's own JW init file: a sphere (on_a_sphere = YES) with pentagons
add_test(NAME history_JW COMMAND ${PY} ${T}/check_history.py ${EXE} ${INIT_JW} ${CPPMPAS_MPAS_REGISTRY} ${W}/ref_JW "${HIST_FIELDS},qv" qv)
# the namelist reader on variants of the MW reference namelist
add_test(NAME namelist_MW COMMAND ${PY} ${T}/check_namelist.py ${EXE} ${REF_MW}/steps/namelist.atmosphere ${W}/namelist)
# --unverified: every array of the init block is proven for an even nVertLevels (initdump tests), and every reference
# has one; the synthetic sc_init.nc has nVertLevels = 3, so exner (the pow tail) needs --unverified there
cppmpas_error_test(unverified "field exner is computed by a port not yet verified against MPAS output; --unverified writes it anyway"
  init ${W}/sc_init.nc ${W}/x.nc --namelist ${REF_SC}/steps/namelist.atmosphere --scalars qv,qc,qr --fields pv_edge,exner)
add_test(NAME init_unverified_SC COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -f '${W}/init_unverified_SC.nc' && '${EXE}' init '${W}/sc_init.nc' '${W}/init_unverified_SC.nc' --namelist '${REF_SC}/steps/namelist.atmosphere' --scalars qv,qc,qr --fields pv_edge,exner --unverified && test -s '${W}/init_unverified_SC.nc' && rm '${W}/init_unverified_SC.nc'")
set_tests_properties(init_unverified_SC PROPERTIES FIXTURES_REQUIRED inputs)
# meshDensity = 2 in one cell and config_h_ScaleWithMesh = true (the MW namelist with JW's value): the meshScaling*
# arrays need --unverified (mesh_scaling_unverified in tools/cppmpas.cpp)
add_test(NAME err_scaled_mesh COMMAND ${CMAKE_COMMAND} -E env sh -c
  "sed 's/config_h_ScaleWithMesh = false/config_h_ScaleWithMesh = true/' '${REF_MW}/steps/namelist.atmosphere' > '${W}/namelist_scaled' && grep -q 'config_h_ScaleWithMesh = true' '${W}/namelist_scaled' && ${Python3_EXECUTABLE} -I '${T}/expect_error.py' 'field meshScalingDel2 is computed by a port not yet verified against MPAS output; --unverified writes it anyway' '${EXE}' init '${W}/mw_init_dense.nc' '${W}/x.nc' --namelist '${W}/namelist_scaled' --scalars qv --fields meshScalingDel2")
set_tests_properties(err_scaled_mesh PROPERTIES FIXTURES_REQUIRED inputs)
# guards and atm_core_init steps on the MW reference init file
add_executable(test_init_guards test_init_guards.cpp)
target_link_libraries(test_init_guards PRIVATE cppmpas)
add_test(NAME unit_init_guards COMMAND test_init_guards ${INIT_MW} ${REF_MW}/steps/namelist.atmosphere)
# the namelist reader's values and refusals (bit patterns of reals, 32-bit integers, groups, StrKIND), the scalars the
# namelist activates, and the table generator's refusals
add_executable(test_namelist test_namelist.cpp)
target_link_libraries(test_namelist PRIVATE cppmpas)
add_test(NAME unit_namelist COMMAND ${CMAKE_COMMAND} -E env sh -c
  "mkdir -p '${W}' && '$<TARGET_FILE:test_namelist>' '${REF_MW}/steps/namelist.atmosphere' '${W}'")
# reals in a process locale whose decimal point is ',' (its own scratch directory: unit_namelist writes the same file
# name, and CTest may run the two at once) (needs a de_DE, fr_FR or tr_TR locale: reported as skipped, not passed, without one)
add_test(NAME unit_namelist_locale COMMAND ${CMAKE_COMMAND} -E env sh -c
  "mkdir -p '${W}/namelist_locale' && '$<TARGET_FILE:test_namelist>' '${REF_MW}/steps/namelist.atmosphere' '${W}/namelist_locale' locale")
set_tests_properties(unit_namelist_locale PROPERTIES SKIP_RETURN_CODE 77)
add_test(NAME namelist_gen_errors COMMAND ${PY} ${T}/check_gen_namelist.py ${CMAKE_SOURCE_DIR}/tools/gen_namelist_table.py
  ${CPPMPAS_MPAS_REGISTRY} ${REF_MW}/steps/namelist.atmosphere ${W}/gen_namelist)
# --scalars must be the scalars the namelist's packages activate
cppmpas_error_test(scalars_extra "--scalars must be qv: the scalars the namelist activates" init ${INIT_MW} ${W}/x.nc
  --namelist ${REF_MW}/steps/namelist.atmosphere --scalars qv,qc,qr --fields dss)
cppmpas_error_test(scalars_missing "--scalars must be qv,qc,qr: the scalars the namelist activates" init ${INIT_SC}
  ${W}/x.nc --namelist ${REF_SC}/steps/namelist.atmosphere --fields dss)

# ---------------------------------------------------------------------------------------------------------------
# Every field atm_mpas_init_block computes, at the initial time, against the initdump stream of the references
# mpas842-dbl-mw-003, -sc-002 and -jw-002 (the short run writes initdump.nc, output_interval="initial_only", at the
# first stream write of atm_core_run, mpas_atm_core.F:701, after the init block; field list INITDUMP_FIELDS of
# cppmpas-data ci/cases.py). One test per case and field, tolerance 0 (tests/compare_arrays.py: type, dimensions,
# attributes, raw bytes).
set(INITDUMP_FIELDS
  xtime
  invAreaCell invDvEdge invDcEdge invAreaTriangle edgesOnCell_sign edgesOnVertex_sign
  zb_cell zb3_cell kiteForCell advCellsForEdge nAdvCellsForEdge adv_coefs adv_coefs_3rd
  deformation_coef_c2 deformation_coef_s2 deformation_coef_cs deformation_coef_c deformation_coef_s
  meshScalingDel2 meshScalingDel4 meshScalingRegionalCell meshScalingRegionalEdge
  etp etm ewp ewm specZoneMaskCell specZoneMaskEdge specZoneMaskVertex nearestRelaxationCell
  theta_m rho_zz ru rw rho_p rtheta_base rtheta_p exner exner_base pressure_p pressure_base
  v rho_edge vorticity divergence ke pv_edge pv_vertex pv_cell gradPVt gradPVn
  uReconstructX uReconstructY uReconstructZ uReconstructZonal uReconstructMeridional)
string(REPLACE ";" "," INITDUMP_LIST "${INITDUMP_FIELDS}")
set(DUMPREF_MW ${CPPMPAS_REFS_DIR}/mpas842-dbl-mw-003/MW)
set(DUMPREF_SC ${CPPMPAS_REFS_DIR}/mpas842-dbl-sc-002/SC)
set(DUMPREF_JW ${CPPMPAS_REFS_DIR}/mpas842-dbl-jw-002/JW)
set(DUMPINIT_MW ${DUMPREF_MW}/init/mtn_wave_init.nc)
set(DUMPINIT_SC ${DUMPREF_SC}/init/supercell_init.nc)
set(DUMPINIT_JW ${DUMPREF_JW}/init/x1.40962.init.nc)
foreach(case MW SC JW)
  add_test(NAME initdump_run_${case} COMMAND ${CMAKE_COMMAND} -E env sh -c
    "mkdir -p '${W}' && rm -f '${W}/initdump_${case}.nc' && '${EXE}' init '${DUMPINIT_${case}}' '${W}/initdump_${case}.nc' --namelist '${DUMPREF_${case}}/steps/namelist.atmosphere' --scalars ${SCAL_${case}} --fields ${INITDUMP_LIST}")
  set_tests_properties(initdump_run_${case} PROPERTIES FIXTURES_SETUP initdump_${case})
  add_test(NAME initdump_vars_${case} COMMAND ${PY} ${T}/initdump_vars.py ${DUMPREF_${case}}/steps/initdump.nc
    ${INITDUMP_LIST})
  foreach(f IN LISTS INITDUMP_FIELDS)
    add_test(NAME initdump_${case}_${f} COMMAND ${PY} ${T}/compare_arrays.py ${DUMPREF_${case}}/steps/initdump.nc
      ${W}/initdump_${case}.nc ${f})
    set_tests_properties(initdump_${case}_${f} PROPERTIES FIXTURES_REQUIRED initdump_${case})
  endforeach()
endforeach()

# atm_compute_solve_diagnostics and mpas_reconstruct_2d as the dycore calls them (on arrays that hold values,
# rk_step 1, 2, 3), against the init state of the initdump references and an independent implementation
add_executable(test_solve_diag test_solve_diag.cpp)
target_link_libraries(test_solve_diag PRIVATE cppmpas m)
add_test(NAME unit_solve_diag COMMAND test_solve_diag
  ${DUMPINIT_MW} ${DUMPREF_MW}/steps/namelist.atmosphere ${SCAL_MW}
  ${DUMPINIT_SC} ${DUMPREF_SC}/steps/namelist.atmosphere ${SCAL_SC}
  ${DUMPINIT_JW} ${DUMPREF_JW}/steps/namelist.atmosphere ${SCAL_JW})

# tests/mutation_check.py's reading of CTest's output: a test that was not run does not kill a mutant
add_test(NAME mutation_check_unit COMMAND ${PY} ${T}/test_mutation_check.py)

# ---------------------------------------------------------------------------------------------------------------
# The dycore (src/core_atmosphere/atm_dycore.cpp, atm_dissipation.cpp): ten 6-s steps of the mountain wave from the
# reference init file with the reference namelist, against the every-step "bisection" stream of
# mpas842-dbl-mw-001 (step.*.nc, written by MPAS at t = 0 and after every step). One test per step and field,
# tolerance 0 (tests/compare_arrays.py: type, dimensions, attributes, raw bytes), and the whole directory through
# tools/mpas_compare.py.
set(DYCORE_MW_FIELDS xtime Time u w rho_zz theta_m rho theta pressure_p qv)
set(DYCORE_MW_STAMPS 00.00.00 00.00.06 00.00.12 00.00.18 00.00.24 00.00.30 00.00.36 00.00.42 00.00.48 00.00.54 00.01.00)
add_test(NAME dycore_run_MW COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/dycore_MW' && mkdir -p '${W}/dycore_MW' && '${EXE}' run '${INIT_MW}' '${W}/dycore_MW' --namelist '${REF_MW}/steps/namelist.atmosphere' --steps 10")
set_tests_properties(dycore_run_MW PROPERTIES FIXTURES_SETUP dycore_MW)
set(step 0)
foreach(stamp IN LISTS DYCORE_MW_STAMPS)
  foreach(f IN LISTS DYCORE_MW_FIELDS)
    add_test(NAME dycore_MW_step${step}_${f} COMMAND ${PY} ${T}/compare_arrays.py
      ${REF_MW}/steps/step.0000-01-01_${stamp}.nc ${W}/dycore_MW/step.0000-01-01_${stamp}.nc ${f})
    set_tests_properties(dycore_MW_step${step}_${f} PROPERTIES FIXTURES_REQUIRED dycore_MW)
  endforeach()
  math(EXPR step "${step} + 1")
endforeach()
# tools/mpas_compare.py on the two directories: every variable of every step identical; the only difference allowed
# is the global attributes (the namelist options MPAS writes as attributes, not written yet, mpas_writer.hpp)
add_test(NAME dycore_MW_mpas_compare COMMAND ${PY} ${T}/dycore_compare.py ${CPPMPAS_MPAS_COMPARE} ${REF_MW}/steps
  ${W}/dycore_MW ${W}/dycore_MW_compare.json 110)
set_tests_properties(dycore_MW_mpas_compare PROPERTIES FIXTURES_REQUIRED dycore_MW)
# the configuration checks and the clock arithmetic
add_executable(test_dycore test_dycore.cpp)
target_link_libraries(test_dycore PRIVATE cppmpas)
add_test(NAME unit_dycore COMMAND test_dycore ${INIT_MW} ${REF_MW}/steps/namelist.atmosphere ${W})
set_tests_properties(unit_dycore PROPERTIES FIXTURES_REQUIRED work_dir TIMEOUT 300)
# configurations of the other references, and meshes, that the dycore refuses before it writes anything
# the supercell (src/core_atmosphere/atm_advance_scalars.cpp, atm_dissipation_perturbation.cpp, atm_microphysics.cpp,
# src/physics/mp_kessler.hpp): ten 3-s steps with scalar transport, perturbation mixing and Kessler microphysics
# against the reference's every-step history (mpas842-dbl-sc-001 steps/), every field at tolerance 0
set(DYCORE_SC_FIELDS xtime Time u w rho_zz theta_m rho theta pressure_p qv qc qr)
set(DYCORE_SC_STAMPS 00.00.00 00.00.03 00.00.06 00.00.09 00.00.12 00.00.15 00.00.18 00.00.21 00.00.24 00.00.27 00.00.30)
add_test(NAME dycore_run_SC COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/dycore_SC' && mkdir -p '${W}/dycore_SC' && '${EXE}' run '${INIT_SC}' '${W}/dycore_SC' --namelist '${REF_SC}/steps/namelist.atmosphere' --scalars qv,qc,qr --steps 10")
set_tests_properties(dycore_run_SC PROPERTIES FIXTURES_SETUP dycore_SC TIMEOUT 900)
set(step 0)
foreach(stamp IN LISTS DYCORE_SC_STAMPS)
  foreach(f IN LISTS DYCORE_SC_FIELDS)
    add_test(NAME dycore_SC_step${step}_${f} COMMAND ${PY} ${T}/compare_arrays.py
      ${REF_SC}/steps/step.0000-01-01_${stamp}.nc ${W}/dycore_SC/step.0000-01-01_${stamp}.nc ${f})
    set_tests_properties(dycore_SC_step${step}_${f} PROPERTIES FIXTURES_REQUIRED dycore_SC)
  endforeach()
  math(EXPR step "${step} + 1")
endforeach()
add_executable(test_microphysics test_microphysics.cpp)
target_link_libraries(test_microphysics PRIVATE cppmpas)
add_test(NAME unit_microphysics COMMAND test_microphysics ${INIT_SC} ${REF_SC}/steps/namelist.atmosphere
  ${REF_MW}/steps/namelist.atmosphere)
add_test(NAME dycore_SC_mpas_compare COMMAND ${PY} ${T}/dycore_compare.py ${CPPMPAS_MPAS_COMPARE} ${REF_SC}/steps
  ${W}/dycore_SC ${W}/dycore_SC_compare.json 132)
set_tests_properties(dycore_SC_mpas_compare PROPERTIES FIXTURES_REQUIRED dycore_SC)
cppmpas_error_test(dycore_odd_levels "nVertLevels = 3 (< 4)" run ${W}/sc_init.nc
  ${W}/x_dir --namelist ${REF_MW}/steps/namelist.atmosphere --steps 1)
cppmpas_error_test(dycore_steps "--steps needs a whole number 0..1000000, not '1x'" run ${INIT_MW} ${W}/x_dir
  --namelist ${REF_MW}/steps/namelist.atmosphere --steps 1x)
cppmpas_error_test(dycore_no_steps "usage: cppmpas run" run ${INIT_MW} ${W}/x_dir
  --namelist ${REF_MW}/steps/namelist.atmosphere)
# 600 steps (label "long", about 5 minutes): the forecast run of the same reference writes output.nc at 30 and 60 min
# (records 1 and 2, after steps 300 and 600) and restart.0000-01-01_00.30.06.nc (after step 301); every field of these
# streams the port holds is compared at tolerance 0
set(DYCORE_RESTART_FIELDS xtime Time u w rho_zz theta_m pressure_p rho theta exner ru rw rtheta_p rho_p qv)
set(DYCORE_OUTPUT_FIELDS xtime u w pressure rho theta divergence vorticity ke uReconstructZonal uReconstructMeridional qv)
set(DYCORE_LONG_FIELDS ${DYCORE_RESTART_FIELDS} pressure divergence vorticity ke uReconstructZonal uReconstructMeridional)
list(REMOVE_DUPLICATES DYCORE_LONG_FIELDS)
string(REPLACE ";" "," DYCORE_LONG_LIST "${DYCORE_LONG_FIELDS}")
add_test(NAME dycore_run_long_MW COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/dycore_long_MW' && mkdir -p '${W}/dycore_long_MW' && '${EXE}' run '${INIT_MW}' '${W}/dycore_long_MW' --namelist '${REF_MW}/steps/namelist.atmosphere' --steps 600 --every 0 --write-at 300,301,600 --write-fields ${DYCORE_LONG_LIST}")
set_tests_properties(dycore_run_long_MW PROPERTIES FIXTURES_SETUP dycore_long_MW LABELS long TIMEOUT 900)
foreach(f IN LISTS DYCORE_RESTART_FIELDS)
  add_test(NAME dycore_MW_restart_${f} COMMAND ${PY} ${T}/compare_arrays.py ${RST_MW}
    ${W}/dycore_long_MW/state.0000-01-01_00.30.06.nc ${f})
  set_tests_properties(dycore_MW_restart_${f} PROPERTIES FIXTURES_REQUIRED dycore_long_MW LABELS long)
endforeach()
foreach(f IN LISTS DYCORE_OUTPUT_FIELDS)
  add_test(NAME dycore_MW_output30_${f} COMMAND ${PY} ${T}/compare_arrays.py ${REF_MW}/forecast/output.nc
    ${W}/dycore_long_MW/state.0000-01-01_00.30.00.nc ${f} --record 1)
  set_tests_properties(dycore_MW_output30_${f} PROPERTIES FIXTURES_REQUIRED dycore_long_MW LABELS long)
  add_test(NAME dycore_MW_output60_${f} COMMAND ${PY} ${T}/compare_arrays.py ${REF_MW}/forecast/output.nc
    ${W}/dycore_long_MW/state.0000-01-01_01.00.00.nc ${f} --record 2)
  set_tests_properties(dycore_MW_output60_${f} PROPERTIES FIXTURES_REQUIRED dycore_long_MW LABELS long)
endforeach()
# The JW baroclinic wave (x1.40962 sphere, 26 levels, dt = 450 s): ten steps against the every-step stream of
# mpas842-dbl-jw-001, one test per step and field at tolerance 0, and the directory through tools/mpas_compare.py.
# JW takes 2d_smagorinsky with the 4th-order filters, no vertical mixing, fEdge != 0 with u_init = v_init = 0, the
# sphere path of mpas_reconstruct and qv (all zero) as its only scalar, not advected.
set(DYCORE_JW_STAMPS 00.00.00 00.07.30 00.15.00 00.22.30 00.30.00 00.37.30 00.45.00 00.52.30 01.00.00 01.07.30 01.15.00)
add_test(NAME dycore_run_JW COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/dycore_JW' && mkdir -p '${W}/dycore_JW' && '${EXE}' run '${INIT_JW}' '${W}/dycore_JW' --namelist '${REF_JW}/steps/namelist.atmosphere' --steps 10")
set_tests_properties(dycore_run_JW PROPERTIES FIXTURES_SETUP dycore_JW TIMEOUT 900)
set(step 0)
foreach(stamp IN LISTS DYCORE_JW_STAMPS)
  foreach(f IN LISTS DYCORE_MW_FIELDS)
    add_test(NAME dycore_JW_step${step}_${f} COMMAND ${PY} ${T}/compare_arrays.py
      ${REF_JW}/steps/step.0000-01-01_${stamp}.nc ${W}/dycore_JW/step.0000-01-01_${stamp}.nc ${f})
    set_tests_properties(dycore_JW_step${step}_${f} PROPERTIES FIXTURES_REQUIRED dycore_JW)
  endforeach()
  math(EXPR step "${step} + 1")
endforeach()
add_test(NAME dycore_JW_mpas_compare COMMAND ${PY} ${T}/dycore_compare.py ${CPPMPAS_MPAS_COMPARE} ${REF_JW}/steps
  ${W}/dycore_JW ${W}/dycore_JW_compare.json 110)
set_tests_properties(dycore_JW_mpas_compare PROPERTIES FIXTURES_REQUIRED dycore_JW)
# 193 steps (label "long", about 10 minutes): output.nc record 1 (24 h, after step 192) and
# restart.0000-01-02_00.07.30.nc (after step 193) of the JW forecast run
add_test(NAME dycore_run_long_JW COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/dycore_long_JW' && mkdir -p '${W}/dycore_long_JW' && '${EXE}' run '${INIT_JW}' '${W}/dycore_long_JW' --namelist '${REF_JW}/steps/namelist.atmosphere' --steps 193 --every 0 --write-at 192,193 --write-fields ${DYCORE_LONG_LIST}")
set_tests_properties(dycore_run_long_JW PROPERTIES FIXTURES_SETUP dycore_long_JW LABELS long TIMEOUT 3600)
foreach(f IN LISTS DYCORE_RESTART_FIELDS)
  add_test(NAME dycore_JW_restart_${f} COMMAND ${PY} ${T}/compare_arrays.py ${RST_JW}
    ${W}/dycore_long_JW/state.0000-01-02_00.07.30.nc ${f})
  set_tests_properties(dycore_JW_restart_${f} PROPERTIES FIXTURES_REQUIRED dycore_long_JW LABELS long)
endforeach()
foreach(f IN LISTS DYCORE_OUTPUT_FIELDS)
  add_test(NAME dycore_JW_output24h_${f} COMMAND ${PY} ${T}/compare_arrays.py ${REF_JW}/forecast/output.nc
    ${W}/dycore_long_JW/state.0000-01-02_00.00.00.nc ${f} --record 1)
  set_tests_properties(dycore_JW_output24h_${f} PROPERTIES FIXTURES_REQUIRED dycore_long_JW LABELS long)
endforeach()
# SC, 601 steps (label "long", about 45 minutes): the forecast run writes output.nc at 30 min (record 1, after step
# 600) and restart.0000-01-01_00.30.03.nc (after step 601); surface_pressure is microphysics_to_MPAS's (PI:950-964)
# output.nc from mpas842-dbl-sc-002 (the same run of the same recipe; its fields are identical to sc-001's at
# tolerance 0, and the workstation copy of sc-001 holds no output.nc)
set(OUT_SC ${CPPMPAS_REFS_DIR}/mpas842-dbl-sc-002/SC/forecast/output.nc)
# the recipes the tests mix are one: sc-001 = sc-002 (init, namelists), and in every reference the steps namelist the
# runs use differs from the forecast namelist of the compared files in config_run_duration only
add_test(NAME ref_recipes COMMAND ${PY} ${T}/check_ref_recipes.py ${CPPMPAS_REFS_DIR})
set(DYCORE_SC_RESTART_FIELDS xtime Time u w rho_zz theta_m pressure_p rho theta exner ru rw rtheta_p rho_p qv qc qr
  surface_pressure)
set(DYCORE_SC_OUTPUT_FIELDS xtime u w rho theta divergence vorticity ke uReconstructZonal uReconstructMeridional
  surface_pressure qv qc qr)
set(DYCORE_SC_LONG_FIELDS ${DYCORE_SC_RESTART_FIELDS} ${DYCORE_SC_OUTPUT_FIELDS})
list(REMOVE_DUPLICATES DYCORE_SC_LONG_FIELDS)
string(REPLACE ";" "," DYCORE_SC_LONG_LIST "${DYCORE_SC_LONG_FIELDS}")
add_test(NAME dycore_run_long_SC COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/dycore_long_SC' && mkdir -p '${W}/dycore_long_SC' && '${EXE}' run '${INIT_SC}' '${W}/dycore_long_SC' --namelist '${REF_SC}/steps/namelist.atmosphere' --scalars qv,qc,qr --steps 601 --every 0 --write-at 600,601 --write-fields ${DYCORE_SC_LONG_LIST}")
set_tests_properties(dycore_run_long_SC PROPERTIES FIXTURES_SETUP dycore_long_SC LABELS long TIMEOUT 7200)
foreach(f IN LISTS DYCORE_SC_RESTART_FIELDS)
  add_test(NAME dycore_SC_restart_${f} COMMAND ${PY} ${T}/compare_arrays.py ${RST_SC}
    ${W}/dycore_long_SC/state.0000-01-01_00.30.03.nc ${f})
  set_tests_properties(dycore_SC_restart_${f} PROPERTIES FIXTURES_REQUIRED dycore_long_SC LABELS long)
endforeach()
foreach(f IN LISTS DYCORE_SC_OUTPUT_FIELDS)
  add_test(NAME dycore_SC_output30_${f} COMMAND ${PY} ${T}/compare_arrays.py ${OUT_SC}
    ${W}/dycore_long_SC/state.0000-01-01_00.30.00.nc ${f} --record 1)
  set_tests_properties(dycore_SC_output30_${f} PROPERTIES FIXTURES_REQUIRED dycore_long_SC LABELS long)
endforeach()
# A restart run (config_do_restart = true; tests/namelists/MW_restart is the MW forecast namelist started at 00:30:06)
# from the reference's restart.0000-01-01_00.30.06.nc. At step 0 the state written is what the restart stream holds
# (no atm_init_coupled_diagnostics, mpas_atm_core.F:522); after 299 steps (01:00:00) every field of output.nc record 2
# the port holds is compared at tolerance 0 (label "long", about 3 minutes).
set(NML_RST_MW ${CMAKE_SOURCE_DIR}/tests/namelists/MW_restart/namelist.atmosphere)
set(DYCORE_RESTART0_FIELDS ${DYCORE_RESTART_FIELDS} rtheta_base exner_base pressure_base)
string(REPLACE ";" "," DYCORE_RESTART0_LIST "${DYCORE_RESTART0_FIELDS}")
add_test(NAME dycore_restart0_run_MW COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/dycore_restart0_MW' && mkdir -p '${W}/dycore_restart0_MW' && '${EXE}' run '${RST_MW}' '${W}/dycore_restart0_MW' --namelist '${NML_RST_MW}' --steps 0 --every 0 --write-at 0 --write-fields ${DYCORE_RESTART0_LIST}")
set_tests_properties(dycore_restart0_run_MW PROPERTIES FIXTURES_SETUP dycore_restart0_MW)
foreach(f IN LISTS DYCORE_RESTART0_FIELDS)
  add_test(NAME dycore_restart0_MW_${f} COMMAND ${PY} ${T}/compare_arrays.py ${RST_MW}
    ${W}/dycore_restart0_MW/state.0000-01-01_00.30.06.nc ${f})
  set_tests_properties(dycore_restart0_MW_${f} PROPERTIES FIXTURES_REQUIRED dycore_restart0_MW)
endforeach()
string(REPLACE ";" "," DYCORE_OUTPUT_LIST "${DYCORE_OUTPUT_FIELDS}")
add_test(NAME dycore_restart_run_long_MW COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/dycore_restart_long_MW' && mkdir -p '${W}/dycore_restart_long_MW' && '${EXE}' run '${RST_MW}' '${W}/dycore_restart_long_MW' --namelist '${NML_RST_MW}' --steps 299 --every 0 --write-at 299 --write-fields ${DYCORE_OUTPUT_LIST}")
set_tests_properties(dycore_restart_run_long_MW PROPERTIES FIXTURES_SETUP dycore_restart_long_MW LABELS long TIMEOUT 900)
foreach(f IN LISTS DYCORE_OUTPUT_FIELDS)
  add_test(NAME dycore_restart_MW_output60_${f} COMMAND ${PY} ${T}/compare_arrays.py ${REF_MW}/forecast/output.nc
    ${W}/dycore_restart_long_MW/state.0000-01-01_01.00.00.nc ${f} --record 2)
  set_tests_properties(dycore_restart_MW_output60_${f} PROPERTIES FIXTURES_REQUIRED dycore_restart_long_MW LABELS long)
endforeach()
# the restart reader's refusals (missing restart-stream fields, a field of an active package) and the start time
add_test(NAME restart_refusals_MW COMMAND ${PY} ${T}/check_restart.py ${EXE} ${RST_MW} ${NML_RST_MW} ${INIT_MW}
  ${W}/restart_refusals)
set_tests_properties(restart_refusals_MW PROPERTIES FIXTURES_REQUIRED work_dir)
cppmpas_error_test(restart_time "does not contain the time 0000-01-01_00:30:06" run ${INIT_MW} ${W}/x_dir
  --namelist ${NML_RST_MW} --steps 1)
cppmpas_error_test(dycore_write_at "--write-at 12 is after the last step" run ${INIT_MW} ${W}/x_dir
  --namelist ${REF_MW}/steps/namelist.atmosphere --steps 10 --write-at 12 --write-fields u)
cppmpas_error_test(dycore_write_fields "variable 'relhum' has no CPPMPAS port" run ${INIT_MW} ${W}/x_dir
  --namelist ${REF_MW}/steps/namelist.atmosphere --steps 10 --write-at 1 --write-fields u,relhum)

# WSM6 and its MPAS driver / interface against MPAS's own code built with the reference compiler and flags
# (tests/wsm6_column/run_check.sh). Needs gfortran 14.2, the MPAS v8.4.2 source tree and MMM-physics a4baf7f
# (cppmpas-data `inputs` release, mmm-physics-a4baf7f.tar.gz, extracted); without them the test fails.
set(CPPMPAS_MPAS_SRC "" CACHE PATH "MPAS-Model v8.4.2 source tree")
set(CPPMPAS_MMM_PHYSICS "" CACHE PATH "extracted MMM-physics a4baf7f (mmm-physics-a4baf7f.tar.gz)")
add_executable(test_wsm6_column wsm6_column/test_wsm6_column.cpp)
target_link_libraries(test_wsm6_column PRIVATE cppmpas)
add_test(NAME wsm6_column COMMAND ${T}/wsm6_column/run_check.sh ${CPPMPAS_MPAS_SRC} ${CPPMPAS_MMM_PHYSICS}
  ${W}/wsm6_column $<TARGET_FILE:test_wsm6_column>)
set_tests_properties(wsm6_column PROPERTIES FIXTURES_REQUIRED work_dir)

# Thompson microphysics (thompson_init with MPAS's .DBL tables, the scheme, its MPAS driver / interface and
# init_thompson_clouddroplets_forMPAS) against MPAS's own code built with the reference compiler and flags
# (tests/thompson_column/run_check.sh). Needs gfortran 14.2, the MPAS v8.4.2 source tree, MMM-physics a4baf7f and the
# Thompson tables of the convection_permitting reference (mpas842-dbl-rgcp-001 RGCP/tables); without them it fails.
# The four MP_THOMPSON_*_DATA.DBL tables are the reference run's own: the tables/ of release mpas842-dbl-rgcp-001
# (also its separate assets), under CPPMPAS_REFS_DIR; run_check.sh checks them against tests/thompson_column/tables.sha256
# (the release's mpas842-dbl-rgcp-001.tables.sha256).
set(CPPMPAS_THOMPSON_TABLES ${CPPMPAS_REFS_DIR}/mpas842-dbl-rgcp-001/RGCP/tables)
add_executable(test_thompson_column thompson_column/test_thompson_column.cpp)
target_link_libraries(test_thompson_column PRIVATE cppmpas)
add_test(NAME thompson_column COMMAND ${T}/thompson_column/run_check.sh ${CPPMPAS_MPAS_SRC} ${CPPMPAS_MMM_PHYSICS}
  ${CPPMPAS_THOMPSON_TABLES} ${W}/thompson_column $<TARGET_FILE:test_thompson_column>)
set_tests_properties(thompson_column PROPERTIES FIXTURES_REQUIRED work_dir)

# The revised MM5 surface layer, YSU and their MPAS drivers, MPAS_to_physics and physics_get_tend against MPAS's own
# code built with the reference compiler and flags (tests/pbl_column/run_check.sh); same requirements as wsm6_column.
add_executable(test_pbl_column pbl_column/test_pbl_column.cpp)
target_link_libraries(test_pbl_column PRIVATE cppmpas)
add_test(NAME pbl_column COMMAND ${T}/pbl_column/run_check.sh ${CPPMPAS_MPAS_SRC} ${CPPMPAS_MMM_PHYSICS}
  ${W}/pbl_column $<TARGET_FILE:test_pbl_column>)
set_tests_properties(pbl_column PROPERTIES FIXTURES_REQUIRED work_dir)

# The cloudiness, RRTMG LW / SW radiation drivers, the RRTMG table initialisation and the ozone climatology against
# MPAS's own code built with the reference compiler and flags (tests/radiation_column/run_check.sh). Needs gfortran
# 14.2, the MPAS v8.4.2 source tree and the reference build's tables directory (RRTMG_*_DATA.DBL, OZONE_*.TBL of the
# cppmpas-data mpas-build asset) and CPPMPAS_CPPWRF_RRTMG (CPPWRF's data/rrtmg: the tables MPAS's initialisation
# computes must be those, byte for byte). Without either the test fails.
set(CPPMPAS_MPAS_TABLES "" CACHE PATH "the MPAS reference build's tables directory (RRTMG_*_DATA.DBL, OZONE_*.TBL)")
set(CPPMPAS_CPPWRF_RRTMG "" CACHE PATH "CPPWRF's data/rrtmg (rrtmg_lw_tables.bin, rrtmg_sw_tables.bin)")
add_executable(test_radiation_column radiation_column/test_radiation_column.cpp)
target_link_libraries(test_radiation_column PRIVATE cppmpas)
add_test(NAME radiation_column COMMAND ${CMAKE_COMMAND} -E env CPPMPAS_CPPWRF_RRTMG=${CPPMPAS_CPPWRF_RRTMG}
  ${T}/radiation_column/run_check.sh ${CPPMPAS_MPAS_SRC} ${CPPMPAS_MPAS_TABLES} ${W}/radiation_column
  $<TARGET_FILE:test_radiation_column>)
set_tests_properties(radiation_column PROPERTIES FIXTURES_REQUIRED work_dir)
# The Noah land surface (driver_lsm, driver_seaice, init_lsm and the tables, atmphys_sfc_diagnostics) against MPAS's
# own code built with the reference compiler and flags (tests/noah_column/run_check.sh). Needs gfortran 14.2, the
# MPAS v8.4.2 source tree and the reference build's tables (VEGPARM.TBL, SOILPARM.TBL, GENPARM.TBL: the mpas-build
# asset's tables/); without them the test fails.
set(CPPMPAS_NOAH_TABLES "" CACHE PATH "directory with the reference build's VEGPARM.TBL, SOILPARM.TBL, GENPARM.TBL")
add_executable(test_noah_column noah_column/test_noah_column.cpp)
target_link_libraries(test_noah_column PRIVATE cppmpas)
add_test(NAME noah_column COMMAND ${T}/noah_column/run_check.sh ${CPPMPAS_MPAS_SRC} ${CPPMPAS_NOAH_TABLES}
  ${W}/noah_column $<TARGET_FILE:test_noah_column>)
set_tests_properties(noah_column PROPERTIES FIXTURES_REQUIRED work_dir)
# GWDO (bl_ysu_gwdo), New Tiedtke (cu_ntiedtke) and their MPAS drivers, with MPAS_to_physics and physics_get_tend,
# against MPAS's own code built with the reference compiler and flags (tests/cu_gwdo_column/run_check.sh); same
# requirements as wsm6_column.
add_executable(test_cu_gwdo_column cu_gwdo_column/test_cu_gwdo_column.cpp)
target_link_libraries(test_cu_gwdo_column PRIVATE cppmpas)
add_test(NAME cu_gwdo_column COMMAND ${T}/cu_gwdo_column/run_check.sh ${CPPMPAS_MPAS_SRC} ${CPPMPAS_MMM_PHYSICS}
  ${W}/cu_gwdo_column $<TARGET_FILE:test_cu_gwdo_column>)
set_tests_properties(cu_gwdo_column PROPERTIES FIXTURES_REQUIRED work_dir)
# Grell-Freitas (cu_grell_freitas) and its MPAS driver, with MPAS_to_physics and physics_get_tend, against MPAS's own
# code built with the reference compiler and flags (tests/gf_column/run_check.sh); same requirements as wsm6_column.
add_executable(test_gf_column gf_column/test_gf_column.cpp)
target_link_libraries(test_gf_column PRIVATE cppmpas)
add_test(NAME gf_column COMMAND ${T}/gf_column/run_check.sh ${CPPMPAS_MPAS_SRC} ${CPPMPAS_MMM_PHYSICS}
  ${W}/gf_column $<TARGET_FILE:test_gf_column>)
set_tests_properties(gf_column PROPERTIES FIXTURES_REQUIRED work_dir)
# The MYNN surface layer (sf_mynn: init_sfclayer, driver_sfclayer with sfclay_mynn, sf_mynn_pre_run and sf_mynn_run)
# with MPAS_to_physics against MPAS's own code built with the reference compiler and flags
# (tests/sfmynn_column/run_check.sh); same requirements as wsm6_column (an empty path fails the test).
add_executable(test_sfmynn_column sfmynn_column/test_sfmynn_column.cpp)
target_link_libraries(test_sfmynn_column PRIVATE cppmpas)
add_test(NAME sfmynn_column COMMAND ${T}/sfmynn_column/run_check.sh "${CPPMPAS_MPAS_SRC}" "${CPPMPAS_MMM_PHYSICS}"
  ${W}/sfmynn_column $<TARGET_FILE:test_sfmynn_column>)
set_tests_properties(sfmynn_column PROPERTIES FIXTURES_REQUIRED work_dir)

# ---------------------------------------------------------------------------------------------------------------
# The real-data case RG (mpas842-dbl-rg-004: the x1.10242 mesh, 55 levels, dt = 600 s, the mesoscale_reference physics
# suite with RRTMG every 30 min, split dynamics-transport with three dynamics substeps, WSM6 scalars). One run of 12
# steps writes the reference's every-step stream (steps 0-10 at tolerance 0, one test per step and field, and the
# directory through tools/mpas_compare.py) and, after steps 6, 7 and 12, the fields of the forecast run's history
# (1 h, 2 h) and restart (01:10) streams. The steps run and the forecast run of the reference hold the same state (the
# forecast run's output alarms only add diagnostics); relhum is left out: MPAS recomputes it in the WSM6 driver when
# the diagnostics alarm rings (compute_relhum, l_diags), which a steps-configuration run never does (the forecast
# configuration is the --streams test below). Needs CPPMPAS_MPAS_TABLES (the reference build's tables) and
# CPPMPAS_CPPWRF_RRTMG (the RRTMG table images).
set(REF_RG ${CPPMPAS_REFS_DIR}/mpas842-dbl-rg-004/RG)
set(DYCORE_RG_FIELDS xtime Time u w rho_zz theta_m rho theta pressure_p qv qc qr qi qs qg)
set(DYCORE_RG_RESTART_FIELDS ru_p rw_p rt_diabatic_tend rthdynten isltyp ivgtyp mminlu isice_lu iswater_lu landmask shdmin shdmax snoalb albedo12m greenfrac var2d con oa1 oa2 oa3 oa4 ol1 ol2 ol3 ol4 qv qc qr qi qs qg u_init v_init t_init xtime Time u w rho_zz theta_m pressure_p rho theta exner exner_base rtheta_base pressure_base rho_base theta_base ru rw rtheta_p rho_p surface_pressure rainprod evapprod re_cloud re_ice re_snow tend_sfc_pressure nsteps_accum ndays_accum tlag tday_accum tyear_mean tyear_accum i_rainnc rainncv rainnc snowncv snownc graupelncv graupelnc sr i_rainc cuprec rainc raincv qcg wstar delta ch kpbl hpbl hfx mavail mol qfx qsfc ust ustm zol znt br cd cda chs chs2 cqs2 ck cka cpm flhc flqc gz1oz0 lh psim psih qgh regime rmol wspd fh fm u10 v10 q2 t2m th2m dusfcg dvsfcg dtaux3d dtauy3d rubldiff rvbldiff i_acswdnb i_acswdnbc i_acswdnt i_acswdntc i_acswupb i_acswupbc i_acswupt i_acswuptc acswdnb acswdnbc acswdnt acswdntc acswupb acswupbc acswupt acswuptc gsw i_aclwdnb i_aclwdnbc i_aclwdnt i_aclwdntc i_aclwupb i_aclwupbc i_aclwupt i_aclwuptc aclwdnb aclwdnbc aclwdnt aclwdntc aclwupb aclwupbc aclwupt aclwuptc glw o3clim o3vmr acsnom acsnow canwat chklowq grdflx lai noahres potevp sfc_albedo sfc_emiss sfc_emibck sfcrunoff smstav smstot snopcx snotime sstsk sstsk_dtc sstsk_dtw thc udrunoff xicem z0 zs rthmpten rqvmpten rqcmpten rqrmpten rqimpten rqsmpten rqgmpten rthcuten rqvcuten rqccuten rqicuten rqvdynten rucuten rvcuten rublten rvblten rthblten rqvblten rqcblten rqiblten rthratensw rthratenlw sfc_albbck skintemp snow snowc snowh sst tmn vegfra seaice xice xland dzs smcrel sh2o smois tslb t_oml t_oml_initial t_oml_200m_initial h_oml h_oml_initial hu_oml hv_oml)
set(DYCORE_RG_HISTORY_FIELDS qv qc qr qi qs qg xtime Time u w pressure surface_pressure rho theta divergence vorticity ke uReconstructZonal uReconstructMeridional i_rainnc rainnc precipw cuprec i_rainc rainc kpbl hpbl hfx qfx cd cda ck cka lh u10 v10 q2 t2m th2m gsw glw acsnow pin ozmixm mminlu isice_lu iswater_lu xland skintemp snow snowh sst vegfra xice sh2o smois tslb)
set(DYCORE_RG_LONG_FIELDS ${DYCORE_RG_RESTART_FIELDS} ${DYCORE_RG_HISTORY_FIELDS})
list(REMOVE_DUPLICATES DYCORE_RG_LONG_FIELDS)
string(REPLACE ";" "," DYCORE_RG_LONG_LIST "${DYCORE_RG_LONG_FIELDS}")
add_test(NAME dycore_run_RG COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/dycore_RG' && mkdir -p '${W}/dycore_RG' && '${EXE}' run '${REF_RG}/init/x1.10242.init.nc' '${W}/dycore_RG' --namelist '${REF_RG}/steps/namelist.atmosphere' --scalars qv,qc,qr,qi,qs,qg --steps 12 --write-at 6,7,12 --write-fields ${DYCORE_RG_LONG_LIST} --tables '${CPPMPAS_MPAS_TABLES}' --rrtmg '${CPPMPAS_CPPWRF_RRTMG}'")
set_tests_properties(dycore_run_RG PROPERTIES FIXTURES_SETUP dycore_RG TIMEOUT 1800)
set(DYCORE_RG_STAMPS 00.00.00 00.10.00 00.20.00 00.30.00 00.40.00 00.50.00 01.00.00 01.10.00 01.20.00 01.30.00 01.40.00)
set(step 0)
foreach(stamp IN LISTS DYCORE_RG_STAMPS)
  foreach(f IN LISTS DYCORE_RG_FIELDS)
    add_test(NAME dycore_RG_step${step}_${f} COMMAND ${PY} ${T}/compare_arrays.py
      ${REF_RG}/steps/step.2026-10-12_${stamp}.nc ${W}/dycore_RG/step.2026-10-12_${stamp}.nc ${f})
    set_tests_properties(dycore_RG_step${step}_${f} PROPERTIES FIXTURES_REQUIRED dycore_RG)
  endforeach()
  math(EXPR step "${step} + 1")
endforeach()
# the reference's 11 step files only (the run also writes steps 11, 12 and the state files)
string(REPLACE ";" " " DYCORE_RG_STAMP_WORDS "${DYCORE_RG_STAMPS}")
add_test(NAME dycore_RG_mpas_compare COMMAND ${CMAKE_COMMAND} -E env sh -c
  "set -e; rm -rf '${W}/dycore_RG_steps' && mkdir -p '${W}/dycore_RG_steps' && for s in ${DYCORE_RG_STAMP_WORDS}; do ln -s '${W}/dycore_RG/step.2026-10-12_'$s.nc '${W}/dycore_RG_steps/'; done && ${Python3_EXECUTABLE} -I '${T}/dycore_compare.py' '${CPPMPAS_MPAS_COMPARE}' '${REF_RG}/steps' '${W}/dycore_RG_steps' '${W}/dycore_RG_compare.json' 165")
set_tests_properties(dycore_RG_mpas_compare PROPERTIES FIXTURES_REQUIRED dycore_RG)
foreach(f IN LISTS DYCORE_RG_RESTART_FIELDS)
  add_test(NAME dycore_RG_restart_${f} COMMAND ${PY} ${T}/compare_arrays.py
    ${REF_RG}/forecast/restart.2026-10-12_01.10.00.nc ${W}/dycore_RG/state.2026-10-12_01.10.00.nc ${f})
  set_tests_properties(dycore_RG_restart_${f} PROPERTIES FIXTURES_REQUIRED dycore_RG)
endforeach()
foreach(f IN LISTS DYCORE_RG_HISTORY_FIELDS)
  add_test(NAME dycore_RG_history1h_${f} COMMAND ${PY} ${T}/compare_arrays.py
    ${REF_RG}/forecast/history.2026-10-12_01.00.00.nc ${W}/dycore_RG/state.2026-10-12_01.00.00.nc ${f})
  set_tests_properties(dycore_RG_history1h_${f} PROPERTIES FIXTURES_REQUIRED dycore_RG)
  add_test(NAME dycore_RG_history2h_${f} COMMAND ${PY} ${T}/compare_arrays.py
    ${REF_RG}/forecast/history.2026-10-12_02.00.00.nc ${W}/dycore_RG/state.2026-10-12_02.00.00.nc ${f})
  set_tests_properties(dycore_RG_history2h_${f} PROPERTIES FIXTURES_REQUIRED dycore_RG)
endforeach()

# RG as the reference's forecast run writes it (cppmpas run --streams with the forecast's namelist.atmosphere and
# streams.atmosphere: history every 1 h, diagnostics every 3 h, restart at 01:10; mpas_atm_core.F:625-919): every
# file of mpas842-dbl-rg-004's forecast directory (diag 00:00, history 00:00 / 01:00 / 02:00, restart 01:10), every
# variable, its attributes and the global attributes (the namelist options included), the restart_timestamp file;
# one test per file and one for the whole directory (about 50 s).
set(RG_STREAM_FILES diag.2026-10-12_00.00.00 history.2026-10-12_00.00.00 history.2026-10-12_01.00.00
  history.2026-10-12_02.00.00 restart.2026-10-12_01.10.00)
add_test(NAME dycore_run_RG_streams COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/dycore_RG_streams' && mkdir -p '${W}/dycore_RG_streams' && '${EXE}' run '${REF_RG}/init/x1.10242.init.nc' '${W}/dycore_RG_streams' --namelist '${REF_RG}/forecast/namelist.atmosphere' --streams '${REF_RG}/forecast/streams.atmosphere' --scalars qv,qc,qr,qi,qs,qg --every 0 --tables '${CPPMPAS_MPAS_TABLES}' --rrtmg '${CPPMPAS_CPPWRF_RRTMG}'")
set_tests_properties(dycore_run_RG_streams PROPERTIES FIXTURES_SETUP dycore_RG_streams TIMEOUT 1800)
foreach(file IN LISTS RG_STREAM_FILES)
  add_test(NAME dycore_RG_streams_${file} COMMAND ${PY} ${CPPMPAS_MPAS_COMPARE}
    ${REF_RG}/forecast/${file}.nc ${W}/dycore_RG_streams/${file}.nc)
  set_tests_properties(dycore_RG_streams_${file} PROPERTIES FIXTURES_REQUIRED dycore_RG_streams)
endforeach()
add_test(NAME dycore_RG_streams_all COMMAND ${PY} ${T}/streams_compare.py ${CPPMPAS_MPAS_COMPARE}
  ${REF_RG}/forecast ${W}/dycore_RG_streams ${W}/dycore_RG_streams.json 5 702)
set_tests_properties(dycore_RG_streams_all PROPERTIES FIXTURES_REQUIRED dycore_RG_streams)

# RG's radiation calls replayed (tests/radiation_replay): cppmpas run of RG's steps configuration for 4 steps with
# CPPMPAS_RADIATION_DUMP set writes the radiation calls of steps 1 and 4 (00:00 and 00:30: the inputs and the outputs
# the run computed); test_radiation_replay runs each call again from its inputs and compares every radiation output
# (o3clim, cldfrac, o3vmr, rre_*, the LW and SW fluxes, coszr, rthratenlw / rthratensw) with the run's at tolerance 0.
# This compares CPPMPAS with CPPMPAS, not with MPAS: the replay proves the step is reproducible from its dumped inputs
# (and gpu_radiation_RG_* that the device reproduces the CPU). The MPAS proof of the radiation comes from the RG tests
# (dycore_RG_*, dycore_RG_streams_*: the run is bit-identical to rg-004) and from radiation_column.
add_executable(test_radiation_replay radiation_replay/test_radiation_replay.cpp)
target_link_libraries(test_radiation_replay PRIVATE cppmpas)
target_compile_options(test_radiation_replay PRIVATE -Wall -Wextra)
set(RADIATION_DUMP_RG ${W}/radiation_dump_RG)
add_test(NAME radiation_dump_RG COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${RADIATION_DUMP_RG}' && mkdir -p '${RADIATION_DUMP_RG}/run' && CPPMPAS_RADIATION_DUMP='${RADIATION_DUMP_RG}' '${EXE}' run '${REF_RG}/init/x1.10242.init.nc' '${RADIATION_DUMP_RG}/run' --namelist '${REF_RG}/steps/namelist.atmosphere' --scalars qv,qc,qr,qi,qs,qg --steps 4 --every 0 --tables '${CPPMPAS_MPAS_TABLES}' --rrtmg '${CPPMPAS_CPPWRF_RRTMG}' && ls -l '${RADIATION_DUMP_RG}' && test -s '${RADIATION_DUMP_RG}/radiation.0.bin' && test -s '${RADIATION_DUMP_RG}/radiation.1800.bin'")
set_tests_properties(radiation_dump_RG PROPERTIES FIXTURES_SETUP radiation_dump_RG TIMEOUT 1800)
set(RADIATION_RG_CALLS 0 1800)
foreach(t IN LISTS RADIATION_RG_CALLS)
  add_test(NAME radiation_replay_RG_${t} COMMAND test_radiation_replay ${RADIATION_DUMP_RG}/radiation.${t}.bin
    ${CPPMPAS_CPPWRF_RRTMG}/rrtmg_lw_tables.bin ${CPPMPAS_CPPWRF_RRTMG}/rrtmg_sw_tables.bin)
  set_tests_properties(radiation_replay_RG_${t} PROPERTIES FIXTURES_REQUIRED radiation_dump_RG)
endforeach()

# RG24 (mpas842-dbl-rg24-001, the v0.1 target: the RG case for 24 h; label "long", about 15 minutes): the forecast
# run as the reference made it (cppmpas run --streams: 144 steps, history and diagnostics every 3 h from 00:00 to
# 2026-10-13_00, restart at 12:10), every file of the reference's forecast directory compared at tolerance 0 with
# tools/mpas_compare.py (all variables, variable and global attributes) and the restart_timestamp file.
set(REF_RG24 ${CPPMPAS_REFS_DIR}/mpas842-dbl-rg24-001/RG24)
set(RG24_STREAM_STAMPS 2026-10-12_00.00.00 2026-10-12_03.00.00 2026-10-12_06.00.00 2026-10-12_09.00.00
  2026-10-12_12.00.00 2026-10-12_15.00.00 2026-10-12_18.00.00 2026-10-12_21.00.00 2026-10-13_00.00.00)
set(RG24_STREAM_FILES restart.2026-10-12_12.10.00)
foreach(stamp IN LISTS RG24_STREAM_STAMPS)
  list(APPEND RG24_STREAM_FILES history.${stamp} diag.${stamp})
endforeach()
add_test(NAME dycore_run_long_RG24 COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/dycore_long_RG24' && mkdir -p '${W}/dycore_long_RG24' && '${EXE}' run '${REF_RG24}/init/x1.10242.init.nc' '${W}/dycore_long_RG24' --namelist '${REF_RG24}/forecast/namelist.atmosphere' --streams '${REF_RG24}/forecast/streams.atmosphere' --scalars qv,qc,qr,qi,qs,qg --every 0 --tables '${CPPMPAS_MPAS_TABLES}' --rrtmg '${CPPMPAS_CPPWRF_RRTMG}'")
set_tests_properties(dycore_run_long_RG24 PROPERTIES FIXTURES_SETUP dycore_long_RG24 LABELS long TIMEOUT 7200)
foreach(file IN LISTS RG24_STREAM_FILES)
  add_test(NAME dycore_RG24_${file} COMMAND ${PY} ${CPPMPAS_MPAS_COMPARE}
    ${REF_RG24}/forecast/${file}.nc ${W}/dycore_long_RG24/${file}.nc)
  set_tests_properties(dycore_RG24_${file} PROPERTIES FIXTURES_REQUIRED dycore_long_RG24 LABELS long)
endforeach()
add_test(NAME dycore_RG24_all COMMAND ${PY} ${T}/streams_compare.py ${CPPMPAS_MPAS_COMPARE}
  ${REF_RG24}/forecast ${W}/dycore_long_RG24 ${W}/dycore_long_RG24.json 19 2168)
set_tests_properties(dycore_RG24_all PROPERTIES FIXTURES_REQUIRED dycore_long_RG24 LABELS long)

# ---------------------------------------------------------------------------------------------------------------
# RGCP (mpas842-dbl-rgcp-001): RG with config_physics_suite = 'convection_permitting' (Thompson with n_microp = 7
# sub-steps, Grell-Freitas, MYNN PBL and surface layer, YSU's GWDO, Noah, RRTMG with cld_fraction; the scalars qv, qc,
# qr, qi, qs, qg, ni, nr). The Thompson tables are the reference run's (RGCP/tables, read through --thompson as
# thompson_init reads them from MPAS's run directory). As for RG: one run of 12 steps writes the every-step stream
# (steps 0-10 at tolerance 0, one test per step and field, and the directory through tools/mpas_compare.py) and,
# after steps 6, 7 and 12, the fields of the forecast run's history (1 h, 2 h) and restart (01:10) streams; relhum
# and refl10cm_max are left out there (the microphysics driver computes them only when the diagnostics alarm rings,
# l_diags, which a steps-configuration run never does; the --streams test below has them). Needs
# CPPMPAS_MPAS_TABLES and CPPMPAS_CPPWRF_RRTMG as RG.
set(REF_RGCP ${CPPMPAS_REFS_DIR}/mpas842-dbl-rgcp-001/RGCP)
set(DYCORE_RGCP_FIELDS xtime Time u w rho_zz theta_m rho theta pressure_p qv qc qr qi qs qg ni nr)
set(DYCORE_RGCP_RESTART_FIELDS ru_p rw_p rt_diabatic_tend rthdynten isltyp ivgtyp mminlu isice_lu iswater_lu landmask shdmin shdmax snoalb albedo12m greenfrac var2d con oa1 oa2 oa3 oa4 ol1 ol2 ol3 ol4 qv qc qr qi qs qg ni nr u_init v_init t_init xtime Time u w rho_zz theta_m pressure_p rho theta exner exner_base rtheta_base pressure_base rho_base theta_base ru rw rtheta_p rho_p surface_pressure rainprod evapprod re_cloud re_ice re_snow tend_sfc_pressure nsteps_accum ndays_accum tlag tday_accum tyear_mean tyear_accum i_rainnc rainncv rainnc snowncv snownc graupelncv graupelnc sr i_rainc cuprec rainc raincv qcg ch kpbl hpbl hfx mavail mol qfx qsfc ust ustm zol znt br cd cda chs chs2 cqs2 ck cka cpm flhc flqc gz1oz0 lh psim psih qgh regime rmol wspd fh fm u10 v10 q2 t2m th2m dusfcg dvsfcg dtaux3d dtauy3d rubldiff rvbldiff i_acswdnb i_acswdnbc i_acswdnt i_acswdntc i_acswupb i_acswupbc i_acswupt i_acswuptc acswdnb acswdnbc acswdnt acswdntc acswupb acswupbc acswupt acswuptc gsw i_aclwdnb i_aclwdnbc i_aclwdnt i_aclwdntc i_aclwupb i_aclwupbc i_aclwupt i_aclwuptc aclwdnb aclwdnbc aclwdnt aclwdntc aclwupb aclwupbc aclwupt aclwuptc glw o3clim o3vmr acsnom acsnow canwat chklowq grdflx lai noahres potevp sfc_albedo sfc_emiss sfc_emibck sfcrunoff smstav smstot snopcx snotime sstsk sstsk_dtc sstsk_dtw thc udrunoff xicem z0 zs rthmpten rqvmpten rqcmpten rqrmpten rqimpten rqsmpten rqgmpten rthcuten rqvcuten rqccuten rqicuten rqvdynten rucuten rvcuten rublten rvblten rthblten rqvblten rqcblten rqiblten rthratensw rthratenlw sfc_albbck skintemp snow snowc snowh sst tmn vegfra seaice xice xland dzs smcrel sh2o smois tslb t_oml t_oml_initial t_oml_200m_initial h_oml h_oml_initial hu_oml hv_oml nt_c mu_c cubot cutop sh3d sm3d qke qke_adv cov qsq tsq tke_pbl el_pbl edmf_a edmf_w edmf_qt edmf_qc edmf_thl edmf_ent sub_thl sub_qv det_thl det_qv cldfrac_bl qc_bl qi_bl k22_shallow kbcon_shallow ktop_shallow ktop_deep xmb_total xmb_shallow qc_cu qi_cu rnimpten rnrmpten rqsblten rncblten rniblten rnifablten rnwfablten qv_init initial_time circulation)
set(DYCORE_RGCP_HISTORY_FIELDS qv qc qr qi qs qg ni nr xtime Time u w pressure surface_pressure rho theta divergence vorticity ke uReconstructZonal uReconstructMeridional i_rainnc rainnc precipw cuprec i_rainc rainc kpbl hpbl hfx qfx cd cda ck cka lh u10 v10 q2 t2m th2m gsw glw acsnow pin ozmixm mminlu isice_lu iswater_lu xland skintemp snow snowh sst vegfra xice sh2o smois tslb initial_time)
set(DYCORE_RGCP_LONG_FIELDS ${DYCORE_RGCP_RESTART_FIELDS} ${DYCORE_RGCP_HISTORY_FIELDS})
list(REMOVE_DUPLICATES DYCORE_RGCP_LONG_FIELDS)
string(REPLACE ";" "," DYCORE_RGCP_LONG_LIST "${DYCORE_RGCP_LONG_FIELDS}")
set(RGCP_SCALARS qv,qc,qr,qi,qs,qg,ni,nr)
add_test(NAME dycore_run_RGCP COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/dycore_RGCP' && mkdir -p '${W}/dycore_RGCP' && '${EXE}' run '${REF_RGCP}/init/x1.10242.init.nc' '${W}/dycore_RGCP' --namelist '${REF_RGCP}/steps/namelist.atmosphere' --scalars ${RGCP_SCALARS} --steps 12 --write-at 6,7,12 --write-fields ${DYCORE_RGCP_LONG_LIST} --tables '${CPPMPAS_MPAS_TABLES}' --rrtmg '${CPPMPAS_CPPWRF_RRTMG}' --thompson '${REF_RGCP}/tables'")
set_tests_properties(dycore_run_RGCP PROPERTIES FIXTURES_SETUP dycore_RGCP TIMEOUT 1800)
set(step 0)
foreach(stamp IN LISTS DYCORE_RG_STAMPS)
  foreach(f IN LISTS DYCORE_RGCP_FIELDS)
    add_test(NAME dycore_RGCP_step${step}_${f} COMMAND ${PY} ${T}/compare_arrays.py
      ${REF_RGCP}/steps/step.2026-10-12_${stamp}.nc ${W}/dycore_RGCP/step.2026-10-12_${stamp}.nc ${f})
    set_tests_properties(dycore_RGCP_step${step}_${f} PROPERTIES FIXTURES_REQUIRED dycore_RGCP)
  endforeach()
  math(EXPR step "${step} + 1")
endforeach()
add_test(NAME dycore_RGCP_mpas_compare COMMAND ${CMAKE_COMMAND} -E env sh -c
  "set -e; rm -rf '${W}/dycore_RGCP_steps' && mkdir -p '${W}/dycore_RGCP_steps' && for s in ${DYCORE_RG_STAMP_WORDS}; do ln -s '${W}/dycore_RGCP/step.2026-10-12_'$s.nc '${W}/dycore_RGCP_steps/'; done && ${Python3_EXECUTABLE} -I '${T}/dycore_compare.py' '${CPPMPAS_MPAS_COMPARE}' '${REF_RGCP}/steps' '${W}/dycore_RGCP_steps' '${W}/dycore_RGCP_compare.json' 187")
set_tests_properties(dycore_RGCP_mpas_compare PROPERTIES FIXTURES_REQUIRED dycore_RGCP)
foreach(f IN LISTS DYCORE_RGCP_RESTART_FIELDS)
  add_test(NAME dycore_RGCP_restart_${f} COMMAND ${PY} ${T}/compare_arrays.py
    ${REF_RGCP}/forecast/restart.2026-10-12_01.10.00.nc ${W}/dycore_RGCP/state.2026-10-12_01.10.00.nc ${f})
  set_tests_properties(dycore_RGCP_restart_${f} PROPERTIES FIXTURES_REQUIRED dycore_RGCP)
endforeach()
foreach(f IN LISTS DYCORE_RGCP_HISTORY_FIELDS)
  add_test(NAME dycore_RGCP_history1h_${f} COMMAND ${PY} ${T}/compare_arrays.py
    ${REF_RGCP}/forecast/history.2026-10-12_01.00.00.nc ${W}/dycore_RGCP/state.2026-10-12_01.00.00.nc ${f})
  set_tests_properties(dycore_RGCP_history1h_${f} PROPERTIES FIXTURES_REQUIRED dycore_RGCP)
  add_test(NAME dycore_RGCP_history2h_${f} COMMAND ${PY} ${T}/compare_arrays.py
    ${REF_RGCP}/forecast/history.2026-10-12_02.00.00.nc ${W}/dycore_RGCP/state.2026-10-12_02.00.00.nc ${f})
  set_tests_properties(dycore_RGCP_history2h_${f} PROPERTIES FIXTURES_REQUIRED dycore_RGCP)
endforeach()

# RGCP as the reference's forecast run writes it (cppmpas run --streams with the forecast's namelist.atmosphere and
# streams.atmosphere, as RG): every file of mpas842-dbl-rgcp-001's forecast directory (diag 00:00, history 00:00 /
# 01:00 / 02:00, restart 01:10), every variable, its attributes and the global attributes, the restart_timestamp
# file; one test per file and one for the whole directory.
add_test(NAME dycore_run_RGCP_streams COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/dycore_RGCP_streams' && mkdir -p '${W}/dycore_RGCP_streams' && '${EXE}' run '${REF_RGCP}/init/x1.10242.init.nc' '${W}/dycore_RGCP_streams' --namelist '${REF_RGCP}/forecast/namelist.atmosphere' --streams '${REF_RGCP}/forecast/streams.atmosphere' --scalars ${RGCP_SCALARS} --every 0 --tables '${CPPMPAS_MPAS_TABLES}' --rrtmg '${CPPMPAS_CPPWRF_RRTMG}' --thompson '${REF_RGCP}/tables'")
set_tests_properties(dycore_run_RGCP_streams PROPERTIES FIXTURES_SETUP dycore_RGCP_streams TIMEOUT 1800)
foreach(file IN LISTS RG_STREAM_FILES)
  add_test(NAME dycore_RGCP_streams_${file} COMMAND ${PY} ${CPPMPAS_MPAS_COMPARE}
    ${REF_RGCP}/forecast/${file}.nc ${W}/dycore_RGCP_streams/${file}.nc)
  set_tests_properties(dycore_RGCP_streams_${file} PROPERTIES FIXTURES_REQUIRED dycore_RGCP_streams)
endforeach()
add_test(NAME dycore_RGCP_streams_all COMMAND ${PY} ${T}/streams_compare.py ${CPPMPAS_MPAS_COMPARE}
  ${REF_RGCP}/forecast ${W}/dycore_RGCP_streams ${W}/dycore_RGCP_streams.json 5 749)
set_tests_properties(dycore_RGCP_streams_all PROPERTIES FIXTURES_REQUIRED dycore_RGCP_streams)

# RGCP24 (mpas842-dbl-rgcp24-001: RGCP for 24 h; label "long"): the forecast run as the reference made it (144 steps,
# history and diagnostics every 3 h from 00:00 to 2026-10-13_00, restart at 12:10), with its own run's Thompson
# tables, every file of the reference's forecast directory compared at tolerance 0 with tools/mpas_compare.py and the
# restart_timestamp file.
set(REF_RGCP24 ${CPPMPAS_REFS_DIR}/mpas842-dbl-rgcp24-001/RGCP24)
add_test(NAME dycore_run_long_RGCP24 COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/dycore_long_RGCP24' && mkdir -p '${W}/dycore_long_RGCP24' && '${EXE}' run '${REF_RGCP24}/init/x1.10242.init.nc' '${W}/dycore_long_RGCP24' --namelist '${REF_RGCP24}/forecast/namelist.atmosphere' --streams '${REF_RGCP24}/forecast/streams.atmosphere' --scalars ${RGCP_SCALARS} --every 0 --tables '${CPPMPAS_MPAS_TABLES}' --rrtmg '${CPPMPAS_CPPWRF_RRTMG}' --thompson '${REF_RGCP24}/tables'")
set_tests_properties(dycore_run_long_RGCP24 PROPERTIES FIXTURES_SETUP dycore_long_RGCP24 LABELS long TIMEOUT 7200)
foreach(file IN LISTS RG24_STREAM_FILES)
  add_test(NAME dycore_RGCP24_${file} COMMAND ${PY} ${CPPMPAS_MPAS_COMPARE}
    ${REF_RGCP24}/forecast/${file}.nc ${W}/dycore_long_RGCP24/${file}.nc)
  set_tests_properties(dycore_RGCP24_${file} PROPERTIES FIXTURES_REQUIRED dycore_long_RGCP24 LABELS long)
endforeach()
add_test(NAME dycore_RGCP24_all COMMAND ${PY} ${T}/streams_compare.py ${CPPMPAS_MPAS_COMPARE}
  ${REF_RGCP24}/forecast ${W}/dycore_long_RGCP24 ${W}/dycore_long_RGCP24.json 19 2227)
set_tests_properties(dycore_RGCP24_all PROPERTIES FIXTURES_REQUIRED dycore_long_RGCP24 LABELS long)

# ---------------------------------------------------------------------------------------------------------------
# GPU proof: the CUDA build of the same sources (cppmpas_gpu, CMake option CPPMPAS_ENABLE_CUDA; every parallel_for a
# kernel on the GPU, core/backend.hpp) against the same references at tolerance 0, and against the CPU build's own
# output (cppmpas, the tests above). Not all of cppmpas_gpu runs on the device: physics_core_init (init_dirs_forphys,
# landuse_init_forMPAS, init_lsm's soil-water loop, the accumulators), physics_update_surface /
# monthly_interp_to_date, radconst, the isobaric level arrays and the refl10cm_1km_max reset are host code in the CUDA
# build too, so a gpu test proves those parts on the CPU only (tests/gpu_test.py lists them). Each test runs through
# tests/gpu_test.py: without a CUDA build or without a CUDA device it exits 77 and CTest reports it as "Skipped"
# (label gpu), never as passed; any other CUDA error fails.
if(CPPMPAS_ENABLE_CUDA)
  set(GEXE $<TARGET_FILE:cppmpas_gpu>)
else()
  set(GEXE NONE)  # tests/gpu_test.py: no CUDA build, every GPU test is skipped
endif()
set(GPU ${PY} ${T}/gpu_test.py ${GEXE})
# ON: a GPU test that cannot run (no CUDA build, no device) fails instead of being skipped (for a machine that must
# prove the GPU)
option(CPPMPAS_GPU_REQUIRED "GPU tests fail instead of skipping when they cannot run" OFF)
function(cppmpas_gpu_props)  # cppmpas_gpu_props(test... [FIXTURES_REQUIRED f] [FIXTURES_SETUP f] [LABELS l])
  cmake_parse_arguments(G "" "" "FIXTURES_REQUIRED;FIXTURES_SETUP;LABELS" ${ARGN})
  set_tests_properties(${G_UNPARSED_ARGUMENTS} PROPERTIES LABELS "gpu;${G_LABELS}")
  if(CPPMPAS_GPU_REQUIRED)
    set_tests_properties(${G_UNPARSED_ARGUMENTS} PROPERTIES ENVIRONMENT "CPPMPAS_GPU_REQUIRED=1")
  else()
    # a skipped run test (FIXTURES_SETUP) still lets CTest run the tests that need it: gpu_test.py --clean has removed
    # its outputs before the device probe, so those tests skip or fail, never compare an earlier run's files
    set_tests_properties(${G_UNPARSED_ARGUMENTS} PROPERTIES SKIP_RETURN_CODE 77)
  endif()
  if(G_FIXTURES_REQUIRED)
    set_tests_properties(${G_UNPARSED_ARGUMENTS} PROPERTIES FIXTURES_REQUIRED "${G_FIXTURES_REQUIRED}")
  endif()
  if(G_FIXTURES_SETUP)
    set_tests_properties(${G_UNPARSED_ARGUMENTS} PROPERTIES FIXTURES_SETUP "${G_FIXTURES_SETUP}")
  endif()
endfunction()
# the device and that the kernels run on it (backend, compute capability, driver and runtime are printed)
add_test(NAME gpu_device COMMAND ${GPU} device)
# no visible device: cppmpas_gpu exits 77 with a SKIP message (a skip, not a pass)
add_test(NAME gpu_skip_without_device COMMAND ${GPU} no-device)
cppmpas_gpu_props(gpu_device gpu_skip_without_device)
# the libm / libmvec wrappers on the device against the host, plain and negated arguments, special values
# included: non-NaN results bit-identical, NaN-ness identical, NaN signs / payloads as recorded for the card
if(CPPMPAS_ENABLE_CUDA)
  add_executable(test_gpu_libm test_gpu_libm.cu)
  target_link_libraries(test_gpu_libm PRIVATE cppmpas_cuda)
  set_target_properties(test_gpu_libm PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  add_test(NAME gpu_libm COMMAND test_gpu_libm 200000)
else()
  add_test(NAME gpu_libm COMMAND ${GPU} device)
endif()
cppmpas_gpu_props(gpu_libm)
# the launch wrappers on the device: element ranges, kernel_fail (message, the thread ends at the call, the first
# launch of a process), work-array capacity, the next launch after a failure, the stack limit never lowered
if(CPPMPAS_ENABLE_CUDA)
  add_executable(test_gpu_backend test_gpu_backend.cu)
  target_link_libraries(test_gpu_backend PRIVATE cppmpas_cuda)
  set_target_properties(test_gpu_backend PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  add_test(NAME gpu_backend COMMAND test_gpu_backend)
  # the CUDA flags that change results are refused at configure time, in CMAKE_CUDA_FLAGS and
  # CMAKE_CUDA_FLAGS_RELEASE, in nvcc's spellings and as host-compiler flags in any form (CMake may wrap the message
  # after "floating-point")
  set(_cuda_bad "--use_fast_math" "-fmad=true" "--fmad true" "--ftz=true" "--prec-div=false" "--prec-sqrt=false"
      "-Xcompiler -ffast-math" "--compiler-options=-O2,-ffp-contract=fast" "-Xcompiler=-funsafe-math-optimizations"
      "-Xcompiler -mfpmath=387")
  set(_i 0)
  foreach(_flag IN LISTS _cuda_bad)
    foreach(_var CMAKE_CUDA_FLAGS CMAKE_CUDA_FLAGS_RELEASE)
      add_test(NAME cmake_rejects_cuda_${_i}_${_var} COMMAND ${CMAKE_COMMAND} -E env sh -c
        "rm -rf '${W}/cmake_cuda_${_i}_${_var}'; out=$(${CMAKE_COMMAND} -S '${CMAKE_SOURCE_DIR}' -B '${W}/cmake_cuda_${_i}_${_var}' -DCMAKE_BUILD_TYPE=Release -DCPPMPAS_ENABLE_CUDA=ON -DCMAKE_CUDA_COMPILER='${CMAKE_CUDA_COMPILER}' -DCMAKE_CUDA_HOST_COMPILER='${CMAKE_CUDA_HOST_COMPILER}' -DCMAKE_CUDA_ARCHITECTURES='${CMAKE_CUDA_ARCHITECTURES}' '-D${_var}=${_flag}' 2>&1); rc=$?; rm -rf '${W}/cmake_cuda_${_i}_${_var}'; echo \"$out\" | tail -5; [ $rc -ne 0 ] && echo \"$out\" | grep -q 'changes the floating-point'")
      set_tests_properties(cmake_rejects_cuda_${_i}_${_var} PROPERTIES FIXTURES_REQUIRED work_dir LABELS gpu)
    endforeach()
    math(EXPR _i "${_i} + 1")
  endforeach()
  # and the strict flags follow the user's: a CUDA source compiles with them last on nvcc's command line
  # (and no flag that changes results comes after them; flags.make of the Makefile generators)
  add_test(NAME cmake_cuda_strict_flags_last COMMAND ${CMAKE_COMMAND} -E env sh -c
    "f='${CMAKE_BINARY_DIR}/CMakeFiles/cppmpas_cuda.dir/flags.make'; test -s \"$f\" && grep -h -- '--fmad=false' \"$f\" | grep -q -- '-O3 -DNDEBUG.*--fmad=false --ftz=false --prec-div=true --prec-sqrt=true -Xcompiler=-ffp-contract=off' && ! grep -E -- '-ffp-contract=off.*(fmad[= ]true|use_fast_math|ffast-math|ftz[= ]true|prec-div[= ]false|prec-sqrt[= ]false|ffp-contract=fast|funsafe-math|Ofast)' \"$f\"")
  set_tests_properties(cmake_cuda_strict_flags_last PROPERTIES LABELS gpu)
else()
  add_test(NAME gpu_backend COMMAND ${GPU} device)
endif()
cppmpas_gpu_props(gpu_backend)
# RG's radiation calls at 00:00 and 00:30 on the device (tests/radiation_replay, the CUDA build of the same test):
# every loop of the cloudiness, ozone and radiation drivers and every RRTMG stage a kernel, against the CPU run's
# outputs at tolerance 0
if(CPPMPAS_ENABLE_CUDA)
  set(_rr_cu ${CMAKE_BINARY_DIR}/cuda_src/tests/radiation_replay/test_radiation_replay.cpp.cu)
  file(CONFIGURE OUTPUT ${_rr_cu} CONTENT "// generated: the CUDA build of tests/radiation_replay/test_radiation_replay.cpp\n#include \"${CMAKE_CURRENT_SOURCE_DIR}/radiation_replay/test_radiation_replay.cpp\"\n")
  add_executable(test_gpu_radiation_replay ${_rr_cu})
  target_link_libraries(test_gpu_radiation_replay PRIVATE cppmpas_cuda)
  set_target_properties(test_gpu_radiation_replay PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  foreach(t IN LISTS RADIATION_RG_CALLS)
    add_test(NAME gpu_radiation_RG_${t} COMMAND test_gpu_radiation_replay ${RADIATION_DUMP_RG}/radiation.${t}.bin
      ${CPPMPAS_CPPWRF_RRTMG}/rrtmg_lw_tables.bin ${CPPMPAS_CPPWRF_RRTMG}/rrtmg_sw_tables.bin)
    cppmpas_gpu_props(gpu_radiation_RG_${t} FIXTURES_REQUIRED radiation_dump_RG)
  endforeach()
else()
  foreach(t IN LISTS RADIATION_RG_CALLS)
    add_test(NAME gpu_radiation_RG_${t} COMMAND ${GPU} device)
    cppmpas_gpu_props(gpu_radiation_RG_${t})
  endforeach()
endif()
# atm_mpas_init_block on the GPU: every initdump field of MW, SC and JW against MPAS, and against the CPU build
foreach(case MW SC JW)
  add_test(NAME gpu_initdump_run_${case} COMMAND ${GPU} --clean ${W}/gpu_initdump_${case}.nc -- ${CMAKE_COMMAND} -E env sh -c
    "mkdir -p '${W}' && rm -f '${W}/gpu_initdump_${case}.nc' && '${GEXE}' init '${DUMPINIT_${case}}' '${W}/gpu_initdump_${case}.nc' --namelist '${DUMPREF_${case}}/steps/namelist.atmosphere' --scalars ${SCAL_${case}} --fields ${INITDUMP_LIST}")
  cppmpas_gpu_props(gpu_initdump_run_${case} FIXTURES_SETUP gpu_initdump_${case})
  foreach(f IN LISTS INITDUMP_FIELDS)
    add_test(NAME gpu_initdump_${case}_${f} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py
      ${DUMPREF_${case}}/steps/initdump.nc ${W}/gpu_initdump_${case}.nc ${f})
    cppmpas_gpu_props(gpu_initdump_${case}_${f} FIXTURES_REQUIRED gpu_initdump_${case})
  endforeach()
  add_test(NAME gpu_initdump_vs_cpu_${case} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py
    ${W}/initdump_${case}.nc ${W}/gpu_initdump_${case}.nc ${INITDUMP_LIST})
  cppmpas_gpu_props(gpu_initdump_vs_cpu_${case} FIXTURES_REQUIRED "gpu_initdump_${case};initdump_${case}")
  # the init-block fields of the restart stream the initdump set leaves out (mpas_initialize_vectors, the RBF
  # coefficients, dss, nominalMinDc, the base state) and the t = 0 state, as init_restart / init_step0 on the CPU
  add_test(NAME gpu_init_run_${case} COMMAND ${GPU} --clean ${W}/gpu_init_${case}.nc -- ${CMAKE_COMMAND} -E env sh -c
    "mkdir -p '${W}' && rm -f '${W}/gpu_init_${case}.nc' && '${GEXE}' init '${INIT_${case}}' '${W}/gpu_init_${case}.nc' --namelist '${REF_${case}}/steps/namelist.atmosphere' --scalars ${SCAL_${case}} --fields ${INIT_RESTART_FIELDS},${INIT_STEP0_FIELDS}")
  cppmpas_gpu_props(gpu_init_run_${case} FIXTURES_SETUP gpu_init_${case})
  add_test(NAME gpu_init_restart_${case} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py ${RST_${case}}
    ${W}/gpu_init_${case}.nc ${INIT_RESTART_FIELDS})
  add_test(NAME gpu_init_step0_${case} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py
    ${REF_${case}}/steps/step.0000-01-01_00.00.00.nc ${W}/gpu_init_${case}.nc ${INIT_STEP0_FIELDS})
  cppmpas_gpu_props(gpu_init_restart_${case} gpu_init_step0_${case} FIXTURES_REQUIRED gpu_init_${case})
endforeach()
# the dycore on the GPU: ten MW steps, every step file against MPAS (per field, and the whole directory through
# tools/mpas_compare.py) and the whole directory against the CPU build's
add_test(NAME gpu_dycore_run_MW COMMAND ${GPU} --clean ${W}/gpu_dycore_MW -- ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/gpu_dycore_MW' && mkdir -p '${W}/gpu_dycore_MW' && '${GEXE}' run '${INIT_MW}' '${W}/gpu_dycore_MW' --namelist '${REF_MW}/steps/namelist.atmosphere' --steps 10")
cppmpas_gpu_props(gpu_dycore_run_MW FIXTURES_SETUP gpu_dycore_MW)
set(step 0)
foreach(stamp IN LISTS DYCORE_MW_STAMPS)
  foreach(f IN LISTS DYCORE_MW_FIELDS)
    add_test(NAME gpu_dycore_MW_step${step}_${f} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py
      ${REF_MW}/steps/step.0000-01-01_${stamp}.nc ${W}/gpu_dycore_MW/step.0000-01-01_${stamp}.nc ${f})
    cppmpas_gpu_props(gpu_dycore_MW_step${step}_${f} FIXTURES_REQUIRED gpu_dycore_MW)
  endforeach()
  math(EXPR step "${step} + 1")
endforeach()
add_test(NAME gpu_dycore_MW_mpas_compare COMMAND ${GPU} -- ${PY} ${T}/dycore_compare.py ${CPPMPAS_MPAS_COMPARE}
  ${REF_MW}/steps ${W}/gpu_dycore_MW ${W}/gpu_dycore_MW_compare.json 110)
add_test(NAME gpu_dycore_MW_vs_cpu COMMAND ${GPU} -- ${PY} ${T}/dycore_compare.py ${CPPMPAS_MPAS_COMPARE}
  ${W}/dycore_MW ${W}/gpu_dycore_MW ${W}/gpu_dycore_MW_vs_cpu.json 110)
cppmpas_gpu_props(gpu_dycore_MW_mpas_compare FIXTURES_REQUIRED gpu_dycore_MW)
cppmpas_gpu_props(gpu_dycore_MW_vs_cpu FIXTURES_REQUIRED "gpu_dycore_MW;dycore_MW")
# 301 steps (labels gpu, long): the 30-min output record and the 30:06 restart of the forecast run, as on the CPU
add_test(NAME gpu_dycore_run_long_MW COMMAND ${GPU} --clean ${W}/gpu_dycore_long_MW -- ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/gpu_dycore_long_MW' && mkdir -p '${W}/gpu_dycore_long_MW' && '${GEXE}' run '${INIT_MW}' '${W}/gpu_dycore_long_MW' --namelist '${REF_MW}/steps/namelist.atmosphere' --steps 301 --every 0 --write-at 300,301 --write-fields ${DYCORE_LONG_LIST}")
cppmpas_gpu_props(gpu_dycore_run_long_MW FIXTURES_SETUP gpu_dycore_long_MW LABELS long)
set_tests_properties(gpu_dycore_run_long_MW PROPERTIES TIMEOUT 1800)
foreach(f IN LISTS DYCORE_RESTART_FIELDS)
  add_test(NAME gpu_dycore_MW_restart_${f} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py ${RST_MW}
    ${W}/gpu_dycore_long_MW/state.0000-01-01_00.30.06.nc ${f})
  cppmpas_gpu_props(gpu_dycore_MW_restart_${f} FIXTURES_REQUIRED gpu_dycore_long_MW LABELS long)
endforeach()
foreach(f IN LISTS DYCORE_OUTPUT_FIELDS)
  add_test(NAME gpu_dycore_MW_output30_${f} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py ${REF_MW}/forecast/output.nc
    ${W}/gpu_dycore_long_MW/state.0000-01-01_00.30.00.nc ${f} --record 1)
  cppmpas_gpu_props(gpu_dycore_MW_output30_${f} FIXTURES_REQUIRED gpu_dycore_long_MW LABELS long)
endforeach()
# 600 steps (labels gpu, long): the 60-min output record (record 2 of the forecast run's output.nc)
add_test(NAME gpu_dycore_run_60min_MW COMMAND ${GPU} --clean ${W}/gpu_dycore_60min_MW -- ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/gpu_dycore_60min_MW' && mkdir -p '${W}/gpu_dycore_60min_MW' && '${GEXE}' run '${INIT_MW}' '${W}/gpu_dycore_60min_MW' --namelist '${REF_MW}/steps/namelist.atmosphere' --steps 600 --every 0 --write-at 600 --write-fields ${DYCORE_LONG_LIST}")
cppmpas_gpu_props(gpu_dycore_run_60min_MW FIXTURES_SETUP gpu_dycore_60min_MW LABELS long)
set_tests_properties(gpu_dycore_run_60min_MW PROPERTIES TIMEOUT 3600)
foreach(f IN LISTS DYCORE_OUTPUT_FIELDS)
  add_test(NAME gpu_dycore_MW_output60_${f} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py ${REF_MW}/forecast/output.nc
    ${W}/gpu_dycore_60min_MW/state.0000-01-01_01.00.00.nc ${f} --record 2)
  cppmpas_gpu_props(gpu_dycore_MW_output60_${f} FIXTURES_REQUIRED gpu_dycore_60min_MW LABELS long)
endforeach()
# the JW baroclinic wave on the GPU: ten steps against MPAS and the CPU build, and (label long) the 24 h output and
# the 24:07:30 restart
add_test(NAME gpu_dycore_run_JW COMMAND ${GPU} --clean ${W}/gpu_dycore_JW -- ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/gpu_dycore_JW' && mkdir -p '${W}/gpu_dycore_JW' && '${GEXE}' run '${INIT_JW}' '${W}/gpu_dycore_JW' --namelist '${REF_JW}/steps/namelist.atmosphere' --steps 10")
cppmpas_gpu_props(gpu_dycore_run_JW FIXTURES_SETUP gpu_dycore_JW)
set_tests_properties(gpu_dycore_run_JW PROPERTIES TIMEOUT 900)
set(step 0)
foreach(stamp IN LISTS DYCORE_JW_STAMPS)
  foreach(f IN LISTS DYCORE_MW_FIELDS)
    add_test(NAME gpu_dycore_JW_step${step}_${f} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py
      ${REF_JW}/steps/step.0000-01-01_${stamp}.nc ${W}/gpu_dycore_JW/step.0000-01-01_${stamp}.nc ${f})
    cppmpas_gpu_props(gpu_dycore_JW_step${step}_${f} FIXTURES_REQUIRED gpu_dycore_JW)
  endforeach()
  math(EXPR step "${step} + 1")
endforeach()
add_test(NAME gpu_dycore_JW_mpas_compare COMMAND ${GPU} -- ${PY} ${T}/dycore_compare.py ${CPPMPAS_MPAS_COMPARE}
  ${REF_JW}/steps ${W}/gpu_dycore_JW ${W}/gpu_dycore_JW_compare.json 110)
add_test(NAME gpu_dycore_JW_vs_cpu COMMAND ${GPU} -- ${PY} ${T}/dycore_compare.py ${CPPMPAS_MPAS_COMPARE}
  ${W}/dycore_JW ${W}/gpu_dycore_JW ${W}/gpu_dycore_JW_vs_cpu.json 110)
cppmpas_gpu_props(gpu_dycore_JW_mpas_compare FIXTURES_REQUIRED gpu_dycore_JW)
cppmpas_gpu_props(gpu_dycore_JW_vs_cpu FIXTURES_REQUIRED "gpu_dycore_JW;dycore_JW")
add_test(NAME gpu_dycore_run_long_JW COMMAND ${GPU} --clean ${W}/gpu_dycore_long_JW -- ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/gpu_dycore_long_JW' && mkdir -p '${W}/gpu_dycore_long_JW' && '${GEXE}' run '${INIT_JW}' '${W}/gpu_dycore_long_JW' --namelist '${REF_JW}/steps/namelist.atmosphere' --steps 193 --every 0 --write-at 192,193 --write-fields ${DYCORE_LONG_LIST}")
cppmpas_gpu_props(gpu_dycore_run_long_JW FIXTURES_SETUP gpu_dycore_long_JW LABELS long)
set_tests_properties(gpu_dycore_run_long_JW PROPERTIES TIMEOUT 3600)
foreach(f IN LISTS DYCORE_RESTART_FIELDS)
  add_test(NAME gpu_dycore_JW_restart_${f} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py ${RST_JW}
    ${W}/gpu_dycore_long_JW/state.0000-01-02_00.07.30.nc ${f})
  cppmpas_gpu_props(gpu_dycore_JW_restart_${f} FIXTURES_REQUIRED gpu_dycore_long_JW LABELS long)
endforeach()
foreach(f IN LISTS DYCORE_OUTPUT_FIELDS)
  add_test(NAME gpu_dycore_JW_output24h_${f} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py ${REF_JW}/forecast/output.nc
    ${W}/gpu_dycore_long_JW/state.0000-01-02_00.00.00.nc ${f} --record 1)
  cppmpas_gpu_props(gpu_dycore_JW_output24h_${f} FIXTURES_REQUIRED gpu_dycore_long_JW LABELS long)
endforeach()
# the MW restart run on the GPU: from the reference's 00:30:06 restart file, step 0 against the restart stream and
# (label long) 01:00:00 against output.nc record 2
add_test(NAME gpu_dycore_restart0_run_MW COMMAND ${GPU} --clean ${W}/gpu_dycore_restart0_MW -- ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/gpu_dycore_restart0_MW' && mkdir -p '${W}/gpu_dycore_restart0_MW' && '${GEXE}' run '${RST_MW}' '${W}/gpu_dycore_restart0_MW' --namelist '${NML_RST_MW}' --steps 0 --every 0 --write-at 0 --write-fields ${DYCORE_RESTART0_LIST}")
cppmpas_gpu_props(gpu_dycore_restart0_run_MW FIXTURES_SETUP gpu_dycore_restart0_MW)
foreach(f IN LISTS DYCORE_RESTART0_FIELDS)
  add_test(NAME gpu_dycore_restart0_MW_${f} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py ${RST_MW}
    ${W}/gpu_dycore_restart0_MW/state.0000-01-01_00.30.06.nc ${f})
  cppmpas_gpu_props(gpu_dycore_restart0_MW_${f} FIXTURES_REQUIRED gpu_dycore_restart0_MW)
endforeach()
add_test(NAME gpu_dycore_restart_run_long_MW COMMAND ${GPU} --clean ${W}/gpu_dycore_restart_long_MW -- ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/gpu_dycore_restart_long_MW' && mkdir -p '${W}/gpu_dycore_restart_long_MW' && '${GEXE}' run '${RST_MW}' '${W}/gpu_dycore_restart_long_MW' --namelist '${NML_RST_MW}' --steps 299 --every 0 --write-at 299 --write-fields ${DYCORE_OUTPUT_LIST}")
cppmpas_gpu_props(gpu_dycore_restart_run_long_MW FIXTURES_SETUP gpu_dycore_restart_long_MW LABELS long)
set_tests_properties(gpu_dycore_restart_run_long_MW PROPERTIES TIMEOUT 1800)
foreach(f IN LISTS DYCORE_OUTPUT_FIELDS)
  add_test(NAME gpu_dycore_restart_MW_output60_${f} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py
    ${REF_MW}/forecast/output.nc ${W}/gpu_dycore_restart_long_MW/state.0000-01-01_01.00.00.nc ${f} --record 2)
  cppmpas_gpu_props(gpu_dycore_restart_MW_output60_${f} FIXTURES_REQUIRED gpu_dycore_restart_long_MW LABELS long)
endforeach()
# the supercell on the GPU (scalar transport, perturbation mixing, Kessler as kernels): ten steps against MPAS and the
# CPU build, and (label long) the 30-min output and the 30:03 restart
add_test(NAME gpu_dycore_run_SC COMMAND ${GPU} --clean ${W}/gpu_dycore_SC -- ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/gpu_dycore_SC' && mkdir -p '${W}/gpu_dycore_SC' && '${GEXE}' run '${INIT_SC}' '${W}/gpu_dycore_SC' --namelist '${REF_SC}/steps/namelist.atmosphere' --scalars qv,qc,qr --steps 10")
cppmpas_gpu_props(gpu_dycore_run_SC FIXTURES_SETUP gpu_dycore_SC)
set_tests_properties(gpu_dycore_run_SC PROPERTIES TIMEOUT 900)
set(step 0)
foreach(stamp IN LISTS DYCORE_SC_STAMPS)
  foreach(f IN LISTS DYCORE_SC_FIELDS)
    add_test(NAME gpu_dycore_SC_step${step}_${f} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py
      ${REF_SC}/steps/step.0000-01-01_${stamp}.nc ${W}/gpu_dycore_SC/step.0000-01-01_${stamp}.nc ${f})
    cppmpas_gpu_props(gpu_dycore_SC_step${step}_${f} FIXTURES_REQUIRED gpu_dycore_SC)
  endforeach()
  math(EXPR step "${step} + 1")
endforeach()
add_test(NAME gpu_dycore_SC_mpas_compare COMMAND ${GPU} -- ${PY} ${T}/dycore_compare.py ${CPPMPAS_MPAS_COMPARE}
  ${REF_SC}/steps ${W}/gpu_dycore_SC ${W}/gpu_dycore_SC_compare.json 132)
add_test(NAME gpu_dycore_SC_vs_cpu COMMAND ${GPU} -- ${PY} ${T}/dycore_compare.py ${CPPMPAS_MPAS_COMPARE}
  ${W}/dycore_SC ${W}/gpu_dycore_SC ${W}/gpu_dycore_SC_vs_cpu.json 132)
cppmpas_gpu_props(gpu_dycore_SC_mpas_compare FIXTURES_REQUIRED gpu_dycore_SC)
cppmpas_gpu_props(gpu_dycore_SC_vs_cpu FIXTURES_REQUIRED "gpu_dycore_SC;dycore_SC")
add_test(NAME gpu_dycore_run_long_SC COMMAND ${GPU} --clean ${W}/gpu_dycore_long_SC -- ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/gpu_dycore_long_SC' && mkdir -p '${W}/gpu_dycore_long_SC' && '${GEXE}' run '${INIT_SC}' '${W}/gpu_dycore_long_SC' --namelist '${REF_SC}/steps/namelist.atmosphere' --scalars qv,qc,qr --steps 601 --every 0 --write-at 600,601 --write-fields ${DYCORE_SC_LONG_LIST}")
cppmpas_gpu_props(gpu_dycore_run_long_SC FIXTURES_SETUP gpu_dycore_long_SC LABELS long)
set_tests_properties(gpu_dycore_run_long_SC PROPERTIES TIMEOUT 3600)
foreach(f IN LISTS DYCORE_SC_RESTART_FIELDS)
  add_test(NAME gpu_dycore_SC_restart_${f} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py ${RST_SC}
    ${W}/gpu_dycore_long_SC/state.0000-01-01_00.30.03.nc ${f})
  cppmpas_gpu_props(gpu_dycore_SC_restart_${f} FIXTURES_REQUIRED gpu_dycore_long_SC LABELS long)
endforeach()
foreach(f IN LISTS DYCORE_SC_OUTPUT_FIELDS)
  add_test(NAME gpu_dycore_SC_output30_${f} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py ${OUT_SC}
    ${W}/gpu_dycore_long_SC/state.0000-01-01_00.30.00.nc ${f} --record 1)
  cppmpas_gpu_props(gpu_dycore_SC_output30_${f} FIXTURES_REQUIRED gpu_dycore_long_SC LABELS long)
endforeach()

# RG with the mesoscale_reference physics on the GPU (the physics drivers and schemes and the radiation as kernels; the
# host parts listed above): the 12 steps of dycore_run_RG against MPAS (step
# files per field and through tools/mpas_compare.py, the 01:10 restart and the 1 h / 2 h histories per field) and the
# step files against the CPU build's; the forecast run's streams of rg-004 (diag 00:00 with l_diags, histories,
# restart); and (label long) RG24's 19 stream files against rg24-001
add_test(NAME gpu_dycore_run_RG COMMAND ${GPU} --clean ${W}/gpu_dycore_RG -- ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/gpu_dycore_RG' && mkdir -p '${W}/gpu_dycore_RG' && '${GEXE}' run '${REF_RG}/init/x1.10242.init.nc' '${W}/gpu_dycore_RG' --namelist '${REF_RG}/steps/namelist.atmosphere' --scalars qv,qc,qr,qi,qs,qg --steps 12 --write-at 6,7,12 --write-fields ${DYCORE_RG_LONG_LIST} --tables '${CPPMPAS_MPAS_TABLES}' --rrtmg '${CPPMPAS_CPPWRF_RRTMG}'")
cppmpas_gpu_props(gpu_dycore_run_RG FIXTURES_SETUP gpu_dycore_RG)
set_tests_properties(gpu_dycore_run_RG PROPERTIES TIMEOUT 1800)
set(step 0)
foreach(stamp IN LISTS DYCORE_RG_STAMPS)
  foreach(f IN LISTS DYCORE_RG_FIELDS)
    add_test(NAME gpu_dycore_RG_step${step}_${f} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py
      ${REF_RG}/steps/step.2026-10-12_${stamp}.nc ${W}/gpu_dycore_RG/step.2026-10-12_${stamp}.nc ${f})
    cppmpas_gpu_props(gpu_dycore_RG_step${step}_${f} FIXTURES_REQUIRED gpu_dycore_RG)
  endforeach()
  math(EXPR step "${step} + 1")
endforeach()
add_test(NAME gpu_dycore_RG_mpas_compare COMMAND ${GPU} --clean ${W}/gpu_dycore_RG_steps -- ${CMAKE_COMMAND} -E env sh -c
  "set -e; rm -rf '${W}/gpu_dycore_RG_steps' && mkdir -p '${W}/gpu_dycore_RG_steps' && for s in ${DYCORE_RG_STAMP_WORDS}; do ln -s '${W}/gpu_dycore_RG/step.2026-10-12_'$s.nc '${W}/gpu_dycore_RG_steps/'; done && ${Python3_EXECUTABLE} -I '${T}/dycore_compare.py' '${CPPMPAS_MPAS_COMPARE}' '${REF_RG}/steps' '${W}/gpu_dycore_RG_steps' '${W}/gpu_dycore_RG_compare.json' 165")
cppmpas_gpu_props(gpu_dycore_RG_mpas_compare FIXTURES_REQUIRED gpu_dycore_RG)
# the step files of steps 0-10 only (both directories also hold the steps 11-12 and the --write-at state files)
add_test(NAME gpu_dycore_RG_vs_cpu COMMAND ${GPU} --clean ${W}/gpu_dycore_RG_vs_cpu_steps -- ${CMAKE_COMMAND} -E env sh -c
  "set -e; rm -rf '${W}/gpu_dycore_RG_vs_cpu_steps' '${W}/gpu_dycore_RG_vs_cpu_ref' && mkdir -p '${W}/gpu_dycore_RG_vs_cpu_steps' '${W}/gpu_dycore_RG_vs_cpu_ref' && for s in ${DYCORE_RG_STAMP_WORDS}; do ln -s '${W}/gpu_dycore_RG/step.2026-10-12_'$s.nc '${W}/gpu_dycore_RG_vs_cpu_steps/' && ln -s '${W}/dycore_RG/step.2026-10-12_'$s.nc '${W}/gpu_dycore_RG_vs_cpu_ref/'; done && ${Python3_EXECUTABLE} -I '${T}/dycore_compare.py' '${CPPMPAS_MPAS_COMPARE}' '${W}/gpu_dycore_RG_vs_cpu_ref' '${W}/gpu_dycore_RG_vs_cpu_steps' '${W}/gpu_dycore_RG_vs_cpu.json' 165")
cppmpas_gpu_props(gpu_dycore_RG_vs_cpu FIXTURES_REQUIRED "gpu_dycore_RG;dycore_RG")
foreach(f IN LISTS DYCORE_RG_RESTART_FIELDS)
  add_test(NAME gpu_dycore_RG_restart_${f} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py
    ${REF_RG}/forecast/restart.2026-10-12_01.10.00.nc ${W}/gpu_dycore_RG/state.2026-10-12_01.10.00.nc ${f})
  cppmpas_gpu_props(gpu_dycore_RG_restart_${f} FIXTURES_REQUIRED gpu_dycore_RG)
endforeach()
foreach(f IN LISTS DYCORE_RG_HISTORY_FIELDS)
  add_test(NAME gpu_dycore_RG_history1h_${f} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py
    ${REF_RG}/forecast/history.2026-10-12_01.00.00.nc ${W}/gpu_dycore_RG/state.2026-10-12_01.00.00.nc ${f})
  add_test(NAME gpu_dycore_RG_history2h_${f} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py
    ${REF_RG}/forecast/history.2026-10-12_02.00.00.nc ${W}/gpu_dycore_RG/state.2026-10-12_02.00.00.nc ${f})
  cppmpas_gpu_props(gpu_dycore_RG_history1h_${f} gpu_dycore_RG_history2h_${f} FIXTURES_REQUIRED gpu_dycore_RG)
endforeach()
add_test(NAME gpu_dycore_run_RG_streams COMMAND ${GPU} --clean ${W}/gpu_dycore_RG_streams -- ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/gpu_dycore_RG_streams' && mkdir -p '${W}/gpu_dycore_RG_streams' && '${GEXE}' run '${REF_RG}/init/x1.10242.init.nc' '${W}/gpu_dycore_RG_streams' --namelist '${REF_RG}/forecast/namelist.atmosphere' --streams '${REF_RG}/forecast/streams.atmosphere' --scalars qv,qc,qr,qi,qs,qg --every 0 --tables '${CPPMPAS_MPAS_TABLES}' --rrtmg '${CPPMPAS_CPPWRF_RRTMG}'")
cppmpas_gpu_props(gpu_dycore_run_RG_streams FIXTURES_SETUP gpu_dycore_RG_streams)
set_tests_properties(gpu_dycore_run_RG_streams PROPERTIES TIMEOUT 1800)
foreach(file IN LISTS RG_STREAM_FILES)
  add_test(NAME gpu_dycore_RG_streams_${file} COMMAND ${GPU} -- ${PY} ${CPPMPAS_MPAS_COMPARE}
    ${REF_RG}/forecast/${file}.nc ${W}/gpu_dycore_RG_streams/${file}.nc)
  cppmpas_gpu_props(gpu_dycore_RG_streams_${file} FIXTURES_REQUIRED gpu_dycore_RG_streams)
endforeach()
add_test(NAME gpu_dycore_RG_streams_all COMMAND ${GPU} -- ${PY} ${T}/streams_compare.py ${CPPMPAS_MPAS_COMPARE}
  ${REF_RG}/forecast ${W}/gpu_dycore_RG_streams ${W}/gpu_dycore_RG_streams.json 5 702)
cppmpas_gpu_props(gpu_dycore_RG_streams_all FIXTURES_REQUIRED gpu_dycore_RG_streams)
add_test(NAME gpu_dycore_run_long_RG24 COMMAND ${GPU} --clean ${W}/gpu_dycore_long_RG24 -- ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/gpu_dycore_long_RG24' && mkdir -p '${W}/gpu_dycore_long_RG24' && '${GEXE}' run '${REF_RG24}/init/x1.10242.init.nc' '${W}/gpu_dycore_long_RG24' --namelist '${REF_RG24}/forecast/namelist.atmosphere' --streams '${REF_RG24}/forecast/streams.atmosphere' --scalars qv,qc,qr,qi,qs,qg --every 0 --tables '${CPPMPAS_MPAS_TABLES}' --rrtmg '${CPPMPAS_CPPWRF_RRTMG}'")
cppmpas_gpu_props(gpu_dycore_run_long_RG24 FIXTURES_SETUP gpu_dycore_long_RG24 LABELS long)
set_tests_properties(gpu_dycore_run_long_RG24 PROPERTIES TIMEOUT 7200)
foreach(file IN LISTS RG24_STREAM_FILES)
  add_test(NAME gpu_dycore_RG24_${file} COMMAND ${GPU} -- ${PY} ${CPPMPAS_MPAS_COMPARE}
    ${REF_RG24}/forecast/${file}.nc ${W}/gpu_dycore_long_RG24/${file}.nc)
  cppmpas_gpu_props(gpu_dycore_RG24_${file} FIXTURES_REQUIRED gpu_dycore_long_RG24 LABELS long)
endforeach()
add_test(NAME gpu_dycore_RG24_all COMMAND ${GPU} -- ${PY} ${T}/streams_compare.py ${CPPMPAS_MPAS_COMPARE}
  ${REF_RG24}/forecast ${W}/gpu_dycore_long_RG24 ${W}/gpu_dycore_long_RG24.json 19 2168)
cppmpas_gpu_props(gpu_dycore_RG24_all FIXTURES_REQUIRED gpu_dycore_long_RG24 LABELS long)

# Isobaric diagnostics (core_atmosphere/diagnostics, mpas_isobaric_diagnostics.F) against the reference's
# diagnostics stream: at t = 0 against mpas842-dbl-rg-004's diag.2026-10-12_00.00.00.nc (the input state; relhum is
# the input file's) and after 18 steps against mpas842-dbl-rg24-001's diag.2026-10-12_03.00.00.nc (about 70 s). A
# --write-at file stands for the ringing diagnostics stream (cppmpas run), so the fields below are exactly those
# computed. relhum_*hPa at 03:00 interpolates the relhum the WSM6 driver computes when the diagnostics alarm rings
# (compute_relhum, l_diags), which only a --streams run has: it is left out here and compared by dycore_RG24_diag.*
# (at t = 0 it is the input file's relhum).
set(DIAG_ISO_FIELDS mslp t_isobaric t_iso_levels z_isobaric z_iso_levels meanT_500_300)
set(DIAG_ISO_LEVEL_FIELDS)
foreach(v relhum dewpoint temperature height uzonal umeridional w vorticity)
  foreach(l 50 100 200 250 500 700 850 925)
    list(APPEND DIAG_ISO_LEVEL_FIELDS ${v}_${l}hPa)
  endforeach()
endforeach()
set(DIAG_ISO_ALL ${DIAG_ISO_FIELDS} ${DIAG_ISO_LEVEL_FIELDS})
string(REPLACE ";" "," DIAG_ISO_LIST "${DIAG_ISO_ALL}")
add_test(NAME diag_run_RG_iso0 COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/diag_RG_iso0' && mkdir -p '${W}/diag_RG_iso0' && '${EXE}' run '${REF_RG}/init/x1.10242.init.nc' '${W}/diag_RG_iso0' --namelist '${REF_RG}/steps/namelist.atmosphere' --scalars qv,qc,qr,qi,qs,qg --steps 0 --every 0 --write-at 0 --write-fields ${DIAG_ISO_LIST} --tables '${CPPMPAS_MPAS_TABLES}' --rrtmg '${CPPMPAS_CPPWRF_RRTMG}'")
set_tests_properties(diag_run_RG_iso0 PROPERTIES FIXTURES_SETUP diag_RG_iso0 TIMEOUT 600)
foreach(f IN LISTS DIAG_ISO_ALL)
  add_test(NAME diag_RG_iso0_${f} COMMAND ${PY} ${T}/compare_arrays.py
    ${REF_RG}/forecast/diag.2026-10-12_00.00.00.nc ${W}/diag_RG_iso0/state.2026-10-12_00.00.00.nc ${f})
  set_tests_properties(diag_RG_iso0_${f} PROPERTIES FIXTURES_REQUIRED diag_RG_iso0)
endforeach()
add_test(NAME diag_run_RG24_iso3h COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/diag_RG24_iso3h' && mkdir -p '${W}/diag_RG24_iso3h' && '${EXE}' run '${REF_RG24}/init/x1.10242.init.nc' '${W}/diag_RG24_iso3h' --namelist '${REF_RG24}/steps/namelist.atmosphere' --scalars qv,qc,qr,qi,qs,qg --steps 18 --every 0 --write-at 18 --write-fields ${DIAG_ISO_LIST} --tables '${CPPMPAS_MPAS_TABLES}' --rrtmg '${CPPMPAS_CPPWRF_RRTMG}'")
set_tests_properties(diag_run_RG24_iso3h PROPERTIES FIXTURES_SETUP diag_RG24_iso3h TIMEOUT 1800)
foreach(f IN LISTS DIAG_ISO_ALL)
  if(f MATCHES "^relhum_")
    continue()
  endif()
  add_test(NAME diag_RG24_iso3h_${f} COMMAND ${PY} ${T}/compare_arrays.py
    ${REF_RG24}/forecast/diag.2026-10-12_03.00.00.nc ${W}/diag_RG24_iso3h/state.2026-10-12_03.00.00.nc ${f})
  set_tests_properties(diag_RG24_iso3h_${f} PROPERTIES FIXTURES_REQUIRED diag_RG24_iso3h)
endforeach()
# ---------------------------------------------------------------------------------------------------------------
# Convective diagnostics (core_atmosphere/diagnostics, mpas_convective_diagnostics.F) against the reference's
# diagnostics stream: at t = 0 (mpas842-dbl-rg-004, diag.2026-10-12_00.00.00.nc; reset, update and compute before the
# first write, mpas_atm_core.F:692-694) and at 3 h and 6 h of RG24 (mpas842-dbl-rg24-001; the maxima since the
# previous diagnostics output, the 00:00 write resetting them). A --write-at file stands for the ringing diagnostics
# stream, so 0 is in the RG24 list: its write resets the t = 0 update as MPAS's does.
set(DIAG_CONV_FIELDS cape cin lcl lfc srh_0_1km srh_0_3km uzonal_surface uzonal_1km uzonal_6km umeridional_surface
    umeridional_1km umeridional_6km temperature_surface dewpoint_surface updraft_helicity_max w_velocity_max
    wind_speed_level1_max)
string(REPLACE ";" "," DIAG_CONV_LIST "${DIAG_CONV_FIELDS}")
add_test(NAME diag_run_RG_convective COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/diag_RG_conv' && mkdir -p '${W}/diag_RG_conv' && '${EXE}' run '${REF_RG}/init/x1.10242.init.nc' '${W}/diag_RG_conv' --namelist '${REF_RG}/steps/namelist.atmosphere' --scalars qv,qc,qr,qi,qs,qg --steps 0 --every 0 --write-at 0 --write-fields ${DIAG_CONV_LIST} --tables '${CPPMPAS_MPAS_TABLES}' --rrtmg '${CPPMPAS_CPPWRF_RRTMG}'")
set_tests_properties(diag_run_RG_convective PROPERTIES FIXTURES_SETUP diag_RG_conv TIMEOUT 600)
foreach(f IN LISTS DIAG_CONV_FIELDS)
  add_test(NAME diag_RG_conv0_${f} COMMAND ${PY} ${T}/compare_arrays.py
    ${REF_RG}/forecast/diag.2026-10-12_00.00.00.nc ${W}/diag_RG_conv/state.2026-10-12_00.00.00.nc ${f})
  set_tests_properties(diag_RG_conv0_${f} PROPERTIES FIXTURES_REQUIRED diag_RG_conv)
endforeach()
add_test(NAME diag_run_long_RG24_convective COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/diag_RG24_conv' && mkdir -p '${W}/diag_RG24_conv' && '${EXE}' run '${REF_RG24}/init/x1.10242.init.nc' '${W}/diag_RG24_conv' --namelist '${REF_RG24}/steps/namelist.atmosphere' --scalars qv,qc,qr,qi,qs,qg --steps 36 --every 0 --write-at 0,18,36 --write-fields ${DIAG_CONV_LIST} --tables '${CPPMPAS_MPAS_TABLES}' --rrtmg '${CPPMPAS_CPPWRF_RRTMG}'")
set_tests_properties(diag_run_long_RG24_convective PROPERTIES FIXTURES_SETUP diag_RG24_conv LABELS long TIMEOUT 3600)
foreach(stamp IN ITEMS 2026-10-12_03.00.00 2026-10-12_06.00.00)
  foreach(f IN LISTS DIAG_CONV_FIELDS)
    add_test(NAME diag_RG24_conv_${stamp}_${f} COMMAND ${PY} ${T}/compare_arrays.py
      ${REF_RG24}/forecast/diag.${stamp}.nc ${W}/diag_RG24_conv/state.${stamp}.nc ${f})
    set_tests_properties(diag_RG24_conv_${stamp}_${f} PROPERTIES FIXTURES_REQUIRED diag_RG24_conv LABELS long)
  endforeach()
endforeach()

# The MYNN PBL (bl_mynn with MPAS's default config_mynn_* options) and its MPAS driver, with MPAS_to_physics and
# physics_get_tend, against MPAS's own code built with the reference compiler and flags (tests/mynn_column/
# run_check.sh): gfortran 14.2, the MPAS v8.4.2 source tree (CPPMPAS_MPAS_SRC) and MMM-physics a4baf7f
# (CPPMPAS_MMM_PHYSICS); without them the reference build fails and so does the test.
add_executable(test_mynn_column mynn_column/test_mynn_column.cpp)
target_link_libraries(test_mynn_column PRIVATE cppmpas)
add_test(NAME mynn_column COMMAND ${T}/mynn_column/run_check.sh ${CPPMPAS_MPAS_SRC} ${CPPMPAS_MMM_PHYSICS}
  ${W}/mynn_column $<TARGET_FILE:test_mynn_column>)
set_tests_properties(mynn_column PROPERTIES FIXTURES_REQUIRED work_dir)
