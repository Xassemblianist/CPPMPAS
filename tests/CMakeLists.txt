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
    "'${EXE}' dump '${${case}}' '${W}/dump_${case}.nc' && ${Python3_EXECUTABLE} -I '${T}/check_order.py' '${${case}}' '${W}/dump_${case}.nc' ${case}")
endforeach()
add_test(NAME order_init_MW COMMAND ${CMAKE_COMMAND} -E env sh -c
  "'${EXE}' dump '${W}/mw_init.nc' '${W}/dump_mw_init.nc' --init && ${Python3_EXECUTABLE} -I '${T}/check_order.py' '${W}/mw_init.nc' '${W}/dump_mw_init.nc' MW --init qv")
add_test(NAME order_init_SC COMMAND ${CMAKE_COMMAND} -E env sh -c
  "'${EXE}' dump '${W}/sc_init.nc' '${W}/dump_sc_init.nc' --init --scalars qv,qc,qr && ${Python3_EXECUTABLE} -I '${T}/check_order.py' '${W}/sc_init.nc' '${W}/dump_sc_init.nc' SC --init qv,qc,qr")
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
cppmpas_error_test(missing_scalar "active scalar qr is missing" info --init --scalars qv,qc,qr ${W}/err_missing_qr.nc)
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
  "${Python3_EXECUTABLE} -I '${CMAKE_SOURCE_DIR}/tools/gen_namelist_table.py' '${CPPMPAS_MPAS_REGISTRY}' ${CPPMPAS_REFS_DIR}/mpas842-dbl-mw-001/MW/steps/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-mw-001/MW/forecast/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-sc-001/SC/steps/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-sc-001/SC/forecast/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-jw-001/JW/steps/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-jw-001/JW/forecast/namelist.atmosphere '${CMAKE_SOURCE_DIR}/tests/namelists/MW_restart/namelist.atmosphere' > '${W}/namelist_table.inc' && ${CMAKE_COMMAND} -E compare_files '${W}/namelist_table.inc' '${CMAKE_SOURCE_DIR}/src/core/namelist_table.inc'")
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
    "mkdir -p '${W}' && '${EXE}' dump '${INIT_${case}}' '${W}/dump_refinit_${case}.nc' --init --scalars ${SCAL_${case}} && ${Python3_EXECUTABLE} -I '${T}/check_order.py' '${INIT_${case}}' '${W}/dump_refinit_${case}.nc' ${case} --init ${SCAL_${case}}")
  # the reader on the reference restart file: its input-stream subset, values and in-memory order
  add_test(NAME readrestart_${case} COMMAND ${CMAKE_COMMAND} -E env sh -c
    "mkdir -p '${W}' && '${EXE}' dump '${RST_${case}}' '${W}/dump_rst_${case}.nc' --init --scalars ${SCAL_${case}} && ${Python3_EXECUTABLE} -I '${T}/check_order.py' '${RST_${case}}' '${W}/dump_rst_${case}.nc' ${case} --init ${SCAL_${case}}")
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
cppmpas_error_test(dycore_SC "config_scalar_advection = true" run ${INIT_SC} ${W}/x_dir
  --namelist ${REF_SC}/steps/namelist.atmosphere --scalars qv,qc,qr --steps 1)
cppmpas_error_test(dycore_odd_levels "nVertLevels = 3 (only an even nVertLevels >= 4 is verified)" run ${W}/sc_init.nc
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

# ---------------------------------------------------------------------------------------------------------------
# GPU proof: the CUDA build of the same sources (cppmpas_gpu, CMake option CPPMPAS_ENABLE_CUDA; every kernel on the
# GPU, core/backend.hpp) against the same references at tolerance 0, and against the CPU build's own output
# (cppmpas, the tests above). Each test runs through tests/gpu_test.py: without a CUDA build or without a usable
# device it exits 77 and CTest reports it as "Skipped" (label gpu), never as passed.
if(CPPMPAS_ENABLE_CUDA)
  set(GEXE $<TARGET_FILE:cppmpas_gpu>)
else()
  set(GEXE NONE)  # tests/gpu_test.py: no CUDA build, every GPU test is skipped
endif()
set(GPU ${PY} ${T}/gpu_test.py ${GEXE})
function(cppmpas_gpu_props)  # cppmpas_gpu_props(test... [FIXTURES_REQUIRED f] [FIXTURES_SETUP f] [LABELS l])
  cmake_parse_arguments(G "" "" "FIXTURES_REQUIRED;FIXTURES_SETUP;LABELS" ${ARGN})
  set_tests_properties(${G_UNPARSED_ARGUMENTS} PROPERTIES SKIP_RETURN_CODE 77 LABELS "gpu;${G_LABELS}")
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
  # CMAKE_CUDA_FLAGS_RELEASE, in nvcc's spellings and as host-compiler flags in any form
  set(_cuda_bad "--use_fast_math" "-fmad=true" "--fmad true" "--ftz=true" "--prec-div=false" "--prec-sqrt=false"
      "-Xcompiler -ffast-math" "--compiler-options=-O2,-ffp-contract=fast" "-Xcompiler=-funsafe-math-optimizations"
      "-Xcompiler -mfpmath=387")
  set(_i 0)
  foreach(_flag IN LISTS _cuda_bad)
    foreach(_var CMAKE_CUDA_FLAGS CMAKE_CUDA_FLAGS_RELEASE)
      add_test(NAME cmake_rejects_cuda_${_i}_${_var} COMMAND ${CMAKE_COMMAND} -E env sh -c
        "rm -rf '${W}/cmake_cuda_${_i}_${_var}'; out=$(${CMAKE_COMMAND} -S '${CMAKE_SOURCE_DIR}' -B '${W}/cmake_cuda_${_i}_${_var}' -DCMAKE_BUILD_TYPE=Release -DCPPMPAS_ENABLE_CUDA=ON -DCMAKE_CUDA_COMPILER='${CMAKE_CUDA_COMPILER}' -DCMAKE_CUDA_HOST_COMPILER='${CMAKE_CUDA_HOST_COMPILER}' -DCMAKE_CUDA_ARCHITECTURES='${CMAKE_CUDA_ARCHITECTURES}' '-D${_var}=${_flag}' 2>&1); rc=$?; rm -rf '${W}/cmake_cuda_${_i}_${_var}'; echo \"$out\" | tail -5; [ $rc -ne 0 ] && echo \"$out\" | grep -q 'changes the floating-point results'")
      set_tests_properties(cmake_rejects_cuda_${_i}_${_var} PROPERTIES FIXTURES_REQUIRED work_dir LABELS gpu)
    endforeach()
    math(EXPR _i "${_i} + 1")
  endforeach()
  # and the strict flags follow the user's: a CUDA source compiles with them last on nvcc's command line
  add_test(NAME cmake_cuda_strict_flags_last COMMAND ${CMAKE_COMMAND} -E env sh -c
    "grep -h -- '--fmad=false' '${CMAKE_BINARY_DIR}/CMakeFiles/cppmpas_cuda.dir/flags.make' | grep -q -- '-O3 -DNDEBUG.*--fmad=false --ftz=false --prec-div=true --prec-sqrt=true -Xcompiler=-ffp-contract=off'")
  set_tests_properties(cmake_cuda_strict_flags_last PROPERTIES LABELS gpu)
else()
  add_test(NAME gpu_backend COMMAND ${GPU} device)
endif()
cppmpas_gpu_props(gpu_backend)
# atm_mpas_init_block on the GPU: every initdump field of MW, SC and JW against MPAS, and against the CPU build
foreach(case MW SC JW)
  add_test(NAME gpu_initdump_run_${case} COMMAND ${GPU} -- ${CMAKE_COMMAND} -E env sh -c
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
  add_test(NAME gpu_init_run_${case} COMMAND ${GPU} -- ${CMAKE_COMMAND} -E env sh -c
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
add_test(NAME gpu_dycore_run_MW COMMAND ${GPU} -- ${CMAKE_COMMAND} -E env sh -c
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
add_test(NAME gpu_dycore_run_long_MW COMMAND ${GPU} -- ${CMAKE_COMMAND} -E env sh -c
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
add_test(NAME gpu_dycore_run_60min_MW COMMAND ${GPU} -- ${CMAKE_COMMAND} -E env sh -c
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
add_test(NAME gpu_dycore_run_JW COMMAND ${GPU} -- ${CMAKE_COMMAND} -E env sh -c
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
add_test(NAME gpu_dycore_run_long_JW COMMAND ${GPU} -- ${CMAKE_COMMAND} -E env sh -c
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
add_test(NAME gpu_dycore_restart0_run_MW COMMAND ${GPU} -- ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/gpu_dycore_restart0_MW' && mkdir -p '${W}/gpu_dycore_restart0_MW' && '${GEXE}' run '${RST_MW}' '${W}/gpu_dycore_restart0_MW' --namelist '${NML_RST_MW}' --steps 0 --every 0 --write-at 0 --write-fields ${DYCORE_RESTART0_LIST}")
cppmpas_gpu_props(gpu_dycore_restart0_run_MW FIXTURES_SETUP gpu_dycore_restart0_MW)
foreach(f IN LISTS DYCORE_RESTART0_FIELDS)
  add_test(NAME gpu_dycore_restart0_MW_${f} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py ${RST_MW}
    ${W}/gpu_dycore_restart0_MW/state.0000-01-01_00.30.06.nc ${f})
  cppmpas_gpu_props(gpu_dycore_restart0_MW_${f} FIXTURES_REQUIRED gpu_dycore_restart0_MW)
endforeach()
add_test(NAME gpu_dycore_restart_run_long_MW COMMAND ${GPU} -- ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/gpu_dycore_restart_long_MW' && mkdir -p '${W}/gpu_dycore_restart_long_MW' && '${GEXE}' run '${RST_MW}' '${W}/gpu_dycore_restart_long_MW' --namelist '${NML_RST_MW}' --steps 299 --every 0 --write-at 299 --write-fields ${DYCORE_OUTPUT_LIST}")
cppmpas_gpu_props(gpu_dycore_restart_run_long_MW FIXTURES_SETUP gpu_dycore_restart_long_MW LABELS long)
set_tests_properties(gpu_dycore_restart_run_long_MW PROPERTIES TIMEOUT 1800)
foreach(f IN LISTS DYCORE_OUTPUT_FIELDS)
  add_test(NAME gpu_dycore_restart_MW_output60_${f} COMMAND ${GPU} -- ${PY} ${T}/compare_arrays.py
    ${REF_MW}/forecast/output.nc ${W}/gpu_dycore_restart_long_MW/state.0000-01-01_01.00.00.nc ${f} --record 2)
  cppmpas_gpu_props(gpu_dycore_restart_MW_output60_${f} FIXTURES_REQUIRED gpu_dycore_restart_long_MW LABELS long)
endforeach()
