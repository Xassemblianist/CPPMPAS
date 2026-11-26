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
  "mkdir -p '${W}' && '$<TARGET_FILE:test_io_core>' '${W}' '${MW}'")

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
cppmpas_error_test(unported_field "variable 'pressure' has no CPPMPAS port" history ${W}/mw_init.nc ${W}/x.nc --fields u,pressure)
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

# src/core/namelist_table.inc is what tools/gen_namelist_table.py makes of the Registry and the reference namelists
add_test(NAME namelist_table_current COMMAND ${CMAKE_COMMAND} -E env sh -c
  "${Python3_EXECUTABLE} -I '${CMAKE_SOURCE_DIR}/tools/gen_namelist_table.py' '${CPPMPAS_MPAS_REGISTRY}' ${CPPMPAS_REFS_DIR}/mpas842-dbl-mw-001/MW/steps/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-mw-001/MW/forecast/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-sc-001/SC/steps/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-sc-001/SC/forecast/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-jw-001/JW/steps/namelist.atmosphere ${CPPMPAS_REFS_DIR}/mpas842-dbl-jw-001/JW/forecast/namelist.atmosphere > '${W}/namelist_table.inc' && ${CMAKE_COMMAND} -E compare_files '${W}/namelist_table.inc' '${CMAKE_SOURCE_DIR}/src/core/namelist_table.inc'")
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
add_test(NAME unit_dycore COMMAND test_dycore ${INIT_MW} ${REF_MW}/steps/namelist.atmosphere)
# configurations of the other references, and meshes, that the dycore refuses before it writes anything
cppmpas_error_test(dycore_JW "config_horiz_mixing = '2d_smagorinsky' has no verified port" run ${INIT_JW} ${W}/x_dir
  --namelist ${REF_JW}/steps/namelist.atmosphere --steps 1)
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
cppmpas_error_test(dycore_fEdge "fEdge != 0" run ${INIT_JW} ${W}/x_dir --namelist ${REF_MW}/steps/namelist.atmosphere
  --steps 1)
cppmpas_error_test(dycore_odd_levels "nVertLevels = 3 (only an even nVertLevels >= 4 is verified)" run ${W}/sc_init.nc
  ${W}/x_dir --namelist ${REF_MW}/steps/namelist.atmosphere --steps 1)
cppmpas_error_test(dycore_steps "--steps needs a whole number 0..1000000, not '1x'" run ${INIT_MW} ${W}/x_dir
  --namelist ${REF_MW}/steps/namelist.atmosphere --steps 1x)
cppmpas_error_test(dycore_no_steps "usage: cppmpas run" run ${INIT_MW} ${W}/x_dir
  --namelist ${REF_MW}/steps/namelist.atmosphere)
# 301 steps (label "long", about 2 minutes): the forecast run of the same reference writes output.nc at 30 min
# (record 1, after step 300) and restart.0000-01-01_00.30.06.nc (after step 301); every field of these streams the
# port holds is compared at tolerance 0
set(DYCORE_RESTART_FIELDS xtime Time u w rho_zz theta_m pressure_p rho theta exner ru rw rtheta_p rho_p qv)
set(DYCORE_OUTPUT_FIELDS xtime u w rho theta divergence vorticity ke uReconstructZonal uReconstructMeridional qv)
set(DYCORE_LONG_FIELDS ${DYCORE_RESTART_FIELDS} divergence vorticity ke uReconstructZonal uReconstructMeridional)
list(REMOVE_DUPLICATES DYCORE_LONG_FIELDS)
string(REPLACE ";" "," DYCORE_LONG_LIST "${DYCORE_LONG_FIELDS}")
add_test(NAME dycore_run_long_MW COMMAND ${CMAKE_COMMAND} -E env sh -c
  "rm -rf '${W}/dycore_long_MW' && mkdir -p '${W}/dycore_long_MW' && '${EXE}' run '${INIT_MW}' '${W}/dycore_long_MW' --namelist '${REF_MW}/steps/namelist.atmosphere' --steps 301 --every 0 --write-at 300,301 --write-fields ${DYCORE_LONG_LIST}")
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
endforeach()
# SC, 601 steps (label "long", about 45 minutes): the forecast run writes output.nc at 30 min (record 1, after step
# 600) and restart.0000-01-01_00.30.03.nc (after step 601); surface_pressure is microphysics_to_MPAS's (PI:950-964)
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
  add_test(NAME dycore_SC_output30_${f} COMMAND ${PY} ${T}/compare_arrays.py ${REF_SC}/forecast/output.nc
    ${W}/dycore_long_SC/state.0000-01-01_00.30.00.nc ${f} --record 1)
  set_tests_properties(dycore_SC_output30_${f} PROPERTIES FIXTURES_REQUIRED dycore_long_SC LABELS long)
endforeach()
cppmpas_error_test(dycore_write_at "--write-at 12 is after the last step" run ${INIT_MW} ${W}/x_dir
  --namelist ${REF_MW}/steps/namelist.atmosphere --steps 10 --write-at 12 --write-fields u)
cppmpas_error_test(dycore_write_fields "variable 'pressure' has no CPPMPAS port" run ${INIT_MW} ${W}/x_dir
  --namelist ${REF_MW}/steps/namelist.atmosphere --steps 10 --write-at 1 --write-fields u,pressure)
