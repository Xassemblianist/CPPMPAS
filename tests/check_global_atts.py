#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Global attributes as MPAS's bootstrapping reads them (mpas_bootstrapping.F:313-339):

* a grid file without file_id / parent_id / mesh_spec gives MPAS's three warnings;
* parent_id = trim(file_id) // NEW_LINE // trim(parent_id), assigned to domain % parent_id, a
  character(len=StrKIND = 512) variable (mpas_domain_types.inc:33): the joined value is cut at 512 characters, and the
  history file carries that value (mpas_subdriver.F:529);
* a blank file_id passes the input's parent_id on unchanged (mpas_bootstrapping.F:327-331);
* NUL bytes: SMIOLf_inquire_att keeps them (smiolf.F90:1632-1633) and Fortran's trim does not remove them, so a
  file_id of NULs is not blank; the history attribute is written as strlen of the C string (smiol.c:1731), which
  ends at the first NUL.

usage: check_global_atts.py <cppmpas> <grid without ids> <init with a long parent_id> <workdir>
"""
import subprocess
import sys
from pathlib import Path

import ctypes
import ctypes.util
import shutil

import netCDF4 as nc

# raw text attributes through the NetCDF C library (netCDF4-python drops NUL bytes from text attributes)
_lib = ctypes.CDLL(ctypes.util.find_library("netcdf"))


def put_text_att(path, name, value):
    ncid = ctypes.c_int()
    if _lib.nc_open(str(path).encode(), 1, ctypes.byref(ncid)) or _lib.nc_redef(ncid):  # NC_WRITE
        raise SystemExit(f"cannot open {path} for writing")
    v = value.encode("latin-1")
    rc = _lib.nc_put_att_text(ncid, -1, name.encode(), ctypes.c_size_t(len(v)), v)  # NC_GLOBAL
    if _lib.nc_close(ncid) or rc:
        raise SystemExit(f"cannot write {name} to {path}")


def get_text_att(path, name):
    ncid, n = ctypes.c_int(), ctypes.c_size_t()
    if _lib.nc_open(str(path).encode(), 0, ctypes.byref(ncid)):
        raise SystemExit(f"cannot open {path}")
    rc = _lib.nc_inq_attlen(ncid, -1, name.encode(), ctypes.byref(n))  # an absent attribute is an error, not ""
    buf = ctypes.create_string_buffer(n.value + 1)
    if rc == 0:
        rc = _lib.nc_get_att_text(ncid, -1, name.encode(), buf)
    if _lib.nc_close(ncid) or rc:
        raise SystemExit(f"cannot read the text attribute {name} of {path} (NetCDF status {rc})")
    return buf.raw[: n.value].decode("latin-1")


def main():
    exe, grid, init, work = sys.argv[1:5]
    fails = []
    r = subprocess.run([exe, "info", grid], capture_output=True, text=True)
    if r.returncode != 0:
        fails.append(f"info {grid} failed: {r.stderr.strip()}")
    for a in ("file_id", "parent_id", "mesh_spec"):
        if f"MPAS warning: Attribute {a} not found" not in r.stderr:
            fails.append(f"no warning for the missing {a}")
    with nc.Dataset(grid) as g:
        if any(a in g.ncattrs() for a in ("file_id", "parent_id", "mesh_spec")):
            fails.append(f"{grid} has file_id / parent_id / mesh_spec: not a test of the warnings")

    out = Path(work) / "long_parent.history.nc"
    out.unlink(missing_ok=True)
    subprocess.run([exe, "history", init, str(out), "--fields", "latCell", "--file-id", "0123456789"], check=True)
    with nc.Dataset(init) as s, nc.Dataset(out) as h:
        joined = s.getncattr("file_id").rstrip(" ") + "\n" + s.getncattr("parent_id").rstrip(" ")
        if len(joined) <= 512:
            fails.append(f"{init}: joined parent_id has {len(joined)} characters, not a test of the cut")
        want = joined[:512].rstrip(" ")
        got = h.getncattr("parent_id")
        if got != want:
            fails.append(f"parent_id has {len(got)} characters, MPAS writes {len(want)}")
    # MPAS_io_put_att_text (mpas_io.F:6036-6045, 6208-6210): trim, cut at StrKIND = 512, trim again: a file_id longer
    # than 512 is cut, and blanks the cut leaves at the end are trimmed
    for fid, want in (("a" * 600, "a" * 512), ("a" * 500 + " " * 20 + "x", "a" * 500), ("ab   ", "ab")):
        out = Path(work) / "long_id.history.nc"
        out.unlink(missing_ok=True)
        subprocess.run([exe, "history", init, str(out), "--fields", "latCell", "--file-id", fid], check=True)
        got = get_text_att(out, "file_id")
        if got != want:
            fails.append(f"file_id of {len(fid)} characters is written as {len(got)} characters {got[-3:]!r}, "
                         f"MPAS writes {len(want)}")
        out.unlink()
    # file_id with NULs: (file_id, parent_id) of the input -> parent_id of the history file. A blank file_id
    # (trim(c_file_id) == '', mpas_bootstrapping.F:327-331) passes the input's parent_id on, trimmed when written
    for fid, pid, want in (("abc\0def", "p", "abc"), ("\0", "p", ""), ("ab ", "p", "ab\np"), ("", "p  ", "p"),
                           ("   ", "p", "p")):
        src = Path(work) / "nul_id.nc"
        out = Path(work) / "nul_id.history.nc"
        out.unlink(missing_ok=True)
        shutil.copyfile(init, src)
        put_text_att(src, "file_id", fid)
        put_text_att(src, "parent_id", pid)
        if get_text_att(src, "file_id") != fid:
            fails.append(f"file_id {fid!r} was stored as {get_text_att(src, 'file_id')!r}: not a test of NULs")
        subprocess.run([exe, "history", str(src), str(out), "--fields", "latCell", "--file-id", "0123456789"], check=True)
        got = get_text_att(out, "parent_id")
        if got != want:
            fails.append(f"file_id {fid!r}: parent_id {got!r}, MPAS writes {want!r}")
        src.unlink()
        out.unlink()
    for f in fails:
        print("FAIL", f)
    if fails:
        return 1
    print("OK global attributes: warnings for missing ids, parent_id cut at StrKIND = 512, NULs as SMIOL")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
