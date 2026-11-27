# SPDX-License-Identifier: GPL-3.0-or-later
"""Compares the MYNN surface layer's object code in two linked executables, function by function:

  objcmp.py <reference bin/atmosphere_model> <sfmynn_column_ref of build_ref.sh>

Each function is disassembled (objdump -d) and normalised: addresses and jump targets dropped, call targets by name
(log@plt, _ZGVbN2vv_pow@plt, a local function), and every rip-relative operand replaced by what it reads: the bytes of
a constant (8 or 16, from the file) or the symbol of a variable. The schemes' own functions must be identical. The
driver's routines (they call the MPAS framework's pools, which the column program stubs) are compared by their
floating-point opcode counts and their calls of libm / libmvec. Also reports FMA (vfmadd...) and ymm registers.
"""
import re
import subprocess
import sys
from collections import Counter

SCHEME = ["__sf_mynn_MOD_sf_mynn_run", "__sf_mynn_MOD_sf_mynn_init", "__sf_mynn_MOD_sf_mynn_finalize",
          "__sf_mynn_MOD_zolrib.constprop.0.isra.0", "__sf_mynn_MOD_li_etal_2010.isra.0",
          "__sf_mynn_MOD_psim_stable.constprop.0.isra.0", "__sf_mynn_MOD_psih_stable.constprop.0.isra.0",
          "__sf_mynn_MOD_psim_unstable.constprop.0.isra.0", "__sf_mynn_MOD_psih_unstable.constprop.0.isra.0",
          "__sf_mynn_MOD_psim_unstable_full.isra.0", "__sf_mynn_MOD_zilitinkevich_1995.isra.0",
          "__sf_mynn_pre_MOD_sf_mynn_pre_run", "__module_sf_mynn_MOD_sfclay_mynn"]
DRIVER = ["__mpas_atmphys_driver_sfclayer_MOD_driver_sfclayer", "__mpas_atmphys_driver_sfclayer_MOD_init_sfclayer",
          "__mpas_atmphys_driver_sfclayer_MOD_allocate_sfclayer",
          "__mpas_atmphys_driver_sfclayer_MOD_deallocate_sfclayer", "sfclayer_from_mpas"]
FP = re.compile(r"^(v?(add|sub|mul|div|sqrt|max|min|and|andn|or|xor|cvt\w*|ucomi|comi|unpck\w*)(s|p)[sd]|vfm\w+)$")


def sections(exe):
    out = subprocess.run(["readelf", "-SW", exe], capture_output=True, text=True, check=True).stdout
    secs = []
    for l in out.splitlines():
        m = re.match(r"\s*\[\s*\d+\]\s+(\S+)\s+\S+\s+([0-9a-f]+)\s+([0-9a-f]+)\s+([0-9a-f]+)", l)
        if m:
            secs.append((int(m.group(2), 16), int(m.group(3), 16), int(m.group(4), 16), m.group(1)))
    return secs


def section_of(secs, addr):
    for a, off, size, name in secs:
        if a and a <= addr < a + size:
            return name
    return None


def read_at(exe, secs, addr, n):
    for a, off, size, name in secs:
        if a and a <= addr < a + size:
            with open(exe, "rb") as f:
                f.seek(off + addr - a)
                return f.read(n).hex()
    return None


def functions(exe):
    syms = subprocess.run(["nm", exe], capture_output=True, text=True, check=True).stdout
    names = {}
    for l in syms.splitlines():
        p = l.split()
        if len(p) == 3 and p[1] in "tT":
            names[p[2]] = int(p[0], 16)
    return names


def dis(exe, secs, sym):
    out = subprocess.run(["objdump", "-d", "--no-show-raw-insn", "--disassemble=" + sym, exe], capture_output=True,
                         text=True).stdout
    ins, on = [], False
    for l in out.splitlines():
        if re.match(r"^[0-9a-f]+ <" + re.escape(sym) + ">:", l):
            on = True
            continue
        if not on:
            continue
        if not l.strip():
            if ins:
                break
            continue
        m = re.match(r"\s*[0-9a-f]+:\s*(.*)", l)
        if not m:
            continue
        s = m.group(1)
        if re.match(r"(nop|xchg %ax,%ax|cs nopw|data16)", s):
            continue
        rip = re.search(r"# ([0-9a-f]+) <([^>]*)>", s)
        s = re.sub(r"\s+#.*", "", s)
        call = re.match(r"(call|jmp)\s+[0-9a-f]+ <([^>]*)>", s)
        if call and not re.match(re.escape(sym) + r"\+", call.group(2)):
            s = call.group(1) + " " + re.sub(r"\+0x[0-9a-f]+$", "", call.group(2))
        else:
            s = re.sub(r"\b(call|jmp|j[a-z]+)\s+[0-9a-f]+( <[^>]*>)?", r"\1 T", s)
        if rip:
            op = s.split()[0]
            target = rip.group(2)
            if (section_of(secs, int(rip.group(1), 16)) or "").startswith(".rodata"):
                n = 16 if re.match(r"(movap|movdqa|andp|andnp|orp|xorp|mulp|addp|subp|divp|maxp|minp|movups|movdqu)", op) else 8
                val = read_at(exe, secs, int(rip.group(1), 16), n)
                s = re.sub(r"-?0x[0-9a-f]+\(%rip\)", f"[{val}]", s)
            else:
                s = re.sub(r"-?0x[0-9a-f]+\(%rip\)", "[" + re.sub(r"\+0x[0-9a-f]+$", "", target) + "]", s)
        ins.append(re.sub(r"\s+", " ", s).strip())
    return ins


def summary(ins):
    fp = Counter(i.split()[0] for i in ins if FP.match(i.split()[0]))
    calls = Counter(i.split()[1] for i in ins
                    if i.startswith("call ") and re.match(r"(log|exp|pow|atan|sqrt|_ZGV\w+)@plt$", i.split()[1]))
    return fp, calls


def main():
    ref, col = sys.argv[1], sys.argv[2]
    rs, cs = sections(ref), sections(col)
    rf, cf = functions(ref), functions(col)
    bad = 0
    for sym in SCHEME + DRIVER:
        rn = [n for n in rf if n == sym or (sym == "sfclayer_from_mpas" and "sfclayer_from_mpas" in n)]
        cn = [n for n in cf if n == sym or (sym == "sfclayer_from_mpas" and "sfclayer_from_mpas" in n)]
        if not rn or not cn:
            print(f"{sym}: missing (reference {rn}, column {cn})")
            bad += 1
            continue
        a, b = dis(ref, rs, rn[0]), dis(col, cs, cn[0])
        fma = sum(1 for i in a if i.startswith("vfm")) + sum(1 for i in b if i.startswith("vfm"))
        ymm = sum(1 for i in a + b if "%ymm" in i)
        fa, ca = summary(a)
        fb, cb = summary(b)
        if sym in SCHEME:
            d = sum(1 for x, y in zip(a, b) if x != y) + abs(len(a) - len(b))
            print(f"{rn[0]}: {len(a)} / {len(b)} instructions, {d} differ, FMA {fma}, ymm {ymm}, calls {dict(ca)}")
            if d:
                bad += 1
                for x, y in zip(a, b):
                    if x != y:
                        print(f"   ref {x}\n   col {y}")
                        break
        else:
            same = fa == fb and ca == cb
            print(f"{rn[0]}: {len(a)} / {len(b)} instructions, FP opcodes {'equal' if fa == fb else 'DIFFER'} "
                  f"({sum(fa.values())}), libm calls {'equal' if ca == cb else 'DIFFER'} {dict(ca)}, FMA {fma}, "
                  f"ymm {ymm}")
            if not same:
                bad += 1
                print(f"   ref {dict(fa)} {dict(ca)}\n   col {dict(fb)} {dict(cb)}")
    print("object code:", "IDENTICAL" if not bad else f"{bad} functions differ")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
