# SPDX-License-Identifier: GPL-3.0-or-later
# The Debian 13 image of the CPU suite (.github/workflows/cpu_suite.yml): the base image and the APT snapshot of the
# MPAS reference builds (Xassemblianist/cppmpas-data ci/Dockerfile: the same digest, the same dated snapshot with its
# InRelease pinned by SHA-256), so glibc 2.41 (libm, libmvec), gcc / g++ / gfortran 14.2 and netCDF are the
# reference's. cmake, pkgconf and locales come from the same snapshot. The de_DE.UTF-8 locale is for
# unit_namelist_locale (a locale whose decimal point is ',').
FROM docker.io/library/debian@sha256:ab3b3c8c5159eebfca4376430132d53c4d425f869c765ecbe47177a030da575d
COPY debian-packages.txt /opt/debian-packages.txt
RUN set -eu; \
    export DEBIAN_FRONTEND=noninteractive; \
    rm -f /etc/apt/sources.list /etc/apt/sources.list.d/*; \
    echo 'deb [check-valid-until=no signed-by=/usr/share/keyrings/debian-archive-keyring.gpg] http://snapshot.debian.org/archive/debian/20261123T000000Z/ trixie main' \
        > /etc/apt/sources.list; \
    printf '%s\n' 'Acquire::Languages "none";' 'Acquire::Retries "5";' > /etc/apt/apt.conf.d/99snapshot; \
    apt-get -o APT::Update::Error-Mode=any update; \
    echo 'bfacad4ca61899ce20d4a0d439f29339f66541561b67506364efb88e27dee870  /var/lib/apt/lists/snapshot.debian.org_archive_debian_20261123T000000Z_dists_trixie_InRelease' \
        | sha256sum -c -; \
    xargs apt-get install -y --no-install-recommends < /opt/debian-packages.txt; \
    sed -i 's/^# *de_DE.UTF-8 UTF-8/de_DE.UTF-8 UTF-8/' /etc/locale.gen; \
    locale-gen; \
    dpkg-query -W -f='${binary:Package}\t${Version}\n' | LC_ALL=C sort > /opt/packages.tsv; \
    apt-get clean; \
    rm -rf /var/lib/apt/lists/*
ENV LC_ALL=C TZ=UTC OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
