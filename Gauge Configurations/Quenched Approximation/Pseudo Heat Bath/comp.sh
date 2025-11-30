#!/bin/sh
# ファイル名: comp.sh

# gfortranを使用。最適化オプション -O3 とレガシー構文許容 -std=legacy を追加
gfortran -O3 -std=legacy -o qcd update.f \
staple.f phbupd.f \
libv_u.f lib_ran.f libv_sort.f siteindex2v.f lib_vec.f libv_reunit2.f
