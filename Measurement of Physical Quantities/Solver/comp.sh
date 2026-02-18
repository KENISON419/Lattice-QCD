#

f77 -o slvm2 solver_llm2.f \
mr2.f opr2.f llmcorr2.f \
libv_ux.f lib_vec.f libv_u.f libv_sort.f siteindex.f libv_gxD.f libv_bfv.f

#f77 -o slvs2 solver_llm2.f \
#bs2.f opr2.f llmcorr2.f \
#libv_ux.f lib_vec.f libv_u.f libv_sort.f siteindex.f libv_gxD.f libv_bfv.f

