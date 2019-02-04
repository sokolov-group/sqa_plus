def compute_S_ip(mr_adc, rdm_ca_so, rdm_ccaa_so, order = None):

# Total S = S(0) + S(1)
    ni_cas = mr_adc.nion_cas
    ni_so = mr_adc.ni_so

    ncore_so = mr_adc.ncore_so
    ncas_so = mr_adc.ncas_so
    nextern_so = mr_adc.nextern_so

    nxyi_so = ncas_so * ncas_so * ncore_so

    naji_so = nextern_so * ncore_so * ncore_so
    naxi_so = nextern_so * ncas_so * ncore_so
    naxy_so = nextern_so * ncas_so * ncas_so
    nxji_so = ncas_so * ncore_so * ncore_so

##################  Only for {h+^1 | h^0} +  {h+^0 | h^1}
#    si = 0
#    fi = ni_so
#    sxyi = fi
#    fxyi = sxyi + nxyi_so
#
#    dim = ni_so + nxyi_so
#
## S(1): {h+^1 | h^0} +  {h+^0 | h^1}
#    S = np.zeros((dim, dim))
#
##{h+^1 | h^0} Overlap matrix
#    # xyi - j
#    S[sxyi:fxyi, si:fi] = np.einsum('xy,ij->xyij', rdm_ca_so, np.identity(ncore_so)).reshape((nxyi_so, ni_so))
#
##{h+^0 | h^1} Overlap matrix
#    # j - xyi
#    S[si:fi, sxyi:fxyi] = S[sxyi:fxyi, si:fi].T.copy()
#
#    print np.linalg.norm(S - S.T)

#    return S
#######################################################################

    sI = mr_adc.sI
    fI = mr_adc.fI
    si = mr_adc.si
    fi = mr_adc.fi
    sxyi = fi
    fxyi = sxyi + nxyi_so

    saji = fxyi
    faji = saji + naji_so
    saxi = faji
    faxi = saxi + naxi_so
    saxy = faxi
    faxy = saxy + naxy_so
    sxji = faxy
    fxji = sxji + nxji_so

#    dim = ni_cas + ni_so + nxyi_so
#    dim = ni_cas + ni_so + nxyi_so + naji_so + naxi_so + naxy_so + nxji_so

#   S : {h+^0 | h^0} 
    if (order == 0):
       dim = ni_cas + ni_so
#   S : {h+^0 | h^0} + {h+^1 | h^0}
    elif (order == 1):
       dim = ni_cas + ni_so + nxyi_so
#   S : {h+^0 | h^0} + {h+^1 | h^0} + {h+^0 | h^1} + {h+^1 | h^1}
    else:
       dim = ni_cas + ni_so + nxyi_so + naji_so + naxi_so + naxy_so + nxji_so

    S = np.zeros((dim, dim))

#   {h+^0 | h^0} Overlap matrix
    # I - J
    S[sI:fI, sI:fI] = np.identity(ni_cas)
    # i - j
    S[si:fi, si:fi] = np.identity(ni_so)

    if (order == 0):
#       print S
#       print S.shape
       return S

#   {h+^1 | h^0}
    # xyi - j
    S[sxyi:fxyi, si:fi] = np.einsum('xy,ij->xyij', rdm_ca_so, np.identity(ncore_so)).reshape((nxyi_so, ni_so))

#   {h+^0 | h^1}
    # j - xyi
    S[si:fi, sxyi:fxyi] = S[sxyi:fxyi, si:fi].T.copy()

#   S = {h+^0 | h^0} + {h+^1 | h^0} +  {h+^0 | h^1}
    if (order == 1):
#       print S[sxyi:fxyi, si:fi]
#       print S.shape
       return S
#
#   S = {h+^0 | h^0} + {h+^1 | h^0} +  {h+^0 | h^1} + {h+^1 | h^1}
#
#   {h+^1 | h^1}
    # xyi - xyi
    S[sxyi:fxyi, sxyi:fxyi] = np.einsum('IJ, XZ, WY->XYIZWJ', np.identity(ncore_so), np.identity(ncas_so), rdm_ca_so, optimize = True).reshape((nxyi_so, nxyi_so))
    S[sxyi:fxyi, sxyi:fxyi] += np.einsum('IJ, WXYZ->XYIZWJ', np.identity(ncore_so), rdm_ccaa_so, optimize = True).reshape((nxyi_so, nxyi_so))
    # aji - aji
    S[saji:faji, saji:faji] = np.einsum('AB, IK, JL->AJIBLK', np.identity(nextern_so), np.identity(ncore_so), np.identity(ncore_so), optimize = True).reshape((naji_so, naji_so))
    S[saji:faji, saji:faji ] -= np.einsum('AB, IL, JK->AJIBLK', np.identity(nextern_so), np.identity(ncore_so), np.identity(ncore_so), optimize = True).reshape((naji_so, naji_so))
    # axi - axi
    S[saxi:faxi, saxi:faxi] = np.einsum('AB, IJ, YX->AXIBYJ', np.identity(nextern_so), np.identity(ncore_so), rdm_ca_so, optimize = True).reshape((naxi_so, naxi_so))
    # axy - axy
    S[saxy:faxy, saxy:faxy] = np.einsum('AB, WZXY->AXYBZW', np.identity(nextern_so), rdm_ccaa_so, optimize = True).reshape((naxy_so, naxy_so))
    # xji - xji
    S[sxji:fxji, sxji:fxji] = np.einsum('IK, JL, XY->XJIYLK', np.identity(ncore_so), np.identity(ncore_so), np.identity(ncas_so), optimize = True).reshape((nxji_so, nxji_so))
    S[sxji:fxji, sxji:fxji] -= np.einsum('IL, JK, XY->XJIYLK', np.identity(ncore_so), np.identity(ncore_so), np.identity(ncas_so), optimize = True).reshape((nxji_so, nxji_so))
    S[sxji:fxji, sxji:fxji] -= np.einsum('IK, JL, XY->XJIYLK', np.identity(ncore_so), np.identity(ncore_so), rdm_ca_so, optimize = True).reshape((nxji_so, nxji_so))
    S[sxji:fxji, sxji:fxji] += np.einsum('IL, JK, XY->XJIYLK', np.identity(ncore_so), np.identity(ncore_so), rdm_ca_so, optimize = True).reshape((nxji_so, nxji_so))

#
#    print S[sI:fI, sI:fI]
#    print S[si:fi, si:fi]
#    print S[sxyi:fxyi, si:fi]
#    print S[sxji:fxji, sxji:fxji]

#    print S
#    print S.shape

#    print np.linalg.norm(S - S.T)

    return S
