#include "energy.h"
#include "structure.h"

/// Energy_BiasingPotential calculates the biasing potential used in the
/// thermalization/equilibration The function calculates (T_current - T_final)
/// and if < 0.001, sets nBiasPotential_Mode_glb = 0.
/// \param beadID
/// \return The energy, given the nBiasPotential_Mode_glb
float Energy_BiasingPotential(const int beadID)
{
    float totEn                = 0.f;
    const int centPos[POS_MAX] = {naBoxSize_glb[0] / 2, naBoxSize_glb[1] / 2, naBoxSize_glb[2] / 2};
    const int tmpPos[POS_MAX]  = {bead_info_glb[beadID][0], bead_info_glb[beadID][1], bead_info_glb[beadID][2]};
    const int posDiff[POS_MAX] = {tmpPos[0] - centPos[0], tmpPos[1] - centPos[1], tmpPos[2] - centPos[2]};

    const int thisType         = bead_info_glb[beadID][BEAD_TYPE];
    switch (nBiasPotential_Mode_glb)
        {
            case 1:
                /*
                 * The potential has the form (r-r_{cen})^2
                 */
                totEn = (float) Dist_VecMagSq(posDiff);
                totEn = fCuTemp_glb * totEn;
                break;

            case 2:
                /*
                 * The potential has the form (r-r_{cen})^2 if the bead is outside the indent radius
                 */
                totEn = (float) Dist_VecMagSq(posDiff);
                if (totEn > fSquishRad_Sq_glb)
                    {
                        totEn = fCuTemp_glb * totEn * fSquish_Stiff_glb;
                    }
                else
                    {
                        totEn = 0.f;
                    }
                break;

            case 3:
                /*
                 * Same as case 2, BUT, ignores the last sticker-type. We can think of the last sticker-type
                 * as a crowder.
                 */
                totEn = (float) Dist_VecMagSq(posDiff);
                if (thisType == ((int) nBeadTypes_glb - 1))
                    {
                        if (totEn < fSquishRad_Sq_glb)
                            {
                                totEn = fCuTemp_glb * fSquish_Stiff_glb;
                            }
                        else
                            {
                                totEn = 0.f;
                            }
                    }
                else
                    {
                        if (totEn > fSquishRad_Sq_glb)
                            {
                                totEn = fCuTemp_glb * totEn * fSquish_Stiff_glb;
                            }
                        else
                            {
                                totEn = 0.f;
                            }
                    }
                break;

            default:
                totEn = 0.f;
                break;
        }
    return totEn;
}

/// Energy_Anisotropic - returns the bond energy for beadID if it is bonded.
/// \param beadID
/// \return
inline float Energy_Anisotropic(const int beadID)
{ // Calculates the SC-SC energy of the bead in question.
    const int bP   = bead_info_glb[beadID][BEAD_FACE];
    const int resi = bead_info_glb[beadID][BEAD_TYPE];
    int resj;
    if (bP != -1)
        {
            resj = bead_info_glb[bP][BEAD_TYPE];
            return faEnergy_glb[resi][resj][E_SC_SC];
        }
    else
        {
            return 0.f;
        }
}

/// Energy_Anisotropic_Self - returns the bond energy for beadID only if it is
/// bonded to a bead in the same chain. Note that this does not account for
/// double counting \param beadID \return
float Energy_Anisotropic_Self(const int beadID)
{                      // Calculates the SC-SC energy of the bead in question.
    float totEn = 0.0; // Storing total overlap energy
    int bP      = bead_info_glb[beadID][BEAD_FACE];
    if (bP != -1)
        { // This bead does have a bond
            if (bead_info_glb[beadID][BEAD_CHAINID] == bead_info_glb[bP][BEAD_CHAINID])
                {
                    totEn += faEnergy_glb[bead_info_glb[beadID][BEAD_TYPE]][bead_info_glb[bP][BEAD_TYPE]][E_SC_SC];
                }
        }
    return totEn;
}

/// Energy_Anisotropic_For_Chain - returns the bond energy for beadID if it is
/// bonded If it is bonded to another in the same chain, divide the energy by 2.
/// This should only be used when calculating the energy of an entire chain!
/// \param beadID
/// \return
float Energy_Anisotropic_For_Chain(const int beadID)
{                       // Calculates the SC-SC energy of the bead in question.
    float totEn  = 0.f; // Storing total overlap energy
    const int bP = bead_info_glb[beadID][BEAD_FACE];
    if (bP != -1)
        { // This bead does have a bond
            if (bead_info_glb[beadID][BEAD_CHAINID] == bead_info_glb[bP][BEAD_CHAINID])
                {
                    totEn +=
                        faEnergy_glb[bead_info_glb[beadID][BEAD_TYPE]][bead_info_glb[bP][BEAD_TYPE]][E_SC_SC] * 0.5f;
                }
            else
                {
                    totEn += faEnergy_glb[bead_info_glb[beadID][BEAD_TYPE]][bead_info_glb[bP][BEAD_TYPE]][E_SC_SC];
                }
        }
    return totEn;
}

/// Energy_Anisotropic_For_Range - returns the bond energy for beadID if it is
/// bonded. This function assumes that we get a beadID range [smallest_bead,
/// largest_bead], and that we will loop over every bead in this range. So if
/// the bead partner is within this range, we divide by 2 to account for double
/// counting. \param beadID \param smBead \param lgBead \return
float Energy_Anisotropic_For_Range(const int beadID, const int smBead, const int lgBead)
{
    float totEn  = 0.f;
    const int bP = bead_info_glb[beadID][BEAD_FACE];
    if (bP != -1)
        {
            if ((beadID <= lgBead) && (beadID >= smBead))
                {
                    totEn +=
                        faEnergy_glb[bead_info_glb[beadID][BEAD_TYPE]][bead_info_glb[bP][BEAD_TYPE]][E_SC_SC] * 0.5f;
                }
            else
                {
                    totEn += faEnergy_glb[bead_info_glb[beadID][BEAD_TYPE]][bead_info_glb[bP][BEAD_TYPE]][E_SC_SC];
                }
        }
    return totEn;
}

/// Energy_Anisotropic_Contiguous_Range - returns the bond energy for beadID if
/// it is bonded. This function assumes that we get a beadID range
/// [smallest_bead, largest_bead], and that we will loop over every bead in this
/// range. So if the bead partner is within this range, we divide by 2 to
/// account for double counting. \param beadID \return
float Energy_Anisotropic_Contiguous_Range(const int beadID, const int smallest_bead, const int largest_bead)
{
    float totEn = 0.0; // Storing total overlap energy
    int bP      = bead_info_glb[beadID][BEAD_FACE];
    if (bP != -1)
        { // This bead does have a bond
            if (bead_info_glb[beadID][BEAD_CHAINID] == bead_info_glb[bP][BEAD_CHAINID])
                { // Intra-molecular
                    if (bP >= smallest_bead && bP <= largest_bead)
                        { // Within subset
                            totEn +=
                                faEnergy_glb[bead_info_glb[beadID][BEAD_TYPE]][bead_info_glb[bP][BEAD_TYPE]][E_SC_SC] /
                                2.;
                        }
                    else
                        {
                            totEn +=
                                faEnergy_glb[bead_info_glb[beadID][BEAD_TYPE]][bead_info_glb[bP][BEAD_TYPE]][E_SC_SC];
                        }
                }
            else
                {
                    totEn += faEnergy_glb[bead_info_glb[beadID][BEAD_TYPE]][bead_info_glb[bP][BEAD_TYPE]][E_SC_SC];
                }
        }
    return totEn;
}

/// Energy_Anisotropic_With_List - returns the bond energy for beadID if it is
/// bonded. This function assumes that we get a list of beadID's of size
/// list_size, but that the last bead is at list_size-1. This also assumes that
/// we will loop over the entire list So if the bead partner is within this
/// list, we divide by 2 to account for double counting. \param beadID \return
float Energy_Anisotropic_With_List(const int beadID, const int* bead_list, const int list_size)
{
    float totEn    = 0.0; // Storing total overlap energy
    int bP         = bead_info_glb[beadID][BEAD_FACE];
    int check_bead = 0;
    int i;
    if (bP != -1)
        { // This bead does have a bond
            for (i = 0; i < list_size; i++)
                { // Checking if we should worry about double counting
                    if (bP == bead_list[i])
                        {
                            check_bead = 1;
                            break;
                        }
                }
            if (check_bead == 1)
                {
                    totEn += faEnergy_glb[bead_info_glb[beadID][BEAD_TYPE]][bead_info_glb[bP][BEAD_TYPE]][E_SC_SC] / 2.;
                }
            else
                {
                    totEn += faEnergy_glb[bead_info_glb[beadID][BEAD_TYPE]][bead_info_glb[bP][BEAD_TYPE]][E_SC_SC];
                }
        }
    return totEn;
}

float Energy_Anisotropic_For_List(const int beadID, const int listSize, const int beadList[MAX_BONDS + 1])
{
    const int bP   = bead_info_glb[beadID][BEAD_FACE];
    const int resi = bead_info_glb[beadID][BEAD_TYPE];
    int resj;
    if (bP != -1)
        {
            resj = bead_info_glb[bP][BEAD_TYPE];
            if (Check_BeadID_InList(bP, listSize, beadList))
                {
                    return faEnergy_glb[resi][resj][E_SC_SC] * 0.5f;
                }
            else
                {
                    return faEnergy_glb[resi][resj][E_SC_SC];
                }
        }
    else
        {
            return 0.f;
        }
}

/// Energy_Iso_Ovlp - calculate the OVLP interaction between two beads of types
/// beadType1 and beadType2, respectively. \param beadType1: Type of first bead.
/// \param beadType2: Type of second bead.
/// \param xDis: Distance between the two beads.
/// \return faEnergy_glb[beadType1][beadType2][E_OVLP].
float Energy_Iso_Ovlp(int const beadType1, int const beadType2, float const xDis)
{
    return faEnergy_glb[beadType1][beadType2][E_OVLP];
}

/// Energy_Iso_Cont - calculate the CONT interaction between two beads of types
/// beadType1 and beadType2, respectively. \param beadType1: Type of first bead.
/// \param beadType2: Type of second bead.
/// \param xDis: Distance between the two beads.
/// \return faEnergy_glb[beadType1][beadType2][E_CONT]/xDis.
float Energy_Iso_Cont(int const beadType1, int const beadType2, float const xDis)
{
    return faEnergy_glb[beadType1][beadType2][E_CONT] / xDis;
}

/// Energy_Iso_fSol: Returns the solvation energy of this bead.
/// \param beadType
/// \param xDis
/// \return
inline float Energy_Iso_fSol(int const beadType)
{
    return faEnergy_glb[beadType][beadType][E_F_SOL];
}

float Energy_Topo_Angle(int const beadID)
{
    const int resi      = bead_info_glb[beadID][BEAD_TYPE];
    const int frontBead = topo_info_glb[beadID][1];
    if (frontBead == -1)
        {
            return 0.f;
        }
    const int backBead = topo_info_glb[beadID][0];

    const int r_pos0[POS_MAX] = {bead_info_glb[beadID][0], bead_info_glb[beadID][1], bead_info_glb[beadID][2]};
    const int r_posB[POS_MAX] = {bead_info_glb[backBead][0], bead_info_glb[backBead][1], bead_info_glb[backBead][2]};
    const int r_posF[POS_MAX] = {bead_info_glb[frontBead][0], bead_info_glb[frontBead][1], bead_info_glb[frontBead][2]};
    int vec1[POS_MAX];
    int vec2[POS_MAX];

    BeadPos_sub_wPBC(vec1, r_pos0, r_posB);
    BeadPos_sub_wPBC(vec2, r_posF, r_pos0);

    const float dumCosTheta = Vec3n_CosTheta(vec1, vec2);

    return -faEnergy_glb[resi][resi][E_STIFF] * (1.f - dumCosTheta * dumCosTheta);
}

/// Energy_OfOvlp_wNeighList: Given this bead, and a supplied list of neighbors
/// and number of neighbors, we loop over all the neighbors and add the
/// energies. This function _only_ calculates the Ovlp energies!
/// \param beadID
/// \param neighList
/// \param neighNum
/// \return
float Energy_OfOvlp_wNeighList(int const beadID, const int* neighList, int const neighNum)
{
    float totEn = 0.f;
    int i;
    const int resi = bead_info_glb[beadID][BEAD_TYPE];
    int resj, tmpID;
    float xDis;

    for (i = 0; i < neighNum; i++)
        {
            tmpID = neighList[i];
            resj  = bead_info_glb[tmpID][BEAD_TYPE];
            xDis  = Dist_BeadToBead(beadID, tmpID);
            totEn += Energy_Iso_Ovlp(resi, resj, xDis);
        }

    return totEn;
}

/// Energy_OfCont_wNeighList: Given this bead, and a supplied list of neighbors
/// and number of neighbors, we loop over all the neighbors and add the
/// energies. This function _only_ calculates the Cont energies! \param beadID
/// \param neighList
/// \param neighNum
/// \return
float Energy_OfCont_wNeighList(int const beadID, const int* neighList, int const neighNum)
{
    float totEn = 0.f;
    int i;
    const int resi = bead_info_glb[beadID][BEAD_TYPE];
    int resj, tmpID;
    float xDis;

    for (i = 0; i < neighNum; i++)
        {
            tmpID = neighList[i];
            resj  = bead_info_glb[tmpID][BEAD_TYPE];
            xDis  = Dist_BeadToBead(beadID, tmpID);
            totEn += Energy_Iso_Cont(resi, resj, xDis);
        }

    return totEn;
}

/// Energy_ofPairs_wNeighList: Given this bead, and a supplied list of neighbors
/// and number of neighbors, we loop over all the neighbors and add the
/// energies. This function _only_ calculate the pairwise energies so the
/// solvation energy has to be handled somewhere else. \param beadID \param
/// neighList \param neighNum \return
float Energy_ofPairs_wNeighList(int const beadID, const int* neighList, int const neighNum)
{
    float totEn = 0.f;
    int i;
    const int resi = bead_info_glb[beadID][BEAD_TYPE];
    int resj, tmpID;
    float xDis;

    for (i = 0; i < neighNum; i++)
        {
            tmpID = neighList[i];
            resj  = bead_info_glb[tmpID][BEAD_TYPE];
            xDis  = Dist_BeadToBead(beadID, tmpID);
            totEn += Energy_Iso_Ovlp(resi, resj, xDis);
            totEn += Energy_Iso_Cont(resi, resj, xDis);
        }

    return totEn;
}

/// Energy_ofSol_wNeighList: Given the neighbor list, we calculate the total
/// solvation energy at the lattice-site where beadID is, but if the site was
/// empty. \param neighList \param neighNum \return
float Energy_ofSol_wNeighList(const int* neighList, int const neighNum)
{
    float totEn = 0.f;
    int i;
    int resj, tmpID;

    for (i = 0; i < neighNum; i++)
        {
            tmpID = neighList[i];
            resj  = bead_info_glb[tmpID][BEAD_TYPE];
            totEn += Energy_Iso_fSol(resj);
        }

    return totEn;
}

/// Energy_OfOvlp_wNeighList_ForChains: Given this bead, and a supplied list of
/// neighbors and number of neighbors, we loop over all the neighbors and add
/// the energies, taking care of intra-chain double-counting. This function
/// _only_ calculates the Ovlp energies! \param beadID \param neighList \param
/// neighNum \return
float Energy_OfOvlp_wNeighList_ForChains(int const beadID, const int* neighList, int const neighNum)
{
    float totEn = 0.f;
    int i;
    const int resi    = bead_info_glb[beadID][BEAD_TYPE];
    const int chain_i = bead_info_glb[beadID][BEAD_CHAINID];
    int resj, chain_j, tmpID;
    const float xDis = 1.f;

    for (i = 0; i < neighNum; i++)
        {
            tmpID   = neighList[i];
            resj    = bead_info_glb[tmpID][BEAD_TYPE];
            chain_j = bead_info_glb[tmpID][BEAD_CHAINID];
            if (chain_i == chain_j)
                {
                    totEn += Energy_Iso_Ovlp(resi, resj, xDis) * 0.5f;
                }
            else
                {
                    totEn += Energy_Iso_Ovlp(resi, resj, xDis);
                }
        }

    return totEn;
}

/// Energy_OfCont_wNeighList_ForChains: Given this bead, and a supplied list of
/// neighbors and number of neighbors,
///// we loop over all the neighbors and add the energies, taking care of
/// intra-chain double-counting.
///// This function _only_ calculates the Contenergies!
/// \param beadID
/// \param neighList
/// \param neighNum
/// \return
float Energy_OfCont_wNeighList_ForChains(int const beadID, const int* neighList, int const neighNum)
{
    float totEn = 0.f;
    int i;
    const int resi    = bead_info_glb[beadID][BEAD_TYPE];
    const int chain_i = bead_info_glb[beadID][BEAD_CHAINID];
    int resj, chain_j, tmpID;
    float xDis;

    for (i = 0; i < neighNum; i++)
        {
            tmpID   = neighList[i];
            resj    = bead_info_glb[tmpID][BEAD_TYPE];
            chain_j = bead_info_glb[tmpID][BEAD_CHAINID];
            xDis    = Dist_BeadToBead(beadID, tmpID);
            if (chain_i == chain_j)
                {
                    totEn += Energy_Iso_Cont(resi, resj, xDis) * 0.5f;
                }
            else
                {
                    totEn += Energy_Iso_Cont(resi, resj, xDis);
                }
        }

    return totEn;
}

/// Energy_OfOvlp_wNeighList_ForRange: Given this bead, and a supplied list of
/// neighbors and number of neighbors, we loop over all the neighbors and add
/// the energies, taking care of intra-range double-counting. Note that the
/// range is [loBead, hiBead]. This function _only_ calculates the Ovlp
/// energies! \param beadID \param neighList \param neighNum \return
float Energy_OfOvlp_wNeighList_ForRange(int const beadID, const int loBead, const int hiBead, const int* neighList,
                                        int const neighNum)
{
    float totEn = 0.f;
    int i;
    const int resi = bead_info_glb[beadID][BEAD_TYPE];
    int resj, tmpID;
    const float xDis = 1.f;

    for (i = 0; i < neighNum; i++)
        {
            tmpID = neighList[i];
            resj  = bead_info_glb[tmpID][BEAD_TYPE];
            //            xDis  = Dist_BeadToBead(beadID, tmpID);
            if ((i >= loBead) && (i <= hiBead))
                {
                    totEn += Energy_Iso_Ovlp(resi, resj, xDis) * 0.5f;
                }
            else
                {
                    totEn += Energy_Iso_Ovlp(resi, resj, xDis);
                }
        }

    return totEn;
}

/// Energy_OfCont_wNeighList_ForRange: Given this bead, and a supplied list of
/// neighbors and number of neighbors, we loop over all the neighbors and add
/// the energies, taking care of intra-range double-counting. Note that the
/// range is [loBead, hiBead]. This function _only_ calculates the Cont
/// energies! \param beadID \param neighList \param neighNum \return
float Energy_OfCont_wNeighList_ForRange(int const beadID, const int loBead, const int hiBead, const int* neighList,
                                        int const neighNum)
{
    float totEn = 0.f;
    int i;
    const int resi = bead_info_glb[beadID][BEAD_TYPE];
    int resj, tmpID;
    float xDis;

    for (i = 0; i < neighNum; i++)
        {
            tmpID = neighList[i];
            resj  = bead_info_glb[tmpID][BEAD_TYPE];
            xDis  = Dist_BeadToBead(beadID, tmpID);
            if ((i >= loBead) && (i <= hiBead))
                {
                    totEn += Energy_Iso_Cont(resi, resj, xDis) * 0.5f;
                }
            else
                {
                    totEn += Energy_Iso_Cont(resi, resj, xDis);
                }
        }

    return totEn;
}

/// Energy_OfOvlp_wNeighList_ForLists: Given the bead, neighbor-list, and
/// bead-list, we loop over all the neighbors and add the energies, taking care
/// of intra-bead-list double-counting. This function _only_ calculates the Ovlp
/// energies! \param beadID \param listSize \param beadList \param neighList
/// \param neighNum
/// \return
float Energy_OfOvlp_wNeighList_ForLists(const int beadID, const int listSize, const int beadList[MAX_BONDS + 1],
                                        const int* neighList, const int neighNum)
{
    float totEn = 0.f;
    int i;
    const int resi = bead_info_glb[beadID][BEAD_TYPE];
    int resj, tmpID;
    const float xDis = 1.f;

    for (i = 0; i < neighNum; i++)
        {
            tmpID = neighList[i];
            resj  = bead_info_glb[tmpID][BEAD_TYPE];
            //            xDis  = Dist_BeadToBead(beadID, tmpID);
            if (Check_BeadID_InList(tmpID, listSize, beadList))
                {
                    totEn += Energy_Iso_Ovlp(resi, resj, xDis) * 0.5f;
                }
            else
                {
                    totEn += Energy_Iso_Ovlp(resi, resj, xDis);
                }
        }

    return totEn;
}

/// Energy_OfOvlp_wNeighList_ForLists: Given the bead, neighbor-list, and
/// bead-list, we loop over all the neighbors and add the energies, taking care
/// of intra-bead-list double-counting. This function _only_ calculates the Cont
/// energies! \param beadID \param listSize \param beadList \param neighList
/// \param neighNum
/// \return
float Energy_OfCont_wNeighList_ForLists(const int beadID, const int listSize, const int beadList[MAX_BONDS + 1],
                                        const int* neighList, const int neighNum)
{
    float totEn = 0.f;
    int i;
    const int resi = bead_info_glb[beadID][BEAD_TYPE];
    int resj, tmpID;
    float xDis;

    for (i = 0; i < neighNum; i++)
        {
            tmpID = neighList[i];
            resj  = bead_info_glb[tmpID][BEAD_TYPE];
            xDis  = Dist_BeadToBead(beadID, tmpID);
            if (Check_BeadID_InList(tmpID, listSize, beadList))
                {
                    totEn += Energy_Iso_Cont(resi, resj, xDis) * 0.5f;
                }
            else
                {
                    totEn += Energy_Iso_Cont(resi, resj, xDis);
                }
        }

    return totEn;
}

/// Energy_Isotroptic_Old calculates the isotropic contribution to the energy by
/// searching the 3^3-1 = 26 'neighbors' The energy function  is like the BFM,
/// where $epislon$ represents the overlap cost for total overlap, which is
/// forbidden explicitly in LASSI, so we have $epislon/2$,$epislon/4$ and
/// $epislon/8$ with increasing distance. \param beadID \return
float Energy_Isotropic_Old(const int beadID)
{                      // Calculate Contact and Overlap energy of bead beadID
    float totEn = 0.f; // Storing total overlap energy
    int i, j;          // Indecies
    int r_pos_0[POS_MAX], r_chck[POS_MAX], r_disp[POS_MAX];
    int x, y, z;      // Lattice indecies
    int secBi, resj;  // Second bead index
    float xDis = 0.f; // Distance between beads.
    int resi   = bead_info_glb[beadID][BEAD_TYPE];
    totEn += nBiasPotential_Mode_glb == -1 ? 0.f : Energy_BiasingPotential(beadID);

    if (nBeadTypeCanOvlp_glb[resi] == 0 && nBeadTypeCanCont_glb[resi] == 0 && nBeadTypeCanFSol_glb[resi] == 0 &&
        nBeadTypeCanTInd_glb[resi] == 0)
        {
            return totEn;
        } // No need to do anything if there's no isotropic interactions.

    int BoxRad = nBeadTypeCanCont_glb[resi] == 0 ? 1 : LARGEST_RADIUS; // No need to search if no cont interactions

    for (j = 0; j < POS_MAX; j++)
        {
            r_pos_0[j] = bead_info_glb[beadID][j];
        }
    for (x = -BoxRad; x <= BoxRad; x++)
        {
            r_disp[0] = x;
            for (y = -BoxRad; y <= BoxRad; y++)
                {
                    r_disp[1] = y;
                    for (z = -BoxRad; z <= BoxRad; z++)
                        {
                            r_disp[2] = z;
                            for (j = 0; j < POS_MAX; j++)
                                {
                                    r_chck[j] = r_pos_0[j] + r_disp[j];
                                    r_chck[j] = r_chck[j] < 0 ? r_chck[j] + naBoxSize_glb[j] : r_chck[j];
                                    r_chck[j] =
                                        r_chck[j] >= naBoxSize_glb[j] ? r_chck[j] - naBoxSize_glb[j] : r_chck[j];
                                }
                            secBi = naTotLattice_glb[Lat_Ind_FromVec(r_chck)];
                            if (secBi != -1 && secBi != beadID)
                                {
                                    resj = bead_info_glb[secBi][BEAD_TYPE];
                                    totEn += Energy_Iso_Ovlp(resi, resj, xDis);
                                    xDis = sqrtf((float) (x * x + y * y + z * z));
                                    // xDis = Dist_BeadToBead(beadID, secBi);
                                    if (xDis <= 1.0)
                                        {
                                            totEn += Energy_Iso_Ovlp(resi, resj, xDis) / 2.0f;
                                        }
                                    else if (xDis <= 1.42)
                                        { // sqrt(2)
                                            totEn += Energy_Iso_Ovlp(resi, resj, xDis) / 4.0f;
                                        }
                                    else if (xDis <= 1.74)
                                        { // sqrt(3)
                                            totEn += Energy_Iso_Ovlp(resi, resj, xDis) / 8.0f;
                                        }
                                    else
                                        { // This way, contact is only outside ovlp
                                            totEn += Energy_Iso_Cont(resi, resj, xDis);
                                        }
                                }
                            else if (secBi == -1)
                                {
                                    if (abs(x) <= 1 && abs(y) <= 1 && abs(z) <= 1)
                                        { // Want solvation radius to be 1
                                            totEn += faEnergy_glb[resi][resi][E_F_SOL];
                                            totEn += faEnergy_glb[resi][resi][E_T_IND] * fCuTemp_glb;
                                        }
                                }
                        }
                }
        }

    return totEn;
}

/// Energy_Isotroptic calculates the isotropic contribution to the energy by
/// searching the 3^3-1 = 26 'neighbors' The energy function  is like the BFM,
/// where $x$ represents the overlap cost for total overlap, which is forbidden
/// explicitly in LASSI, so we have $x/2$,$x/4$ and $x/8$ with increasing
/// distance.
/// \param beadID
/// \return
float Energy_Isotropic(const int beadID)
{                      // Calculate Contact and Overlap energy of bead beadID
    float totEn = 0.0; // Storing total overlap energy
    int i, j;          // Indecies
    int r_pos_0[POS_MAX], r_chck[POS_MAX], r_disp[POS_MAX];
    int x, y, z;     // Lattice indecies
    int secBi, resj; // Second bead index
    float xDis = 0.; // Distance between beads.
    int resi   = bead_info_glb[beadID][BEAD_TYPE];
    totEn += nBiasPotential_Mode_glb == -1 ? 0.f : Energy_BiasingPotential(beadID);

    if (nBeadTypeCanOvlp_glb[resi] == 0 && nBeadTypeCanCont_glb[resi] == 0 && nBeadTypeCanFSol_glb[resi] == 0 &&
        nBeadTypeCanTInd_glb[resi] == 0)
        {
            return totEn;
        } // No need to do anything if there's no isotropic interactions.

    int BoxRad = nBeadTypeCanCont_glb[resi] == 0 ? 1 : LARGEST_RADIUS; // No need to search if no cont interactions

    for (j = 0; j < POS_MAX; j++)
        {
            r_pos_0[j] = bead_info_glb[beadID][j];
        }
    for (x = -BoxRad; x <= BoxRad; x++)
        {
            r_disp[0] = x;
            for (y = -BoxRad; y <= BoxRad; y++)
                {
                    r_disp[1] = y;
                    for (z = -BoxRad; z <= BoxRad; z++)
                        {
                            r_disp[2] = z;
                            LatPos_add_wPBC(r_chck, r_pos_0, r_disp);
                            secBi = naTotLattice_glb[Lat_Ind_FromVec(r_chck)];
                            if (secBi != -1 && secBi != beadID)
                                {
                                    resj = bead_info_glb[secBi][BEAD_TYPE];
                                    xDis = sqrtf((float) (x * x + y * y + z * z));
                                    totEn += Energy_Iso_Ovlp(resi, resj,
                                                             xDis); /// xDis / xDis / xDis;

                                    totEn += Energy_Iso_Cont(resi, resj, xDis);
                                }
                            else if (secBi == -1)
                                {
                                    if (abs(x) <= 1 && abs(y) <= 1 && abs(z) <= 1)
                                        { // Want solvation radius to be 1
                                            totEn += faEnergy_glb[resi][resi][E_F_SOL];
                                        }
                                }
                        }
                }
        }

    return totEn;
}

/// Energy_Isotropic_Self - returns the bond energy for beadID only if it is
/// interacting with a bead within the same molecule. Note that it does not
/// account for double counting \param beadID \return
float Energy_Isotropic_Self(const int beadID)
{ // Calculate Contact and Overlap energy of beadID but
  // only intra-molecular
    // interactions
    float totEn = 0.0; // Storing total overlap energy
    int i, j;          // Indecies
    int r_pos_0[POS_MAX], r_chck[POS_MAX], r_disp[POS_MAX];
    int x, y, z;     // Lattice indecies
    int secBi, resj; // Second bead index
    float xDis = 0.; // Distance between beads.
    int resi   = bead_info_glb[beadID][BEAD_TYPE];
    // totEn += nBiasPotential_Mode_glb == -1 ? 0. : Energy_BiasingPotential(beadID);

    if (nBeadTypeCanOvlp_glb[resi] == 0 && nBeadTypeCanCont_glb[resi] == 0 && nBeadTypeCanFSol_glb[resi] == 0 &&
        nBeadTypeCanTInd_glb[resi] == 0)
        {
            return totEn;
        } // No need to do anything if there's no isotropic interactions.

    int BoxRad = nBeadTypeCanCont_glb[resi] == 0 ? 1 : LARGEST_RADIUS; // No need to search if no cont interactions

    for (j = 0; j < POS_MAX; j++)
        {
            r_pos_0[j] = bead_info_glb[beadID][j];
        }
    for (x = -BoxRad; x <= BoxRad; x++)
        {
            r_disp[0] = x;
            for (y = -BoxRad; y <= BoxRad; y++)
                {
                    r_disp[1] = y;
                    for (z = -BoxRad; z <= BoxRad; z++)
                        {
                            r_disp[2] = z;
                            for (j = 0; j < POS_MAX; j++)
                                {
                                    r_chck[j] = r_pos_0[j] + r_disp[j];
                                    r_chck[j] = r_chck[j] < 0 ? r_chck[j] + naBoxSize_glb[j] : r_chck[j];
                                    r_chck[j] =
                                        r_chck[j] >= naBoxSize_glb[j] ? r_chck[j] - naBoxSize_glb[j] : r_chck[j];
                                }
                            secBi = naTotLattice_glb[Lat_Ind_FromVec(r_chck)];
                            if (secBi != -1 && secBi != beadID)
                                {
                                    if (bead_info_glb[secBi][BEAD_CHAINID] == bead_info_glb[beadID][BEAD_CHAINID])
                                        {
                                            resj = bead_info_glb[secBi][BEAD_TYPE];
                                            xDis = sqrtf((float) (x * x + y * y + z * z));
                                            totEn += Energy_Iso_Ovlp(resi, resj,
                                                                     xDis); /// xDis / xDis / xDis;

                                            totEn += Energy_Iso_Cont(resi, resj, xDis);
                                        }
                                }
                            else if (secBi == -1)
                                {
                                    if (abs(x) <= 1 && abs(y) <= 1 && abs(z) <= 1)
                                        { // Want solvation radius to be 1
                                            totEn += faEnergy_glb[resi][resi][E_F_SOL];
                                        }
                                }
                        }
                }
        }

    return totEn;
}

/// Energy_Isotropic_For_Chain - returns the isotropic interaction energy for
/// beadID If we have an intra-molecular interaction, divide the energy by 2.
/// This should only be used when calculating the energy of an entire chain!
/// \param beadID
/// \return
float Energy_Isotropic_For_Chain(const int beadID)
{ // Calculate Contact and Overlap energy of beadID
    // Takes care of intra-molecular double counting
    float totEn = 0.0; // Storing total overlap energy
    int i, j;          // Indecies
    int r_pos_0[POS_MAX], r_chck[POS_MAX], r_disp[POS_MAX];
    int x, y, z;     // Lattice indecies
    int secBi, resj; // Second bead index
    float xDis = 0.; // Distance between beads.
    int resi   = bead_info_glb[beadID][BEAD_TYPE];
    totEn += nBiasPotential_Mode_glb == -1 ? 0. : Energy_BiasingPotential(beadID);

    if (nBeadTypeCanOvlp_glb[resi] == 0 && nBeadTypeCanCont_glb[resi] == 0 && nBeadTypeCanFSol_glb[resi] == 0 &&
        nBeadTypeCanTInd_glb[resi] == 0)
        {
            return totEn;
        } // No need to do anything if there's no isotropic interactions.

    int BoxRad = nBeadTypeCanCont_glb[resi] == 0 ? 1 : LARGEST_RADIUS; // No need to search if no cont interactions

    for (j = 0; j < POS_MAX; j++)
        {
            r_pos_0[j] = bead_info_glb[beadID][j];
        }
    for (x = -BoxRad; x <= BoxRad; x++)
        {
            r_disp[0] = x;
            for (y = -BoxRad; y <= BoxRad; y++)
                {
                    r_disp[1] = y;
                    for (z = -BoxRad; z <= BoxRad; z++)
                        {
                            r_disp[2] = z;
                            for (j = 0; j < POS_MAX; j++)
                                {
                                    r_chck[j] = r_pos_0[j] + r_disp[j];
                                    r_chck[j] = r_chck[j] < 0 ? r_chck[j] + naBoxSize_glb[j] : r_chck[j];
                                    r_chck[j] =
                                        r_chck[j] >= naBoxSize_glb[j] ? r_chck[j] - naBoxSize_glb[j] : r_chck[j];
                                }
                            secBi = naTotLattice_glb[Lat_Ind_FromVec(r_chck)];
                            if (secBi != -1 && secBi != beadID)
                                {
                                    resj = bead_info_glb[secBi][BEAD_TYPE];
                                    xDis = sqrtf((float) (x * x + y * y + z * z));
                                    if (bead_info_glb[secBi][BEAD_CHAINID] == bead_info_glb[beadID][BEAD_CHAINID])
                                        { // Intra-molecular
                                            totEn += Energy_Iso_Ovlp(resi, resj, xDis) / 2.f;

                                            totEn += Energy_Iso_Cont(resi, resj, xDis) / 2.f;
                                        }
                                    else
                                        { // Inter-molecular
                                            totEn += Energy_Iso_Ovlp(resi, resj, xDis);

                                            totEn += Energy_Iso_Cont(resi, resj, xDis);
                                        }
                                }
                            else if (secBi == -1)
                                {
                                    if (abs(x) <= 1 && abs(y) <= 1 && abs(z) <= 1)
                                        { // Want solvation radius to be 1
                                            totEn += faEnergy_glb[resi][resi][E_F_SOL];
                                        }
                                }
                        }
                }
        }

    return totEn;
}

/// Energy_Isotropic_Contiguous_Range - returns the isotropic interactions for
/// beadID This function assumes that we get a beadID range [smallest_bead,
/// largest_bead], and that we will loop over every bead in this range. So if
/// the interactor is within this range, we divide by 2 to account for double
/// counting. \param beadID \return
float Energy_Isotropic_Contiguous_Range(const int beadID, const int smallest_bead, const int largest_bead)
{ // Calculate Contact and Overlap energy of beadID
    // Takes care of intra-molecular double counting
    // If the bead is between smallest_bead and largest_bead, we divide by two.
    // This assumes that every bead between smallest_bead and largest_bead will
    // be looped over
    float totEn = 0.0; // Storing total overlap energy
    int i, j;          // Indecies
    int r_pos_0[POS_MAX], r_chck[POS_MAX], r_disp[POS_MAX];
    int x, y, z;     // Lattice indecies
    int secBi, resj; // Second bead index
    float xDis = 0.; // Distance between beads.
    int resi   = bead_info_glb[beadID][BEAD_TYPE];
    totEn += nBiasPotential_Mode_glb == -1 ? 0. : Energy_BiasingPotential(beadID);

    if (nBeadTypeCanOvlp_glb[resi] == 0 && nBeadTypeCanCont_glb[resi] == 0 && nBeadTypeCanFSol_glb[resi] == 0 &&
        nBeadTypeCanTInd_glb[resi] == 0)
        {
            return totEn;
        } // No need to do anything if there's no isotropic interactions.

    int BoxRad = nBeadTypeCanCont_glb[resi] == 0 ? 1 : LARGEST_RADIUS; // No need to search if no cont interactions

    for (j = 0; j < POS_MAX; j++)
        {
            r_pos_0[j] = bead_info_glb[beadID][j];
        }
    for (x = -BoxRad; x <= BoxRad; x++)
        {
            r_disp[0] = x;
            for (y = -BoxRad; y <= BoxRad; y++)
                {
                    r_disp[1] = y;
                    for (z = -BoxRad; z <= BoxRad; z++)
                        {
                            r_disp[2] = z;
                            for (j = 0; j < POS_MAX; j++)
                                {
                                    r_chck[j] = r_pos_0[j] + r_disp[j];
                                    r_chck[j] = r_chck[j] < 0 ? r_chck[j] + naBoxSize_glb[j] : r_chck[j];
                                    r_chck[j] =
                                        r_chck[j] >= naBoxSize_glb[j] ? r_chck[j] - naBoxSize_glb[j] : r_chck[j];
                                }
                            secBi = naTotLattice_glb[Lat_Ind_FromVec(r_chck)];
                            if (secBi != -1 && secBi != beadID)
                                {
                                    resj = bead_info_glb[secBi][BEAD_TYPE];
                                    xDis = sqrtf((float) (x * x + y * y + z * z));
                                    if (bead_info_glb[secBi][BEAD_CHAINID] == bead_info_glb[beadID][BEAD_CHAINID])
                                        { // Intra-molecular
                                            if (secBi >= smallest_bead && secBi <= largest_bead)
                                                { // Within subset
                                                    totEn += Energy_Iso_Ovlp(resi, resj, xDis) /
                                                             2.; // xDis / xDis / xDis / 2.;

                                                    totEn += Energy_Iso_Cont(resi, resj, xDis) / 2.;
                                                }
                                            else
                                                {
                                                    totEn += Energy_Iso_Ovlp(resi, resj, xDis); // xDis / xDis / xDis;

                                                    totEn += Energy_Iso_Cont(resi, resj, xDis);
                                                }
                                        }
                                    else
                                        {                                               // Inter-molecular
                                            totEn += Energy_Iso_Ovlp(resi, resj, xDis); // / xDis / xDis / xDis;

                                            totEn += Energy_Iso_Cont(resi, resj, xDis);
                                        }
                                }
                            else if (secBi == -1)
                                {
                                    if (abs(x) <= 1 && abs(y) <= 1 && abs(z) <= 1)
                                        { // Want solvation radius to be 1
                                            totEn += faEnergy_glb[resi][resi][E_F_SOL];
                                        }
                                }
                        }
                }
        }

    return totEn;
}

/// Energy_Isotropic_With_List - returns the isotropic interactions for beadID.
/// This function assumes that we get a list of beadID's of size list_size, but
/// that the last bead is at list_size-1. This also assumes that we will loop
/// over the entire list So if the interactor is within this list, we divide by
/// 2 to account for double counting. \param beadID \return
float Energy_Isotropic_With_List(const int beadID, const int* bead_list, const int list_size)
{ // Calculate Contact and Overlap energy of beadID
    // Takes care of intra-molecular double counting
    float totEn = 0.0; // Storing total overlap energy
    int i, j;          // Indecies
    int bead_check = 0;
    int r_pos_0[POS_MAX], r_chck[POS_MAX], r_disp[POS_MAX];
    int x, y, z;     // Lattice indecies
    int secBi, resj; // Second bead index
    float xDis = 0.; // Distance between beads.
    int resi   = bead_info_glb[beadID][BEAD_TYPE];
    totEn += nBiasPotential_Mode_glb == -1 ? 0.f : Energy_BiasingPotential(beadID);

    if (nBeadTypeCanOvlp_glb[resi] == 0 && nBeadTypeCanCont_glb[resi] == 0 && nBeadTypeCanFSol_glb[resi] == 0 &&
        nBeadTypeCanTInd_glb[resi] == 0)
        {
            return totEn;
        } // No need to do anything if there's no isotropic interactions.

    int BoxRad = nBeadTypeCanCont_glb[resi] == 0 ? 1 : LARGEST_RADIUS; // No need to search if no cont interactions

    for (j = 0; j < POS_MAX; j++)
        {
            r_pos_0[j] = bead_info_glb[beadID][j];
        }
    for (x = -BoxRad; x <= BoxRad; x++)
        {
            r_disp[0] = x;
            for (y = -BoxRad; y <= BoxRad; y++)
                {
                    r_disp[1] = y;
                    for (z = -BoxRad; z <= BoxRad; z++)
                        {
                            r_disp[2] = z;
                            LatPos_add_wPBC(r_chck, r_pos_0, r_disp);
                            //                for (j = 0; j < POS_MAX; j++){
                            //                    r_chck[j] = r_pos_0[j] + r_disp[j];
                            //                    r_chck[j] = r_chck[j] < 0 ? r_chck[j] +
                            //                    naBoxSize_glb[j] : r_chck[j]; r_chck[j] =
                            //                    r_chck[j] >= naBoxSize_glb[j] ? r_chck[j] -
                            //                    naBoxSize_glb[j] : r_chck[j];
                            //                }
                            secBi = naTotLattice_glb[Lat_Ind_FromVec(r_chck)];
                            if (secBi != -1 && secBi != beadID)
                                {
                                    resj = bead_info_glb[secBi][BEAD_TYPE];
                                    xDis = sqrtf((float) (x * x + y * y + z * z));

                                    bead_check = 0;
                                    for (i = 0; i < list_size; i++)
                                        {
                                            if (secBi == bead_list[i])
                                                {
                                                    bead_check = 1;
                                                    break;
                                                }
                                        }
                                    if (bead_check == 1)
                                        { // Intra-list
                                            totEn +=
                                                Energy_Iso_Ovlp(resi, resj, xDis) / 2.; // / xDis / xDis / xDis /2.;

                                            totEn += Energy_Iso_Cont(resi, resj, xDis) / 2.;
                                        }
                                    else
                                        {                                               // Inter-list
                                            totEn += Energy_Iso_Ovlp(resi, resj, xDis); // / xDis / xDis / xDis;

                                            totEn += Energy_Iso_Cont(resi, resj, xDis);
                                        }
                                }
                            else if (secBi == -1)
                                {
                                    if (abs(x) <= 1 && abs(y) <= 1 && abs(z) <= 1)
                                        { // Want solvation radius to be 1
                                            totEn += faEnergy_glb[resi][resi][E_F_SOL];
                                        }
                                }
                        }
                }
        }

    return totEn;
}

/// Energy_Total_System: Calculate the total energy of the system.
void Energy_Total_System(void)
{
    int i; // Indecies
    // initialization
    for (i = 0; i < MAX_E; i++)
        {
            faCurrEn_glb[i] = 0.0f;
        }
    // printf("We have %d beads and %d chains", tot_beads_glb, tot_chains_glb);

    for (i = 0; i < tot_beads_glb; i++)
        {
            faCurrEn_glb[E_SC_SC] += Energy_Anisotropic(i);
        }

    int resi, ovlp_num, cont_num;

    for (i = 0; i < tot_beads_glb; i++)
        {
            resi = bead_info_glb[i][BEAD_TYPE];
            if (nBeadTypeCanCont_glb[resi])
                {
                    cont_num = NeighborSearch_ForCont(i, bead_info_glb[i], naOldContNeighs_glb, naOldOvlpNeighs_glb,
                                                      &ovlp_num);
                    faCurrEn_glb[E_CONT] += Energy_OfCont_wNeighList(i, naOldContNeighs_glb, cont_num);
                }
            else if (nBeadTypeCanOvlp_glb[resi] || nBeadTypeCanFSol_glb[resi])
                {
                    ovlp_num = NeighborSearch_ForOvlp(i, bead_info_glb[i], naOldOvlpNeighs_glb);
                    faCurrEn_glb[E_OVLP] += Energy_OfOvlp_wNeighList(i, naOldOvlpNeighs_glb, ovlp_num);
                    faCurrEn_glb[E_F_SOL] += (float) (26 - ovlp_num) * faEnergy_glb[resi][resi][E_F_SOL];
                }

        }

    if (bSystemHasTopo_glb)
        {
            for (i = 0; i < tot_beads_glb; i++)
                {
                    faCurrEn_glb[E_STIFF] += Energy_Topo_Angle(i);
                }
        }

    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = 0; i < tot_beads_glb; i++)
                {
                    faCurrEn_glb[E_BIAS] += Energy_BiasingPotential(i);
                }
        }

    // Taking care of double-counting energies.
    faCurrEn_glb[E_SC_SC] *= 0.5f;
    faCurrEn_glb[E_OVLP] *= 0.5f;
    faCurrEn_glb[E_CONT] *= 0.5f;
    for (i = 1; i < MAX_E; i++)
        {
            faCurrEn_glb[E_TOT] += faCurrEn_glb[i];
        }
}

/// Energy_Of_Chain - calculates the total energy of a molecule. Takes care of
/// double-counting intra interactions \param chainID - ID of the molecule to
/// calculate the energy. \return The total aniso + isotropic energy of this
/// chain.
float Energy_Of_Chain(const int chainID)
{ // Calculates the energy of the given chain.
    float totEn = 0.0;
    int i; // Looping index
    int fB = chain_info_glb[chainID][CHAIN_START];
    int lB = chain_info_glb[chainID][CHAIN_START] + chain_info_glb[chainID][CHAIN_LENGTH];
    for (i = fB; i < lB; i++)
        {
            totEn += Energy_Anisotropic_For_Chain(i) + Energy_Isotropic_For_Chain(i);
        }

    return totEn;
}

/// Energy_Of_Chain_Self - only calculates the intra-molecular energy of this
/// molecule \param chainID - ID of the molecule to calculate the energy.
/// \return The total aniso + isotropic energy of this chain.
float Energy_Of_Chain_Self(const int chainID)
{ // Calculates only intra-molecular interactions
    float totEn = 0.0;
    int i; // Looping index
    int fB = chain_info_glb[chainID][CHAIN_START];
    int lB = chain_info_glb[chainID][CHAIN_START] + chain_info_glb[chainID][CHAIN_LENGTH];
    for (i = fB; i < lB; i++)
        {
            totEn += Energy_Anisotropic_Self(i) + Energy_Isotropic_Self(i);
        }

    return totEn / 2.;
}

/// Energy_Of_Chain_OLD - calculates the total energy of a molecule. Double
/// counts intra-molecular energies. \param chainID - ID of the molecule to
/// calculate the energy. \return The total aniso + isotropic energy of this
/// chain. Note that this sub-routine is dumb and does not account for double
/// counting within the same chain. As such, this function should be used in
/// both the old and new energy calculations so that the energy difference is
/// still correct.
float Energy_Of_Chain_OLD(const int chainID)
{ // Calculates the energy of the given chain
    float totEn = 0.0;
    int i; // Looping index

    for (i = chain_info_glb[chainID][CHAIN_START];
         i < chain_info_glb[chainID][CHAIN_START] + chain_info_glb[chainID][CHAIN_LENGTH]; i++)
        {
            totEn += Energy_Anisotropic(i) + Energy_Isotropic(i);
        }

    return totEn;
}

void Energy_Iso_ForLocal(const int beadID, const int resi, const int* r_pos0, long double* oldEn, long double* newEn,
                         int* ovlp_num, int* cont_num, int* ovlp_neighs, int* cont_neighs)
{

    *ovlp_num = 0;
    *cont_num = 0;
    if (nBeadTypeCanCont_glb[resi])
        { // CONT neighbors.
            *cont_num = NeighborSearch_ForCont(beadID, r_pos0, cont_neighs, ovlp_neighs, ovlp_num);
        }
    else if (nBeadTypeIsSticker_glb[resi] || nBeadTypeCanOvlp_glb[resi] || nBeadTypeCanFSol_glb[resi])
        {
            // OVLP neighbors.
            *ovlp_num = NeighborSearch_ForOvlp(beadID, r_pos0, ovlp_neighs);
        }

    if (nBeadTypeCanCont_glb[resi])
        { // CONT energy.
            *oldEn = *oldEn + Energy_OfCont_wNeighList(beadID, cont_neighs, *cont_num);
        }

    if (nBeadTypeCanOvlp_glb[resi])
        { // OVLP energy.
            *oldEn = *oldEn + Energy_OfOvlp_wNeighList(beadID, ovlp_neighs, *ovlp_num);
        }

    if (nBeadTypeCanFSol_glb[resi])
        { // Solvation energy.
            *oldEn = *oldEn + (float) (26 - *ovlp_num) * faEnergy_glb[resi][resi][E_F_SOL];
            *newEn = *newEn + Energy_ofSol_wNeighList(ovlp_neighs, *ovlp_num);
        }
}

void Energy_Iso_ForLocalEquil(const int beadID, const int resi, const int* r_pos0, long double* oldEn,
                              long double* newEn, int* ovlp_num, int* cont_num, int* ovlp_neighs, int* cont_neighs)
{

    *ovlp_num = 0;
    *cont_num = 0;
    if (nBeadTypeCanCont_glb[resi])
        { // CONT neighbors.
            *cont_num = NeighborSearch_ForCont(beadID, r_pos0, cont_neighs, ovlp_neighs, ovlp_num);
        }
    else if (nBeadTypeCanOvlp_glb[resi] || nBeadTypeCanFSol_glb[resi])
        {
            // OVLP neighbors.
            *ovlp_num = NeighborSearch_ForOvlp(beadID, r_pos0, ovlp_neighs);
        }

    if (nBeadTypeCanCont_glb[resi])
        { // CONT energy.
            *oldEn = *oldEn + Energy_OfCont_wNeighList(beadID, cont_neighs, *cont_num);
        }

    if (nBeadTypeCanOvlp_glb[resi])
        { // OVLP energy.
            *oldEn = *oldEn + Energy_OfOvlp_wNeighList(beadID, ovlp_neighs, *ovlp_num);
        }

    if (nBeadTypeCanFSol_glb[resi])
        { // Solvation energy.
            *oldEn = *oldEn + (float) (26 - *ovlp_num) * faEnergy_glb[resi][resi][E_F_SOL];
            *newEn = *newEn + Energy_ofSol_wNeighList(ovlp_neighs, *ovlp_num);
        }
}

void Energy_Iso_ForChains(const int beadID, long double* oldEn, long double* newEn, int* ovlp_num, int* cont_num,
                          int* ovlp_neighs, int* cont_neighs)
{
    const int resi            = bead_info_glb[beadID][BEAD_TYPE];
    const int r_pos0[POS_MAX] = {bead_info_glb[beadID][0], bead_info_glb[beadID][1], bead_info_glb[beadID][2]};
    *ovlp_num                 = 0;
    *cont_num                 = 0;
    if (nBeadTypeCanCont_glb[resi])
        { // CONT neighbors.
            *cont_num = NeighborSearch_ForCont(beadID, r_pos0, cont_neighs, ovlp_neighs, ovlp_num);
        }
    else if (nBeadTypeIsSticker_glb[resi] || nBeadTypeCanOvlp_glb[resi] || nBeadTypeCanFSol_glb[resi])
        {
            // OVLP neighbors.
            *ovlp_num = NeighborSearch_ForOvlp(beadID, r_pos0, ovlp_neighs);
        }

    if (nBeadTypeCanCont_glb[resi])
        { // CONT energy.
            *oldEn = *oldEn + Energy_OfCont_wNeighList_ForChains(beadID, cont_neighs, *cont_num);
        }

    if (nBeadTypeCanOvlp_glb[resi])
        { // OVLP energy.
            *oldEn = *oldEn + Energy_OfOvlp_wNeighList_ForChains(beadID, ovlp_neighs, *ovlp_num);
        }

    if (nBeadTypeCanFSol_glb[resi])
        { // Solvation energy.
            *oldEn = *oldEn + (float) (26 - *ovlp_num) * faEnergy_glb[resi][resi][E_F_SOL];
            *newEn = *newEn + Energy_ofSol_wNeighList(ovlp_neighs, *ovlp_num);
        }
}

void Energy_Iso_ForChainsEquil(const int beadID, long double* oldEn, long double* newEn, int* ovlp_num, int* cont_num,
                               int* ovlp_neighs, int* cont_neighs)
{
    const int resi            = bead_info_glb[beadID][BEAD_TYPE];
    const int r_pos0[POS_MAX] = {bead_info_glb[beadID][0], bead_info_glb[beadID][1], bead_info_glb[beadID][2]};
    *ovlp_num                 = 0;
    *cont_num                 = 0;
    if (nBeadTypeCanCont_glb[resi])
        { // CONT neighbors.
            *cont_num = NeighborSearch_ForCont(beadID, r_pos0, cont_neighs, ovlp_neighs, ovlp_num);
        }
    else if (nBeadTypeCanOvlp_glb[resi] || nBeadTypeCanFSol_glb[resi])
        {
            // OVLP neighbors.
            *ovlp_num = NeighborSearch_ForOvlp(beadID, r_pos0, ovlp_neighs);
        }

    if (nBeadTypeCanCont_glb[resi])
        { // CONT energy.
            *oldEn = *oldEn + Energy_OfCont_wNeighList_ForChains(beadID, cont_neighs, *cont_num);
        }

    if (nBeadTypeCanOvlp_glb[resi])
        { // OVLP energy.
            *oldEn = *oldEn + Energy_OfOvlp_wNeighList_ForChains(beadID, ovlp_neighs, *ovlp_num);
        }

    if (nBeadTypeCanFSol_glb[resi])
        { // Solvation energy.
            *oldEn = *oldEn + (float) (26 - *ovlp_num) * faEnergy_glb[resi][resi][E_F_SOL];
            *newEn = *newEn + Energy_ofSol_wNeighList(ovlp_neighs, *ovlp_num);
        }
}

void Energy_Iso_ForCoLocal(const int thisBeadID, const int otherBeadID, const int* r_pos0, long double* oldEn,
                           long double* newEn, int* ovlp_num, int* cont_num, int* ovlp_neighs, int* cont_neighs)
{
    *ovlp_num      = 0;
    *cont_num      = 0;
    const int resi = bead_info_glb[thisBeadID][BEAD_TYPE];
    if (nBeadTypeCanCont_glb[resi])
        { // CONT neighbors.
            *cont_num = NeighborSearch_ForCont(thisBeadID, r_pos0, cont_neighs, ovlp_neighs, ovlp_num);
        }
    else
        { // Since CoLocal assumes bonded beads, it is a sticker.
            // OVLP neighbors.
            *ovlp_num = NeighborSearch_ForOvlp(thisBeadID, r_pos0, ovlp_neighs);
        }

    if (nBeadTypeCanCont_glb[resi])
        { // CONT energy.
            *oldEn = *oldEn +
                     Energy_OfCont_wNeighList_ForRange(thisBeadID, otherBeadID, otherBeadID, cont_neighs, *cont_num);
        }

    if (nBeadTypeCanOvlp_glb[resi])
        { // OVLP energy.
            *oldEn = *oldEn +
                     Energy_OfOvlp_wNeighList_ForRange(thisBeadID, otherBeadID, otherBeadID, ovlp_neighs, *ovlp_num);
        }

    if (nBeadTypeCanFSol_glb[resi])
        { // Solvation energy.
            *oldEn = *oldEn + (float) (26 - *ovlp_num) * faEnergy_glb[resi][resi][E_F_SOL];
            *newEn = *newEn + Energy_ofSol_wNeighList(ovlp_neighs, *ovlp_num);
        }
}

void Energy_Iso_ForLists(const int beadIdx, int const listSize, const int beadList[MAX_BONDS + 1],
                         const int beadPos[MAX_BONDS + 1][POS_MAX], long double* oldEn, long double* newEn,
                         int* ovlp_num, int* cont_num, int* ovlp_neighs, int* cont_neighs)
{

    *ovlp_num                 = 0;
    *cont_num                 = 0;
    const int beadID          = beadList[beadIdx];
    const int resi            = bead_info_glb[beadID][BEAD_TYPE];
    const int r_pos0[POS_MAX] = {beadPos[beadIdx][0], beadPos[beadIdx][1], beadPos[beadIdx][2]};

    if (nBeadTypeCanCont_glb[resi])
        { // CONT neighbors.
            *cont_num = NeighborSearch_ForCont(beadID, r_pos0, cont_neighs, ovlp_neighs, ovlp_num);
        }
    else if (nBeadTypeIsSticker_glb[resi] || nBeadTypeCanOvlp_glb[resi] || nBeadTypeCanFSol_glb[resi])
        {
            // OVLP neighbors.
            *ovlp_num = NeighborSearch_ForOvlp(beadID, r_pos0, ovlp_neighs);
        }

    if (nBeadTypeCanFSol_glb[resi])
        { // Solvation energy.
            *oldEn = *oldEn + (float) (26 - *ovlp_num) * faEnergy_glb[resi][resi][E_F_SOL];
            *newEn = *newEn + Energy_ofSol_wNeighList(ovlp_neighs, *ovlp_num);
        }

    if (nBeadTypeCanCont_glb[resi])
        { // CONT energy.
            *oldEn = *oldEn + Energy_OfCont_wNeighList_ForLists(beadID, listSize, beadList, cont_neighs, *cont_num);
        }

    if (nBeadTypeCanOvlp_glb[resi])
        { // OVLP energy.
            *oldEn = *oldEn + Energy_OfOvlp_wNeighList_ForLists(beadID, listSize, beadList, ovlp_neighs, *ovlp_num);
        }
}

void Energy_Iso_ForRange(const int beadID, const int smallestBead, const int largestBead, long double* oldEn,
                         long double* newEn, int* ovlp_num, int* cont_num, int* ovlp_neighs, int* cont_neighs)
{
    *ovlp_num                 = 0;
    *cont_num                 = 0;
    const int resi            = bead_info_glb[beadID][BEAD_TYPE];
    const int r_pos0[POS_MAX] = {bead_info_glb[beadID][0], bead_info_glb[beadID][1], bead_info_glb[beadID][2]};

    if (nBeadTypeCanCont_glb[resi])
        { // CONT neighbors.
            *cont_num = NeighborSearch_ForCont(beadID, r_pos0, cont_neighs, ovlp_neighs, ovlp_num);
        }
    else if (nBeadTypeIsSticker_glb[resi] || nBeadTypeCanOvlp_glb[resi] || nBeadTypeCanFSol_glb[resi])
        { // Since CoLocal assumes bonded
          // beads, it is a sticker.
            // OVLP neighbors.
            *ovlp_num = NeighborSearch_ForOvlp(beadID, r_pos0, ovlp_neighs);
        }

    if (nBeadTypeCanCont_glb[resi])
        { // CONT energy.
            *oldEn =
                *oldEn + Energy_OfCont_wNeighList_ForRange(beadID, smallestBead, largestBead, cont_neighs, *cont_num);
        }

    if (nBeadTypeCanOvlp_glb[resi])
        { // OVLP energy.
            *oldEn =
                *oldEn + Energy_OfOvlp_wNeighList_ForRange(beadID, smallestBead, largestBead, ovlp_neighs, *ovlp_num);
        }

    if (nBeadTypeCanFSol_glb[resi])
        { // Solvation energy.
            *oldEn = *oldEn + (float) (26 - *ovlp_num) * faEnergy_glb[resi][resi][E_F_SOL];
            *newEn = *newEn + Energy_ofSol_wNeighList(ovlp_neighs, *ovlp_num);
        }
}

void Energy_Iso_ForRangeEquil(const int beadID, const int smallestBead, const int largestBead, long double* oldEn,
                              long double* newEn, int* ovlp_num, int* cont_num, int* ovlp_neighs, int* cont_neighs)
{
    *ovlp_num                 = 0;
    *cont_num                 = 0;
    const int resi            = bead_info_glb[beadID][BEAD_TYPE];
    const int r_pos0[POS_MAX] = {bead_info_glb[beadID][0], bead_info_glb[beadID][1], bead_info_glb[beadID][2]};

    if (nBeadTypeCanCont_glb[resi])
        { // CONT neighbors.
            *cont_num = NeighborSearch_ForCont(beadID, r_pos0, cont_neighs, ovlp_neighs, ovlp_num);
        }
    else if (nBeadTypeCanOvlp_glb[resi] || nBeadTypeCanFSol_glb[resi])
        { // Since CoLocal assumes bonded
          // beads, it is a sticker.
            // OVLP neighbors.
            *ovlp_num = NeighborSearch_ForOvlp(beadID, r_pos0, ovlp_neighs);
        }

    if (nBeadTypeCanCont_glb[resi])
        { // CONT energy.
            *oldEn =
                *oldEn + Energy_OfCont_wNeighList_ForRange(beadID, smallestBead, largestBead, cont_neighs, *cont_num);
        }

    if (nBeadTypeCanOvlp_glb[resi])
        { // OVLP energy.
            *oldEn =
                *oldEn + Energy_OfOvlp_wNeighList_ForRange(beadID, smallestBead, largestBead, ovlp_neighs, *ovlp_num);
        }

    if (nBeadTypeCanFSol_glb[resi])
        { // Solvation energy.
            *oldEn = *oldEn + (float) (26 - *ovlp_num) * faEnergy_glb[resi][resi][E_F_SOL];
            *newEn = *newEn + Energy_ofSol_wNeighList(ovlp_neighs, *ovlp_num);
        }
}

float Energy_Topo_Angle_ForList(const int bondNum, const int* bondList)
{
    float totEn = 0.f;
    int i, tmpBead;

    for (i = 0; i < bondNum; ++i)
        {
            tmpBead = bondList[i];
            totEn += Energy_Topo_Angle(tmpBead);
        }

    return totEn;
}
