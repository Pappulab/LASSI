#include "mcmove.h"
#include "cluster.h"
#include "energy.h"
#include "structure.h"

/// MC_Step - given the MC move frequencies, picks which move to try.
/// \param fMCTemp
/// \return bAccept - the acceptance state where MCMove=bAccept/12 and
/// MCAccept=bAccept%2.
int MC_Step(float fMCTemp)
{
    int mode = 0; // Which move to do
    int nAccept;  // Used in MC steps
    float prob = (float) rand() / (float) RAND_MAX;
    int i; // Internal iterator, and int.
    for (i = 0; i < MAX_MV; i++)
        {
            if (prob < faMCFreq_glb[i])
                { // Pick this move!
                    mode = i;
                    break;
                }
        }
    nAccept = 0;

    // Performing the appropriate move.
    switch (mode)
        {
            // Translation
            case MV_TRANS:
                i       = rand() % tot_chains_glb; // Pick random chain
                nAccept = Move_Trans(i, fMCTemp);
                break;

                // Cluster translation
            case MV_CLSTR:
                nAccept = Move_Clus_Network(fMCTemp);
                break;

                // Cluster translation iff ClusSize <= 5
            case MV_SMCLSTR:
                i       = rand() % tot_chains_glb; // Pick a random chain
                nAccept = Move_SmallClus_Network(i, fMCTemp);
                break;

                // Change Rotational State
            case MV_STROT:
                i       = rand() % tot_beads_glb;
                nAccept = Move_Rot(i, fMCTemp);
                break;

                // Local Move
            case MV_LOCAL:
                i       = rand() % tot_beads_glb;
                nAccept = Move_Local(i, fMCTemp);
                break;

                // Slithering Snake
            case MV_SNAKE:
                i       = rand() % tot_chains_glb;
                nAccept = Move_Snake(i, fMCTemp);
                break;

                // Double Pivot
            case MV_DBPVT:
                i       = rand() % tot_beads_glb; // Pick a random bead
                nAccept = Move_DbPvt(i, fMCTemp);
                break;

                // Co-Local
            case MV_COLOCAL:
                i       = rand() % tot_beads_glb; // Pick a random bead
                nAccept = Move_CoLocal(i, fMCTemp);
                break;

                // Shake Move
            case MV_MTLOCAL:
                i       = rand() % tot_beads_glb;
                nAccept = Move_MultiLocal(i, fMCTemp);
                break;

                // Pivot
            case MV_PIVOT:
                i       = rand() % tot_chains_glb;
                nAccept = Move_Pivot(i, fMCTemp);
                break;

                // Branched Rotation
            case MV_BRROT:
                i       = rand() % tot_chains_glb;
                nAccept = Move_BranchedRot(i, fMCTemp);
                break;
            case MV_PR_SMCLSTR:
                i       = rand() % tot_chains_glb;
                nAccept = Move_SmallClus_Proximity(i, fMCTemp);
                break;

            case MV_PR_CLSTR:
                nAccept = Move_Clus_Proximity(fMCTemp);
            default:
                nAccept = 0;
                break;
        }

    naMCAccepMat_glb[nAccept][mode]++; // Just adding which move got accepted/rejected.
    return mode * 12 + nAccept;
}

/// MC_Step_Therm - given the MC move frequencies, picks non-anisotropic
/// variants of the moves. \param fMCTemp \return bAccept - the acceptance state
/// where MCMove=bAccept/12 and MCAccept=bAccept%2.
int MC_Step_Therm(float fMCTemp)
{
    int mode = 0; // Which move to do
    int nAccept;  // Used in MC steps
    float prob = (float) rand() / (float) RAND_MAX;
    int i; // Internal iterator, and int.
    for (i = 0; i < MAX_MV; i++)
        {
            if (prob < faMCFreq_glb[i])
                { // Pick this move!
                    mode = i;
                    break;
                }
        }
    // Start by assuming failure.
    nAccept = 0;

    // Performing the appropriate move.
    switch (mode)
        {
            // translation of a random chain
            case MV_TRANS:
                i       = rand() % tot_chains_glb; // Pick random chain
                nAccept = Move_Trans_Equil(i, fMCTemp);
                break;

                // cluster translation: moves largest cluster to another spot.
            case MV_CLSTR:                         // In equil, just a translation of the chain
                i       = rand() % tot_chains_glb; // Pick random chain
                nAccept = Move_Trans_Equil(i, fMCTemp);
                break;

            case MV_SMCLSTR:                       // In equil, just a translation of the chain
                i       = rand() % tot_chains_glb; // Pick random chain
                nAccept = Move_Trans_Equil(i, fMCTemp);
                break;

                // face change
            case MV_STROT: // In equil, just a local move
                i       = rand() % tot_beads_glb;
                nAccept = Move_Local_Equil(i, fMCTemp);
                break;

                // local moves
            case MV_LOCAL:
                i       = rand() % tot_beads_glb;
                nAccept = Move_Local_Equil(i, fMCTemp);
                break;

                // slithering snake
            case MV_SNAKE:
                i       = rand() % tot_chains_glb;
                nAccept = Move_Snake_Equil(i, fMCTemp);
                break;

                // double pivot
            case MV_DBPVT:
                i       = rand() % tot_beads_glb; // Pick a random bead
                nAccept = Move_DbPvt(i, fMCTemp);
                break;

                // co-local
            case MV_COLOCAL: // In equil, just a local move
                i       = rand() % tot_beads_glb;
                nAccept = Move_Local_Equil(i, fMCTemp);
                break;

                // shake
            case MV_MTLOCAL:
                i       = rand() % tot_beads_glb;
                nAccept = Move_MultiLocal_Equil(i, fMCTemp);
                break;

                // pivot
            case MV_PIVOT:
                i       = rand() % tot_chains_glb;
                nAccept = Move_Pivot_Equil(i, fMCTemp);
                break;

                // branched rotate
            case MV_BRROT:
                i       = rand() % tot_chains_glb;
                nAccept = Move_BranchedRot_Equil(i, fMCTemp);
                break;

            case MV_PR_SMCLSTR:
                // In equil, just a translation of the chain
                i       = rand() % tot_chains_glb; // Pick random chain
                nAccept = Move_Trans_Equil(i, fMCTemp);
                break;

            default:
                nAccept = 0;
                break;
        }

    naMCAccepMat_glb[nAccept][mode]++; // Just adding which move got accepted/rejected.
    return mode * 12 + nAccept;
}

/*
 * The following functions are sub-routines that perform the various Monte-Carlo
 * (MC) moves.
 */

/// Move_Rot - searches the 3^3-1=26 possible sites for forming a bond, and
/// performs Metropolis-Hastings. Rejects move if beadID cannot forms bonds.
/// \param beadID
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_Rot(int beadID, float MyTemp)
{
    // Performs a rotational MC-Move on beadID
    int bAccept;    // Used in MC steps
    int resi, resj; // To track bead types.
    // Firstly -- make sure that bead i can even have rotational states.
    resi = bead_info_glb[beadID][BEAD_TYPE];
    // printf("Beginning ROT\n");
    if (nBeadTypeIsSticker_glb[resi] == 0)
        { // Skip beads that cannot rotate!
            bAccept = 0;
            return bAccept;
        }

    lLDub MCProb, oldEn, newEn; // For Metropolis Hastings
    oldEn = 0.;
    newEn = 0.;
    int i; // General looping iterators
    int yTemp;
    int FWWeight;

    if (bead_info_glb[beadID][BEAD_FACE] != -1)
        {
            resj  = bead_info_glb[bead_info_glb[beadID][BEAD_FACE]][BEAD_TYPE]; // This is who I'm currently bonded to
            oldEn = (lLDub) faEnergy_glb[resi][resj][E_SC_SC];
        }
    OP_ShuffleRotIndecies();
    FWWeight = Check_RotStatesOld(beadID, resi, MyTemp);
    OP_NormalizeRotState(0, FWWeight);

    yTemp = OP_PickRotState(FWWeight);

    if (yTemp != -1)
        { // There is a bead here
            resj  = bead_info_glb[yTemp][BEAD_TYPE];
            newEn = (lLDub) faEnergy_glb[resi][resj][E_SC_SC];
        }
    // See if we can accept this move
    MCProb      = (lLDub) rand() / (lLDub) RAND_MAX;
    lLDub MHAcc = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub) MyTemp);
    if (MCProb < MHAcc)
        { // Accept this state
            if (bead_info_glb[beadID][BEAD_FACE] != -1)
                {                                                                    // Break old bond
                    bead_info_glb[bead_info_glb[beadID][BEAD_FACE]][BEAD_FACE] = -1; // Breaking bond with old partner
                }
            bead_info_glb[beadID][BEAD_FACE] = yTemp;
            if (yTemp != -1)
                { // New bond
                    bead_info_glb[yTemp][BEAD_FACE] = beadID;
                }
            bAccept = 1;
            return bAccept;
        }
    else
        {
            bAccept = 0;
            return bAccept;
        }
}

/// Move_Local - performs a biased local move on beadID by:
/// 1. Seeing if there is a free spot to move to in a <+-2,+-2,+-2> random
/// location. Rejects if no space found/exists.
/// 2. Calculating the Rosenbluth weight by searching the 3^3-1=26 possible
/// bonding locations, and calculate the rest of the energy.
/// 3. Move the bead. Calculate the Rosenbluth weight, and energy. Propose a
/// bond and perform Metropolis-Hastings
/// \param beadID
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_Local(int beadID, float MyTemp)
{ // Performs a local translation MC-move on beadID

    int bAccept = 0;                                         // Used in MC steps
    int yTemp;                                               // Random numbers to store things
    int r_pos0[POS_MAX], r_posNew[POS_MAX], r_disp[POS_MAX]; // Vectors to stores coordinates.

    // Attempt to find an empty lattice point.
    LatPos_copy(r_pos0, bead_info_glb[beadID]);
    LatPos_gen_rand_wRad(r_disp, 2);
    //    LatPos_gen_rand_wRad(r_disp, linker_len_glb[beadID][0]);
    LatPos_add_wPBC(r_posNew, r_pos0, r_disp);

    // Checking to see validity of new point.
    yTemp = Check_MoveBeadTo(r_posNew);
    if (yTemp == 1)
        { // This means we found an empty lattice site. So let's
          // check if the linkers are okay.
            yTemp = Check_LinkerConstraint(beadID, r_posNew);
        }
    if (yTemp != 1)
        { // No space -- reject.
            bAccept = 0;
            return bAccept;
        }
    // Have successfully found a good lattice spot. Let's perform the usual
    // Metropolis-Hastings shenanigans.

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;

    const int resi = bead_info_glb[beadID][BEAD_TYPE];

    lLDub oldEn = nBiasPotential_Mode_glb == -1 ? 0. : Energy_BiasingPotential(beadID);
    lLDub newEn = 0.;

    Energy_Iso_ForLocal(beadID, resi, r_pos0, &oldEn, &newEn, &old_ovlp_num, &old_cont_num, naOldOvlpNeighs_glb,
                        naOldContNeighs_glb);

    int bondList[MAX_BONDS + 1];
    int bondNum = bSystemHasTopo_glb ? OP_GetTopoBonds(beadID, bondList) : 0;
    bondNum     = BeadList_CanTopoAngle(bondNum, bondList);

    oldEn += bondNum ? Energy_Topo_Angle_ForList(bondNum, bondList) : 0.;

    lLDub BWRos = MC_RosenbluthSampling_ForLocal_AtOld(beadID, resi, &oldEn, old_ovlp_num);

    OP_System_MoveBeadTo(beadID, r_posNew);

    newEn += nBiasPotential_Mode_glb == -1 ? 0.f : Energy_BiasingPotential(beadID);

    Energy_Iso_ForLocal(beadID, resi, r_posNew, &newEn, &oldEn, &new_ovlp_num, &new_cont_num, naNewOvlpNeighs_glb,
                        naNewContNeighs_glb);

    newEn += bondNum ? Energy_Topo_Angle_ForList(bondNum, bondList) : 0.;

    lLDub FWRos = MC_RosenbluthSampling_ForLocal_AtNew(beadID, resi, &yTemp, &newEn, new_ovlp_num);
    int resj;
    if (yTemp != -1)
        {
            resj = bead_info_glb[yTemp][BEAD_TYPE];
            newEn += faEnergy_glb[resi][resj][E_SC_SC];
        }

    lLDub MCProb = (lLDub) rand() / (lLDub) RAND_MAX;
    lLDub MHAcc  = OP_GenMHValue(FWRos, BWRos, oldEn - newEn, (lLDub) MyTemp);
    // Accept or reject this state
    if (MCProb < MHAcc)
        { // Accept this state
            if (bead_info_glb[beadID][BEAD_FACE] != -1)
                { // Breaking old bond
                    bead_info_glb[bead_info_glb[beadID][BEAD_FACE]][BEAD_FACE] = -1;
                }
            bead_info_glb[beadID][BEAD_FACE] = yTemp;
            if (yTemp != -1)
                { // Making new bond
                    bead_info_glb[yTemp][BEAD_FACE] = beadID;
                }
            bAccept = 1;
            return bAccept;
        }
    else
        {
            OP_System_MoveBeadTo_Inv(beadID);
            bAccept = 0;
            return bAccept;
        }
}

/// Move_Snake - performs a biased reptation move on chainID by:
/// 1. Randomly picking which end to reptate.
/// 2. Performing sort of a local move by finding an empty spot to move to,
/// while recording
///    the Rosenbluth weights, but for the whole chain. Calculate the total
///    energy of the chain.
/// Break all the physical bonds.
/// 3. Moving the chain along and recalculating the total weights. While
/// calculating the total weights, also randomly assign a physical bond from the
/// boltzmann distribution. Calculate the new energy. Note that rather than
/// \f\$prod_{i=1}^{N}W_{i}\f$, I calculate the \f$\log_{10}\f$ so I calculate
/// the sum of smallish numbers rather than the product of large numbers.
/// 4. Perform the Metropolis-Hastings step.
/// The move is automatically rejected for molecules that are not chains.
/// The move is NOT smart enough to detect that the given chain has all the same
/// linker lengths within the chain.
/// TODO: In initialization, add a sub-routine that checks if the molecules have
/// the same linker lengths throughout.
/// \param chainID
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_Snake(int chainID, float MyTemp)
{ // Performs a slither MC-move on chainID

    int bAccept = 0; // Used in MC steps
    // Finding the bounds for looping over the molecule/chain
    const int firstB = chain_info_glb[chainID][CHAIN_START];
    const int lastB  = firstB + chain_info_glb[chainID][CHAIN_LENGTH];
    if (lastB - firstB == 1)
        { // This means we have a monomer. Reject the
          // move, because Local or Trans
            // moves should be the ones that move monomers.
            bAccept = 0;
            return bAccept;
        }
    else
        {
            if (nChainTypeIsLinear_glb[chain_info_glb[chainID][CHAIN_TYPE]] != 1)
                {
                    // If chain is not linear. Reject move because slithering will not
                    // work!
                    bAccept = 0;
                    return bAccept;
                }
        }

    // This chain is suitable to perform a slithering-snake move on.

    int i; // Loop iterator
    int resi;
    int yTemp; // Random numbers to store things
    int r_posNew[POS_MAX], r_posTmp1[POS_MAX],
        r_posTmp2[POS_MAX]; // Vectors to store positions.
    int angBead_old, angBead_new;
    lLDub MCProb        = (lLDub) rand() / (lLDub) RAND_MAX; // To decide if we slither forwards or backwards
    const char SnakeFwd = MCProb < 0.5f ? 1 : 0;             // 1 for Fwds, and 0 for Bkwds.

    if (SnakeFwd)
        { // Forwards slither, so lastB-1 (last bead) is anchor
            LatPos_gen_rand_wRad(r_posTmp1,
                                 linker_len_glb[lastB - 1][0]); // lastB-1 will be replaced by lastB-2
            LatPos_copy(r_posTmp2, bead_info_glb[lastB - 1]);
            LatPos_add_wPBC(r_posNew, r_posTmp1, r_posTmp2);
            yTemp       = Check_MoveBeadTo(r_posNew); // 0: there is no space, 1: there is space
            angBead_old = firstB + 1;
            angBead_new = lastB - 2;
        }
    else
        { // Backwards slither, so firstB is anchor
            LatPos_gen_rand_wRad(r_posTmp1,
                                 linker_len_glb[firstB][0]); // firstB will be replaced by firstB+1
            LatPos_copy(r_posTmp2, bead_info_glb[firstB]);
            LatPos_add_wPBC(r_posNew, r_posTmp1, r_posTmp2);
            yTemp       = Check_MoveBeadTo(r_posNew); // 0: there is no space, 1: there is space
            angBead_new = firstB + 1;
            angBead_old = lastB - 2;
        }

    if (yTemp == 0)
        { // Couldn't find a spot, so reject the damn move
            bAccept = 0;
            return bAccept;
        }

    // Let's remember where this chain exists.
    OP_CopyBeadsToOld(firstB, lastB);

    lLDub newEn = 0.;
    lLDub oldEn = 0.;
    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = firstB; i < lastB; i++)
                {
                    oldEn += Energy_BiasingPotential(i);
                }
        }

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;

    lLDub BSum = 0.;
    for (i = firstB; i < lastB; i++)
        {
            Energy_Iso_ForChains(i, &oldEn, &newEn, &old_ovlp_num, &old_cont_num, naOldOvlpNeighs_glb,
                                 naOldContNeighs_glb);
            resi = bead_info_glb[i][BEAD_TYPE];
            BSum += MC_RosenbluthSampling_ForChains_AtOld(i, resi, &oldEn, old_ovlp_num);
        }

    resi = bead_info_glb[angBead_old][BEAD_TYPE];
    oldEn += faEnergy_glb[resi][resi][E_STIFF] ? Energy_Topo_Angle(angBead_old) : 0.;

    if (SnakeFwd)
        { // Slithering the chain forwards in ID-space
            OP_System_Snake_SlitherFwd(firstB, lastB, r_posNew);
        }
    else
        { // Slithering backwards in ID-space
            OP_System_Snake_SlitherBck(firstB, lastB, r_posNew);
        }

    // Have slithered the chain whichever way
    for (i = firstB; i < lastB; i++)
        { // Break the old bonds!
            yTemp = bead_info_glb[i][BEAD_FACE];
            if (yTemp != -1)
                { // Need to break this bond
                    bead_info_glb[yTemp][BEAD_FACE] = -1;
                    bead_info_glb[i][BEAD_FACE]     = -1;
                }
        }

    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = firstB; i < lastB; i++)
                {
                    newEn += Energy_BiasingPotential(i);
                }
        }

    lLDub FSum = 0.;
    for (i = firstB; i < lastB; i++)
        {
            Energy_Iso_ForChains(i, &newEn, &oldEn, &new_ovlp_num, &new_cont_num, naNewOvlpNeighs_glb,
                                 naNewContNeighs_glb);
            resi = bead_info_glb[i][BEAD_TYPE];
            FSum += MC_RosenbluthSampling_ForChains_AtNew(i, resi, &yTemp, &newEn, new_ovlp_num);
            if (yTemp != -1)
                { // An appropriate partner has been selected. Form
                  // the bonds and add the energy
                    bead_info_glb[i][BEAD_FACE]     = yTemp;
                    bead_info_glb[yTemp][BEAD_FACE] = i;
                    newEn += Energy_Anisotropic_For_Chain(i);
                }
        }

    resi = bead_info_glb[angBead_new][BEAD_TYPE];
    newEn += faEnergy_glb[resi][resi][E_STIFF] ? Energy_Topo_Angle(angBead_new) : 0.;

    // Doing the Metropolis-Hastings thing
    MCProb      = (lLDub) rand() / (lLDub) RAND_MAX;
    lLDub MHAcc = OP_GenMHValue(FSum, BSum, oldEn - newEn, (lLDub) MyTemp);
    if (MCProb < MHAcc)
        { // Accept this state. Bonds have already been formed!
            bAccept = 1;
            return bAccept;
        }
    else
        {
            OP_RestoreChain_ForSnake(firstB, lastB);
            bAccept = 0;
            return bAccept;
        }
}

/// Move_Trans - performs a biased translation of chainID by:
/// 1. Seeing if there is a spot to move the entire chain in a
/// <+-L/2,+-L/2,+-L/2> random location. Move fails if no spot available.
/// 2. Calculate the Rosenbluth weights of the whole chain like Move_Snake().
/// 3. Move the chain and recalculate the total weights.
/// 4. Perform the Metropolis-Hastings step.
/// \param chainID
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_Trans(int chainID, float MyTemp)
{                        // Performs a translation move with orientational bias
    int bAccept = 0;     // Used in MC steps
    int yTemp;           // Random numbers to store things
    int r_disp[POS_MAX]; // Vectors to store coordinates.

    // All moves are L/2 radius
    LatPos_gen_rand_wRad(r_disp, naBoxSize_glb[0] / 2);

    yTemp = Check_ChainDisp(chainID, r_disp); // yTemp=0 means clash
    if (yTemp == 0)
        { // We have failed to find a good spot for this chain.
            bAccept = 0;
            return bAccept;
        }

    // We now have a chain which when moved does not overlap.

    int i, j; // Loop iterators

    int resi, resj, firstB, lastB;
    // Finding the bounds for looping over the molecule/chain
    firstB = chain_info_glb[chainID][CHAIN_START];
    lastB  = firstB + chain_info_glb[chainID][CHAIN_LENGTH];

    lLDub newEn = 0.;
    lLDub oldEn = 0.;
    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = firstB; i < lastB; i++)
                {
                    oldEn += Energy_BiasingPotential(i);
                }
        }

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;

    lLDub BSum = 0.;
    for (i = firstB; i < lastB; i++)
        {
            Energy_Iso_ForChains(i, &oldEn, &newEn, &old_ovlp_num, &old_cont_num, naOldOvlpNeighs_glb,
                                 naOldContNeighs_glb);
            resi = bead_info_glb[i][BEAD_TYPE];
            BSum += MC_RosenbluthSampling_ForChains_AtOld(i, resi, &oldEn, old_ovlp_num);
        }

    OP_System_DispChain_ForTrans(chainID, r_disp); // Moved the chain, broke bonds, and remembered stuff.

    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = firstB; i < lastB; i++)
                {
                    newEn += Energy_BiasingPotential(i);
                }
        }

    lLDub FSum = 0.;
    for (i = firstB; i < lastB; i++)
        {
            Energy_Iso_ForChains(i, &newEn, &oldEn, &new_ovlp_num, &new_cont_num, naNewOvlpNeighs_glb,
                                 naNewContNeighs_glb);
            resi = bead_info_glb[i][BEAD_TYPE];
            FSum += MC_RosenbluthSampling_ForChains_AtNew(i, resi, &yTemp, &newEn, new_ovlp_num);
            if (yTemp != -1)
                { // An appropriate partner has been selected. Form
                  // the bonds and add the energy
                    bead_info_glb[i][BEAD_FACE]     = yTemp;
                    bead_info_glb[yTemp][BEAD_FACE] = i;
                    newEn += Energy_Anisotropic_For_Chain(i);
                }
        }

    lLDub MCProb = (lLDub) rand() / (lLDub) RAND_MAX;
    lLDub MHAcc  = OP_GenMHValue(FSum, BSum, oldEn - newEn, (lLDub) MyTemp);
    if (MCProb < MHAcc)
        { // Accept the move. Remember that the bonds were
          // assigned above!
            bAccept = 1;
            return bAccept;
        }
    else
        {
            OP_System_RestoreChain_ForTrans(chainID);
            bAccept = 0;
            return bAccept;
        }
}

/// Move_Clus - translates the second largest cluster of the system by:
/// 1. Performs a total clustering analysis of the system and finds the second
/// largest cluster. If there are more than 1 second largest clusters, randomly
/// pick 1.
/// 2. Try to find an empty spot in a <+-L/2,+-L/2,+-L/2> random location.
/// Calculate old energy.
/// 3. Move the cluster over.
/// Note that since the definition of a cluster is based on the network of
/// existing physical bonds, the system energy, the clusters do not change post
/// translations.
/// 4. Calculate the new energy and perform the Metropolis-Hastings step.
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_Clus_Network(float MyTemp)
{
    // Attempts to move the second largest cluster

    // Attempts to move the second largest cluster

    int bAccept     = 0; // Used in MC steps, assume that move fails initially.
    int* naClusList = (int*) malloc(sizeof(int) * (tot_chains_glb + 1));

    const int ClusSize = ClusUtil_AnisoCluster_OfSystem_SecondLargest_ForMCMove(naClusList);

    if (ClusSize < 1)
        {
            bAccept = 0;
            free(naClusList);
            return bAccept;
        }

    int r_Disp[POS_MAX];
    int yTemp;
    int i, j;
    // Radii for translation moves. All moves are L/4 radius
    // I guess moving single chains around as well is not a bad idea

    LatPos_gen_rand_wRad(r_Disp, naBoxSize_glb[0] / 2);
    for (i = 0; i < ClusSize; i++)
        {
            yTemp = Check_ChainDisp(naClusList[i], r_Disp); // Checking for steric clash
            if (yTemp == 0)
                {
                    bAccept = 0;
                    free(naClusList);

                    // printf("End CLUS - No space\n");
                    return bAccept;
                }
        }
    // This means that there is no steric clash when the cluster is moved.

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;
    int thisChain, firstB, lastB;

    lLDub oldEn = 0.;
    lLDub newEn = 0.;
    if (nBiasPotential_Mode_glb != -1)
        {
            for (j = 0; j < ClusSize; j++)
                {
                    thisChain = naClusList[j];
                    firstB    = chain_info_glb[thisChain][CHAIN_START];
                    lastB     = firstB + chain_info_glb[thisChain][CHAIN_LENGTH];
                    for (i = firstB; i < lastB; i++)
                        {
                            oldEn += Energy_BiasingPotential(i);
                        }
                }
        }

    if (bSystemHasCont_glb || bSystemHasOvlp_glb || bSystemHasFSol_glb)
        {
            for (j = 0; j < ClusSize; j++)
                {
                    thisChain = naClusList[j];
                    firstB    = chain_info_glb[thisChain][CHAIN_START];
                    lastB     = firstB + chain_info_glb[thisChain][CHAIN_LENGTH];
                    for (i = firstB; i < lastB; i++)
                        {
                            Energy_Iso_ForChains(i, &oldEn, &newEn, &old_ovlp_num, &old_cont_num, naOldOvlpNeighs_glb,
                                                 naOldContNeighs_glb);
                        }
                }
        }

    for (j = 0; j < ClusSize; j++)
        {
            OP_System_DispChain(naClusList[j], r_Disp); // Moving the cluster properly
        }

    if (nBiasPotential_Mode_glb != -1)
        {
            for (j = 0; j < ClusSize; j++)
                {
                    thisChain = naClusList[j];
                    firstB    = chain_info_glb[thisChain][CHAIN_START];
                    lastB     = firstB + chain_info_glb[thisChain][CHAIN_LENGTH];
                    for (i = firstB; i < lastB; i++)
                        {
                            newEn += Energy_BiasingPotential(i);
                        }
                }
        }

    if (bSystemHasCont_glb || bSystemHasOvlp_glb || bSystemHasFSol_glb)
        {
            for (j = 0; j < ClusSize; j++)
                {
                    thisChain = naClusList[j];
                    firstB    = chain_info_glb[thisChain][CHAIN_START];
                    lastB     = firstB + chain_info_glb[thisChain][CHAIN_LENGTH];
                    for (i = firstB; i < lastB; i++)
                        {
                            Energy_Iso_ForChains(i, &newEn, &oldEn, &new_ovlp_num, &new_cont_num, naNewOvlpNeighs_glb,
                                                 naNewContNeighs_glb);
                        }
                }
        }

    const lLDub MCProb = (lLDub) rand() / (lLDub) RAND_MAX;
    const lLDub MHAcc  = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub) MyTemp);
    if (MCProb < MHAcc)
        {                // Accept this state
            bAccept = 1; // Accept the move
            free(naClusList);
            return bAccept;
            // printf("End CLUS - Yes\n");
        }
    else
        {
            for (i = 0; i < ClusSize; i++)
                {
                    OP_System_RestoreChain(naClusList[i]); // Placing  the cluster back properly
                }
            bAccept = 0;
            free(naClusList);
            return bAccept;
            // printf("End CLUS - Failed.\n");
        }
}

/// Move_SmallClus - translates this cluster only if it is smaller than 5 total
/// molecules by
/// 1. Performs a clustering analysis on chainID. Move fails if the cluster is
/// larger than 5.
/// 2. Try to find an empty spot in a <+-L/2,+-L/2,+-L/2> random location.
/// Calculate old energy.
/// 3. Move the cluster over.
/// Note that since the definition of a cluster is based on the network of
/// existing physical bonds, the system energy, the clusters do not change post
/// translations.
/// 4. Calculate the new energy and perform the Metropolis-Hastings step.
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_SmallClus_Network(int chainID, float MyTemp)
{
    // Performs a cluster move where a given chain and it's cluster are moved.
    // No new 'bonds' are made so the move is reversible....

    int bAccept = 0; // Used in MC steps, assume that move fails initially.

    //    const int ClusSize = Clus_Network_LimitedCluster(chainID); // Looking at everything that is connected to
    //    chainID

    char* oldHashTab = calloc(tot_chains_glb, sizeof(char));
    int* clusList    = calloc(tot_chains_glb, sizeof(int));

    const int ClusSize = ClusUtil_AnisoCluster_OfChain_wMaxSize(chainID, oldHashTab, clusList, nLimitedClusterSize_glb);
    if ((ClusSize < 2))
        {
            bAccept = 0;
            free(oldHashTab);
            free(clusList);
            return bAccept;
        }

    int r_Disp[POS_MAX];
    int yTemp;
    int i, j;
    // Radii for translation moves. All moves are L/2 radius
    // I guess moving single chains around as well is not a bad idea

    LatPos_gen_rand_wRad(r_Disp, naBoxSize_glb[0] / 2);
    for (i = 0; i < ClusSize; i++)
        {
            // printf("%d\n", naList_glb[i]);
            yTemp = Check_ChainDisp(clusList[i], r_Disp); // Checking for steric clash
            if (yTemp == 0)
                {
                    bAccept = 0;
                    // printf("End CLUS - No space\n");
                    free(oldHashTab);
                    free(clusList);
                    return bAccept;
                }
        }
    // This means that there is no steric clash when the cluster is moved.

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;
    int thisChain, firstB, lastB;

    lLDub oldEn = 0.;
    lLDub newEn = 0.;
    if (nBiasPotential_Mode_glb != -1)
        {
            for (j = 0; j < ClusSize; j++)
                {
                    thisChain = clusList[j];
                    firstB    = chain_info_glb[thisChain][CHAIN_START];
                    lastB     = firstB + chain_info_glb[thisChain][CHAIN_LENGTH];
                    for (i = firstB; i < lastB; i++)
                        {
                            oldEn += Energy_BiasingPotential(i);
                        }
                }
        }

    if (bSystemHasCont_glb || bSystemHasOvlp_glb || bSystemHasFSol_glb)
        {
            for (j = 0; j < ClusSize; j++)
                {
                    thisChain = clusList[j];
                    firstB    = chain_info_glb[thisChain][CHAIN_START];
                    lastB     = firstB + chain_info_glb[thisChain][CHAIN_LENGTH];
                    for (i = firstB; i < lastB; i++)
                        {
                            Energy_Iso_ForChains(i, &oldEn, &newEn, &old_ovlp_num, &old_cont_num, naOldOvlpNeighs_glb,
                                                 naOldContNeighs_glb);
                        }
                }
        }

    for (j = 0; j < ClusSize; j++)
        {
            OP_System_DispChain(clusList[j], r_Disp); // Moving the cluster properly
        }

    if (nBiasPotential_Mode_glb != -1)
        {
            for (j = 0; j < ClusSize; j++)
                {
                    thisChain = clusList[j];
                    firstB    = chain_info_glb[thisChain][CHAIN_START];
                    lastB     = firstB + chain_info_glb[thisChain][CHAIN_LENGTH];
                    for (i = firstB; i < lastB; i++)
                        {
                            newEn += Energy_BiasingPotential(i);
                        }
                }
        }

    if (bSystemHasCont_glb || bSystemHasOvlp_glb || bSystemHasFSol_glb)
        {
            for (j = 0; j < ClusSize; j++)
                {
                    thisChain = clusList[j];
                    firstB    = chain_info_glb[thisChain][CHAIN_START];
                    lastB     = firstB + chain_info_glb[thisChain][CHAIN_LENGTH];
                    for (i = firstB; i < lastB; i++)
                        {
                            Energy_Iso_ForChains(i, &newEn, &oldEn, &new_ovlp_num, &new_cont_num, naNewOvlpNeighs_glb,
                                                 naNewContNeighs_glb);
                        }
                }
        }
    const lLDub MCProb = (lLDub) rand() / (lLDub) RAND_MAX;
    const lLDub MHAcc  = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub) MyTemp);
    if (MCProb < MHAcc)
        {                // Accept this state
            bAccept = 1; // Accept the move
            free(oldHashTab);
            free(clusList);
            return bAccept;
            // printf("End CLUS - Yes\n");
        }
    else
        {
            for (i = 0; i < ClusSize; i++)
                {
                    OP_System_RestoreChain(clusList[i]); // Placing  the cluster back properly
                }
            bAccept = 0;
            free(oldHashTab);
            free(clusList);
            return bAccept;
            // printf("End CLUS - Failed.\n");
        }
}

/// Move_DbPvt - performs a Double-Pivot move on beadID by:
/// 1. Calculate the number of possible bridges in a <+-2,+-2,+-2> around
/// beadID. The move is rejected if the molecule is not a chain.
/// 2. Pick one of these bridges, and remember which beadID'+1 will be beadID+1
/// after the move.
/// 3. Calculate the number of bridges possible around beadID'+1.
/// 4. Calculate the ratio between the number of possible bridges, and perform
/// the Metropolis-Hastings step. Note that bridging only changes the local
/// connectivity of the chains, and thus not the energy.
/// \param beadID
/// \return 1 if accepted, 0 if rejected.
int Move_DbPvt(const int beadID, const float myTemp)
{ // Performs a double-pivot move.
    /* Molecule MUST be LINEAR
  The move requires selecting a random bead, which is beadID. Then, we'll search
  the lattice in +-min(l,3) sites around beadID. Let i be the position of beadID along
  it's chain. Let i' denote same position along another chain of the same type.
  We want dist(i,i'+1) < linker_len_glb[i] && Dist_BeadToBead (i',i+1) <
  linker_len_glb[i']. We'll count however many candidates there are and select one
  randomly. Then we just swap beads from i to then end of the chain. Do a
  Metropolis thing, and decide. In other words, i'+1 becomes i+1, i'+2 becomes
  i+2 until N, and i+1 become i'+1 and so on.
   */
    int bAccept        = 0;                                    // Move acceptance and such LEL
    const int PChainID = bead_info_glb[beadID][BEAD_CHAINID];  // The proposed chainID
    const int PType    = chain_info_glb[PChainID][CHAIN_TYPE]; // Type of chain.

    if (nChainTypeIsLinear_glb[PType] == 0)
        {
            // Reject the move because the chain is not linear.
            return bAccept;
        }
    // Now make sure that the proposed bead is neither the start or end of a
    // chain.
    const int PStart = chain_info_glb[PChainID][CHAIN_START];           // Start of this chain.
    const int PEnd   = PStart + chain_info_glb[PChainID][CHAIN_LENGTH]; // LastBead+1 for this chain.
    if (PStart == beadID || PEnd - 1 == beadID)
        {
            return bAccept;
        }

    int i, j; // Loop iterators
    int candList[345] = {0};
    int filtList[345] = {0};

    const int thisBead         = beadID;
    const int thisBeadChainPos = thisBead - PStart; // 0-indexed chain position.
    const int nextBead         = beadID + 1;
    const int nextBeadChainPos = thisBeadChainPos + 1;

    // The length of the box is min(LARGEST_RADIUS, linker_len_glb[beadID][1])
    const int n_rad = LARGEST_RADIUS < linker_len_glb[thisBead][1] ?
                          LARGEST_RADIUS :
                          linker_len_glb[thisBead][1]; // The length of the box we'll be searching.

    const int r_pos0[POS_MAX] = {bead_info_glb[thisBead][0], bead_info_glb[thisBead][1], bead_info_glb[thisBead][2]};
    int neigh_num             = NeighborSearch_AroundPoint_wRad_IgnBead(thisBead, r_pos0, n_rad, candList);

    BeadListOP_GetChainIDs(neigh_num, candList, filtList);
    int other_chains_num = BeadListOP_InvFilter_wrt_SecondList(neigh_num, candList, filtList, PChainID);

    BeadListOP_GetChainTypes(other_chains_num, candList, filtList);
    int same_chain_types = BeadListOP_Filter_wrt_SecondList(other_chains_num, candList, filtList, PType);

    BeadListOP_GetIntraChainID(same_chain_types, candList, filtList);
    int corr_chain_pos = BeadListOP_Filter_wrt_SecondList(same_chain_types, candList, filtList, nextBeadChainPos);

    const int fwd_cand_num = BeadListOP_Filter_DbPvtLinkerConFwd(corr_chain_pos, candList, nextBead);

    if (fwd_cand_num == 0) // Rejecting if no forward moves
        {
            bAccept = 0;
            return bAccept;
        }

    // For detailed balance, we need the reverse move. So we look at i', and find suitable i+1 candidates.

    const int nextBead_P = candList[rand() % fwd_cand_num];
    const int thisBead_P = nextBead_P - 1;
    const int PChainID_P = bead_info_glb[thisBead_P][BEAD_CHAINID];

    const int r_pos0_P[POS_MAX] = {bead_info_glb[thisBead_P][0], bead_info_glb[thisBead_P][1],
                                   bead_info_glb[thisBead_P][2]};

    neigh_num = NeighborSearch_AroundPoint_wRad_IgnBead(thisBead_P, r_pos0_P, n_rad, candList);

    BeadListOP_GetChainIDs(neigh_num, candList, filtList);
    other_chains_num = BeadListOP_InvFilter_wrt_SecondList(neigh_num, candList, filtList, PChainID_P);

    BeadListOP_GetChainTypes(other_chains_num, candList, filtList);
    same_chain_types = BeadListOP_Filter_wrt_SecondList(other_chains_num, candList, filtList, PType);

    BeadListOP_GetIntraChainID(same_chain_types, candList, filtList);
    corr_chain_pos = BeadListOP_Filter_wrt_SecondList(same_chain_types, candList, filtList, nextBeadChainPos);

    const int bck_cand_num = BeadListOP_Filter_DbPvtLinkerConFwd(corr_chain_pos, candList, nextBead_P);

    if (bck_cand_num == 0)
        {          // We swap the chains.
            j = 0; // Tracks the beads
            for (i = nextBead; i < PEnd; i++)
                { // Swapping from beadID+1 onwards
                    OP_SwapBeads(nextBead + j,
                                 nextBead_P + j); // This is pretty dumb, but easier to read UGH
                    j++;
                }
            bAccept = 1;
            return bAccept;
        }
    lLDub oldEn, newEn;

    const int resi = bead_info_glb[thisBead][BEAD_TYPE];
    oldEn          = 0.;
    if (faEnergy_glb[resi][resi][E_STIFF])
        {
            oldEn = Energy_Topo_Angle(thisBead) + Energy_Topo_Angle(nextBead);
            oldEn += Energy_Topo_Angle(thisBead_P) + Energy_Topo_Angle(nextBead_P);
        }

    // We swap the chains.
    j = 0; // Tracks the beads
    for (i = nextBead; i < PEnd; i++)
        { // Swapping from beadID+1 onwards
            OP_SwapBeads(nextBead + j,
                         nextBead_P + j); // This is pretty dumb, but easier to read UGH
            j++;
        }

    newEn = 0.;
    if (faEnergy_glb[resi][resi][E_STIFF])
        {
            newEn = Energy_Topo_Angle(thisBead) + Energy_Topo_Angle(nextBead);
            newEn += Energy_Topo_Angle(thisBead_P) + Energy_Topo_Angle(nextBead_P);
        }

    const lLDub MCProb  = (lLDub) rand() / (lLDub) RAND_MAX;
    const lLDub MHAcc_W = (lLDub) (fwd_cand_num) / (lLDub) (bck_cand_num);
    if (MCProb < MHAcc_W * expl((oldEn - newEn) / myTemp))
        {
            bAccept = 1;
            return bAccept;
        }
    else
        { // Rejecting, so we swap back.
            // We swap the chains.
            j = 0; // Tracks the beads
            for (i = nextBead; i < PEnd; i++)
                { // Swapping from beadID+1 onwards
                    OP_SwapBeads(nextBead + j,
                                 nextBead_P + j); // This is pretty dumb, but easier to read UGH
                    j++;
                }
            bAccept = 0;
            return bAccept;
        } //
}

/// Move_CoLocal - move thisBeadID and it's physical bond partner to a new
/// location. Move fails if thisBeadID has no partner, and if no space is found.
/// A standard unbiased translation of two beads where a Metropolis-Hastings
/// step is performed at the end.
/// \param thisBeadID
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_CoLocal(int thisBeadID, float MyTemp)
{
    /*
  Translate a bead and its partner in tandem. If no partner, reject move.
  */

    int bAccept = 0; // Used in MC steps

    if (bead_info_glb[thisBeadID][BEAD_FACE] == -1)
        {
            bAccept = 0;
            return bAccept;
        }
    int r_disp[POS_MAX]; // Random translation vector.
    int r_posNew1[POS_MAX], r_posNew2[POS_MAX];
    int r_pos1[POS_MAX], r_pos2[POS_MAX];
    const int otherBeadID = bead_info_glb[thisBeadID][BEAD_FACE];
    int yTemp;

    LatPos_copy(r_pos1, bead_info_glb[thisBeadID]);
    LatPos_copy(r_pos2, bead_info_glb[otherBeadID]);

    LatPos_gen_rand_wRad(r_disp, 2);
    LatPos_add_wPBC(r_posNew1, r_disp, r_pos1);
    LatPos_add_wPBC(r_posNew2, r_disp, r_pos2);

    yTemp = Check_MoveBeadTo(r_posNew1) * Check_MoveBeadTo(r_posNew2);
    if (yTemp == 1)
        { // This means we found an empty lattice site. So let's
          // check if the linkers are okay.
            yTemp = Check_LinkerConstraint(thisBeadID, r_posNew1) * Check_LinkerConstraint(otherBeadID, r_posNew2);
        }
    if (yTemp == 0)
        { // Steric-clash or bad location for the beads.
            bAccept = 0;
            return bAccept;
        }

    lLDub newEn = 0.;
    lLDub oldEn = 0.;
    if (nBiasPotential_Mode_glb != -1)
        {
            oldEn += Energy_BiasingPotential(thisBeadID);
            oldEn += Energy_BiasingPotential(otherBeadID);
        }

    int ovlp_num, cont_num;
    Energy_Iso_ForCoLocal(thisBeadID, otherBeadID, r_pos1, &oldEn, &newEn, &ovlp_num, &cont_num, naOldOvlpNeighs_glb,
                          naOldContNeighs_glb);
    Energy_Iso_ForCoLocal(otherBeadID, thisBeadID, r_pos2, &oldEn, &newEn, &ovlp_num, &cont_num, naOldOvlpNeighs_glb,
                          naOldContNeighs_glb);

    int bondList[2 * (MAX_BONDS + 1)];
    int bond_tmpList[MAX_BONDS + 1];
    int bondNum, bond_tmpNum;

    if (bSystemHasTopo_glb)
        {
            bond_tmpNum = OP_GetTopoBonds(thisBeadID, bond_tmpList);
            bondNum     = BeadList_AppendBeads(0, bondList, bond_tmpList, bond_tmpNum);
            bond_tmpNum = OP_GetTopoBonds(otherBeadID, bond_tmpList);
            bondNum     = BeadList_AppendBeads(bondNum, bondList, bond_tmpList, bond_tmpNum);
            qsort(bondList, bondNum, sizeof(int), UtilFunc_CompareInts);
            bondNum = ListOP_UniqueElementsOfSortedList_Int(bondNum, bondList);
            bondNum = BeadList_CanTopoAngle(bondNum, bondList);
        }
    else
        {
            bondNum = 0;
        }

    oldEn += bondNum ? Energy_Topo_Angle_ForList(bondNum, bondList) : 0.;

    OP_System_MoveBeadTo(thisBeadID, r_posNew1);
    OP_System_MoveBeadTo(otherBeadID, r_posNew2);

    if (nBiasPotential_Mode_glb != -1)
        {
            newEn += Energy_BiasingPotential(thisBeadID);
            newEn += Energy_BiasingPotential(otherBeadID);
        }

    Energy_Iso_ForCoLocal(thisBeadID, otherBeadID, r_posNew1, &newEn, &oldEn, &ovlp_num, &cont_num, naOldOvlpNeighs_glb,
                          naOldContNeighs_glb);
    Energy_Iso_ForCoLocal(otherBeadID, thisBeadID, r_posNew2, &newEn, &oldEn, &ovlp_num, &cont_num, naOldOvlpNeighs_glb,
                          naOldContNeighs_glb);

    newEn += bondNum ? Energy_Topo_Angle_ForList(bondNum, bondList) : 0.;

    lLDub MCProb = (lLDub) rand() / (lLDub) RAND_MAX;
    lLDub MHAcc  = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub) MyTemp);
    if (MCProb < MHAcc)
        { // Accept this state
            bAccept = 1;
            // printf("Accepted!\n");
            return bAccept;
        }
    else
        {
            OP_System_MoveBeadTo_Inv(thisBeadID);
            OP_System_MoveBeadTo_Inv(otherBeadID);
            bAccept = 0;
            // printf("Rejected move; undoing!\n");
            return bAccept;
        }
}

/// Move_MultiLocal - performs a biased set of local moves where beadID and
/// every covalently bonded bead of beadID is moved by:
/// 1. Finding an empty spot for all the beads involved (in a <+-2,+-2,+-2>).
/// This is a little crude because I put the beads on the lattice into the new
/// proposed spot because the structure constraint sub-routine is designed as
/// such.
/// 2. If there is space for all of the beads. Calculate the energies and the
/// total Rosenbluth weights.
/// 3. Move the beads to new location and recalculate the energies and weights.
/// Perform a Metropolis-Hastings step.
/// \param beadID
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_MultiLocal(int beadID, float MyTemp)
{

    int beadsList[MAX_BONDS + 1];
    const int beadNum = OP_GetTopoBonds(beadID, beadsList);
    int bAccept;

    if (beadNum < 2)
        { // Bead must be bonded to at least 2 things.
            bAccept = 0;
            return bAccept;
        }

    int beads_pos0[MAX_BONDS + 1][POS_MAX];
    int beads_posNew[MAX_BONDS + 1][POS_MAX];
    int r_posTmp1[POS_MAX];
    int tmpBead;
    int i;
    int yTemp;

    OP_Beads_CopyBeadsInListToOld(beadNum, beadsList);
    OP_Beads_CopyBeadsInListToPosList(beadNum, beadsList, beads_pos0);

    OP_Lattice_EmptySitesForListOfPos(beadNum, beads_pos0);

    for (i = 0; i < beadNum; i++)
        {
            LatPos_gen_rand_wRad(r_posTmp1, 2);
            LatPos_add_wPBC(beads_posNew[i], r_posTmp1, beads_pos0[i]);
            yTemp = Check_MoveBeadTo(beads_posNew[i]);
            if (yTemp == 0)
                {
                    break;
                }
            tmpBead                                            = beadsList[i];
            naTotLattice_glb[Lat_Ind_FromVec(beads_posNew[i])] = tmpBead;
        }

    OP_Lattice_EmptySitesForListOfPos(i, beads_posNew);

    if (yTemp == 1)
        { // No steric clash so check for topology constraint.
            OP_Beads_MoveBeadsInListToPos(beadNum, beadsList, beads_posNew);
            yTemp = Check_LinkerConstraints_ForBeadList(beadNum, beadsList);
            OP_Beads_MoveBeadsInListToPos(beadNum, beadsList, beads_pos0);
        }

    OP_Lattice_PlaceBeadsInList(beadNum, beadsList);

    if (yTemp == 0)
        {
            bAccept = 0;
            return bAccept;
        }

    lLDub oldEn = 0.;

    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = 0; i < beadNum; i++)
                {
                    tmpBead = beadsList[i];
                    oldEn += Energy_BiasingPotential(tmpBead);
                }
        }

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;

    lLDub BSum  = 0.;
    lLDub newEn = 0.;
    for (i = 0; i < beadNum; i++)
        {
            Energy_Iso_ForLists(i, beadNum, beadsList, beads_pos0, &oldEn, &newEn, &old_ovlp_num, &old_cont_num,
                                naOldOvlpNeighs_glb, naOldContNeighs_glb);
            BSum += MC_RosenbluthSampling_ForLists_AtOld(i, beadNum, beadsList, &oldEn, old_ovlp_num);
        }

    int bondAngNum = 0;
    int tmpAngList[MAX_BONDS + 1];
    int bondAngList[(MAX_BONDS + 1) * (MAX_BONDS + 1)];
    if (bSystemHasTopo_glb)
        {
            bondAngNum = BeadList_AppendBeads(0, bondAngList, beadsList, beadNum);
            for (i = 1; i < beadNum; i++)
                {
                    tmpBead    = beadsList[i];
                    yTemp      = OP_GetTopoBonds(tmpBead, tmpAngList);
                    bondAngNum = BeadList_AppendBeads(bondAngNum, bondAngList, tmpAngList, yTemp);
                }
            qsort(bondAngList, bondAngNum, sizeof(int), UtilFunc_CompareInts);
            bondAngNum = ListOP_UniqueElementsOfSortedList_Int(bondAngNum, bondAngList);
        }
    else
        {
            bondAngNum = 0;
        }

    oldEn += bondAngNum ? Energy_Topo_Angle_ForList(bondAngNum, bondAngList) : 0.;

    OP_System_MoveBeadsInListToPos(beadNum, beadsList, beads_posNew);
    OP_Beads_BreakBondsInList(beadNum, beadsList);

    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = 0; i < beadNum; i++)
                {
                    tmpBead = beadsList[i];
                    newEn += Energy_BiasingPotential(tmpBead);
                }
        }

    lLDub FSum = 0.;
    for (i = 0; i < beadNum; i++)
        {
            Energy_Iso_ForLists(i, beadNum, beadsList, beads_posNew, &newEn, &oldEn, &new_ovlp_num, &new_cont_num,
                                naNewOvlpNeighs_glb, naNewContNeighs_glb);
            FSum += MC_RosenbluthSampling_ForLists_AtNew(i, beadNum, beadsList, &yTemp, &newEn, new_ovlp_num);
            if (yTemp != -1)
                { // An appropriate partner has been selected. Form
                  // the bonds and add the energy
                    tmpBead = beadsList[i];
#if DEBUG_BUILD
                    if (bead_info_glb[tmpBead][BEAD_FACE] != -1)
                        {
                            printf("HOW?!\n");
                            exit(1);
                        }
#endif
                    bead_info_glb[tmpBead][BEAD_FACE] = yTemp;
                    bead_info_glb[yTemp][BEAD_FACE]   = tmpBead;
                    newEn += Energy_Anisotropic_For_List(tmpBead, beadNum, beadsList);
                }
        }

    newEn += bondAngNum ? Energy_Topo_Angle_ForList(bondAngNum, bondAngList) : 0.;

    const lLDub MCProb = (lLDub) rand() / (lLDub) RAND_MAX;
    const lLDub MHAcc  = OP_GenMHValue(FSum, BSum, oldEn - newEn, (lLDub) MyTemp);

    if (MCProb < MHAcc)
        { // Accept the move. Remember that the bonds were
          // assigned above
            bAccept = 1;
            return bAccept;
        }
    else
        {
            OP_Beads_BreakBondsInList(beadNum, beadsList);
            OP_Lattice_EmptySitesForListOfBeads(beadNum, beadsList);
            OP_Beads_CopyBeadsInListFromOld(beadNum, beadsList);
            OP_Beads_RestoreBondsInList(beadNum, beadsList);
            OP_Lattice_PlaceBeadsInList(beadNum, beadsList);
            bAccept = 0;
            return bAccept;
        }
}

/// Move_Pivot - performs a biased pivot move on chainID.
/// Move fails if chainID is branched. The process is:
/// 1. Pick a random bead within chainID that shall act as the pivot around
/// which a rotation operation shall be performed. The operation is performed on
/// the shorter part of the chain. If the selected bead is an end, the move is
/// rejected.
/// 2. Pick a rotation operation. The list of operation is a rotation by 90,
/// 180, 270 or 360 (null) degrees around the X, Y, or Z axis.
/// 3. Calculate the Rosenbluth weights and energies, and check if there is
/// steric clash after the proposed move. If there is a clash, reject the move.
/// 4. Rotate the beads, and recalculate the weights, propose new bonds and
/// recalculate the energies.
/// 5. Perform the Metropolis-Hastings step.
/// \param chainID
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_Pivot(int chainID, float MyTemp)
{
    // Performs a pivot move on chainID
    /* Randomly pick a bead in this chain (anchorBead), and perform a symmetry
     * operation on the beads after the anchorBead. Note that if anchorBead is
     * the second-to-last, or last, bead, the move is rejected. Local and
     * slither moves would be better for those. The protein must be linear and
     * the move is rejected if the chain is branched. It is assumed that a
     * linear chain has been passed!
     */
    int bAccept = 0;
    if (! nChainTypeIsLinear_glb[chain_info_glb[chainID][CHAIN_TYPE]])
        {
            return bAccept;
        }
    // Check if the chain is longer than 3 or if it is linear
    const int chainLength = chain_info_glb[chainID][CHAIN_LENGTH];
    if (chainLength <= 3)
        {
            return bAccept;
        }

    int PivotM;
    PivotM = rand() % 10;
    if (PivotM == 0)
        { // Null move
            bAccept = 1;
            return bAccept;
        }

    // Get the beadID's for the first and last bead
    const int firstB = chain_info_glb[chainID][CHAIN_START];
    const int lastB  = firstB + chainLength;

    // Randomly select a bead that is neither the first nor last
    int anchorBead = chainLength - 2;
    anchorBead     = 1 + (rand() % anchorBead);

    //-1 Means backwards, +1 means forwards. Always Pivot the smaller portion
    const int PivotDir = anchorBead > chainLength / 2 ? 1 : -1;
    anchorBead         = firstB + anchorBead;
    // printf("Bead: %d\n", anchorBead);

    int i, j;
    int xTemp, yTemp;
    int anchorPos[POS_MAX];
    int beadsList[MAX_CHAINLEN];
    int beadNum = 0;

    if (PivotDir == 1)
        {
            for (i = anchorBead + 1; i < lastB; i++)
                {
                    beadsList[beadNum++] = i;
                }
        }
    else
        {
            for (i = firstB; i < anchorBead; i++)
                {
                    beadsList[beadNum++] = i;
                }
        }

    const int smBead = beadsList[0];
    const int lgBead = beadsList[beadNum - 1];

    LatPos_copy(anchorPos, bead_info_glb[anchorBead]);

    int tmpBead;

    for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            // Use anchorbead as origin, and check what happens after the rotation
            // to every bead
            OP_Rotation(PivotM, tmpBead, anchorPos);
            yTemp = Check_MoveBeadTo(naTempR_glb);
            if (yTemp == 0)
                {
                    break;
                }
        }
    if (yTemp == 0)
        {
            bAccept = 0;
            return bAccept;
        }

    lLDub oldEn = 0.;
    lLDub newEn = 0.;
    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = 0; i < beadNum; i++)
                {
                    tmpBead = beadsList[i];
                    oldEn += Energy_BiasingPotential(tmpBead);
                }
        }

    int old_ovlp_num, old_cont_num;
    int resi;

    lLDub BSum = 0.;
    for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            Energy_Iso_ForRange(tmpBead, smBead, lgBead, &oldEn, &newEn, &old_ovlp_num, &old_cont_num,
                                naOldOvlpNeighs_glb, naOldContNeighs_glb);
            resi = bead_info_glb[tmpBead][BEAD_TYPE];
            BSum += MC_RosenbluthSampling_ForRange_AtOld(tmpBead, resi, smBead, lgBead, &oldEn, old_ovlp_num);
        }

    resi = bead_info_glb[anchorBead][BEAD_TYPE];
    oldEn += faEnergy_glb[resi][resi][E_STIFF] ? Energy_Topo_Angle(anchorBead) : 0.;

    for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            OP_Rotation(PivotM, tmpBead, anchorPos);
            OP_System_MoveBeadTo(tmpBead, naTempR_glb);
        }

    OP_Beads_BreakBondsInList(beadNum, beadsList);

    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = 0; i < beadNum; i++)
                {
                    tmpBead = beadsList[i];
                    newEn += Energy_BiasingPotential(tmpBead);
                }
        }

    int new_ovlp_num, new_cont_num;
    int resj;

    lLDub FSum = 0.;
    for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            resi    = bead_info_glb[tmpBead][BEAD_TYPE];
            Energy_Iso_ForRange(tmpBead, smBead, lgBead, &newEn, &oldEn, &new_ovlp_num, &new_cont_num,
                                naNewOvlpNeighs_glb, naNewContNeighs_glb);
            FSum += MC_RosenbluthSampling_ForRange_AtNew(tmpBead, resi, &yTemp, &newEn, new_ovlp_num);
            if (yTemp != -1)
                {
                    resj                              = bead_info_glb[yTemp][BEAD_TYPE];
                    bead_info_glb[yTemp][BEAD_FACE]   = tmpBead;
                    bead_info_glb[tmpBead][BEAD_FACE] = yTemp;
                    newEn += Energy_Anisotropic_For_Range(tmpBead, smBead, lgBead);
                }
        }

    resi = bead_info_glb[anchorBead][BEAD_TYPE];
    newEn += faEnergy_glb[resi][resi][E_STIFF] ? Energy_Topo_Angle(anchorBead) : 0.;

    const lLDub MCProb = (lLDub) rand() / (lLDub) RAND_MAX;
    lLDub MHAcc        = OP_GenMHValue(FSum, BSum, oldEn - newEn, (lLDub) MyTemp);
    if (MCProb < MHAcc)
        { // Accept the move. Remember that the bonds were
          // assigned above!
            bAccept = 1;
            return bAccept;
        }
    else
        { // Rejecting move
            OP_Beads_BreakBondsInList(beadNum, beadsList);
            for (i = 0; i < beadNum; i++)
                {
                    tmpBead = beadsList[i];
                    OP_System_MoveBeadTo_Inv(tmpBead);
                }
            bAccept = 0;
            return bAccept;
        }
}

/// Move_BranchedRot - performs the branched analog for the pivot move where the
/// central bead of the branched molecule is considered the anchor, and the
/// whole molecule is rotated with the biased sampling, and re-assigning of
/// bonds. The move is rejected if chainID is linear.
/// \param chainID
/// \param MyTemp
/// \return 1 if accepted, 0 if rejected.
int Move_BranchedRot(int chainID, float MyTemp)
{
    // Rotates a branched molecule about the branching, which is assumed to be
    // the firstB of chainID Performs a Move_Pivot() on molecule where the
    // rotation occurs around firstB
    /*
      Set the first bead as anchorBead, and perform a symmetry operation on the
      beads after the anchorBead (the whole molecule). Note that if the molecule
      is linear, the move is outright rejected. Again, it is assumed that the
      first bead in that molecule is the 'node'.
      */
    int bAccept = 0;
    // Reject if the molecule is linear
    if (nChainTypeIsLinear_glb[chain_info_glb[chainID][CHAIN_TYPE]])
        {
            return bAccept;
        }
    const int chainLength = chain_info_glb[chainID][CHAIN_LENGTH];

    int PivotM;
    PivotM = rand() % 10;
    if (PivotM == 0)
        { // Null move
            bAccept = 1;
            return bAccept;
        }

    // Get the beadID's for the first and last bead
    const int firstB = chain_info_glb[chainID][CHAIN_START];
    const int lastB  = firstB + chainLength;

    // FirstB is the anchor.
    const int anchorBead = firstB;

    int i, j;
    int xTemp, yTemp;
    int anchorPos[POS_MAX];
    int beadsList[MAX_CHAINLEN];
    int beadNum = 0;

    for (i = anchorBead + 1; i < lastB; i++)
        {
            beadsList[beadNum++] = i;
        }

    const int smBead = beadsList[0];
    const int lgBead = beadsList[beadNum - 1];

    LatPos_copy(anchorPos, bead_info_glb[anchorBead]);

    int tmpBead;

    for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            // Use anchorbead as origin, and check what happens after the rotation
            // to every bead
            OP_Rotation(PivotM, tmpBead, anchorPos);
            yTemp = Check_MoveBeadTo(naTempR_glb);
            if (yTemp == 0)
                {
                    break;
                }
        }
    if (yTemp == 0)
        {
            bAccept = 0;
            return bAccept;
        }

    lLDub oldEn = 0.;
    lLDub newEn = 0.;
    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = 0; i < beadNum; i++)
                {
                    tmpBead = beadsList[i];
                    oldEn += Energy_BiasingPotential(tmpBead);
                }
        }

    int old_ovlp_num, old_cont_num;
    int resi;

    lLDub BSum = 0.;
    for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            Energy_Iso_ForRange(tmpBead, smBead, lgBead, &oldEn, &newEn, &old_ovlp_num, &old_cont_num,
                                naOldOvlpNeighs_glb, naOldContNeighs_glb);
            resi = bead_info_glb[tmpBead][BEAD_TYPE];
            BSum += MC_RosenbluthSampling_ForRange_AtOld(tmpBead, resi, smBead, lgBead, &oldEn, old_ovlp_num);
        }

    for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            OP_Rotation(PivotM, tmpBead, anchorPos);
            OP_System_MoveBeadTo(tmpBead, naTempR_glb);
        }

    OP_Beads_BreakBondsInList(beadNum, beadsList);

    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = 0; i < beadNum; i++)
                {
                    tmpBead = beadsList[i];
                    newEn += Energy_BiasingPotential(tmpBead);
                }
        }

    int new_ovlp_num, new_cont_num;
    int resj;

    lLDub FSum = 0.;
    for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            resi    = bead_info_glb[tmpBead][BEAD_TYPE];
            Energy_Iso_ForRange(tmpBead, smBead, lgBead, &newEn, &oldEn, &new_ovlp_num, &new_cont_num,
                                naNewOvlpNeighs_glb, naNewContNeighs_glb);
            FSum += MC_RosenbluthSampling_ForRange_AtNew(tmpBead, resi, &yTemp, &newEn, new_ovlp_num);
            if (yTemp != -1)
                {
                    resj                              = bead_info_glb[yTemp][BEAD_TYPE];
                    bead_info_glb[yTemp][BEAD_FACE]   = tmpBead;
                    bead_info_glb[tmpBead][BEAD_FACE] = yTemp;
                    newEn += Energy_Anisotropic_For_Range(tmpBead, smBead, lgBead);
                }
        }

    lLDub MCProb = (lLDub) rand() / (lLDub) RAND_MAX;
    lLDub MHAcc  = OP_GenMHValue(FSum, BSum, oldEn - newEn, (lLDub) MyTemp);
    if (MCProb < MHAcc)
        { // Accept the move. Remember that the bonds were
          // assigned above!
            bAccept = 1;
            return bAccept;
        }
    else
        { // Rejecting move
            OP_Beads_BreakBondsInList(beadNum, beadsList);
            for (i = 0; i < beadNum; i++)
                {
                    tmpBead = beadsList[i];
                    OP_System_MoveBeadTo_Inv(tmpBead);
                }
            bAccept = 0;
            return bAccept;
        }
}

/// Move_SmallClus_Proximity
/// \param chainID
/// \param myTemp
/// \return
int Move_SmallClus_Proximity(const int chainID, const float myTemp)
{
    // Performs a cluster move where a given chain and it's cluster are moved.

    int bAccept = 0; // Used in MC steps, assume that move fails initially.

    int* oldClusList = calloc(tot_chains_glb, sizeof(int));
    char* oldHashTab = calloc(tot_chains_glb, sizeof(char));
    const int ClusSize =
        ClusUtil_OvlpCluster_OfChain_wMaxSize(chainID, oldHashTab, oldClusList, nLimitedClusterSize_glb);

    if (ClusSize < 2)
        {
            bAccept = 0;
            free(oldClusList);
            free(oldHashTab);
            return bAccept;
        }
    //    printf("\n%d\n", ClusSize);

    // Radii for translation moves. All moves are L/2 radius
    int r_Disp[POS_MAX];
    int yTemp;
    int i, j;
    LatPos_gen_rand_wRad(r_Disp, naBoxSize_glb[0] / 2);

    for (i = 0; i < ClusSize; i++)
        {
            yTemp = Check_ChainDisp(oldClusList[i], r_Disp); // Checking for steric clash
            if (yTemp == 0)
                {
                    bAccept = 0;
                    free(oldClusList);
                    free(oldHashTab);
                    return bAccept;
                }
        }
    // No  clash

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;
    int thisChain, firstB, lastB;

    lLDub oldEn = 0.;
    lLDub newEn = 0.;
    if (nBiasPotential_Mode_glb != -1)
        {
            for (j = 0; j < ClusSize; j++)
                {
                    thisChain = oldClusList[j];
                    firstB    = chain_info_glb[thisChain][CHAIN_START];
                    lastB     = firstB + chain_info_glb[thisChain][CHAIN_LENGTH];
                    for (i = firstB; i < lastB; i++)
                        {
                            oldEn += Energy_BiasingPotential(i);
                        }
                }
        }

    // If no CONT interactions, then just checking the cluster is enough. If we have a different cluster, we immediately
    // reject. Therefore, the clusters are moved as rigid bodies and will have the same OVLP energies before and after.
    // If we have CONT interactions, we could have different energies. Therefore, we should calculate the energies.
    if (bSystemHasCont_glb)
        {
            for (j = 0; j < ClusSize; j++)
                {
                    thisChain = oldClusList[j];
                    firstB    = chain_info_glb[thisChain][CHAIN_START];
                    lastB     = firstB + chain_info_glb[thisChain][CHAIN_LENGTH];
                    for (i = firstB; i < lastB; i++)
                        {
                            Energy_Iso_ForChains(i, &oldEn, &newEn, &old_ovlp_num, &old_cont_num, naOldOvlpNeighs_glb,
                                                 naOldContNeighs_glb);
                        }
                }
        }

    for (i = 0; i < ClusSize; i++)
        {
            OP_System_DispChain(oldClusList[i], r_Disp); // Moving the cluster properly
        }
    // Recalculating cluster to see if we have the same cluster or not. If not, we reject.

    const int ClusCheckN = ClusUtil_OvlpCluster_OfChain_CheckForSame(oldHashTab, oldClusList, ClusSize);
    free(oldHashTab);

    if (ClusCheckN == -1)
        {
            for (i = 0; i < ClusSize; i++)
                {
                    OP_System_RestoreChain(oldClusList[i]); // Placing  the cluster back properly
                }
            bAccept = 0;
            free(oldClusList);
            return bAccept;
        }

    if (nBiasPotential_Mode_glb != -1)
        {
            for (j = 0; j < ClusSize; j++)
                {
                    thisChain = oldClusList[j];
                    firstB    = chain_info_glb[thisChain][CHAIN_START];
                    lastB     = firstB + chain_info_glb[thisChain][CHAIN_LENGTH];
                    for (i = firstB; i < lastB; i++)
                        {
                            newEn += Energy_BiasingPotential(i);
                        }
                }
        }

    if (bSystemHasCont_glb)
        {
            for (j = 0; j < ClusSize; j++)
                {
                    thisChain = oldClusList[j];
                    firstB    = chain_info_glb[thisChain][CHAIN_START];
                    lastB     = firstB + chain_info_glb[thisChain][CHAIN_LENGTH];
                    for (i = firstB; i < lastB; i++)
                        {
                            Energy_Iso_ForChains(i, &newEn, &oldEn, &new_ovlp_num, &new_cont_num, naNewOvlpNeighs_glb,
                                                 naNewContNeighs_glb);
                        }
                }
        }

    const lLDub MCProb = (lLDub) rand() / (lLDub) RAND_MAX;
    const lLDub MHAcc  = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub) myTemp);
    if (MCProb < MHAcc)
        {                // Accept this state
            bAccept = 1; // Accept the move
            free(oldClusList);
            return bAccept;
            // printf("End CLUS - Yes\n");
        }
    else
        {
            for (i = 0; i < ClusSize; i++)
                {
                    OP_System_RestoreChain(oldClusList[i]); // Placing  the cluster back properly
                }
            bAccept = 0;
            free(oldClusList);
            return bAccept;
            // printf("End CLUS - Failed.\n");
        }
}

/// Move_Clus_Proximity - Attempts to move the second largest cluster in the system
///  \param myTemp
///  \return
int Move_Clus_Proximity(const float myTemp)
{
    // Performs a cluster move where a given chain and it's cluster are moved.
    // No new 'bonds' are made so the move is reversible.... So we have to check if we get the same cluster back

    int bAccept = 0; // Used in MC steps, assume that move fails initially.

    int* naClusList    = (int*) malloc(sizeof(int) * (tot_chains_glb + 1));
    const int ClusSize = ClusUtil_OvlpCluster_OfSystem_SecondLargest_ForMCMove(naClusList);
    if (ClusSize < 2)
        {
            bAccept = 0;
            free(naClusList);
            return bAccept;
        }

    // Radii for translation moves. All moves are L/2 radius
    int r_Disp[POS_MAX];
    int yTemp;
    int i, j;
    LatPos_gen_rand_wRad(r_Disp, naBoxSize_glb[0] / 2);

    for (i = 0; i < ClusSize; i++)
        {
            yTemp = Check_ChainDisp(naClusList[i], r_Disp); // Checking for steric clash
            if (yTemp == 0)
                {
                    bAccept = 0;
                    free(naClusList);
                    return bAccept;
                }
        }
    // No  clash

    char* oldHashTab = calloc(tot_chains_glb, sizeof(char));

    for (i = 0; i < ClusSize; i++)
        {
            oldHashTab[naClusList[i]] = 1;
        }

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;
    int thisChain, firstB, lastB;

    lLDub oldEn = 0.;
    lLDub newEn = 0.;
    if (nBiasPotential_Mode_glb != -1)
        {
            for (j = 0; j < ClusSize; j++)
                {
                    thisChain = naClusList[j];
                    firstB    = chain_info_glb[thisChain][CHAIN_START];
                    lastB     = firstB + chain_info_glb[thisChain][CHAIN_LENGTH];
                    for (i = firstB; i < lastB; i++)
                        {
                            oldEn += Energy_BiasingPotential(i);
                        }
                }
        }

    // If no CONT interactions, then just checking the cluster is enough. If we have a different cluster, we immediately
    // reject. Therefore, the clusters are moved as rigid bodies and will have the same OVLP energies before and after.
    // If we have CONT interactions, we could have different energies. Therefore, we should calculate the energies.
    if (bSystemHasCont_glb)
        {
            for (j = 0; j < ClusSize; j++)
                {
                    thisChain = naClusList[j];
                    firstB    = chain_info_glb[thisChain][CHAIN_START];
                    lastB     = firstB + chain_info_glb[thisChain][CHAIN_LENGTH];
                    for (i = firstB; i < lastB; i++)
                        {
                            Energy_Iso_ForChains(i, &oldEn, &newEn, &old_ovlp_num, &old_cont_num, naOldOvlpNeighs_glb,
                                                 naOldContNeighs_glb);
                        }
                }
        }

    for (i = 0; i < ClusSize; i++)
        {
            OP_System_DispChain(naClusList[i], r_Disp); // Moving the cluster properly
        }
    // Recalculating cluster to see if we have the same cluster or not. If not, we reject.
    const int ClusCheckN = ClusUtil_OvlpCluster_OfChain_CheckForSame(oldHashTab, naClusList, ClusSize);
    free(oldHashTab);

    if (ClusCheckN == -1)
        {
            for (i = 0; i < ClusSize; i++)
                {
                    OP_System_RestoreChain(naClusList[i]); // Placing  the cluster back properly
                }
            bAccept = 0;
            free(naClusList);
            return bAccept;
        }

    if (nBiasPotential_Mode_glb != -1)
        {
            for (j = 0; j < ClusSize; j++)
                {
                    thisChain = naClusList[j];
                    firstB    = chain_info_glb[thisChain][CHAIN_START];
                    lastB     = firstB + chain_info_glb[thisChain][CHAIN_LENGTH];
                    for (i = firstB; i < lastB; i++)
                        {
                            newEn += Energy_BiasingPotential(i);
                        }
                }
        }

    if (bSystemHasCont_glb)
        {
            for (j = 0; j < ClusSize; j++)
                {
                    thisChain = naClusList[j];
                    firstB    = chain_info_glb[thisChain][CHAIN_START];
                    lastB     = firstB + chain_info_glb[thisChain][CHAIN_LENGTH];
                    for (i = firstB; i < lastB; i++)
                        {
                            Energy_Iso_ForChains(i, &newEn, &oldEn, &new_ovlp_num, &new_cont_num, naNewOvlpNeighs_glb,
                                                 naNewContNeighs_glb);
                        }
                }
        }

    const lLDub MCProb = (lLDub) rand() / (lLDub) RAND_MAX;
    const lLDub MHAcc  = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub) myTemp);
    if (MCProb < MHAcc)
        {                // Accept this state
            bAccept = 1; // Accept the move
            free(naClusList);
            return bAccept;
            // printf("End CLUS - Yes\n");
        }
    else
        {
            for (i = 0; i < ClusSize; i++)
                {
                    OP_System_RestoreChain(naClusList[i]); // Placing  the cluster back properly
                }
            bAccept = 0;
            free(naClusList);
            return bAccept;
            // printf("End CLUS - Failed.\n");
        }
}

// All the _Equil variants of the moves are spatially the same as their non
// _Equil variants. The only difference is that the anisotropic interaction is
// ignored, and thus no bias is applied.

int Move_Local_Equil(int beadID, float MyTemp)
{ // Performs a local translation MC-move on beadID

    int bAccept = 0; // Used in MC steps
    int yTemp;
    int r_pos0[POS_MAX], r_posNew[POS_MAX],
        r_disp[POS_MAX]; // Vectors to stores coordinates.
    // printf("Beginning LOCAL\n");
    LatPos_copy(r_pos0, bead_info_glb[beadID]);
    //    LatPos_gen_rand_wRad(r_disp, linker_len_glb[beadID][0]);
    LatPos_gen_rand_wRad(r_disp, 2);
    LatPos_add_wPBC(r_posNew, r_pos0, r_disp);

    yTemp = Check_MoveBeadTo(r_posNew);

    if (yTemp == 1)
        { // This means we found an empty lattice site. So let's
          // check if the linkers are okay.
            yTemp = Check_LinkerConstraint(beadID, r_posNew);
        }

    if (yTemp != 1)
        {
            // This means that we have failed to find an appropriate spot for this
            // bead to be moved to.
            //  Therefore, the move is rejected!
            bAccept = 0;
            // printf("End LOCAL - No space\n");
            return bAccept;
        }
    // Have successfully found a good lattice spot.
    lLDub MCProb, oldEn, newEn; // For Metropolis Hastings
    oldEn = 0.;
    newEn = 0.;

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;

    const int resi = bead_info_glb[beadID][BEAD_TYPE];

    oldEn = nBiasPotential_Mode_glb == -1 ? 0.f : Energy_BiasingPotential(beadID);
    newEn = 0.;

    Energy_Iso_ForLocalEquil(beadID, resi, r_pos0, &oldEn, &newEn, &old_ovlp_num, &old_cont_num, naOldOvlpNeighs_glb,
                             naOldContNeighs_glb);

    int bondList[MAX_BONDS + 1];
    int bondNum = bSystemHasTopo_glb ? OP_GetTopoBonds(beadID, bondList) : 0;
    bondNum     = BeadList_CanTopoAngle(bondNum, bondList);

    oldEn += bondNum ? Energy_Topo_Angle_ForList(bondNum, bondList) : 0.;

    OP_System_MoveBeadTo(beadID, r_posNew);

    newEn += nBiasPotential_Mode_glb == -1 ? 0.f : Energy_BiasingPotential(beadID);

    Energy_Iso_ForLocalEquil(beadID, resi, r_posNew, &newEn, &oldEn, &new_ovlp_num, &new_cont_num, naNewOvlpNeighs_glb,
                             naNewContNeighs_glb);

    newEn += bondNum ? Energy_Topo_Angle_ForList(bondNum, bondList) : 0.;

    MCProb      = (lLDub) rand() / (lLDub) RAND_MAX;
    lLDub MHAcc = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub) MyTemp);
    if (MCProb < MHAcc)
        {                // Accept this state
            bAccept = 1; // Accepting!
            // printf("End LOCAL - Accept\n");
            return bAccept;
        }
    else
        {
            OP_System_MoveBeadTo_Inv(beadID);
            bAccept = 0;
            // printf("End LOCAL - Fail\n");
            return bAccept;
        }
}

int Move_Snake_Equil(int chainID, float MyTemp)
{ // Performs a slither MC-move on chainID

    int bAccept = 0; // Used in MC steps
    // Finding the bounds for looping over the molecule/chain
    const int firstB = chain_info_glb[chainID][CHAIN_START];
    const int lastB  = firstB + chain_info_glb[chainID][CHAIN_LENGTH];
    if (lastB - firstB == 1)
        { // This means we have a monomer. Reject the
          // move, because Local or Trans
            // moves should be the ones that move monomers.
            bAccept = 0;
            return bAccept;
        }
    else
        {
            if (nChainTypeIsLinear_glb[chain_info_glb[chainID][CHAIN_TYPE]] != 1)
                { // If chain is not linear. Reject move because
                  // slithering will not werk!
                    bAccept = 0;
                    return bAccept;
                }
        }
    // This chain is suitable to perform a slithering-snake move on.

    int i;     // Loop iterators
    int yTemp; // Random numbers to store things
    int r_posNew[POS_MAX], r_posTmp1[POS_MAX],
        r_posTmp2[POS_MAX];       // Vectors to store positions.
    int angBead_old, angBead_new; // BeadIDs for angular energies.

    lLDub MCProb        = (lLDub) rand() / (lLDub) RAND_MAX; // To decide if we slither forwards or backwards
    const char SnakeFwd = MCProb < 0.5f ? 1 : 0;             // 1 For forwards, 0 for backwards.
    if (SnakeFwd)
        { // Forwards slither, so lastB-1 (last bead) is anchor
            LatPos_gen_rand_wRad(r_posNew,
                                 linker_len_glb[lastB - 1][0]); // lastB-1 will be replaced by lastB-2
            LatPos_copy(r_posTmp1, r_posNew);
            LatPos_add_wPBC(r_posNew, bead_info_glb[lastB - 1], r_posTmp1);
            yTemp       = Check_MoveBeadTo(r_posNew); // 0: there is no space, 1: there is space
            angBead_old = firstB + 1;
            angBead_new = lastB - 2;
        }
    else
        { // Backwards slither, so firstB is anchor
            LatPos_gen_rand_wRad(r_posNew,
                                 linker_len_glb[firstB][0]); // firstB will be replaced by firstB+1
            LatPos_copy(r_posTmp1, r_posNew);
            LatPos_add_wPBC(r_posNew, bead_info_glb[firstB], r_posTmp1);
            yTemp       = Check_MoveBeadTo(r_posNew); // 0: there is no space, 1: there is space
            angBead_new = firstB + 1;
            angBead_old = lastB - 2;
        }

    if (yTemp == 0)
        { // Couldn't find a spot, so reject the damn move
            bAccept = 0;
            return bAccept;
        }
    // We should have a spot to move to! r_posNew has the location

    // Let's remember where this chain exists.
    OP_CopyBeadsToOld(firstB, lastB);

    lLDub newEn = 0.;
    lLDub oldEn = 0.;
    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = firstB; i < lastB; i++)
                {
                    oldEn += Energy_BiasingPotential(i);
                }
        }

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;

    for (i = firstB; i < lastB; i++)
        {
            Energy_Iso_ForChainsEquil(i, &oldEn, &newEn, &old_ovlp_num, &old_cont_num, naOldOvlpNeighs_glb,
                                      naOldContNeighs_glb);
        }
    const int resi = bead_info_glb[angBead_old][BEAD_TYPE];
    oldEn += faEnergy_glb[resi][resi][E_STIFF] ? Energy_Topo_Angle(angBead_old) : 0.;

    if (SnakeFwd)
        { // Slithering the chain forwards in ID-space
            OP_System_Snake_SlitherFwd(firstB, lastB, r_posNew);
        }
    else
        { // Slithering backwards in ID-space
            OP_System_Snake_SlitherBck(firstB, lastB, r_posNew);
        }

    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = firstB; i < lastB; i++)
                {
                    newEn += Energy_BiasingPotential(i);
                }
        }

    lLDub FSum = 0.;
    for (i = firstB; i < lastB; i++)
        {
            Energy_Iso_ForChainsEquil(i, &newEn, &oldEn, &new_ovlp_num, &new_cont_num, naNewOvlpNeighs_glb,
                                      naNewContNeighs_glb);
        }

    newEn += faEnergy_glb[resi][resi][E_STIFF] ? Energy_Topo_Angle(angBead_new) : 0.;

    // Doing the Metropolis-Hastings thing
    MCProb      = (lLDub) rand() / (lLDub) RAND_MAX;
    lLDub MHAcc = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub) MyTemp);
    if (MCProb < MHAcc)
        { // Accept the move. Remember that the bonds were assigned
          // above!Accept. Bonds have been handled before!
            bAccept = 1;
            return bAccept;
        }
    else
        {
            OP_RestoreChain_ForSnake(firstB, lastB);
            bAccept = 0;
            return bAccept;
        }
}

int Move_Trans_Equil(int chainID, float MyTemp)
{                               // Performs a translation move with orientational bias
    int bAccept = 0;            // Used in MC steps
    lLDub MCProb, oldEn, newEn; // For Metropolis Hastings
    oldEn = 0.;
    newEn = 0.;
    int i, j; // Loop iterators
    int firstB, lastB;
    int yTemp;           // Random numbers to store things
    int r_disp[POS_MAX]; // Vectors to store coordinates.
    // Finding the bounds for looping over the molecule/chain
    firstB = chain_info_glb[chainID][CHAIN_START];
    lastB  = firstB + chain_info_glb[chainID][CHAIN_LENGTH];
    // Radii for translation moves. All moves are L/4 radius

    LatPos_gen_rand_wRad(r_disp, naBoxSize_glb[0] / 2);

    yTemp = Check_ChainDisp(chainID, r_disp); // yTemp=0 means clash
    if (yTemp == 0)
        { // We have failed to find a good spot for this chain.
            bAccept = 0;
            // printf("Ending TRANS no space\n");
            return bAccept;
        }
    // We now have a chain which when moved does not overlap.

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;

    oldEn = 0.;
    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = firstB; i < lastB; i++)
                {
                    oldEn += Energy_BiasingPotential(i);
                }
        }

    for (i = firstB; i < lastB; i++)
        {
            Energy_Iso_ForChainsEquil(i, &oldEn, &newEn, &old_ovlp_num, &old_cont_num, naOldOvlpNeighs_glb,
                                      naOldContNeighs_glb);
        }

    OP_System_DispChain_ForTrans(chainID, r_disp); // Moved the chain, broke bonds, and remembered stuff

    newEn = 0.;
    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = firstB; i < lastB; i++)
                {
                    newEn += Energy_BiasingPotential(i);
                }
        }

    for (i = firstB; i < lastB; i++)
        {
            Energy_Iso_ForChainsEquil(i, &newEn, &oldEn, &new_ovlp_num, &new_cont_num, naNewOvlpNeighs_glb,
                                      naNewContNeighs_glb);
        }

    MCProb      = (lLDub) rand() / (lLDub) RAND_MAX;
    lLDub MHAcc = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub) MyTemp);
    if (MCProb < MHAcc)
        { // Accept the move. Remember that the bonds were
          // assigned above!
            bAccept = 1;
            return bAccept;
        }
    else
        {
            OP_System_RestoreChain_ForTrans(chainID);
            bAccept = 0;
            return bAccept;
        }
}

int Move_MultiLocal_Equil(int beadID, float MyTemp)
{

    int beadsList[MAX_BONDS + 1];
    const int beadNum = OP_GetTopoBonds(beadID, beadsList);
    int bAccept;

    if (beadNum < 2)
        { // Bead must be bonded to at least 2 things.
            bAccept = 0;
            return bAccept;
        }

    int beads_pos0[MAX_BONDS + 1][POS_MAX];
    int beads_posNew[MAX_BONDS + 1][POS_MAX];
    int r_posTmp1[POS_MAX];
    int tmpBead;
    int i;
    int yTemp;

    OP_Beads_CopyBeadsInListToOld(beadNum, beadsList);
    OP_Beads_CopyBeadsInListToPosList(beadNum, beadsList, beads_pos0);
    OP_Lattice_EmptySitesForListOfPos(beadNum, beads_pos0);

    for (i = 0; i < beadNum; i++)
        {
            LatPos_gen_rand_wRad(r_posTmp1, 2);
            LatPos_add_wPBC(beads_posNew[i], r_posTmp1, beads_pos0[i]);
            yTemp = Check_MoveBeadTo(beads_posNew[i]);
            if (yTemp == 0)
                {
                    break;
                }
            tmpBead                                            = beadsList[i];
            naTotLattice_glb[Lat_Ind_FromVec(beads_posNew[i])] = tmpBead;
        }

    OP_Lattice_EmptySitesForListOfPos(i, beads_posNew);

    if (yTemp == 1)
        { // No steric clash so check for topology constraint.
            OP_Beads_MoveBeadsInListToPos(beadNum, beadsList, beads_posNew);
            yTemp = Check_LinkerConstraints_ForBeadList(beadNum, beadsList);
            OP_Beads_MoveBeadsInListToPos(beadNum, beadsList, beads_pos0);
        }

    OP_Lattice_PlaceBeadsInList(beadNum, beadsList);

    if (yTemp == 0)
        {
            bAccept = 0;
            return bAccept;
        }

    lLDub oldEn = 0.;

    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = 0; i < beadNum; i++)
                {
                    tmpBead = beadsList[i];
                    oldEn += Energy_BiasingPotential(tmpBead);
                }
        }

    int old_ovlp_num, old_cont_num, new_ovlp_num, new_cont_num;

    lLDub newEn = 0.;
    for (i = 0; i < beadNum; i++)
        {
            Energy_Iso_ForLists(i, beadNum, beadsList, beads_pos0, &oldEn, &newEn, &old_ovlp_num, &old_cont_num,
                                naOldOvlpNeighs_glb, naOldContNeighs_glb);
        }

    int bondAngNum = 0;
    int tmpAngList[MAX_BONDS + 1];
    int bondAngList[(MAX_BONDS + 1) * (MAX_BONDS + 1)];
    if (bSystemHasTopo_glb)
        {
            bondAngNum = BeadList_AppendBeads(0, bondAngList, beadsList, beadNum);
            for (i = 1; i < beadNum; i++)
                {
                    tmpBead    = beadsList[i];
                    yTemp      = OP_GetTopoBonds(tmpBead, tmpAngList);
                    bondAngNum = BeadList_AppendBeads(bondAngNum, bondAngList, tmpAngList, yTemp);
                }
            qsort(bondAngList, bondAngNum, sizeof(int), UtilFunc_CompareInts);
            bondAngNum = ListOP_UniqueElementsOfSortedList_Int(bondAngNum, bondAngList);
        }
    else
        {
            bondAngNum = 0;
        }

    oldEn += bondAngNum ? Energy_Topo_Angle_ForList(bondAngNum, bondAngList) : 0.;

    OP_System_MoveBeadsInListToPos(beadNum, beadsList, beads_posNew);

    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = 0; i < beadNum; i++)
                {
                    tmpBead = beadsList[i];
                    newEn += Energy_BiasingPotential(tmpBead);
                }
        }

    for (i = 0; i < beadNum; i++)
        {
            Energy_Iso_ForLists(i, beadNum, beadsList, beads_posNew, &newEn, &oldEn, &new_ovlp_num, &new_cont_num,
                                naNewOvlpNeighs_glb, naNewContNeighs_glb);
        }

    newEn += bondAngNum ? Energy_Topo_Angle_ForList(bondAngNum, bondAngList) : 0.;

    lLDub MCProb = (lLDub) rand() / (lLDub) RAND_MAX;
    lLDub MHAcc  = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub) MyTemp);

    if (MCProb < MHAcc)
        { // Accept the move. Remember that the bonds were
          // assigned above
            bAccept = 1;
            return bAccept;
        }
    else
        {
            OP_Lattice_EmptySitesForListOfBeads(beadNum, beadsList);
            OP_Beads_CopyBeadsInListFromOld(beadNum, beadsList);
            OP_Lattice_PlaceBeadsInList(beadNum, beadsList);
            bAccept = 0;
            return bAccept;
        }
}

int Move_Pivot_Equil(int chainID, float MyTemp)
{
    // Performs a pivot move on chainID
    /*
    Randomly pick a bead in this chain (anchorBead), and perform a symmetry
    operation on the beads after the anchorBead. Note that if anchorBead is the
    second-to-last, or last, bead, the move is rejected. Local and slither moves
    would be better for those. The protein must be linear and the move is
    rejected if the chain is branched. It is assumed that a linear chain has
    been passed!
    */
    int bAccept = 0;
    if (! nChainTypeIsLinear_glb[chain_info_glb[chainID][CHAIN_TYPE]])
        {
            return bAccept;
        }
    // Check if the chain is longer than 3 or if it is linear
    const int chainLength = chain_info_glb[chainID][CHAIN_LENGTH];
    if (chainLength <= 3)
        {
            return bAccept;
        }

    int PivotM;
    PivotM = rand() % 10;
    if (PivotM == 0)
        { // Null move
            bAccept = 1;
            return bAccept;
        }

    // Get the beadID's for the first and last bead
    const int firstB = chain_info_glb[chainID][CHAIN_START];
    const int lastB  = firstB + chainLength;

    // Randomly select a bead that is neither the first nor last
    int anchorBead = chainLength - 2;
    anchorBead     = 1 + (rand() % anchorBead);

    //-1 Means backwards, +1 means forwards. Always Pivot the smaller portion
    const int PivotDir = anchorBead > chainLength / 2 ? 1 : -1;
    anchorBead         = firstB + anchorBead;
    // printf("Bead: %d\n", anchorBead);

    int i, j;
    int xTemp, yTemp;
    int anchorPos[POS_MAX];
    int beadsList[MAX_CHAINLEN];
    int beadNum = 0;

    if (PivotDir == 1)
        {
            for (i = anchorBead + 1; i < lastB; i++)
                {
                    beadsList[beadNum++] = i;
                }
        }
    else
        {
            for (i = firstB; i < anchorBead; i++)
                {
                    beadsList[beadNum++] = i;
                }
        }

    const int smBead = beadsList[0];
    const int lgBead = beadsList[beadNum - 1];

    LatPos_copy(anchorPos, bead_info_glb[anchorBead]);

    int tmpBead;

    for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            // Use anchorbead as origin, and check what happens after the rotation
            // to every bead
            OP_Rotation(PivotM, tmpBead, anchorPos);
            yTemp = Check_MoveBeadTo(naTempR_glb);
            if (yTemp == 0)
                {
                    break;
                }
        }
    if (yTemp == 0)
        {
            bAccept = 0;
            return bAccept;
        }

    lLDub oldEn = 0.;
    lLDub newEn = 0.;
    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = 0; i < beadNum; i++)
                {
                    tmpBead = beadsList[i];
                    oldEn += Energy_BiasingPotential(tmpBead);
                }
        }

    int old_ovlp_num, old_cont_num;

    lLDub BSum = 0.;
    for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            Energy_Iso_ForRangeEquil(tmpBead, smBead, lgBead, &oldEn, &newEn, &old_ovlp_num, &old_cont_num,
                                     naOldOvlpNeighs_glb, naOldContNeighs_glb);
        }

    const int resi = bead_info_glb[anchorBead][BEAD_TYPE];
    oldEn += faEnergy_glb[resi][resi][E_STIFF] ? Energy_Topo_Angle(anchorBead) : 0.;

    for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            OP_Rotation(PivotM, tmpBead, anchorPos);
            OP_System_MoveBeadTo(tmpBead, naTempR_glb);
        }

    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = 0; i < beadNum; i++)
                {
                    tmpBead = beadsList[i];
                    newEn += Energy_BiasingPotential(tmpBead);
                }
        }

    int new_ovlp_num, new_cont_num;

    lLDub FSum = 0.;
    for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            Energy_Iso_ForRangeEquil(tmpBead, smBead, lgBead, &newEn, &oldEn, &new_ovlp_num, &new_cont_num,
                                     naNewOvlpNeighs_glb, naNewContNeighs_glb);
        }

    newEn += faEnergy_glb[resi][resi][E_STIFF] ? Energy_Topo_Angle(anchorBead) : 0.;

    lLDub MCProb = (lLDub) rand() / (lLDub) RAND_MAX;
    lLDub MHAcc  = OP_GenMHValue(FSum, BSum, oldEn - newEn, (lLDub) MyTemp);
    if (MCProb < MHAcc)
        { // Accept the move. Remember that the bonds were
          // assigned above!
            bAccept = 1;
            return bAccept;
        }
    else
        { // Rejecting move
            for (i = 0; i < beadNum; i++)
                {
                    tmpBead = beadsList[i];
                    OP_System_MoveBeadTo_Inv(tmpBead);
                }
            bAccept = 0;
            return bAccept;
        }
}

int Move_BranchedRot_Equil(int chainID, float MyTemp)
{
    // Rotates a branched molecule about the branching, which is assumed to be
    // the firstB of chainID Performs a Move_Pivot() on molecule where the
    // rotation occurs around firstB
    /*
      Set the first bead as anchorBead, and perform a symmetry operation on the
      beads after the anchorBead (the whole molecule). Note that if the molecule
      is linear, the move is outright rejected. Again, it is assumed that the
      first bead in that molecule is the 'node'.
      */
    int bAccept = 0;
    // Reject if the molecule is linear
    if (nChainTypeIsLinear_glb[chain_info_glb[chainID][CHAIN_TYPE]])
        {
            return bAccept;
        }
    const int chainLength = chain_info_glb[chainID][CHAIN_LENGTH];

    int PivotM;
    PivotM = rand() % 10;
    if (PivotM == 0)
        { // Null move
            bAccept = 1;
            return bAccept;
        }

    // Get the beadID's for the first and last bead
    const int firstB = chain_info_glb[chainID][CHAIN_START];
    const int lastB  = firstB + chainLength;

    // FirstB is the anchor.
    const int anchorBead = firstB;

    int i;
    int yTemp;
    int anchorPos[POS_MAX];
    int beadsList[MAX_CHAINLEN];
    int beadNum = 0;

    for (i = anchorBead + 1; i < lastB; i++)
        {
            beadsList[beadNum++] = i;
        }

    const int smBead = beadsList[0];
    const int lgBead = beadsList[beadNum - 1];

    LatPos_copy(anchorPos, bead_info_glb[anchorBead]);

    int tmpBead;

    for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            // Use anchorbead as origin, and check what happens after the rotation
            // to every bead
            OP_Rotation(PivotM, tmpBead, anchorPos);
            yTemp = Check_MoveBeadTo(naTempR_glb);
            if (yTemp == 0)
                {
                    break;
                }
        }
    if (yTemp == 0)
        {
            bAccept = 0;
            return bAccept;
        }

    lLDub oldEn = 0.;
    lLDub newEn = 0.;
    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = 0; i < beadNum; i++)
                {
                    tmpBead = beadsList[i];
                    oldEn += Energy_BiasingPotential(tmpBead);
                }
        }

    int old_ovlp_num, old_cont_num;

    lLDub BSum = 0.;
    for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            Energy_Iso_ForRangeEquil(tmpBead, smBead, lgBead, &oldEn, &newEn, &old_ovlp_num, &old_cont_num,
                                     naOldOvlpNeighs_glb, naOldContNeighs_glb);
        }

    for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            OP_Rotation(PivotM, tmpBead, anchorPos);
            OP_System_MoveBeadTo(tmpBead, naTempR_glb);
        }

    if (nBiasPotential_Mode_glb != -1)
        {
            for (i = 0; i < beadNum; i++)
                {
                    tmpBead = beadsList[i];
                    newEn += Energy_BiasingPotential(tmpBead);
                }
        }

    int new_ovlp_num, new_cont_num;

    for (i = 0; i < beadNum; i++)
        {
            tmpBead = beadsList[i];
            Energy_Iso_ForRangeEquil(tmpBead, smBead, lgBead, &newEn, &oldEn, &new_ovlp_num, &new_cont_num,
                                     naNewOvlpNeighs_glb, naNewContNeighs_glb);
        }

    lLDub MCProb = (lLDub) rand() / (lLDub) RAND_MAX;
    lLDub MHAcc  = OP_GenMHValue(0., 0., oldEn - newEn, (lLDub) MyTemp);
    if (MCProb < MHAcc)
        { // Accept the move. Remember that the bonds were
          // assigned above!
            bAccept = 1;
            return bAccept;
        }
    else
        { // Rejecting move
            for (i = 0; i < beadNum; i++)
                {
                    tmpBead = beadsList[i];
                    OP_System_MoveBeadTo_Inv(tmpBead);
                }
            bAccept = 0;
            return bAccept;
        }
}

/// Check_ChainDisp - checks if moving chainID by vec_disp[POS_MAX] will lead to
/// sterix clash.
/// \param chainID
/// \param vec_disp
/// \return 0 if clash, 1 if no clash.
int Check_ChainDisp(const int chainID, const int* vec_disp)
{ // Checks if chain can be displaced by vec_disp
    int i;
    int r_chck[POS_MAX];
    const int firstB = chain_info_glb[chainID][CHAIN_START];
    const int lastB  = firstB + chain_info_glb[chainID][CHAIN_LENGTH];
    for (i = firstB; i < lastB; i++)
        {
            LatPos_add_wPBC(r_chck, vec_disp, bead_info_glb[i]);
            if (naTotLattice_glb[Lat_Ind_FromVec(r_chck)] != -1)
                { // Steric clash
                    return 0;
                }
        }

    return 1;
}

/// OP_System_DispChain - displaced chainID by movR[POS_MAX], while remembering the
/// chain in old_beads. Also handles the lattice placement.
/// \param chainID
/// \param movR
void OP_System_DispChain(int chainID, const int* movR)
{
    int i, l;
    const int fB = chain_info_glb[chainID][CHAIN_START];
    const int lB = fB + chain_info_glb[chainID][CHAIN_LENGTH];
    int r_pos0[POS_MAX], r_posNew[POS_MAX];

    // Copy beads to old.
    OP_CopyBeadsToOld(fB, lB);

    for (i = fB; i < lB; i++)
        {
            LatPos_copy(r_pos0, bead_info_glb[i]);
            LatPos_add_wPBC(bead_info_glb[i], r_pos0, movR);
            LatPos_copy(r_posNew, bead_info_glb[i]);

            naTotLattice_glb[Lat_Ind_FromVec(r_pos0)]   = -1; // Removing from old place
            naTotLattice_glb[Lat_Ind_FromVec(r_posNew)] = i;  // Adding to new place
        }
}

/// OP_System_DispChain_ForTrans - displaced chainID by movR[POS_MAX], while
/// remembering the chain in old_beads. Also handles the lattice placement.
/// Move_Trans variant where I break all the physical bonds.
/// \param chainID
/// \param movR
void OP_System_DispChain_ForTrans(const int chainID, const int* movR)
{
    // Displaces current chain by movR and handles lattice
    // Specific for Move_Trans because it breaks old bonds!
    // Also remembers where everything was moved and saves into old_bead_glb

    int i, l;
    const int fB = chain_info_glb[chainID][CHAIN_START];
    const int lB = fB + chain_info_glb[chainID][CHAIN_LENGTH];
    int r_pos0[POS_MAX], r_posNew[POS_MAX];

    // Copy beads to old.
    OP_CopyBeadsToOld(fB, lB);

    for (i = fB; i < lB; i++)
        {
            LatPos_copy(r_pos0, bead_info_glb[i]);
            LatPos_add_wPBC(bead_info_glb[i], r_pos0, movR);
            LatPos_copy(r_posNew, bead_info_glb[i]);

            naTotLattice_glb[Lat_Ind_FromVec(r_pos0)]   = -1; // Removing from old place
            naTotLattice_glb[Lat_Ind_FromVec(r_posNew)] = i;  // Adding to new place
            OP_Beads_BreakBond(i);
        }
}

/// OP_System_RestoreChain - uses old_bead_glb to undo what OP_System_DispChain does.
/// \param chainID
void OP_System_RestoreChain(int chainID)
{ // Uses old_bead_glb to undo what OP_System_DispChain does.
    int i, l;
    const int fB = chain_info_glb[chainID][CHAIN_START];
    const int lB = fB + chain_info_glb[chainID][CHAIN_LENGTH];

    // Removing from 'new' place.
    for (i = fB; i < lB; ++i)
        {
            naTotLattice_glb[Lat_Ind_OfBead(i)] = -1;
        }

    // Remembering where we were.
    for (i = fB; i < lB; ++i)
        {
            OP_CopyBead(bead_info_glb[i], old_bead_glb[i]);
        }

    // Placing back.
    for (i = fB; i < lB; ++i)
        {
            naTotLattice_glb[Lat_Ind_OfBead(i)] = i;
        }
}

/// OP_CopyBead: Copies BEADINFO_MAX elements of orig_arr into copy_arr
/// \param copy_arr
/// \param orig_arr
inline void OP_CopyBead(int* restrict copy_arr, const int* restrict orig_arr)
{
    int j;
    for (j = 0; j < BEADINFO_MAX; j++)
        {
            copy_arr[j] = orig_arr[j];
        }
}

/// OP_CopyBeadsToOld: Copies all beads in [firstB, lastB) to old_beads. This is
/// sort of a wrapper for ease.
/// \param firstB
/// \param lastB
inline void OP_CopyBeadsToOld(const int firstB, const int lastB)
{
    int i;
    for (i = firstB; i < lastB; i++)
        {
            OP_CopyBead(old_bead_glb[i], bead_info_glb[i]);
        }
}

/// OP_RestoreBeadsFromOld: Copies all beads in [firstB, lastB) from old_beads.
/// This is sort of a wrapper for ease.
/// \param firstB
/// \param lastB
inline void OP_RestoreBeadsFromOld(const int firstB, const int lastB)
{
    int i, j;
    for (i = firstB; i < lastB; i++)
        {
            OP_CopyBead(bead_info_glb[i], old_bead_glb[i]);
        }
}

/// OP_System_RestoreChain_ForTrans - translation variant to restore chainID.
/// \param chainID
void OP_System_RestoreChain_ForTrans(int chainID)
{ // Uses old_bead_glb to undo what OP_System_DispChain_ForTrans does.
    // Note that removing and placing separately makes sure that no incorrect
    // bonds are formed, or unformed.
    int i, l;
    const int fB = chain_info_glb[chainID][CHAIN_START];
    const int lB = fB + chain_info_glb[chainID][CHAIN_LENGTH];

    // Breaking all current bonds
    for (i = fB; i < lB; ++i)
        {
            OP_Beads_BreakBond(i);
        }

    // Removing from 'new' place.
    for (i = fB; i < lB; ++i)
        {
            naTotLattice_glb[Lat_Ind_OfBead(i)] = -1;
        }

    // Remembering where we were.
    for (i = fB; i < lB; ++i)
        {
            OP_CopyBead(bead_info_glb[i], old_bead_glb[i]);
        }

    // Placing back.
    for (i = fB; i < lB; ++i)
        {
            naTotLattice_glb[Lat_Ind_OfBead(i)] = i;
        }

    // Restoring older bonds.
    for (i = fB; i < lB; ++i)
        {
            OP_Beads_RestoreBond(i);
        }
}

/// OP_RestoreChain_ForSnake - restores the chain: beadID's betweem fB and lB-1,
/// to before the snake move.
/// \param fB
/// \param lB
void OP_RestoreChain_ForSnake(const int fB, const int lB)
{
    int i, j;

    for (i = fB; i < lB; i++)
        { // Resetting the lattice
            if (bead_info_glb[i][BEAD_FACE] != -1)
                { // Need to break the newly proposed bond
                    bead_info_glb[bead_info_glb[i][BEAD_FACE]][BEAD_FACE] = -1;
                }
            naTotLattice_glb[Lat_Ind_OfBead(i)] = -1;
        }
    for (i = fB; i < lB; i++)
        {
            for (j = 0; j < BEADINFO_MAX; j++)
                {
                    bead_info_glb[i][j] = old_bead_glb[i][j]; // Restoring
                }
            if (bead_info_glb[i][BEAD_FACE] != -1)
                { // I was bonded so restore the bond
                    bead_info_glb[bead_info_glb[i][BEAD_FACE]][BEAD_FACE] = i;
                }
            naTotLattice_glb[Lat_Ind_OfBead(i)] = i;
        }
}

/// Check_MoveBeadTo - checks if newPos[POS_MAX] is empty on the lattice.
/// \param newPos
/// \return
inline int Check_MoveBeadTo(const int* newPos)
{ // Checks if I can move here

    if (naTotLattice_glb[Lat_Ind_FromVec(newPos)] != -1)
        {
            return 0; // There's something here already, brah
        }
    return 1;
}

/// OP_System_MoveBeadTo - move beadID to newPos[POS_MAX], while remembering the
/// bead properties, and handling the lattice.
/// \param beadID
/// \param newPos
void OP_System_MoveBeadTo(const int beadID, const int* newPos)
{ // Updates position to new one and handles lattice
    OP_CopyBead(old_bead_glb[beadID], bead_info_glb[beadID]);
    LatPos_copy(bead_info_glb[beadID], newPos);
    naTotLattice_glb[Lat_Ind_FromVec(old_bead_glb[beadID])] = -1;     // Removing from old place.
    naTotLattice_glb[Lat_Ind_FromVec(newPos)]               = beadID; // Placing in new place.
}

/// OP_System_MoveBeadTo_Inv - performs the inverse of OP_System_MoveBeadTo to
/// restore beadID using old_bead_glb.
/// \param beadID
void OP_System_MoveBeadTo_Inv(int beadID)
{                                                  // Undoes what OP_System_MoveBeadTo does
    naTotLattice_glb[Lat_Ind_OfBead(beadID)] = -1; // Removing from newly proposed place.
    OP_CopyBead(bead_info_glb[beadID], old_bead_glb[beadID]);
    naTotLattice_glb[Lat_Ind_OfBead(beadID)] = beadID; // Placing back where we used to be.
    OP_Beads_RestoreBond(beadID);
}

/// OP_MoveBeadTo_ForMTLocal - moves the bead over to newPos[POS_MAX], handles
/// the lattice and breaks bonds. This function is specifically for the
/// MultiLocal move.
/// \param beadID
/// \param newPos
void OP_MoveBeadTo_ForMTLocal(int beadID, const int* newPos)
{ // Updates position to new one and handles
  // lattice specifically for shake move
    int i;
    int tmpR2[POS_MAX];
    for (i = 0; i < POS_MAX; i++)
        {
            bead_info_glb[beadID][i] = newPos[i];
            tmpR2[i]                 = bead_info_glb[beadID][i];
        }
    naTotLattice_glb[Lat_Ind_FromVec(tmpR2)] = beadID;
    i                                        = bead_info_glb[beadID][BEAD_FACE];
    if (i != -1)
        {
            bead_info_glb[i][BEAD_FACE]      = -1;
            bead_info_glb[beadID][BEAD_FACE] = -1;
        }
}

void OP_System_MoveBeadsInListToPos(const int listSize, const int* beadList, const int (*newPos)[POS_MAX])
{

    OP_Lattice_EmptySitesForListOfBeads(listSize, beadList);
    OP_Beads_MoveBeadsInListToPos(listSize, beadList, newPos);
    OP_Lattice_PlaceBeadsInList(listSize, beadList);
}

void OP_Inv_MoveBeadsTo_InList(const int listSize, const int beadList[MAX_BONDS + 1],
                               const int newPos[MAX_BONDS + 1][POS_MAX])
{
    OP_Lattice_EmptySitesForListOfBeads(listSize, beadList);
    //    OP_Lattice_PlaceBeadsInList()
}

/// OP_SwapBeads - swaps all the properties, including bonding partners between
/// bead1 and bead2. Used in the DbPvt move.
/// \param bead1
/// \param bead2
void OP_SwapBeads(int bead1, int bead2)
{
    // Swaps ALL properties of the beads, and exchanged bonding partners
    int MyF1, MyF2;
    int i;
    int tmpR[POS_MAX], tmpR2[POS_MAX];
    // First the coordinates
    for (i = 0; i < POS_MAX; i++)
        {
            tmpR[i]                 = bead_info_glb[bead1][i];
            tmpR2[i]                = bead_info_glb[bead2][i];
            bead_info_glb[bead1][i] = tmpR2[i];
            bead_info_glb[bead2][i] = tmpR[i];
        }

    // Onto bonds and partners
    MyF1 = bead_info_glb[bead1][BEAD_FACE]; // Bead 1's bonding partner
    MyF2 = bead_info_glb[bead2][BEAD_FACE]; // Bead 2's bonding partner.
    // If they are bonded to each other, don't do anything
    if (MyF1 != bead2)
        {
            bead_info_glb[bead1][BEAD_FACE] = MyF2;
            bead_info_glb[bead2][BEAD_FACE] = MyF1;
            if (MyF1 != -1)
                {                                           // Need to swap partners -- very 2019
                    bead_info_glb[MyF1][BEAD_FACE] = bead2; // It's bonded to bead2 now
                }
            if (MyF2 != -1)
                {
                    bead_info_glb[MyF2][BEAD_FACE] = bead1;
                }
        }
    // Swap them on the lattice
    naTotLattice_glb[Lat_Ind_FromVec(tmpR)]  = bead2;
    naTotLattice_glb[Lat_Ind_FromVec(tmpR2)] = bead1;
}

/// OP_Rotation - given the rotation operation PivotM, rotate beadID using
/// tmpR[POS_MAX] as the origin. Note the global array naTempR_glb[POS_MAX] stores
/// the locations of the post-rotated beads. Furthermore, if PivotM is 10, the
/// null operation is still performed. Mathematically, we have that
/// f$vec{r}^prime = hat{R}_i(theta)vec{r}f$ where f$vec{r}f$ is the new
/// position, f$hat{R}_i(theta)f$ is the standard rotation matrix where i=X,Y,Z
/// and $theta=45^circ,180^circ,270^circ$. For now, only 90 degree rotations are
/// implemented where the Rot_{#1}_{#2} functions/sub-routines below explicitly
/// calculate what the new vector is given the axis ({#1}) and angle ({#2}).
/// Aribitrary rotations are a little harder on lattices because some rotations
/// do not fully overlap with points on the lattice being used. In our case, it
/// is a simple cubic lattice.
/// \param PivotM - this is the operation to perform
/// \param beadID - the ID of the bead to move
/// \param tmpR   - the location of the anchor point, or, origin around which to
/// rotate
void OP_Rotation(int PivotM, int beadID, int* tmpR)
{
    int j;
    switch (PivotM)
        {
            case 1:
                Rot_X_90(beadID, tmpR);
                break;

            case 2:
                Rot_X_180(beadID, tmpR);
                break;

            case 3:
                Rot_X_270(beadID, tmpR);
                break;

            case 4:
                Rot_Y_90(beadID, tmpR);
                break;

            case 5:
                Rot_Y_180(beadID, tmpR);
                break;

            case 6:
                Rot_Y_270(beadID, tmpR);
                break;

            case 7:
                Rot_Z_90(beadID, tmpR);
                break;

            case 8:
                Rot_Z_180(beadID, tmpR);
                break;

            case 9:
                Rot_Z_270(beadID, tmpR);
                break;

            default:
                // printf("How\n");
                for (j = 0; j < POS_MAX; j++)
                    {
                        naTempR_glb[j] = bead_info_glb[beadID][j]; // Does nothing; should make the move fail.
                    }
                break;
        }
}

void Rot_X_90(int beadID, const int tmpR[POS_MAX])
{

    int j;
    // Performs a rotation by 90 degrees around x-axis. tmpR is the new origin.
    for (j = 0; j < POS_MAX; j++)
        {
            naTempR_glb[j] = bead_info_glb[beadID][j];
        }
    naTempR_glb[POS_Y] = tmpR[POS_Y] + tmpR[POS_Z] - bead_info_glb[beadID][POS_Z];
    naTempR_glb[POS_Z] = bead_info_glb[beadID][POS_Y] - tmpR[POS_Y] + tmpR[POS_Z];
    // Adjusting for periodic boundaries
    naTempR_glb[POS_Y] = (naTempR_glb[POS_Y] + naBoxSize_glb[POS_Y]) % naBoxSize_glb[POS_Y];
    naTempR_glb[POS_Z] = (naTempR_glb[POS_Z] + naBoxSize_glb[POS_Z]) % naBoxSize_glb[POS_Z];
}

void Rot_X_180(int beadID, const int tmpR[POS_MAX])
{

    int j;
    // Performs a rotation by 180 degrees around x-axis. tmpR is the new origin
    // for the rotation.
    for (j = 0; j < POS_MAX; j++)
        {
            naTempR_glb[j] = bead_info_glb[beadID][j];
        }
    naTempR_glb[POS_Y] = 2 * tmpR[POS_Y] - bead_info_glb[beadID][POS_Y];
    naTempR_glb[POS_Z] = 2 * tmpR[POS_Z] - bead_info_glb[beadID][POS_Z];
    // Adjusting for periodic boundaries
    naTempR_glb[POS_Y] = (naTempR_glb[POS_Y] + naBoxSize_glb[POS_Y]) % naBoxSize_glb[POS_Y];
    naTempR_glb[POS_Z] = (naTempR_glb[POS_Z] + naBoxSize_glb[POS_Z]) % naBoxSize_glb[POS_Z];
}

void Rot_X_270(int beadID, const int tmpR[POS_MAX])
{

    int j;
    // Performs a rotation by 270 degrees around x-axis. tmpR is the new origin.
    for (j = 0; j < POS_MAX; j++)
        {
            naTempR_glb[j] = bead_info_glb[beadID][j];
        }
    naTempR_glb[POS_Y] = bead_info_glb[beadID][POS_Z] + tmpR[POS_Y] - tmpR[POS_Z];
    naTempR_glb[POS_Z] = tmpR[POS_Y] + tmpR[POS_Z] - bead_info_glb[beadID][POS_Y];
    // Adjusting for periodic boundaries
    naTempR_glb[POS_Y] = (naTempR_glb[POS_Y] + naBoxSize_glb[POS_Y]) % naBoxSize_glb[POS_Y];
    naTempR_glb[POS_Z] = (naTempR_glb[POS_Z] + naBoxSize_glb[POS_Z]) % naBoxSize_glb[POS_Z];
}

void Rot_Y_90(int beadID, const int tmpR[POS_MAX])
{
    int j;
    // Performs a rotation by 90 degrees around x-axis. tmpR is the new origin.
    for (j = 0; j < POS_MAX; j++)
        {
            naTempR_glb[j] = bead_info_glb[beadID][j];
        }
    naTempR_glb[POS_X] = (bead_info_glb[beadID][POS_Z] - tmpR[POS_Z] + tmpR[POS_X]);
    naTempR_glb[POS_Z] = (tmpR[POS_X] + tmpR[POS_Z] - bead_info_glb[beadID][POS_X]);
    // Adjusting for periodic boundaries
    naTempR_glb[POS_X] = (naTempR_glb[POS_X] + naBoxSize_glb[POS_X]) % naBoxSize_glb[POS_X];
    naTempR_glb[POS_Z] = (naTempR_glb[POS_Z] + naBoxSize_glb[POS_Z]) % naBoxSize_glb[POS_Z];
}

void Rot_Y_180(int beadID, const int tmpR[POS_MAX])
{
    int j;
    // Performs a rotation by 270 degrees around x-axis. tmpR is the new origin.
    for (j = 0; j < POS_MAX; j++)
        {
            naTempR_glb[j] = bead_info_glb[beadID][j];
        }
    naTempR_glb[POS_X] = (2 * tmpR[POS_X] - bead_info_glb[beadID][POS_X]);
    naTempR_glb[POS_Z] = (2 * tmpR[POS_Z] - bead_info_glb[beadID][POS_Z]);
    // Adjusting for periodic boundaries
    naTempR_glb[POS_X] = (naTempR_glb[POS_X] + naBoxSize_glb[POS_X]) % naBoxSize_glb[POS_X];
    naTempR_glb[POS_Z] = (naTempR_glb[POS_Z] + naBoxSize_glb[POS_Z]) % naBoxSize_glb[POS_Z];
}

void Rot_Y_270(int beadID, const int tmpR[POS_MAX])
{
    int j;
    // Performs a rotation by 270 degrees around x-axis. tmpR is the new origin.
    for (j = 0; j < POS_MAX; j++)
        {
            naTempR_glb[j] = bead_info_glb[beadID][j];
        }
    naTempR_glb[POS_X] = (tmpR[POS_X] + tmpR[POS_Z] - bead_info_glb[beadID][POS_Z]);
    naTempR_glb[POS_Z] = (bead_info_glb[beadID][POS_X] - tmpR[POS_X] + tmpR[POS_Z]);
    // Adjusting for periodic boundaries
    naTempR_glb[POS_X] = (naTempR_glb[POS_X] + naBoxSize_glb[POS_X]) % naBoxSize_glb[POS_X];
    naTempR_glb[POS_Z] = (naTempR_glb[POS_Z] + naBoxSize_glb[POS_Z]) % naBoxSize_glb[POS_Z];
}

void Rot_Z_90(int beadID, const int tmpR[POS_MAX])
{
    int j;
    // Performs a rotation by 270 degrees around x-axis. tmpR is the new origin.
    for (j = 0; j < POS_MAX; j++)
        {
            naTempR_glb[j] = bead_info_glb[beadID][j];
        }
    naTempR_glb[POS_X] = (tmpR[POS_X] + tmpR[POS_Y] - bead_info_glb[beadID][POS_Y]);
    naTempR_glb[POS_Y] = (bead_info_glb[beadID][POS_X] - tmpR[POS_X] + tmpR[POS_Y]);
    // Adjusting for periodic boundaries
    naTempR_glb[POS_X] = (bead_info_glb[beadID][POS_X] + naBoxSize_glb[POS_X]) % naBoxSize_glb[POS_X];
    naTempR_glb[POS_Y] = (bead_info_glb[beadID][POS_Y] + naBoxSize_glb[POS_Y]) % naBoxSize_glb[POS_Y];
}

void Rot_Z_180(int beadID, const int tmpR[POS_MAX])
{
    int j;
    // Performs a rotation by 270 degrees around x-axis. tmpR is the new origin.
    for (j = 0; j < POS_MAX; j++)
        {
            naTempR_glb[j] = bead_info_glb[beadID][j];
        }
    naTempR_glb[POS_X] = (2 * tmpR[POS_X] - bead_info_glb[beadID][POS_X]);
    naTempR_glb[POS_Y] = (2 * tmpR[POS_Y] - bead_info_glb[beadID][POS_Y]);
    // Adjusting for periodic boundaries
    naTempR_glb[POS_X] = (bead_info_glb[beadID][POS_X] + naBoxSize_glb[POS_X]) % naBoxSize_glb[POS_X];
    naTempR_glb[POS_Y] = (bead_info_glb[beadID][POS_Y] + naBoxSize_glb[POS_Y]) % naBoxSize_glb[POS_Y];
}

void Rot_Z_270(int beadID, const int tmpR[POS_MAX])
{
    int j;
    // Performs a rotation by 270 degrees around x-axis. tmpR is the new origin.
    for (j = 0; j < POS_MAX; j++)
        {
            naTempR_glb[j] = bead_info_glb[beadID][j];
        }
    naTempR_glb[POS_X] = (bead_info_glb[beadID][POS_Y] + tmpR[POS_X] - tmpR[POS_Y]);
    naTempR_glb[POS_Y] = (tmpR[POS_X] + tmpR[POS_Y] - bead_info_glb[beadID][POS_X]);
    // Adjusting for periodic boundaries
    naTempR_glb[POS_X] = (bead_info_glb[beadID][POS_X] + naBoxSize_glb[POS_X]) % naBoxSize_glb[POS_X];
    naTempR_glb[POS_Y] = (bead_info_glb[beadID][POS_Y] + naBoxSize_glb[POS_Y]) % naBoxSize_glb[POS_Y];
}

/// OP_ShuffleRotIndecies - randomly shuffles the incdecies to sample the
/// rotational states around a bead. Ensures that we fully, but randomly, sample
/// the states. The implementation was taken from the provided URL.
void OP_ShuffleRotIndecies(void)
{
    // Algorithm taken from
    // https://www.w3resource.com/c-programming-exercises/array/c-array-exercise-77.php
    // Takes the array and shuffles the numbers in it. Used to randomly, but
    // fully, sample the rotational states.

    int i, j;
    int i_val, j_val;

    for (i = MAX_ROTSTATES - 2; i > 0; i--)
        {
            j                   = rand() % (i + 1);
            i_val               = naRot_IndArr_glb[i];
            j_val               = naRot_IndArr_glb[j];
            naRot_IndArr_glb[i] = j_val;
            naRot_IndArr_glb[j] = i_val;
        }
}

/// OP_ShuffleArray: Given an array of at-most size arr_size, we shuffle the
/// elements in dum_arr. Note that if arr_size is smaller than the actual size
/// of the array, only the elements less than arr_size shall be shuffled. So
/// this _could_ be used to shuffle a starting subset of the array.
/// \param arr_size
/// \param dum_arr
void OP_ShuffleArray(const int arr_size, int* dum_arr)
{
    // Algorithm taken from
    // https://www.w3resource.com/c-programming-exercises/array/c-array-exercise-77.php
    // Takes the array and shuffles the numbers in it. Used to randomly, but
    // fully, sample arrays.
    int i;
    int jVal, j;
    for (i = arr_size - 1; i > 0; i--)
        {
            j          = rand() % (i + 1);
            jVal       = dum_arr[i];
            dum_arr[i] = dum_arr[j];
            dum_arr[j] = jVal;
        }
}

/// Check_RotStates_wNeighList: Given this bead and the supplied neighbor-list
/// and numbers, we see which of the neighbors are viable anisotropic partners.
/// For each viable bead, we also record the Rosenbluth factors, and store which
/// bead it is.
/// \param beadID
/// \param resi
/// \param neighList
/// \param neighNum
/// \return The number of beads that beadID could bond to.
int Check_RotStates_wNeighList(int const beadID, int const resi, const int* neighList, int const neighNum)
{
    int CandNums = 0;
    int tmpBead;
    int resj, i;

    for (i = 0; i < neighNum; i++)
        {
            tmpBead = neighList[i];
            resj    = bead_info_glb[tmpBead][BEAD_TYPE];

            if (faEnergy_glb[resi][resj][E_SC_SC] != 0.f)
                {
                    if (bead_info_glb[tmpBead][BEAD_FACE] == -1 || bead_info_glb[tmpBead][BEAD_FACE] == beadID)
                        {
                            ldaBoltzFac_glb[CandNums]   = ldaBoltzFacNorm_glb[resi][resj];
                            naRotTrial_glb[0][CandNums] = tmpBead;
                            CandNums++;
                        }
                }
        }

    return CandNums;
}

/// Check_RotStatesOld -- goes around beadID, of type resi to built the
/// boltzmann distribution if the rotational states. Used in the old location.
/// Note that this implementation is one where the solvent interaction does not
/// exist for the rotational state. The commented version below is one where the
/// solvent also has a contribution to the boltzmann distribution of the
/// rotational states.
/// \param beadID
/// \param resi
/// \param MyTemp
/// \return The total number of possible candidates
int Check_RotStatesOld(int const beadID, int const resi, float const MyTemp)
{

    int i, j, k, tmpBead;
    int CandNums = 0;
    int r_pos0[POS_MAX], r_search[POS_MAX];
    for (j = 0; j < POS_MAX; j++)
        {
            r_pos0[j] = bead_info_glb[beadID][j];
        }
    for (k = 0; k < MAX_ROTSTATES - 1; k++)
        {
            i = naRot_IndArr_glb[k];

            LatPos_add_wPBC(r_search, r_pos0, naLocalArr_glb[i]);

            tmpBead = Lat_Ind_FromVec(r_search);
            tmpBead = naTotLattice_glb[tmpBead];
            if (tmpBead != -1)
                {
                    j = bead_info_glb[tmpBead][BEAD_TYPE];
                    if (faEnergy_glb[resi][j][E_SC_SC] != 0 &&
                        (bead_info_glb[tmpBead][BEAD_FACE] == -1 || bead_info_glb[tmpBead][BEAD_FACE] == beadID))
                        {
                            ldaBoltzFac_glb[CandNums]   = ldaBoltzFacNorm_glb[resi][j];
                            naRotTrial_glb[0][CandNums] = tmpBead;
                            CandNums++;
                        }
                }
        }

    return CandNums;
}

/// Check_RotStatesNew -- goes around beadID, of type resi to built the
/// boltzmann distribution if the rotational states. Used in the new proposed
/// location. Note that this implementation is one where the solvent interaction
/// does not exist for the rotaitonal state. The commented version below is one
/// where the solvent also has a contribution to the boltzmann distribution of
/// the rotational states.
/// \param beadID
/// \param resi
/// \param MyTemp
/// \return The total number of possible candidates
int Check_RotStatesNew(int const beadID, int const resi, float const MyTemp)
{

    int i, j, k, tmpBead;
    int CandNums = 0;
    int tmpR[POS_MAX], tmpR2[POS_MAX];
    for (j = 0; j < POS_MAX; j++)
        {
            tmpR[j] = bead_info_glb[beadID][j];
        }
    for (k = 0; k < MAX_ROTSTATES - 1; k++)
        {
            i = naRot_IndArr_glb[k];
            for (j = 0; j < POS_MAX; j++)
                {
                    tmpR2[j] = (tmpR[j] + naLocalArr_glb[i][j] + naBoxSize_glb[j]) % naBoxSize_glb[j];
                }
            tmpBead = Lat_Ind_FromVec(tmpR2);
            tmpBead = naTotLattice_glb[tmpBead];
            if (tmpBead != -1)
                {
                    j = bead_info_glb[tmpBead][BEAD_TYPE];
                    if (faEnergy_glb[resi][j][E_SC_SC] != 0 &&
                        (bead_info_glb[tmpBead][BEAD_FACE] == -1 || bead_info_glb[tmpBead][BEAD_FACE] == beadID))
                        {
                            ldaBoltzFac_glb[CandNums]   = ldaBoltzFacNorm_glb[resi][j];
                            naRotTrial_glb[0][CandNums] = tmpBead;
                            CandNums++;
                        }
                }
        }

    return CandNums;
}

/// OP_NormalizeRotState -- normalizes the rotational boltzmann distribution,
/// and generates the cumulative distribution so that we can sample from it.
/// \param beadVal - ldaBoltzNorm_glb[MAX_VALENCY] stores the total Rosenbluth weights
/// for moves that sample multiple stickers. Note that ldaBoltzFac_glb[] is now the
/// cumulative boltzmann distribution, which is sampled to propose bonds.
/// \param CandNums - total possible bonding partners for this particular sticker case.
void OP_NormalizeRotState(const int beadVal, const int CandNums)
{
    int i;

    if (CandNums > 0)
        { // There is a possible candidate, so normalize ldaBoltzFac_glb
            ldaBoltzNorm_glb[beadVal] = 0.;
            for (i = 0; i < CandNums; i++)
                {
                    ldaBoltzNorm_glb[beadVal] += ldaBoltzFac_glb[i];
                }
            for (i = 0; i < CandNums; i++)
                {
                    ldaBoltzFac_glb[i] /= ldaBoltzNorm_glb[beadVal];
                }
            for (i = 1; i < CandNums; i++)
                {
                    ldaBoltzFac_glb[i] += ldaBoltzFac_glb[i - 1];
                }
        }
    else
        {
            ldaBoltzNorm_glb[beadVal] = 1.; // If no candidates, we set it to 1 because
                                            // this will not be used
        }
}

/// OP_PickRotState - propose a new bonding partner.
/// This is the implementation where the solvent is ignored. Thus, there is a
/// 1/(CandNums+1) to break the bond, and (CandNums)/(CandNums+1) to form a
/// bond. If forming a bond, sample from the boltzmann distribution of
/// rotational states already built up and stored in ldaBoltzFac_glb[]
/// \param CandNums - the number of possible bonding candidates.
/// \return The new proposed bonding partner.
int OP_PickRotState(int CandNums)
{
    int newRot = -1;
    int i;
    lLDub fProb;
    int nCheck = CandNums + 1;
    nCheck     = rand() % nCheck;
    if (nCheck == 0)
        {
            newRot = -1;
        }
    else
        {
            fProb = (lLDub) rand() / (lLDub) RAND_MAX;
            for (i = 0; i < CandNums; i++)
                {
                    if (fProb < ldaBoltzFac_glb[i])
                        {
                            break;
                        }
                }
            newRot = naRotTrial_glb[0][i];
        }
    return newRot;
}

/// OP_GenMHValue - calculates the acceptance ratio for a given state.
/// In this implementation, it is assumed that all the inputs to the function
/// correspond to the logl() of the values (to keep the numbers small)
/// If the total value x, in log space, is larger than 0, it means that expl(x)
/// > 1. Therefore the move will be accepted regardless, so the function shall
/// return 2 (a value larger than 1). Similarly, if x < -21.5, expl(x)
/// ~1x10^{-10} which is smaller than 1/RAND_MAX; so it is essentially 0 -- thus
/// the function returns 0. In all other cases, the function returns expl(x),
/// which is then compared to a number between 0 and 1. Remember that
/// log(exp(a)) = a, so log(exp((old_en-new_en)/myTemp)) =
/// (old_en-new_en)/myTemp
/// \param fRos: Forward Rosenbluth weight.
/// \param bRos: Backwards Rosenbluth weight.
/// \param Delta_En: Difference in energy.
/// \param Cur_Temp: Current temperature.
/// \return The new proposed bonding partner.
lLDub OP_GenMHValue(lLDub fRos, lLDub bRos, lLDub Delta_En, lLDub Cur_Temp)
{
    lLDub MH_Value;

    MH_Value = fRos - bRos + Delta_En / Cur_Temp;
    if (MH_Value <= ldLogOfSmallestPossibleProb_glb)
        {
            MH_Value = 0.;
        }
    else if (MH_Value >= 0.)
        {
            MH_Value = 2.0;
        }
    else
        {
            MH_Value = expl(MH_Value);
        }
    return MH_Value;
}

///
/// \param firstB
/// \param lastB
/// \param r_posNew
void OP_System_Snake_SlitherFwd(const int firstB, const int lastB, const int* r_posNew)
{
    int i;
    int r_tmp[POS_MAX];
    naTotLattice_glb[Lat_Ind_OfBead(firstB)] = -1; // Only the firstB's location is empty
    for (i = firstB; i < lastB - 1; i++)
        {
            LatPos_copy(bead_info_glb[i], old_bead_glb[i + 1]); // Hopping over by one bead
            LatPos_copy(r_tmp, bead_info_glb[i]);
            naTotLattice_glb[Lat_Ind_FromVec(r_tmp)] = i;
        }
    // Moving the last bead, and it has to be done independently because lastB-1
    // -> r_posNew
    i = lastB - 1;
    LatPos_copy(bead_info_glb[i], r_posNew);
    naTotLattice_glb[Lat_Ind_FromVec(r_posNew)] = i;
}

void OP_System_Snake_SlitherBck(const int firstB, const int lastB, const int* r_posNew)
{
    int i;
    int r_tmp[POS_MAX];
    naTotLattice_glb[Lat_Ind_OfBead(lastB - 1)] = -1; // Only the lastB-1's location is empty
    for (i = firstB + 1; i < lastB; i++)
        {
            LatPos_copy(bead_info_glb[i], old_bead_glb[i - 1]); // Hopping back by one bead
            LatPos_copy(r_tmp, bead_info_glb[i]);
            naTotLattice_glb[Lat_Ind_FromVec(r_tmp)] = i;
        }

    // Moving the first bead, and it has to be done independently because firstB
    // -> r_posNew
    i = firstB;
    LatPos_copy(bead_info_glb[i], r_posNew);
    naTotLattice_glb[Lat_Ind_FromVec(r_posNew)] = i;
}

void OP_Beads_CopyBeadsInListToOld(const int listSize, const int* beadList)
{

    int i, thisBead;
    for (i = 0; i < listSize; i++)
        {
            thisBead = beadList[i];
            OP_CopyBead(old_bead_glb[thisBead], bead_info_glb[thisBead]);
        }
}

void OP_Beads_CopyBeadsInListFromOld(const int listSize, const int* beadList)
{

    int i, thisBead;
    for (i = 0; i < listSize; i++)
        {
            thisBead = beadList[i];
            OP_CopyBead(bead_info_glb[thisBead], old_bead_glb[thisBead]);
        }
}

void OP_Beads_CopyBeadsInListToPosList(const int listSize, const int* beadList, int (*beadPos)[POS_MAX])
{

    int i, thisBead;
    for (i = 0; i < listSize; i++)
        {
            thisBead = beadList[i];
            LatPos_copy(beadPos[i], bead_info_glb[thisBead]);
        }
}

void OP_Lattice_EmptySitesForListOfBeads(const int listSize, const int* beadList)
{
    int i, thisBead;
    for (i = 0; i < listSize; i++)
        {
            thisBead                                   = beadList[i];
            naTotLattice_glb[Lat_Ind_OfBead(thisBead)] = -1;
        }
}

void OP_Lattice_PlaceBeadsInList(const int listSize, const int* beadList)
{
    int i, thisBead;
    for (i = 0; i < listSize; i++)
        {
            thisBead                                   = beadList[i];
            naTotLattice_glb[Lat_Ind_OfBead(thisBead)] = thisBead;
        }
}

void OP_Lattice_EmptySitesForListOfPos(const int listSize, const int (*beadPos)[POS_MAX])
{
    int i;
    for (i = 0; i < listSize; i++)
        {
            naTotLattice_glb[Lat_Ind_FromVec(beadPos[i])] = -1;
        }
}

void OP_Beads_MoveBeadsInListToPos(const int listSize, const int* beadList, const int (*newPos)[POS_MAX])
{
    int i, tmpBead;
    for (i = 0; i < listSize; i++)
        {
            tmpBead = beadList[i];
            LatPos_copy(bead_info_glb[tmpBead], newPos[i]);
        }
}

void OP_Beads_BreakBondsInList(const int listSize, const int* beadList)
{
    int i, tmpID, bP;
    for (i = 0; i < listSize; i++)
        {
            tmpID = beadList[i];
            OP_Beads_BreakBond(tmpID);
        }
}

void OP_Beads_BreakBond(const int beadID)
{
    const int bP = bead_info_glb[beadID][BEAD_FACE];
    if (bP != -1)
        {
            bead_info_glb[bP][BEAD_FACE]     = -1;
            bead_info_glb[beadID][BEAD_FACE] = -1;
        }
}

void OP_Beads_RestoreBondsInList(const int listSize, const int* beadList)
{
    int i, tmpID, bP;
    for (i = 0; i < listSize; i++)
        {
            tmpID = beadList[i];
            OP_Beads_RestoreBond(tmpID);
        }
}

void OP_Beads_RestoreBond(const int beadID)
{
    const int bP = bead_info_glb[beadID][BEAD_FACE];
    if (bP != -1)
        {
            bead_info_glb[bP][BEAD_FACE] = beadID;
        }
}

void OP_Inv_MoveBeads_InList_ToPos(const int listSize, const int* beadList)
{
    int i, tmpBead;
    for (i = 0; i < listSize; i++)
        {
            tmpBead = beadList[i];
            LatPos_copy(bead_info_glb[tmpBead], old_bead_glb[tmpBead]);
        }
}

lLDub MC_RosenbluthSampling_ForLocal_AtOld(const int beadID, const int resi, long double* oldEn, const int neigh_num)
{
    int ros_num;
    if (nBeadTypeIsSticker_glb[resi])
        {
            *oldEn  = *oldEn + Energy_Anisotropic(beadID);
            ros_num = Check_RotStates_wNeighList(beadID, resi, naOldOvlpNeighs_glb, neigh_num);
            OP_NormalizeRotState(0, ros_num);
            return logl(ldaBoltzNorm_glb[0]);
        }
    else
        {
            return 0.;
        }
}

lLDub MC_RosenbluthSampling_ForLocal_AtNew(const int beadID, const int resi, int* bead_part, long double* newEn,
                                           const int neigh_num)
{
    int ros_num;

    *bead_part = -1;
    if (nBeadTypeIsSticker_glb[resi])
        {
            OP_ShuffleArray(neigh_num, naNewOvlpNeighs_glb);
            ros_num = Check_RotStates_wNeighList(beadID, resi, naNewOvlpNeighs_glb, neigh_num);
            OP_NormalizeRotState(0, ros_num);

            *bead_part = OP_PickRotState(ros_num);

            return logl(ldaBoltzNorm_glb[0]);
        }
    else
        {
            return 0.;
        }
}

lLDub MC_RosenbluthSampling_ForChains_AtOld(const int beadID, const int resi, long double* oldEn, const int neigh_num)
{
    int ros_num;
    if (nBeadTypeIsSticker_glb[resi])
        {
            *oldEn  = *oldEn + Energy_Anisotropic_For_Chain(beadID);
            ros_num = Check_RotStates_wNeighList(beadID, resi, naOldOvlpNeighs_glb, neigh_num);
            OP_NormalizeRotState(0, ros_num);
            return logl(ldaBoltzNorm_glb[0]);
        }
    else
        {
            return 0.;
        }
}

lLDub MC_RosenbluthSampling_ForChains_AtNew(const int beadID, const int resi, int* bead_part, long double* newEn,
                                            const int neigh_num)
{
    int ros_num;
    const int cur_part = bead_info_glb[beadID][BEAD_FACE];
    *bead_part         = -1;
    if (nBeadTypeIsSticker_glb[resi])
        {
            OP_ShuffleArray(neigh_num, naNewOvlpNeighs_glb);
            ros_num = Check_RotStates_wNeighList(beadID, resi, naNewOvlpNeighs_glb, neigh_num);
            OP_NormalizeRotState(0, ros_num);
            if (cur_part == -1)
                {
                    *bead_part = OP_PickRotState(ros_num);
                }
            return logl(ldaBoltzNorm_glb[0]);
        }
    else
        {
            return 0.;
        }
}

lLDub MC_RosenbluthSampling_ForLists_AtOld(const int beadIdx, const int listSize, const int beadList[MAX_BONDS + 1],
                                           lLDub* oldEn, const int neigh_num)
{
    const int beadID = beadList[beadIdx];
    const int resi   = bead_info_glb[beadID][BEAD_TYPE];
    int ros_num;
    if (nBeadTypeIsSticker_glb[resi])
        {
            *oldEn  = *oldEn + Energy_Anisotropic_For_List(beadID, listSize, beadList);
            ros_num = Check_RotStates_wNeighList(beadID, resi, naOldOvlpNeighs_glb, neigh_num);
            OP_NormalizeRotState(0, ros_num);
            return logl(ldaBoltzNorm_glb[0]);
        }
    else
        {
            return 0.;
        }
}

lLDub MC_RosenbluthSampling_ForLists_AtNew(const int beadIdx, const int listSize, const int beadList[MAX_BONDS + 1],
                                           int* bead_part, lLDub* oldEn, const int neigh_num)
{
    *bead_part         = -1;
    const int beadID   = beadList[beadIdx];
    const int resi     = bead_info_glb[beadID][BEAD_TYPE];
    const int cur_part = bead_info_glb[beadID][BEAD_FACE];
    int ros_num;
    if (nBeadTypeIsSticker_glb[resi])
        {
            *oldEn  = *oldEn + Energy_Anisotropic_For_List(beadID, listSize, beadList);
            ros_num = Check_RotStates_wNeighList(beadID, resi, naNewOvlpNeighs_glb, neigh_num);
            if (cur_part == -1)
                {
                    *bead_part = OP_PickRotState(ros_num);
                }
            OP_NormalizeRotState(0, ros_num);
            return logl(ldaBoltzNorm_glb[0]);
        }
    else
        {
            return 0.;
        }
}

lLDub MC_RosenbluthSampling_ForRange_AtOld(const int beadID, const int resi, const int smallestBead,
                                           const int largestBead, lLDub* oldEn, const int neigh_num)
{
    int ros_num;
    if (nBeadTypeIsSticker_glb[resi])
        {
            *oldEn  = *oldEn + Energy_Anisotropic_For_Range(beadID, smallestBead, largestBead);
            ros_num = Check_RotStates_wNeighList(beadID, resi, naOldOvlpNeighs_glb, neigh_num);
            OP_NormalizeRotState(0, ros_num);
            return logl(ldaBoltzNorm_glb[0]);
        }
    else
        {
            return 0.;
        }
}

lLDub MC_RosenbluthSampling_ForRange_AtNew(const int beadID, const int resi, int* bead_part, lLDub* newEn,
                                           const int neigh_num)
{
    int ros_num;
    const int cur_part = bead_info_glb[beadID][BEAD_FACE];
    *bead_part         = -1;
    if (nBeadTypeIsSticker_glb[resi])
        {
            OP_ShuffleArray(neigh_num, naNewOvlpNeighs_glb);
            ros_num = Check_RotStates_wNeighList(beadID, resi, naNewOvlpNeighs_glb, neigh_num);
            OP_NormalizeRotState(0, ros_num);
            if (cur_part == -1)
                {
                    *bead_part = OP_PickRotState(ros_num);
                }
            return logl(ldaBoltzNorm_glb[0]);
        }
    else
        {
            return 0.;
        }
}
