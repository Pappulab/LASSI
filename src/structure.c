#include "structure.h"
#include "cluster.h"
#include "energy.h"
#include "mcmove.h"

/// Lat_Ind_FromCoords - helper function to get the correct 1D index of this position
/// \param i
/// \param j
/// \param k
/// \return the 1D index location for (i,j,k) position
int Lat_Ind_FromCoords(const int i, const int j, const int k)
{ // Lattice index from 3D to 1D array
    return i + naBoxSize_glb[POS_X] * (j + naBoxSize_glb[POS_Y] * k);
}

/// Lat_Ind_FromVec - returns the 1D index given the array xArr
/// \param xArr
/// \return
int Lat_Ind_FromVec(const int* const xArr)
{ // Just the vector form of the above function. Easier to read sometimes.
    return xArr[POS_X] + naBoxSize_glb[POS_X] * (xArr[POS_Y] + naBoxSize_glb[POS_Y] * xArr[POS_Z]);
}

/// Lat_Ind_OfBead - returns the 1D index of this bead's location
/// \param beadID
/// \return
int Lat_Ind_OfBead(const int beadID)
{
    return bead_info_glb[beadID][POS_X] +
           naBoxSize_glb[POS_X] * (bead_info_glb[beadID][POS_Y] + naBoxSize_glb[POS_Y] * bead_info_glb[beadID][POS_Z]);
}

// Note that the distance functions account for periodic boundaries

/// Dist_PointTotPoint_Float - euclidean distance between the vectors (arrays) f1 and f2 where f1 and f2 are floats.
/// \param f1
/// \param f2
/// \return
float Dist_PointToPoint_Float(const float* f1, const float* f2)
{
    float d[POS_MAX];
    int i;
    for (i = 0; i < POS_MAX; i++)
        {
            d[i] = fabsf(f1[i] - f2[i]);
            d[i] = d[i] > (float) naBoxSize_glb[i] / 2. ? (float) naBoxSize_glb[i] - d[i] : d[i];
        }

    return sqrtf(d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]);
}

/// Dist_PointToPoint - euclidean distance between the vectors (arrays) f1 and f2 where f1 and f2 are integer points.
/// \param f1
/// \param f2
/// \return
float Dist_PointToPoint(const int* const restrict f1, const int* const restrict f2)
{
    int d[POS_MAX];
    int i;
    for (i = 0; i < POS_MAX; i++)
        {
            d[i] = abs(f1[i] - f2[i]);
            d[i] = d[i] > naBoxSize_glb[i] / 2 ? naBoxSize_glb[i] - d[i] : d[i];
        }

    return sqrtf((float) (d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]));
}

/// Dist_BeadToPoint - euclidean distance between beadID and the vector f1.
/// \param beadID
/// \param f1
/// \return
float Dist_BeadToPoint(const int beadID, const int* const f1)
{
    int d[POS_MAX];
    int i;
    for (i = 0; i < POS_MAX; i++)
        {
            d[i] = abs(bead_info_glb[beadID][i] - f1[i]);
            d[i] = d[i] > naBoxSize_glb[i] / 2 ? naBoxSize_glb[i] - d[i] : d[i];
        }
    return sqrtf((float) (d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]));
}

/// Dist_BeadToPoint_Double - euclidean distance between beadID and the double array f1.
/// \param beadID
/// \param f1
/// \return
float Dist_BeadToPoint_Double(const int beadID, const lDub* f1)
{
    lDub d[POS_MAX];
    int i;
    for (i = 0; i < POS_MAX; i++)
        {
            d[i] = fabs((lDub) bead_info_glb[beadID][i] - f1[i]);
            d[i] = d[i] > (lDub) naBoxSize_glb[i] / 2. ? (lDub) naBoxSize_glb[i] - d[i] : d[i];
        }
    return sqrtf((float) (d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]));
}

/// Dist_BeadToPoint_Float - euclidean distance between beadID and the float array f1.
/// \param beadID
/// \param f1
/// \return
float Dist_BeadToPoint_Float(const int beadID, const float* f1)
{
    float d[POS_MAX];
    int i;
    for (i = 0; i < POS_MAX; i++)
        {
            d[i] = fabsf((float) bead_info_glb[beadID][i] - f1[i]);
            d[i] = d[i] > (float) naBoxSize_glb[i] / 2. ? (float) naBoxSize_glb[i] - d[i] : d[i];
        }
    return sqrtf((d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]));
}

/// Dist_BeadToBead - euclidean distance between the two beads.
/// \param n1
/// \param n2
/// \return
float Dist_BeadToBead(const int n1, const int n2)
{
    lInt d[POS_MAX];
    lInt i;

    for (i = 0; i < POS_MAX; i++)
        {
            d[i] = abs(bead_info_glb[n1][i] - bead_info_glb[n2][i]);
            d[i] = d[i] > naBoxSize_glb[i] / 2 ? naBoxSize_glb[i] - d[i] : d[i];
        }

    return sqrtf((float) (d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]));
}

/// Check_System_Structure - performs a intensive and extensive sanity check. Check:
/// 1. All distances between beads in a molecule are legal.
/// 2. All bonds are symmetric.
/// 3. Bonds are only between possibly interacting beads.
/// 4. Beads' locations correspond with actual lattice location.
/// 5. Bonds are not between beads that are too far apart.
/// \return 0 if everything is okay, beadID+1 if failed.
int Check_System_Structure(void)
{
    int i, j;          // Looping variables
    int idx;           // Internal iterators for covalent bonds.
    int tmpR[POS_MAX]; // Just a vector to store coordinates
    int bondPart;
    for (i = 0; i < tot_beads_glb; i++)
        {
            idx      = 0;
            bondPart = topo_info_glb[i][idx];
            for (j = 0; j < POS_MAX; j++)
                {
                    tmpR[j] = bead_info_glb[i][j];
                }
            if (naTotLattice_glb[Lat_Ind_FromVec(tmpR)] == -1)
                {
                    printf("Lattice Position for bead %d is empty! Chain: %d\n", i, bead_info_glb[i][BEAD_CHAINID]);
                    return i + 1;
                }
            if (i - naTotLattice_glb[Lat_Ind_FromVec(tmpR)] != 0)
                { // This means there is a mismatch between where the bead is and where the lattice thinks the bead is
                    printf("Bead position and lattice value not the same. Crashing\t\t");
                    printf("B1:%d B2:%d\t C1:%d C2:%d\n", i, naTotLattice_glb[Lat_Ind_FromVec(tmpR)],
                           bead_info_glb[i][BEAD_CHAINID],
                           bead_info_glb[naTotLattice_glb[Lat_Ind_FromVec(tmpR)]][BEAD_CHAINID]);
                    return i + 1;
                }
            while (topo_info_glb[i][idx] != -1 && idx < MAX_BONDS)
                {
                    bondPart = topo_info_glb[i][idx];
                    if (Dist_BeadToBead(i, bondPart) > LINKER_RSCALE * linker_len_glb[i][idx])
                        {
                            printf("Bad beads! %d\t(%d %d %d)\t\tTopo:(%d %d %d)\t\tLinkers:(%.5f\t%.5f\t%.5f)\n", i,
                                   bead_info_glb[i][0], bead_info_glb[i][1], bead_info_glb[i][2], topo_info_glb[i][0],
                                   topo_info_glb[i][1], topo_info_glb[i][2], (float) linker_len_glb[i][0],
                                   (float) linker_len_glb[i][1], (float) linker_len_glb[i][2]);
                            printf("\t\t\t\t\t-------------------->\t\t%f\tSHOULD BE\t%f\n",
                                   Dist_BeadToBead(i, bondPart), LINKER_RSCALE * (float) linker_len_glb[i][idx]);
                            printf("Bad beads! %d\t(%d %d %d)\t\tTopo:(%d %d %d)\t\tLinkers:(%.5f\t%.5f\t%.5f)\n\n",
                                   bondPart, bead_info_glb[bondPart][0], bead_info_glb[bondPart][1],
                                   bead_info_glb[bondPart][2], topo_info_glb[bondPart][0], topo_info_glb[bondPart][1],
                                   topo_info_glb[bondPart][2], (float) linker_len_glb[bondPart][0],
                                   (float) linker_len_glb[bondPart][1], (float) linker_len_glb[bondPart][2]);
                            return i + 1;
                        }
                    idx++;
                }
            if (bead_info_glb[i][BEAD_FACE] != -1)
                {
                    if (nBeadTypeIsSticker_glb[bead_info_glb[i][BEAD_TYPE]] == 0)
                        {
                            printf("This bead -- %d -- should not have a bond. Crashing.\n", i);
                            return i + 1;
                        }
                    if (bead_info_glb[i][BEAD_FACE] == i)
                        {
                            printf("Self bonded.\n");
                            return i + 1;
                        }
                    if (i != bead_info_glb[bead_info_glb[i][BEAD_FACE]][BEAD_FACE])
                        {
                            printf("Bad bond!\n\t%d %d %d %f\nCrashing.\n", i, bead_info_glb[i][BEAD_FACE],
                                   bead_info_glb[bead_info_glb[i][BEAD_FACE]][BEAD_FACE],
                                   faEnergy_glb[bead_info_glb[i][BEAD_TYPE]]
                                               [bead_info_glb[bead_info_glb[i][BEAD_FACE]][BEAD_TYPE]][E_SC_SC]);
                            return i + 1;
                        }
                    if (Dist_BeadToBead(i, bead_info_glb[i][BEAD_FACE]) > LINKER_RSCALE)
                        {
                            printf("Bad bond! Distance is wrong\n\t%d %d %d\n Distance is %f. Crashing.\n", i,
                                   bead_info_glb[i][BEAD_FACE], bead_info_glb[bead_info_glb[i][BEAD_FACE]][BEAD_FACE],
                                   Dist_BeadToBead(i, bead_info_glb[i][BEAD_FACE]));
                            return i + 1;
                        }
                }
        }
    return 0;
}

/// Dist_Vec3n - non periodic boundary euclidean magnitude of vector
/// \param f1: The array where indicies 0,1 and 3 correspond to x y and z.
/// \return
float Dist_Vec3n(const int* f1)
{ // Outputs the magnitude of the vector
    return sqrtf((float) (f1[0] * f1[0] + f1[1] * f1[1] + f1[2] * f1[2]));
}

/// Dist_VecMagSq - non periodic boundary euclidean magnitude of vector
/// \param f1: The array where indicies 0,1 and 3 correspond to x y and z.
/// \return Square of magnitude: x^2 + y^2 + z^2
int Dist_VecMagSq(const int* f1)
{ // Outputs the magnitude of the vector
    return (f1[0] * f1[0] + f1[1] * f1[1] + f1[2] * f1[2]);
}

/// GyrTensor_ClusterSpecific - calculates the total Gyration Tensor for a given cluster
/// \param ClusSize - the total size of the cluster.
/// \param ClusIndex - the index on naClusterMatrix_glb where the cluster is stored.
/// THIS IS VERY OLD AND HASN'T BEEN LOOKED AT IN A WHILE
/// TODO: Update this for the new version
void GyrTensor_ClusterSpecific(int ClusSize, int ClusIndex)
{
    // Calculate the components of the gyration tensor for a given cluster.
    // ClusSize is the size of the cluster -- obviously -- whereas ClusIndex tell us
    // where in naClusterMatrix_glb the chain indecies are located. naClusterMatrix_glb[ClusIndex][0-ClusSize] is all
    // the chainID's I need for the calculation Remember that the Gyration Tensor is a 3x3 symmetric object so we only
    // need 6 numbers.
    int i, k, j, j2; // Basic indecies for loops
    for (i = 0; i < 7; i++)
        {
            faGyrTensor_glb[i] = 0.;
        }                          // Initializing
    int firstB, lastB;             // Tracks the first and last bead of the given chain
    float tot_COM[POS_MAX] = {0.}; // This is where we shall store the COM of the cluster.
    int NumRes             = 0;    // Tracks how many residues are in this cluster

    // The only thing one needs to be careful about is to take PBC into account; the rest is tedium.
    // Use good old differential geometry to map each coordinate to two new coordinates, and unpack at the end.
    float theta[POS_MAX] = {0.};
    float zeta[POS_MAX]  = {0.}; // Extra coordinates for goodness
    float dumArg         = 0.;   // Just a dummy variable to be more efficient
    float dumArg2        = 0.;   // Another one
    // printf("Starting with COM\n");
    // Calculating the COM
    for (i = 0; i < ClusSize; i++)
        {
            firstB = chain_info_glb[naClusterMatrix_glb[ClusIndex][i]][CHAIN_START];
            lastB  = firstB + chain_info_glb[naClusterMatrix_glb[ClusIndex][i]][CHAIN_LENGTH];
            // printf("%d %d\n", firstB, lastB);
            // Just easier to track each chain like this
            for (k = firstB; k < lastB; k++)
                {
                    NumRes++; // Adding a residue to the total
                    for (j = 0; j < POS_MAX; j++)
                        {
                            dumArg = 2. * M_PI * ((float) bead_info_glb[k][j] / (float) naBoxSize_glb[j]);
                            theta[j] += cosf(dumArg);
                            zeta[j] += sinf(dumArg); // Since I am taking the average just keep adding
                        }
                }
            // printf("Next bead\n");
        }
    // printf("Done with COM\n");
    // Calculating average, and then COM
    for (j = 0; j < POS_MAX; j++)
        {
            theta[j]   = theta[j] / (float) NumRes;
            zeta[j]    = zeta[j] / (float) NumRes;
            tot_COM[j] = atan2f(-theta[j], -zeta[j]) + M_PI;
            tot_COM[j] = naBoxSize_glb[j] * (tot_COM[j] / 2. / M_PI);
        }

    // Using the COM to calculate the Gyration Tensor
    //  GyrTen_{ij} = 1/N sum_1^N (r_i-com_i)(r_j-com_j) so just take the sums and divide at the end
    for (i = 0; i < ClusSize; i++)
        {
            firstB = chain_info_glb[naClusterMatrix_glb[ClusIndex][i]][CHAIN_START];
            lastB  = firstB + chain_info_glb[naClusterMatrix_glb[ClusIndex][i]][CHAIN_LENGTH];
            // Just easier to track each chain like this
            for (k = firstB; k < lastB; k++)
                {
                    for (j = 0; j < POS_MAX; j++)
                        {
                            for (j2 = j; j2 < POS_MAX; j2++)
                                {
                                    dumArg =
                                        fabsf((float) bead_info_glb[k][j] - tot_COM[j]) <
                                                (float) naBoxSize_glb[j] -
                                                    fabsf((float) bead_info_glb[k][j] - tot_COM[j]) ?
                                            fabsf((float) bead_info_glb[k][j] - tot_COM[j]) :
                                            (float) naBoxSize_glb[j] - fabsf((float) bead_info_glb[k][j] - tot_COM[j]);
                                    dumArg2 = fabsf((float) bead_info_glb[k][j2] - tot_COM[j2]) <
                                                      (float) naBoxSize_glb[j2] -
                                                          fabsf((float) bead_info_glb[k][j2] - tot_COM[j2]) ?
                                                  fabsf((float) bead_info_glb[k][j2] - tot_COM[j2]) :
                                                  (float) naBoxSize_glb[j2] -
                                                      fabsf((float) bead_info_glb[k][j2] - tot_COM[j2]);
                                    faGyrTensor_glb[j + 3 * (j2 - j)] += dumArg * dumArg2;
                                    // 0 = xx; 1 = yy; 2 = zz; 3 = xy; 6 = xz; 4 = yz; Need smarter indexing
                                }
                        }
                }
        }

    // Calculating the average;
    for (i = 0; i < 7; i++)
        {
            faGyrTensor_glb[i] /= (float) NumRes;
        } // printf("%f\n",faGyrTensor_glb[i]);}printf("\n");

    // exit(1);
}

/// GyrTensor_GyrRad - calculates the total Gyration Tensor of the system.
/// THIS IS VERY OLD AND HASN'T BEEN LOOKED AT IN A WHILE
/// TODO: Update this for the new version
void GyrTensor_GyrRad(void)
{                    // Calculates the gyration tensor for the whole system
    int i, k, j, j2; // Basic indecies for loops
    for (i = 0; i < 7; i++)
        {
            faGyrTensor_glb[i] = 0.f;
        }                           // Initializing
    float tot_COM[POS_MAX] = {0.f}; // This is where we shall store the COM of the cluster.

    // The only thing one needs to be careful about is to take PBC into account; the rest is tedium.
    // Use good old differential geometry to map each coordinate to two new coordinates, and unpack at the end.
    float theta[POS_MAX] = {0.f};
    float zeta[POS_MAX]  = {0.f}; // Extra coordinates for goodness
    float dumArg         = 0.f;   // Just a dummy variable to be more efficient
    float dumArg2        = 0.f;   // Another one

    for (i = 0; i < tot_beads_glb; i++)
        {
            for (j = 0; j < POS_MAX; j++)
                {
                    tot_COM[j] += bead_info_glb[i][j];
                }
        }

    // Calculating average, and then COM
    for (j = 0; j < POS_MAX; j++)
        {
            tot_COM[j] /= (float) tot_beads_glb;
        }
    // printf("\n");
    // Using the COM to calculate the Gyration Tensor
    //  GyrTen_{ij} = 1/N sum_1^N (r_i-com_i)(r_j-com_j) so just take the sums and divide at the end
    for (i = 0; i < tot_beads_glb; i++)
        {
            for (j = 0; j < POS_MAX; j++)
                {
                    dumArg = (float) bead_info_glb[i][j] - tot_COM[j];
                    for (j2 = j; j2 < POS_MAX; j2++)
                        {
                            dumArg2 = (float) bead_info_glb[i][j2] - tot_COM[j2];
                            faGyrTensor_glb[j + 3 * (j2 - j)] += dumArg * dumArg2;
                            // printf("%.2f\t", dumArg2);
                            // 0 = xx; 1 = yy; 2 = zz; 3 = xy; 6 = xz; 4 = yz; Need smarter indexing
                        }
                } // printf("\n");
        }
    // exit(1);
    // Calculating the average;
    for (i = 0; i < 7; i++)
        {
            faGyrTensor_glb[i] /= (float) tot_beads_glb;
        } // printf("%f\n",faGyrTensor_glb[i]);}printf("\n");
}

/// GyrTensor_GyrRad_Avg - calculates the total radius of gyration of the system, while not being smart about the
/// periodic boundaries. This is used as a proxy to detect phase separation, but is a relic of the old formalism. The
/// RDF should be used in general. Although this can be used without the need for a non-interacting prior.
void GyrTensor_GyrRad_Avg(void)
{
    /*
    Only calculates the diagonals of the gyration tensor, and calculates the sum of the
    diagonals. Remember that Rg^2 = Tr(GyrTen) so we only need to calculate the diagonals, and then then sum.
    I shall borrow most of the code from above, and so read GyrTensor_ClusterSpecific for what's happening here.
    */
    int i, j; // Loop indecies
    lDub tot_COM[POS_MAX] = {0.};

    Calc_SystemCenterOfMass(tot_COM);
    int dumArg[POS_MAX] = {0};
    for (i = 0; i < 7; i++)
        { // Initializing to 0
            faGyrTensor_glb[i] = 0.f;
        }

    for (i = 0; i < tot_beads_glb; i++)
        {
            for (j = 0; j < POS_MAX; j++)
                {
                    dumArg[j] = abs(bead_info_glb[i][j] - (int) tot_COM[j]);
                    dumArg[j] = dumArg[j] > naBoxSize_glb[j] / 2 ? naBoxSize_glb[j] - dumArg[j] : dumArg[j];
                    faGyrTensor_glb[j] += (float) (dumArg[j] * dumArg[j]);
                }
        }

    // Adding to the total faSysGyrRad_glb to be averaged at the end.
    faSysGyrRad_glb += sqrtf((faGyrTensor_glb[0] + faGyrTensor_glb[1] + faGyrTensor_glb[2]) / (float) tot_beads_glb);
    nTotGyrRadCounter_glb++; // Remembering that we have calculated the radius; for final averaging.
}

/// RDF_ComponentIndex - 1D index for the symmetric g_{ij} matrix given i and j.
/// \param i
/// \param j
/// \return The index of the array
/// The way it is set up, the indexing goes through the diagonal and then 0-1, 0-2, ... 0-N, 1-2, ... 1-N and so on
int RDF_ComponentIndex(const int i, const int j)
{
    if (i == j)
        {
            return i + 1;
        }
    else if (i < j)
        {
            return nBeadTypes_glb + j - (i * (3 + i - 2 * nBeadTypes_glb)) / 2;
        }
    else
        {
            return nBeadTypes_glb + i - (j * (3 + j - 2 * nBeadTypes_glb)) / 2;
        }
}

/// RDFArr_Index - 1D index for ldaTOTRDF_Arr_glb which is used to globally store the different RDFs
/// \param run_cycle
/// \param rdf_comp
/// \param x_pos
/// \return 1D index for the totalRDFArray
int RDFArr_Index(const int run_cycle, const int rdf_comp, const int x_pos)
{
    return x_pos + nRDF_TotBins_glb * (rdf_comp + nRDF_TotComps_glb * run_cycle);
}

int RadDen_ComponentIndex(const int i, const int j)
{
    if (i < 0)
        {
            return j;
        }
    else
        {
            return tot_chain_types_glb + j + tot_chain_types_glb * i;
        }
}

int RadDenArr_Index(const int run_cycle, const int rad_comp, const int x_pos)
{
    return x_pos + nRDF_TotBins_glb * (rad_comp + nRadDen_TotComps_glb * run_cycle);
}

int MolClusArr_Index(const int run_cycle, const int chain_type, const int clus_size)
{
    return clus_size + (int) tot_chains_glb * (chain_type + (int) tot_chain_types_glb * run_cycle);
}

/// RDF_ComponentWise_Avg - calculates the pair-distribution of the system where every bead acts as the center
/// of a radial histogram of pairs. Note that dr = 1/4 lattice units.
void RDF_ComponentWise_Avg(void)
{
    /*
    Calculates the RDF and adds it all up so that it can be averaged out at the end of the run.
    */
    float x; // For distance
    int i, j, k;
    int resi, resj;
    int myBin = 0;
    int array_pos;

    // Calculating where, and how many, pairs exist
    for (i = 0; i < tot_beads_glb; i++)
        {
            resi = bead_info_glb[i][BEAD_TYPE];
            for (j = i + 1; j < tot_beads_glb; j++)
                {
                    resj = bead_info_glb[j][BEAD_TYPE];
                    x    = Dist_BeadToBead(i, j);
                    // Note that Dist_BeadToBead(i,j) automatically ensures no distance is greater than (L/2)*sqrt(3)
                    myBin = (int) floor(4. * x);                      // I am assuming for now that dr=1/4
                    ldaRDF_Arr_glb[RDFArr_Index(0, 0, myBin)] += 2.0; // Adding a pair to that bin
                    array_pos = RDF_ComponentIndex(resi, resj);
                    ldaRDF_Arr_glb[RDFArr_Index(0, array_pos, myBin)] += 2.0;
                }
        }
    nTotRDFCounter_glb++;
}

/// Check_LinkerConstraint - if I move beadID to tmpR, do I still satisfy the linker lengths for beadID?
/// For the proposed location, we loop over all bonded partners of beadID and check if the distance is within the
/// linker_len_glb distance.
/// \param beadID
/// \param tmpR
/// \return 1 means all is good, 0 means bad.
int Check_LinkerConstraint(const int beadID, const int* tmpR)
{
    // Check if the proposed new location for beadID is such that all the linkers are unbroken.
    int idx;         // Iterator to loop over bond Partners
    int bondPartner; // It is what it is.
    idx         = 0;
    bondPartner = topo_info_glb[beadID][idx]; // Initializing the two.
    while (idx < MAX_BONDS && topo_info_glb[beadID][idx] != -1)
        { // Keep going till we run out of partners
            bondPartner = topo_info_glb[beadID][idx];
            if (Dist_PointToPoint(bead_info_glb[bondPartner], tmpR) >
                LINKER_RSCALE * (float) linker_len_glb[beadID][idx])
                {
                    return 0; // This means that we have broken one of the linkers.
                }
            idx++;
        }
    return 1; // This means that all linker constraints are satisfied.
}

/// Check_MTLinkerConstraint_OLD - if I move beadID and all of it's covalent bonded beads to the locations stored in
/// tmpR[][], would all the linker constraints be satisfied?
/// \param beadID
/// \param tmpR
/// \return
int Check_MTLinkerConstraint_OLD(int beadID, int (*tmpR)[POS_MAX])
{

    int curID = beadID;
    int idx, bPart;
    int j;
    int topIt = 0;
    int canI  = 1;

    while (curID != -1)
        {
            for (j = 0; j < POS_MAX; j++)
                {
                    bead_info_glb[curID][j] = tmpR[topIt][j]; // Moving
                }
            curID = topo_info_glb[beadID][topIt++];
        }

    curID = beadID;
    topIt = 0;
    while (curID != -1 && canI == 1)
        {
            idx   = 0;
            bPart = topo_info_glb[curID][idx];
            while (bPart != -1 && idx < MAX_BONDS)
                {
                    if (Dist_BeadToBead(curID, bPart) > LINKER_RSCALE * (float) linker_len_glb[curID][idx])
                        {
                            canI = 0;
                            break;
                        }
                    bPart = topo_info_glb[curID][++idx];
                }
            curID = topo_info_glb[beadID][topIt++];
        }
    curID = beadID;
    topIt = 0;
    while (curID != -1)
        {
            for (j = 0; j < POS_MAX; j++)
                {
                    bead_info_glb[curID][j] = old_bead_glb[curID][j]; // Moving back
                }
            curID = topo_info_glb[beadID][topIt++];
        }

    return canI;
}

/// Check_LinkerConstraints_ForBeadList - given that all the beads in beadList are at the new
/// positions already, as part of the MTLocal move, we iterate over every bead and check if
/// any linker constraints are broken. If so, return 0. If not, return 1 for success.
/// \param beadID
/// \param tmpR
/// \return
int Check_LinkerConstraints_ForBeadList(const int listSize, const int* beadList)
{

    int i, j, bead1, bead2;
    int dumBonds[MAX_BONDS + 1];
    int dumNum;
    float xDis;
    for (i = 0; i < listSize; i++)
        {
            bead1  = beadList[i];
            dumNum = OP_GetTopoBonds(bead1, dumBonds);
            for (j = 1; j < dumNum; j++)
                { // No need for self-distance = 0.
                    bead2 = dumBonds[j];
                    xDis  = Dist_BeadToBead(bead1, bead2);
                    if (xDis > LINKER_RSCALE * (float) linker_len_glb[bead1][j - 1])
                        {
                            return 0;
                        }
                }
        }

    return 1;
}

/// Check_BeadID_InList: Iterate over the list of beads to see if the bead is in the list.
/// If so, return 1. Otherwise return 0.
/// \param thisBeadID
/// \param listSize
/// \param beadList
/// \return
int Check_BeadID_InList(const int thisBeadID, const int listSize, const int beadList[MAX_BONDS + 1])
{
    int i;
    for (i = 0; i < listSize; i++)
        {
            if (thisBeadID == beadList[i])
                {
                    return 1;
                }
        }
    return 0;
}

///
/// \param tmpR
void Calc_SystemCenterOfMass(lDub* const tmpR)
{

    int i, j;                     // Iterators
    lDub tot_COM[POS_MAX] = {0.}; // This is where the COM will be stored

    lDub zeta[POS_MAX] = {0.};
    lDub xi[POS_MAX]   = {0.};
    lDub dumArg        = 0.;

    lDub dumConst[POS_MAX] = {0.};
    for (j = 0; j < POS_MAX; j++)
        {
            dumConst[j] = 2. * M_PI / (lDub) naBoxSize_glb[j];
        }

    for (i = 0; i < tot_beads_glb; i++)
        {
            for (j = 0; j < POS_MAX; j++)
                {
                    dumArg = dumConst[j] * (lDub) bead_info_glb[i][j];
                    zeta[j] += sin(dumArg);
                    xi[j] += cos(dumArg);
                }
        }

    int nCheck[POS_MAX] = {0};

    for (j = 0; j < POS_MAX; j++)
        {
            if ((zeta[j] == 0.) && (xi[j] == 0.))
                { // If both 0, then undefined, so skip.
                    nCheck[j] = 0;
                }
            else
                {
                    nCheck[j] = 1;
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            if (nCheck[j] == 1)
                {
                    xi[j] /= (lDub) tot_beads_glb;
                    zeta[j] /= (lDub) tot_beads_glb;
                    tot_COM[j] = atan2(-zeta[j], -xi[j]) + M_PI;
                    tot_COM[j] /= dumConst[j];
                }
            else
                {
                    tot_COM[j] = 0.;
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            tmpR[j] = tot_COM[j];
        }
}

///
/// \param naCOM_out
/// \param cluster_size
/// \param naClusList_in
void Calc_CenterOfMass_OfCluster(lDub* const naCOM_out, const int cluster_size, const int* const restrict naClusList_in)
{
    // This version measures the COM of a cluster of size c
    //  cluster size, given the molecule ID's in naList_glb.
    // The COM from this is not necessarily the COM of the system as a whole.
    int thisMol, i, j, k;         // Iterators
    int fB, lB;                   // Keep track of first and last beads of a given molecule
    lDub tot_COM[POS_MAX] = {0.}; // This is where the COM will be stored
    int bead_total_now    = 0;

    lDub zeta[POS_MAX] = {0.};
    lDub xi[POS_MAX]   = {0.};
    lDub dumArg        = 0.;

    lLDub dumConst[POS_MAX] = {0.};
    for (j = 0; j < POS_MAX; j++)
        {
            dumConst[j] = 2. * M_PI / (lLDub) naBoxSize_glb[j];
        }

    for (k = 0; k < cluster_size; k++)
        {
            thisMol = naClusList_in[k];
            fB      = chain_info_glb[thisMol][CHAIN_START];
            lB      = fB + chain_info_glb[thisMol][CHAIN_LENGTH];
            for (i = fB; i < lB; i++)
                {
                    bead_total_now++;
                    for (j = 0; j < POS_MAX; j++)
                        {
                            dumArg = dumConst[j] * (lDub) bead_info_glb[i][j];
                            zeta[j] += sin(dumArg);
                            xi[j] += cos(dumArg);
                        }
                }
        }
    int nCheck[POS_MAX] = {0};

    for (j = 0; j < POS_MAX; j++)
        {
            if ((zeta[j] == 0.) && (xi[j] == 0.))
                { // If both 0, then undefined, so skip.
                    nCheck[j] = 0;
                }
            else
                {
                    nCheck[j] = 1;
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            if (nCheck[j] == 1)
                {
                    xi[j] /= (lLDub) bead_total_now;
                    zeta[j] /= (lLDub) bead_total_now;
                    tot_COM[j] = atan2(-zeta[j], -xi[j]) + M_PI;
                    tot_COM[j] /= dumConst[j];
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            naCOM_out[j] = tot_COM[j];
        }
}

/// \param naCOM_out
/// \param cluster_size
/// \param naClusList_in
void Calc_CenterOfMass_OfChain(float* const naCOM_out, const int chainID)
{
    // This version measures the COM of a cluster of size c
    //  cluster size, given the molecule ID's in naList_glb.
    // The COM from this is not necessarily the COM of the system as a whole.
    int thisMol, i, j, k;         // Iterators
    thisMol=chainID;
    int fB, lB;                   // Keep track of first and last beads of a given molecule
    float tot_COM[POS_MAX] = {0.}; // This is where the COM will be stored
    int bead_total_now    = 0;

    lDub zeta[POS_MAX] = {0.};
    lDub xi[POS_MAX]   = {0.};
    lDub dumArg        = 0.;

    lLDub dumConst[POS_MAX] = {0.};
    for (j = 0; j < POS_MAX; j++)
        {
            dumConst[j] = 2. * M_PI / (lLDub) naBoxSize_glb[j];
        }


    fB      = chain_info_glb[thisMol][CHAIN_START];
    lB      = fB + chain_info_glb[thisMol][CHAIN_LENGTH];
    for (i = fB; i < lB; i++)
        {
            bead_total_now++;
            for (j = 0; j < POS_MAX; j++)
                {
                    dumArg = dumConst[j] * (lDub) bead_info_glb[i][j];
                    zeta[j] += sin(dumArg);
                    xi[j] += cos(dumArg);
                }
        }

    int nCheck[POS_MAX] = {0};

    for (j = 0; j < POS_MAX; j++)
        {
            if ((zeta[j] == 0.) && (xi[j] == 0.))
                { // If both 0, then undefined, so skip.
                    nCheck[j] = 0;
                }
            else
                {
                    nCheck[j] = 1;
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            if (nCheck[j] == 1)
                {
                    xi[j] /= (lLDub) bead_total_now;
                    zeta[j] /= (lLDub) bead_total_now;
                    tot_COM[j] = atan2(-zeta[j], -xi[j]) + M_PI;
                    tot_COM[j] /= dumConst[j];
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            naCOM_out[j] = tot_COM[j];
        }
}

///
/// \param tmpR
/// \param thisType
void Calc_SystemCenterOfMass_OfMolType(lDub* tmpR, const int thisType)
{
    // This version measures the COM of only type thisType
    // The COM from this is not necessarily the COM of the system as a whole.
    int thisMol, i, j, k;         // Iterators
    int fB, lB;                   // Keep track of first and last beads of a given molecule
    lDub tot_COM[POS_MAX] = {0.}; // This is where the COM will be stored
    int bead_total_now    = 0;

    lDub zeta[POS_MAX] = {0.};
    lDub xi[POS_MAX]   = {0.};
    lDub dumArg        = 0.;

    lDub dumConst[POS_MAX] = {0.};
    for (j = 0; j < POS_MAX; j++)
        {
            dumConst[j] = 2. * M_PI / (lDub) naBoxSize_glb[j];
        }

    for (k = 0; k < tot_chains_glb; k++)
        {
            thisMol = chain_info_glb[k][CHAIN_TYPE];
            if (thisMol != thisType)
                {
                    continue;
                }
            fB = chain_info_glb[k][CHAIN_START];
            lB = fB + chain_info_glb[k][CHAIN_LENGTH];
            for (i = fB; i < lB; i++)
                {
                    bead_total_now++;
                    for (j = 0; j < POS_MAX; j++)
                        {
                            dumArg = dumConst[j] * (lDub) bead_info_glb[i][j];
                            zeta[j] += sin(dumArg);
                            xi[j] += cos(dumArg);
                        }
                }
        }
    int nCheck[POS_MAX] = {0};

    for (j = 0; j < POS_MAX; j++)
        {
            if ((zeta[j] == 0.) && (xi[j] == 0.))
                { // If both 0, then undefined, so skip.
                    nCheck[j] = 0;
                }
            else
                {
                    nCheck[j] = 1;
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            if (nCheck[j] == 1)
                {
                    xi[j] /= (lDub) bead_total_now;
                    zeta[j] /= (lDub) bead_total_now;
                    tot_COM[j] = atan2(-zeta[j], -xi[j]) + M_PI;
                    tot_COM[j] /= dumConst[j];
                }
        }

    /*for (j=0; j<POS_MAX; j++){
        tmpR[j] = ((int) tot_COM[j]);
    }*/
    for (j = 0; j < POS_MAX; j++)
        {
            tmpR[j] = tot_COM[j];
        }
}

///
/// \param tmpR
/// \param thisType
void Calc_SystemCenterOfMass_WithoutMolType(lDub* tmpR, const int thisType)
{
    // This version measures the COM of everything except type thisType
    // The COM from this is not necessarily the COM of the system as a whole.
    int thisMol, i, j, k;         // Iterators
    int fB, lB;                   // Keep track of first and last beads of a given molecule
    lDub tot_COM[POS_MAX] = {0.}; // This is where the COM will be stored
    int bead_total_now    = 0;

    lDub zeta[POS_MAX] = {0.};
    lDub xi[POS_MAX]   = {0.};
    lDub dumArg        = 0.;

    lDub dumConst[POS_MAX] = {0.};
    for (j = 0; j < POS_MAX; j++)
        {
            dumConst[j] = 2. * M_PI / (lDub) naBoxSize_glb[j];
        }

    for (k = 0; k < tot_chains_glb; k++)
        {
            thisMol = chain_info_glb[k][CHAIN_TYPE];
            if (thisMol == thisType)
                {
                    continue;
                }
            fB = chain_info_glb[k][CHAIN_START];
            lB = fB + chain_info_glb[k][CHAIN_LENGTH];
            for (i = fB; i < lB; i++)
                {
                    bead_total_now++;
                    for (j = 0; j < POS_MAX; j++)
                        {
                            dumArg = dumConst[j] * (lDub) bead_info_glb[i][j];
                            zeta[j] += sin(dumArg);
                            xi[j] += cos(dumArg);
                        }
                }
        }
    int nCheck[POS_MAX] = {0};

    for (j = 0; j < POS_MAX; j++)
        {
            if ((zeta[j] == 0.) && (xi[j] == 0.))
                { // If both 0, then undefined, so skip.
                    nCheck[j] = 0;
                }
            else
                {
                    nCheck[j] = 1;
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            if (nCheck[j] == 1)
                {
                    xi[j] /= (lDub) bead_total_now;
                    zeta[j] /= (lDub) bead_total_now;
                    tot_COM[j] = atan2(-zeta[j], -xi[j]) + M_PI;
                    tot_COM[j] /= dumConst[j];
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            tmpR[j] = tot_COM[j];
        }
}

void RadDen_Avg_MolTypeWise_FromSysCen(void)
{

    int i;        // Iterator for loop
    int thisType; // Tracks the type of the chain
    lDub sysCOM[POS_MAX] = {0.};
    int myBin;
    float xDis = 0.; // Tracks the distance between the COM and the specific bead

    Calc_SystemCenterOfMass(sysCOM);
    for (i = 0; i < tot_beads_glb; i++)
        {
            thisType = bead_info_glb[i][BEAD_CHAINID];
            thisType = chain_info_glb[thisType][CHAIN_TYPE];
            xDis     = Dist_BeadToPoint_Double(i, sysCOM);
            myBin    = (int) floor(4. * xDis);
            ldaRadDen_Arr_glb[RadDenArr_Index(0, thisType, myBin)] += 1.0;
        }
    nTotRadDenCounter_glb++;
}

void RadDen_Avg_MolTypeWise_FromMolTypeCen_Old_CorrectVersion(void)
{

    int i, j;     // Iterators for loop
    int thisType; // Tracks the type of the chain
    int thisComp; // Tracks which component of ldRadDen
    lDub typeCOM[POS_MAX] = {0};
    int COM_int[POS_MAX]  = {0.};
    int myBin;
    float xDis   = 0.; // Tracks the distance between the COM and the specific bead
    int cur_type = 0;

    Calc_SystemCenterOfMass(typeCOM);
    for (j = 0; j < POS_MAX; j++)
        {
            COM_int[j] = (int) typeCOM[j];
        }
    for (i = 0; i < tot_beads_glb; i++)
        {
            thisType = bead_info_glb[i][BEAD_CHAINID];
            thisType = chain_info_glb[thisType][CHAIN_TYPE];
            thisComp = RadDen_ComponentIndex(-1, thisType);
            xDis     = Dist_BeadToPoint(i, COM_int);
            myBin    = (int) (4. * xDis);
            ldaRadDen_Arr_glb[RadDenArr_Index(0, thisComp, myBin)] += 1.0;
        }

    for (cur_type = 0; cur_type < tot_chain_types_glb; cur_type++)
        { // Go through each molType's center
            Calc_SystemCenterOfMass_OfMolType(typeCOM, cur_type);
            for (j = 0; j < POS_MAX; j++)
                {
                    COM_int[j] = (int) typeCOM[j];
                }
            for (i = 0; i < tot_beads_glb; i++)
                {
                    thisType = bead_info_glb[i][BEAD_CHAINID];
                    thisType = chain_info_glb[thisType][CHAIN_TYPE];
                    thisComp = RadDen_ComponentIndex(cur_type, thisType);
                    xDis     = Dist_BeadToPoint(i, COM_int);
                    myBin    = (int) (4. * xDis);
                    ldaRadDen_Arr_glb[RadDenArr_Index(0, thisComp, myBin)] += 1.0;
                }
        }

    /*int tmpBead;
    int cur_POS[POS_MAX] = {0.};
    int cur_DIS[POS_MAX] = {0.};
    int radRange;
    radRange = naBoxSize_glb[0];
    thisComp = RadDen_ComponentIndex(2,2);
    for(j=0; j<POS_MAX; j++){
        COM_int[j] = 0;
    }
    for (i = 0; i < radRange; i++) {
        //cur_POS[POS_X] = (COM_int[POS_X] + i) % naBoxSize_glb[POS_X];
        cur_DIS[POS_X] = i > naBoxSize_glb[0]/2 ? naBoxSize_glb[0] - i: i;
        cur_DIS[POS_X] *= cur_DIS[POS_X];
        for (j = 0; j < radRange; j++) {
            //cur_POS[POS_Y] = (COM_int[POS_Y] + j ) % naBoxSize_glb[POS_Y];
            cur_DIS[POS_Y] = j > naBoxSize_glb[0]/2 ? naBoxSize_glb[0] - j: j;
            cur_DIS[POS_Y] *= cur_DIS[POS_Y];
            for (k = 0; k < radRange; k++) {
                //cur_POS[POS_Z] = (COM_int[POS_Z] + k) % naBoxSize_glb[POS_Z];
                cur_DIS[POS_Z] = k > naBoxSize_glb[0]/2 ? naBoxSize_glb[0] - k: k;
                cur_DIS[POS_Z] *= cur_DIS[POS_Z];
                xDis = sqrtf((float)(cur_DIS[POS_X]+cur_DIS[POS_Y]+cur_DIS[POS_Z] ));
                myBin = (int) (4.*xDis);
                ldaRadDen_Arr_glb[RadDenArr_Index(0, thisComp, myBin)] += 1.0;
            }
        }
    }
    */

    nTotRadDenCounter_glb++;
}

void RadialDensityAnalysis_Perform_Analysis(void)
{

	RadDenHistUtil_ForSystem_FromLargestClusterOfMolTypes_AllChains();
    RadDenHistUtil_ForSystem_FromCenterOfMassOfMolTypes();

    nTotRadDenCounter_glb++;
}

void RadDenHistUtil_ForSystem_FromCenterOfMassOfMolTypes(void)
{

    int i, j, k;
    int naCOM_r[POS_MAX] = {0};
    int nBin;
    float fDis = 0.f; // Tracks the distance between the COM and the specific bead
    int nCurrentType;
    int thisType, thisComp;
    lDub LaTypeCOM_r[POS_MAX];
    int fB, lB;

    for (nCurrentType = 0; nCurrentType < tot_chain_types_glb + 1; nCurrentType++)
        {
            if (nCurrentType == 0)
                {
                    Calc_SystemCenterOfMass(LaTypeCOM_r);
                }
            else
                {
                    Calc_SystemCenterOfMass_OfMolType(LaTypeCOM_r, nCurrentType - 1);
                }
            for (j = 0; j < POS_MAX; j++)
                {
                    naCOM_r[j] = (int) LaTypeCOM_r[j];
                }

            for (k = 0; k < tot_chains_glb; k++)
                {
                    fB       = chain_info_glb[k][CHAIN_START];
                    lB       = fB + chain_info_glb[k][CHAIN_LENGTH];
                    thisType = chain_info_glb[k][CHAIN_TYPE];
                    thisComp = nRadDen_CompShift_glb + RadDen_ComponentIndex(nCurrentType - 1, thisType);
                    for (i = fB; i < lB; i++)
                        {
                            fDis = Dist_BeadToPoint(i, naCOM_r);
                            nBin = (int) (4.f * fDis);
                            ldaRadDen_Arr_glb[RadDenArr_Index(0, thisComp, nBin)] += 1.0;
                        }
                }
        }
}

void RadDenHistUtil_ForSystem_FromLargestClusterOfMolTypes(void)
{

    int i, j, k;        // Iterators for loop
    int thisType;       // Tracks the type of the chain
    int nRadDenMolComp; // Tracks which component of ldRadDen
    int thisChain;
    lDub LaTypeCOM_r[POS_MAX] = {0.};
    int naCOM_r[POS_MAX]      = {0};
    int nRadBin;
    float fDis       = 0.f; // Tracks the distance between the COM and the specific bead
    int nCurrentType = 0;
    int fB, lB;
    int nMolWiseClusSize = 0;

    int* const naMolWiseClusIDs   = malloc((tot_chain_types_glb + 2) * sizeof(int));
    int* const naTmpClusChainList = malloc((tot_chains_glb + 1) * sizeof(int));
    int* const naFullClusList     = malloc((tot_chains_glb + 1) * sizeof(int));
    int* const naFullClusSizes    = calloc((tot_chains_glb + 1), sizeof(int));
    int* const naFullClusCumSizes = calloc((tot_chains_glb + 1), sizeof(int));

    const int nMolWiseClusNum = ClusUtil_OfSystem_MolWise_GetLargestClusters(naMolWiseClusIDs, naFullClusList,
                                                                             naFullClusCumSizes, naFullClusSizes);

    for (nCurrentType = 0; nCurrentType < tot_chain_types_glb + 1; nCurrentType++)
        {

            nMolWiseClusSize = ClusUtil_GetCluster_FromFullClusAndCumSizes(
                naMolWiseClusIDs[nCurrentType], naTmpClusChainList, naFullClusList, naFullClusCumSizes, naFullClusSizes,
                nMolWiseClusNum);

            Calc_CenterOfMass_OfCluster(LaTypeCOM_r, nMolWiseClusSize, naTmpClusChainList);

            for (j = 0; j < POS_MAX; j++)
                {
                    naCOM_r[j] = (int) LaTypeCOM_r[j];
                }

            for (k = 0; k < nMolWiseClusSize; k++)
                {
                    thisChain      = naTmpClusChainList[k];
                    fB             = chain_info_glb[thisChain][CHAIN_START];
                    lB             = fB + chain_info_glb[thisChain][CHAIN_LENGTH];
                    thisType       = chain_info_glb[thisChain][CHAIN_TYPE];
                    nRadDenMolComp = RadDen_ComponentIndex(nCurrentType - 1, thisType);
                    for (i = fB; i < lB; i++)
                        {
                            fDis    = Dist_BeadToPoint(i, naCOM_r);
                            nRadBin = (int) (4.f * fDis);
                            ldaRadDen_Arr_glb[RadDenArr_Index(0, nRadDenMolComp, nRadBin)] += 1.0;
                        }
                }
        }

    free(naTmpClusChainList);
    free(naFullClusSizes);
    free(naFullClusList);
    free(naFullClusCumSizes);
    free(naMolWiseClusIDs);
}

void RadDenHistUtil_ForSystem_FromLargestClusterOfMolTypes_AllChains(void)
{

    int i, j, k;        // Iterators for loop
    int thisType;       // Tracks the type of the chain
    int nRadDenMolComp; // Tracks which component of ldRadDen
    lDub LaTypeCOM_r[POS_MAX] = {0.};
    int naCOM_r[POS_MAX]      = {0};
    int nRadBin;
    float fDis       = 0.f; // Tracks the distance between the COM and the specific bead
    int nCurrentType = 0;
    int fB, lB;
    int nMolWiseClusSize = 0;

    int* const naMolWiseClusIDs   = malloc((tot_chain_types_glb + 2) * sizeof(int));
    int* const naTmpClusChainList = malloc((tot_chains_glb + 1) * sizeof(int));
    int* const naFullClusList     = malloc((tot_chains_glb + 1) * sizeof(int));
    int* const naFullClusSizes    = calloc((tot_chains_glb + 1), sizeof(int));
    int* const naFullClusCumSizes = calloc((tot_chains_glb + 1), sizeof(int));

    const int nMolWiseClusNum = ClusUtil_OfSystem_MolWise_GetLargestClusters(naMolWiseClusIDs, naFullClusList,
                                                                             naFullClusSizes, naFullClusCumSizes);

    for (nCurrentType = 0; nCurrentType < tot_chain_types_glb + 1; nCurrentType++)
        {

            nMolWiseClusSize = ClusUtil_GetCluster_FromFullClusAndCumSizes(
                naMolWiseClusIDs[nCurrentType], naTmpClusChainList, naFullClusList, naFullClusCumSizes, naFullClusSizes,
                nMolWiseClusNum);

            Calc_CenterOfMass_OfCluster(LaTypeCOM_r, nMolWiseClusSize, naTmpClusChainList);

            for (j = 0; j < POS_MAX; j++)
                {
                    naCOM_r[j] = (int) LaTypeCOM_r[j];
                }

            for (k = 0; k < tot_chains_glb; k++)
                {
                    fB       = chain_info_glb[k][CHAIN_START];
                    lB       = fB + chain_info_glb[k][CHAIN_LENGTH];
                    thisType = chain_info_glb[k][CHAIN_TYPE];
                    nRadDenMolComp = RadDen_ComponentIndex(nCurrentType - 1, thisType);
                    for (i = fB; i < lB; i++)
                        {
                    	    fDis    = Dist_BeadToPoint(i, naCOM_r);
                    	    nRadBin = (int) (4.f * fDis);
                    	    ldaRadDen_Arr_glb[RadDenArr_Index(0, nRadDenMolComp, nRadBin)] += 1.0;
                        }
                }
        }

    free(naTmpClusChainList);
    free(naFullClusSizes);
    free(naFullClusList);
    free(naFullClusCumSizes);
    free(naMolWiseClusIDs);
}

void Calculate_Distances_For_Radius(float* thisList, const int nRad)
{
    int r_disp[POS_MAX] = {0};

    int dum_iter = 0;

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
        {
            for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
                {
                    for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
                        {
                            thisList[dum_iter++] = Dist_Vec3n(r_disp);
                        }
                }
        }
    //    printf("%d", dum_iter);
    //    exit(1);
}

/// NeighborSearch_ForOvlp: Given the beadID, which should be at startVec, we calculate all the neighbors
/// within the +-1 cube. We return the number of neighbors.
/// \param beadID
/// \param startVec
/// \param neighList
/// \return Number of neighbors.
int NeighborSearch_ForOvlp(int const beadID, const int* startVec, int* neighList)
{
    int neigh_num = 0;

    int tmpBead;
    int r_disp[POS_MAX] = {0};
    int r_chck[POS_MAX] = {0};
    const int nRad      = 1;

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
        {
            LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_X);
            for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
                {
                    LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_Y);
                    for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
                        {
                            LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_Z);
                            tmpBead = naTotLattice_glb[Lat_Ind_FromVec(r_chck)];
                            if (tmpBead != -1 && tmpBead != beadID)
                                {
                                    neighList[neigh_num++] = tmpBead;
                                }
                        }
                }
        }
    neighList[neigh_num] = -1;
    return neigh_num;
}

/// NeighborSearch_ForOvlp: Given the beadID, which should be at startVec, we calculate all the neighbors
/// within the +-1 cube. We return the number of neighbors.
/// \param beadID
/// \param startVec
/// \param neighList
/// \return Number of neighbors.
int NeighborSearch_ForCont(int const beadID, const int* startVec, int* contList, int* ovlpList, int* ovlpNum)
{
    int neigh_num = 0;

    int tmpBead;
    int r_disp[POS_MAX] = {0};
    int r_chck[POS_MAX] = {0};
    const int nRad      = LARGEST_RADIUS;

    *ovlpNum = 0;

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
        {
            LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_X);
            for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
                {
                    LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_Y);
                    for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
                        {
                            LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_Z);
                            tmpBead = naTotLattice_glb[Lat_Ind_FromVec(r_chck)];
                            if (tmpBead != -1 && tmpBead != beadID)
                                {
                                    if ((abs(r_disp[0]) <= 1) && (abs(r_disp[1]) <= 1) && (abs(r_disp[2]) <= 1))
                                        {
                                            ovlpList[*ovlpNum] = tmpBead;
                                            *ovlpNum           = *ovlpNum + 1;
                                        }
                                    contList[neigh_num++] = tmpBead;
                                }
                        }
                }
        }
    contList[neigh_num] = -1;
    ovlpList[*ovlpNum]  = -1;
    return neigh_num;
}

/// NeighborSearch_EmptySitesAround_wRad: Given the starting point, and radius, we calculate how many sites
/// around this point are empty.
/// \param startVec
/// \param nRad
/// \return Number of empty lattice-sites around startVec
int NeighborSearch_EmptySitesAround_wRad(const int* startVec, const int nRad)
{
    int empty_num = 0;
    int tmpBead;
    int r_disp[POS_MAX]   = {0};
    int r_search[POS_MAX] = {0};

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
        {
            for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
                {
                    for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
                        {
                            LatPos_add_wPBC(r_search, startVec, r_disp);
                            tmpBead = naTotLattice_glb[Lat_Ind_FromVec(r_search)];
                            if (tmpBead == -1)
                                {
                                    empty_num++;
                                }
                        }
                }
        }
    return empty_num;
}

/// NeighborSearch_SolSitesAround: number of solvation sites around this point. Finds the number of empty lattice
/// sites in the 26 closest sites.
/// \param startVec
/// \return Empty sites around this point.
int NeighborSearch_SolSitesAround(const int* startVec)
{
    return NeighborSearch_EmptySitesAround_wRad(startVec, 1);
}

int NeighborSearch_AroundPoint_UptoIndex(const int beadID, const int* startVec, const int upTO, int* neighList)
{
    int i;
    int neigh_num = 0;
    int tmpBead;

    int r_search[POS_MAX] = {0};

    int r_all[100][POS_MAX] = {0};

    for (i = 0; i < upTO; ++i)
        {
            LatPos_add_wPBC(r_search, startVec, r_all[i]);
            tmpBead = naTotLattice_glb[Lat_Ind_FromVec(r_search)];
            if (tmpBead != -1)
                {
                    neighList[neigh_num++] = tmpBead;
                }
        }

    return neigh_num;
}

int NeighborSearch_AroundPoint_wRad_IgnBead(const int beadID, const int* startVec, const int nRad, int* neighList)
{

    int neigh_num = 0;
    int tmpBead;
    int r_disp[POS_MAX]   = {0};
    int r_search[POS_MAX] = {0};

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
        {
            for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
                {
                    for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
                        {
                            LatPos_add_wPBC(r_search, startVec, r_disp);
                            tmpBead = naTotLattice_glb[Lat_Ind_FromVec(r_search)];
                            if (tmpBead != -1 && tmpBead != beadID)
                                {
                                    neighList[neigh_num++] = tmpBead;
                                }
                        }
                }
        }
    return neigh_num;
}

int NeighborSearch_AroundPoint_wRad_wDists(const int beadID, const int* startVec, const int nRad, int* neighList,
                                           float* distList)
{

    int neigh_num = 0;
    int tmpBead;
    int r_disp[POS_MAX] = {0};

    int r_search[POS_MAX] = {0};

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
        {
            for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
                {
                    for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
                        {
                            LatPos_add_wPBC(r_search, startVec, r_disp);
                            tmpBead = naTotLattice_glb[Lat_Ind_FromVec(r_search)];
                            if (tmpBead != -1 && tmpBead != beadID)
                                {
                                    neighList[neigh_num] = tmpBead;
                                    distList[neigh_num]  = Dist_Vec3n(r_disp);
                                    neigh_num++;
                                }
                        }
                }
        }
    return neigh_num;
}

/// NeighborSearch_ForOvlp: Given the beadID, which should be at startVec, we calculate all the neighbors
/// within the +-1 cube. We return the number of neighbors.
/// \param beadID
/// \param startVec
/// \param neighList
/// \return Number of neighbors.
int NeighborSearch_ForCluster_Ovlp(int const beadID, const int* startVec, int* neighList)
{
    int neigh_num = 0;

    int tmpBead;
    int r_disp[POS_MAX] = {0};
    int r_chck[POS_MAX] = {0};
    const int nRad      = 1;

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
        {
            LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_X);
            for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
                {
                    LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_Y);
                    for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
                        {
                            LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_Z);
                            tmpBead = naTotLattice_glb[Lat_Ind_FromVec(r_chck)];
                            if (tmpBead != -1)
                                {
                                    neighList[neigh_num] = tmpBead;
                                    neigh_num++;
                                }
                        }
                }
        }
    neighList[neigh_num] = -1;
    return neigh_num;
}

/// LatPos_copy: Copies the position elements to copy_vec.
/// \param copy_vec
/// \param in_vec
void LatPos_copy(int* copy_vec, const int* in_vec)
{
    copy_vec[0] = in_vec[0];
    copy_vec[1] = in_vec[1];
    copy_vec[2] = in_vec[2];
}

/// LatPos_gen_rand_wRad. Given a radius, we generate a random set of coordinates between -Radius and Radius
/// for each dimension. The result is stored in outVec
/// \param outVec
/// \param nRadius
void LatPos_gen_rand_wRad(int outVec[POS_MAX], const int nRadius)
{
    int j;
    int radUp = 2 * nRadius + 1;
    for (j = 0; j < POS_MAX; j++)
        {
            outVec[j] = rand() % radUp;
        }
    for (j = 0; j < POS_MAX; j++)
        {
            outVec[j] -= nRadius;
        }
}

/// LatPos_add_wPBC: Given firVec and secVec, we add the two vectors, and take care of PBCs, storing in outVec.
/// The implementation assumes that the sum of any coordinate can not be larger than twice the lattice-size,
/// because there is no explicit modulo operation.
/// \param outVec
/// \param firVec
/// \param secVec
void LatPos_add_wPBC(int* outVec, const int* firVec, const int* secVec)
{
    short j;
    for (j = 0; j < POS_MAX; j++)
        {
            outVec[j] = firVec[j] + secVec[j];
        }
    for (j = 0; j < POS_MAX; j++)
        {
            outVec[j] = outVec[j] < 0 ? outVec[j] + naBoxSize_glb[j] : outVec[j];
        }
    for (j = 0; j < POS_MAX; j++)
        {
            outVec[j] = outVec[j] >= naBoxSize_glb[j] ? outVec[j] - naBoxSize_glb[j] : outVec[j];
        }
}

/// LatPos_add_wPBC_ofComp: Given firVec and secVec, we add the compNum component of the two vectors,
/// and take care of PBCs, storing in outVec.
/// The implementation assumes that the sum of any coordinate can not be larger than twice the lattice-size,
/// because there is no explicit modulo operation.
/// \param outVec
/// \param firVec
/// \param secVec
inline void LatPos_add_wPBC_ofComp(int* outVec, const int* firVec, const int* secVec, const int compNum)
{
    outVec[compNum] = firVec[compNum] + secVec[compNum];
    outVec[compNum] = outVec[compNum] < 0 ? outVec[compNum] + naBoxSize_glb[compNum] : outVec[compNum];
    outVec[compNum] =
        outVec[compNum] >= naBoxSize_glb[compNum] ? outVec[compNum] - naBoxSize_glb[compNum] : outVec[compNum];
}

/// LatPos_add_noPBC: Given the two arrays firVec and secVec, we store the sum of the two arrays in outVec.
/// \param outVec
/// \param firVec
/// \param secVec
void LatPos_add_noPBC(int* outVec, const int* firVec, const int* secVec)
{
    short j;
    for (j = 0; j < POS_MAX; j++)
        {
            outVec[j] = firVec[j] + secVec[j];
        }
}

/// BeadPos_sub_wPBC: Given firVec and secVec, we subtract second from first, and take care of PBCs, storing in outVec.
/// The implementation assumes that the sum of any coordinate can not be larger than the lattice size
/// because there is no explicit modulo operation.
/// This version is specifically for bead-positions, or for inter-bead calculations.
/// \param outVec
/// \param firVec
/// \param secVec
void BeadPos_sub_wPBC(int* outVec, const int* firVec, const int* secVec)
{
    short j;
    for (j = 0; j < POS_MAX; j++)
        {
            outVec[j] = firVec[j] - secVec[j];
        }
    for (j = 0; j < POS_MAX; j++)
        {
            outVec[j] = outVec[j] < -naBoxSize_glb[j] / 2 ? outVec[j] + naBoxSize_glb[j] : outVec[j];
        }
    for (j = 0; j < POS_MAX; j++)
        {
            outVec[j] = outVec[j] >= naBoxSize_glb[j] / 2 ? outVec[j] - naBoxSize_glb[j] : outVec[j];
        }
}

/// OP_GetTopoBonds - including beadID, get all the beads that are covalently bonded to beadID.
/// \param beadID
/// \param dum_list
/// \return Number_of_bonds+1. Can be used to loop over list which includes beadID.
int OP_GetTopoBonds(const int beadID, int* dum_list)
{
    int top_it = 0;
    int curID  = beadID;

    while (curID != -1)
        {
            dum_list[top_it] = curID;
            curID            = topo_info_glb[beadID][top_it++];
        }

    return top_it;
}

int Vec3n_DotProd(const int* vec_1, const int* vec_2)
{
    return vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2];
}

/// Vec3n_CosTheta - calculates the Cos(theta) between the vectors v1 and v2.
/// Cos(theta) = v_1 . v_2 / |v_1| / |v_2|
/// \param v1
/// \param v2
/// \return
float Vec3n_CosTheta(const int* v1, const int* v2)
{
    const int dotProd       = Vec3n_DotProd(v1, v2);
    const int vec1MagSq     = Vec3n_DotProd(v1, v1);
    const int vec2MagSq     = Vec3n_DotProd(v2, v2);
    const float dumDenom    = (float) vec1MagSq * (float) vec2MagSq;
    const float dumCosTheta = (float) dotProd / sqrtf(dumDenom);

    return dumCosTheta;
}

float Vec3n_AngleBetweenVecs(const int* vec_1, const int* vec_2)
{

    const float dumAng = Vec3n_CosTheta(vec_1, vec_2);
    return acosf(dumAng);
}

float BeadPos_CosThetaOfBeads(const int bead1, const int bead2)
{

    const int r_pos1[POS_MAX] = {bead_info_glb[bead1][0], bead_info_glb[bead1][1], bead_info_glb[bead1][2]};

    const int r_pos2[POS_MAX] = {bead_info_glb[bead2][0], bead_info_glb[bead2][1], bead_info_glb[bead2][2]};

    return Vec3n_CosTheta(r_pos1, r_pos2);
}

/// BeadIDListOP_GetChainIDs: Given beadList of size beadNum, we record the chainIDs of every bead into chainList.
/// Make sure chainList is as big as beadList.
/// \param beadNum
/// \param beadList
/// \param chainList
void BeadListOP_GetChainIDs(const int beadNum, const int* restrict beadList, int* restrict chainList)
{
    int i, tmp_bead, dum_chain;

    for (i = 0; i < beadNum; ++i)
        {
            tmp_bead     = beadList[i];
            dum_chain    = bead_info_glb[tmp_bead][BEAD_CHAINID];
            chainList[i] = dum_chain;
        }

    chainList[i] = -1;
}

/// BeadIDListOP_GetChainTypes: Given beadList of size beadNum, we record the chain_types of every bead into
/// chainList.
/// Make sure chainList is as big as beadList.
/// \param beadNum
/// \param beadList
/// \param chainList
void BeadListOP_GetChainTypes(const int beadNum, const int* beadList, int* chainList)
{
    int i, tmp_bead, dum_chain;

    for (i = 0; i < beadNum; ++i)
        {
            tmp_bead     = beadList[i];
            dum_chain    = bead_info_glb[tmp_bead][BEAD_CHAINID];
            chainList[i] = chain_info_glb[dum_chain][CHAIN_TYPE];
        }

    chainList[i] = -1;
}

/// BeadListOP_Filter_wrt_SecondList. We overwrite beadList such that whichever index propList[i] == prop_val, we
/// add to beadList, and ignore all others. Should act like a masking operation over beadList, where all masked
/// elements are removed.
/// \param beadNum
/// \param beadList
/// \param propList
/// \param prop_val
/// \return newSize - size of overwritten beadList.
int BeadListOP_Filter_wrt_SecondList(const int beadNum, int* beadList, const int* propList, const int prop_val)
{
    int newSize = 0;
    int i;
    for (i = 0; i < beadNum; i++)
        {
            if (propList[i] == prop_val)
                {
                    beadList[newSize] = beadList[i];
                    newSize++;
                }
        }
    beadList[newSize] = -1;
    return newSize;
}

/// BeadListOP_InvFilter_wrt_SecondList. We overwrite beadList such that whichever index propList[i] != prop_val, we
/// add to beadList, and ignore all others. Should act like a masking operation over beadList, where all masked
/// elements are removed.
/// Acts as the logical not for Filter_wrt_SecondList
/// \param beadNum
/// \param beadList
/// \param propList
/// \param prop_val
/// \return newSize - size of overwritten beadList.
int BeadListOP_InvFilter_wrt_SecondList(const int beadNum, int* beadList, const int* propList, const int prop_val)
{
    int newSize = 0;
    int i;
    for (i = 0; i < beadNum; i++)
        {
            if (propList[i] != prop_val)
                {
                    beadList[newSize] = beadList[i];
                    newSize++;
                }
        }
    beadList[newSize] = -1;
    return newSize;
}

/// BeadListOP_GetIntraChainID. For the beadList of size beadNum, we store the intra-chain ID of each bead.
/// tmpBeadID - chainStartID
/// \param beadNum
/// \param beadList
/// \param chainList
void BeadListOP_GetIntraChainID(const int beadNum, const int* beadList, int* chainList)
{
    int i, tmp_bead, dum_chain, dum_len;

    for (i = 0; i < beadNum; ++i)
        {
            tmp_bead     = beadList[i];
            dum_chain    = bead_info_glb[tmp_bead][BEAD_CHAINID];
            dum_len      = chain_info_glb[dum_chain][CHAIN_START];
            chainList[i] = tmp_bead - dum_len;
        }
    chainList[i] = -1;
}

/// BeadListOP_Filter_DbPvtLinkerConFwd. For every bead (tmpBead) in beadList, we calculate the
/// distance between thisBead and tmpBead-1. If the distance satisfies linker constraints,
/// we add into beadList.
/// This is for DbPvt move specifically. d(i+1, i') < linker[thisBead][0]
/// \param beadNum
/// \param beadList
/// \param thisBead
int BeadListOP_Filter_DbPvtLinkerConFwd(const int beadNum, int* beadList, const int thisBead)
{
    int newSize = 0;
    int i, tmpBead;
    const float linker_cons = (float) linker_len_glb[thisBead][0] * LINKER_RSCALE;

    float xDis;

    for (i = 0; i < beadNum; ++i)
        {
            tmpBead = beadList[i];
            xDis    = Dist_BeadToBead(thisBead, tmpBead - 1);

            if (xDis < linker_cons)
                {
                    beadList[newSize] = tmpBead;
                    newSize++;
                }
        }

    return newSize;
}

/// BeadListOP_Filter_DbPvtLinkerConBck. For every bead (tmpBead) in beadList, we calculate the
/// distance between thisBead and tmpBead-1. If the distance satisfies linker constraints,
/// we add into beadList.
/// This is for DbPvt move specifically. d(i+1, i') < linker[thisBead][0]
/// \param beadNum
/// \param beadList
/// \param thisBead
int BeadListOP_Filter_DbPvtLinkerConBck(const int beadNum, int* beadList, const int thisBead)
{
    int newSize = 0;
    int i, tmpBead;
    const float linker_cons = (float) linker_len_glb[thisBead][1] * LINKER_RSCALE;

    float xDis = 0.f;

    for (i = 0; i < beadNum; ++i)
        {
            tmpBead = beadList[i];
            xDis    = Dist_BeadToBead(thisBead, tmpBead);

            if (xDis < linker_cons)
                {
                    beadList[newSize] = tmpBead;
                    newSize++;
                }
        }

    return newSize;
}

/// ListOP_UniqueElementsOfSortedList_Int: Overwrite the provided _sorted_ list of integers with only the unique
/// elements, and return the size of the new array. The overwritten list is also sorted. The implementation only works
/// on sorted lists. \param dum_list: Sorted list of beadID's (or integers) \return number of unique elements.
int ListOP_UniqueElementsOfSortedList_Int(const int size, int* sorted_list)
{
    int* first  = sorted_list;
    int* second = first;
    int newSize = 1;
    while (++second != (sorted_list + size))
        {
            if (*second != *first)
                {
                    *(++first) = *second;
                    newSize++;
                }
        }
    return newSize;
}

/// BeadList_AppendBeads
/// Note that old_list should be at least as big as old_size + app_size. No bounds checking is done since
/// this dumb implementation has no idea about the size of the array.
/// \param old_size Number of elements in old_list.
/// \param old_list
/// \param app_list
/// \param app_size Number of elements of app_list to be appended to old_list.
/// \return old_size + app_size. Size of the new array.
int BeadList_AppendBeads(const int old_size, int* old_list, const int* app_list, const int app_size)
{
    int i;
    for (i = 0; i < app_size; i++)
        {
            old_list[i + old_size] = app_list[i];
        }

    return old_size + app_size;
}

/// UtilFunc_CompareInts: A simple function that is used in CSTDLIB qsort for integers.
/// \param a
/// \param b
/// \return
int UtilFunc_CompareInts(const void* a, const void* b)
{
    return (*(int*) a - *(int*) b);
}

/// BeadList_CanTopoAngle: For a list of beadIDs, we only keep the ones that have TopoAngle energies.
/// \param size
/// \param beadList
/// \return An overwritten list containing only beads that have a TopoAngle energy
int BeadList_CanTopoAngle(const int size, int* beadList)
{

    int i, tmpBead, resi;
    int newSize = 0;
    for (i = 0; i < size; i++)
        {
            tmpBead = beadList[i];
            resi    = bead_info_glb[tmpBead][BEAD_TYPE];
            if (faEnergy_glb[resi][resi][E_STIFF])
                {
                    beadList[newSize++] = tmpBead;
                }
        }

    return newSize;
}

/// ChainListOP_AddUniqueChains_wSize: Add any chains that are not in clusList from chainList. This function assumes
/// that chainList contains only unique elements, and that clusList has the smallest chainID as the first element, and
/// the largest chainID the last element (at clusList[clusSize-1]). \param clusSize \param clusList \param newChainNums
/// \param chainList
int ChainListOP_AddUniqueChains_wSize(const int clusSize, int* clusList, const int newChainNums, const int* chainList,
                                      const int maxClus)
{
    int i;
    int curChain;
    int newClusSize = clusSize;

    for (i = 0; i < newChainNums; i++)
        {
            curChain    = chainList[i];
            newClusSize = ClusListOP_AddChainID(newClusSize, clusList, curChain);
            if (newClusSize >= maxClus)
                {
                    return newClusSize;
                }
        }

    return newClusSize;
}

/// ChainListOP_AddUniqueChains_wSize: Add any chains that are not in clusList from chainList. This function assumes
/// that chainList contains only unique elements, and that clusList has the smallest chainID as the first element, and
/// the largest chainID the last element (at clusList[clusSize-1]).
/// \param clusSize
/// \param clusList
/// \param newChainNums
/// \param chainList
int ChainListOP_AddUniqueChains_wSize_Check(const int clusSize, int* clusList, const int newChainNums,
                                            const int* chainList, const int maxClus, const int* oldList)
{
    int i;
    int curChain;
    int newClusSize = clusSize;

    for (i = 0; i < newChainNums; i++)
        {
            curChain    = chainList[i];
            newClusSize = ClusListOP_AddChainID(newClusSize, clusList, curChain);
            if (newClusSize >= maxClus)
                {
                    return newClusSize;
                }
            if (clusList[newClusSize - 1] != oldList[newClusSize - 1])
                {
                    return maxClus;
                }
        }

    return newClusSize;
}

/// ChainListOP_GetChainTypes - Given a list that contains chainIDs, we get the type of each of the chains.
/// \param naTypesList Array where the chain-types will be recorded.
/// \param naChainIDs Array that contains all the chainIDs.
/// \param nListSize Number of chains in the naChainIDs array
void ChainListOP_GetChainTypes(int* const naTypesList, const int* const naChainIDs, const int nListSize)
{
    int i;
    int tmpChain;
    for (i = 0; i < nListSize; i++)
        {
            tmpChain       = naChainIDs[i];
            naTypesList[i] = chain_info_glb[tmpChain][CHAIN_TYPE];
        }
}

/// ClusListOP_AddChainID - Given a cluster-chainIDs list containing clusSize chains, we try to add the proposed
/// chainID. Assumes that clusList[0] is the smallest, and clusList[clusSize-1] is the largest.
/// - If the chain is smaller than the smallest, we add it to become the first element of the chain, because it is now
/// the smallest chainID. The largest chainID moves up one, and the old smallest is placed behind the largest.
/// - Else, if the chain is larger than the largest, we add it infront of the old largest.
/// - Else, we loop over the clusList to see if the proposed chainID is unique.
/// \param clusSize
/// \param clusList
/// \param chainID
/// \return
int ClusListOP_AddChainID(const int clusSize, int* clusList, const int chainID)
{

    int i;
    for (i = 0; i < clusSize; i++)
        {
            if (clusList[i] == chainID)
                {
                    return clusSize;
                }
        }

    clusList[clusSize] = chainID;

    return clusSize + 1;
}

/// ClusListOP_AddIfUniqueChainID - Loop over the clusList. If none of the chains match, we insert the chainID behind
/// last chain in clusList. So we move the last chain up by one.
/// \param clusSize
/// \param clusList
/// \param chainID
/// \return clusSize+1 if the chain was added, clusSize if not.
int ClusListOP_AddIfUniqueChainID(const int clusSize, int* clusList, const int chainID)
{
    int i;
    for (i = 0; i < clusSize; i++)
        {
            if (clusList[i] == chainID)
                {
                    return clusSize;
                }
        }
    clusList[clusSize]     = clusList[clusSize - 1];
    clusList[clusSize - 1] = chainID;
    return clusSize + 1;
}

/// LatticeUtil_GetOvlpLattIndecies: Given the starting location r_pos0, we store all of the naTotLattice_glb indecies
/// corresponding to the +-1 neighbors. We store all the indecies in numList. NOTE: This list _will_ include the
/// point at r_pos0 as well.
/// \param r_pos0
/// \param numList - MUST be at least 27 = 3^3 elements long. Undefined behavior if not that long.
void LatticeUtil_GetOvlpLattIndecies(const int* restrict const r_pos0, int* restrict numList)
{
    int i               = 0;
    int r_disp[POS_MAX] = {0};
    int r_chck[POS_MAX] = {0};
    const int nRad      = 1;
    int tmpBead;

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
        {
            LatPos_add_wPBC_ofComp(r_chck, r_pos0, r_disp, POS_X);
            for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
                {
                    LatPos_add_wPBC_ofComp(r_chck, r_pos0, r_disp, POS_Y);
                    for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
                        {
                            LatPos_add_wPBC_ofComp(r_chck, r_pos0, r_disp, POS_Z);
                            tmpBead      = Lat_Ind_FromVec(r_chck);
                            numList[i++] = tmpBead;
                        }
                }
        }
}

/// LatticeUtil_GetLattVals_FromList - loop over nIndexList and store the naTotLattice_glb values in nLatValsList.
/// Note that this will include -1 values as well, or empty sites.
/// \param nIndexList
/// \param nLatValsList - Must be at least 27 elements long.
void LatticeUtil_GetLattVals_FromList(const int* restrict const nIndexList, int* restrict nLatValsList,
                                      const int listSize)
{
    int i;
    for (i = 0; i < listSize; ++i)
        {
            nLatValsList[i] = naTotLattice_glb[nIndexList[i]];
        }
}

/// LatticeUtil_GetBeadIDs_FromList - Given the list of Lattice values, we get all the values that are not -1.
/// Or, the we get the values that correspond to non-empty lattice sites.
/// \param nLatValsList
/// \param nBeadsList
/// \return
int LatticeUtil_GetBeadIDs_FromList(const int* restrict const nLatValsList, int* restrict nBeadsList,
                                    const int nListSize)
{
    int nBeadsNum = 0;
    int tmpBead;
    int i;

    for (i = 0; i < nListSize; ++i)
        {
            tmpBead = nLatValsList[i];
            if (tmpBead != -1)
                {
                    nBeadsList[nBeadsNum++] = tmpBead;
                }
        }
    return nBeadsNum;
}

/// LatticeUtil_GetNeighBeads_AtPos - given the starting position
/// \param r_pos0
/// \param nBeadsList
/// \return
int LatticeUtil_GetNeighBeads_AtPos(const int* restrict const r_pos0, int* restrict nBeadsList)
{

    int nLatInd[CLUS_CONTACT_NEIGHS];
    LatticeUtil_GetOvlpLattIndecies(r_pos0, nLatInd);

    int nLatVals[CLUS_CONTACT_NEIGHS];
    LatticeUtil_GetLattVals_FromList(nLatInd, nLatVals, CLUS_CONTACT_NEIGHS);

    return LatticeUtil_GetBeadIDs_FromList(nLatVals, nBeadsList, CLUS_CONTACT_NEIGHS);
}

/// ListOP_GetMaxVal_Int - Given an integer array naList, of size nListSize, we return the max val
/// \param nListSize
/// \param naList
/// \return Maximum value in the supplied array.
int ListOP_GetMaxVal_Int(const int nListSize, const int* const naList)
{
    int i;
    int maxVal = naList[0];

    for (i = 0; i < nListSize; i++)
        {
            if (naList[i] > maxVal)
                {
                    maxVal = naList[i];
                }
        }

    return maxVal;
}

///
/// \param nVal
/// \param naOutList
/// \param nListSize
/// \param naInputList
/// \return
int ListOP_IndicesForVal_Int(const int nVal, int* const naOutList, const int nListSize, const int* const naInputList)
{
    int i;
    int nNumberOfInds = 0;
    for (i = 0; i < nListSize; i++)
        {
            if (naInputList[i] == nVal)
                {
                    naOutList[nNumberOfInds] = i;
                    nNumberOfInds++;
                }
        }
    return nNumberOfInds;
}

///
/// \param nVal
/// \param nListSize
/// \param naList
/// \return
int ListOP_GetRandomIndexForVal_Int(const int nVal, const int nListSize, const int* const naList)
{
    int nOutIdx          = -1;
    int* const naIdxList = (int*) malloc(sizeof(int) * nListSize);

    const int nTotIds = ListOP_IndicesForVal_Int(nVal, naIdxList, nListSize, naList);

    if (nTotIds == 1)
        {
            nOutIdx = naIdxList[0];
        }
    else if (nTotIds > 1)
        {
            nOutIdx = naIdxList[rand() % nTotIds];
        }

    free(naIdxList);

    return nOutIdx;
}

/// ListOP_UniqueElementsOfList_Int: Overwrite the provided _unsorted_ list of integers with only the unique elements,
/// and return the nSize of the new array. The list is sorted first, and then ListOP_UniqueElementsOfSortedList_Int is
/// run on that newly sorted list. This is a wrapper-like function for ease. The overwritten list is also sorted. \param
/// naUnsortedList: Array of ints \return number of unique elements.
int ListOP_UniqueElementsOfList_Int(const int nSize, int* const naUnsortedList)
{
    qsort(naUnsortedList, nSize, sizeof(int), UtilFunc_CompareInts);
    return ListOP_UniqueElementsOfSortedList_Int(nSize, naUnsortedList);
}

/// ListOP_Get2ndLargestVal_Int - Given an integer array naInputList, of size nListSize, we return the second largest
/// value in the array. If the array only has one value, we just return that value back. We get the sorted unique
/// elements in the input array and return the second last. \param nListSize \param naInputList \return 2nd largest
/// value in the input array.
int ListOP_Get2ndLargestVal_Int(const int nListSize, const int* const naInputList)
{
    int* naTmpList = (int*) malloc(nListSize * sizeof(int));
    int i;
    for (i = 0; i < nListSize; i++)
        {
            naTmpList[i] = naInputList[i];
        }

    // Extract Unique-sizes (this is sorted)
    const int nUniqueVals = ListOP_UniqueElementsOfList_Int(nListSize, naTmpList);

    int thisVal = naTmpList[0];

    if (nUniqueVals > 1)
        {
            thisVal = naTmpList[nUniqueVals - 2];
        }

    free(naTmpList);

    return thisVal;
}

/// ListOP_GenHistFromCounts_Int - For every occurence in the naInputList, we increment that index by 1 in naOutFreqHist
/// naOutFreqHist[naInputList[i]]++ over all elements
/// \param naOutFreqHist
/// \param naInputList
/// \param nSize
void ListOP_GenHistFromCounts_Int(int* const restrict naOutFreqHist, const int* const restrict naInputList,
                                  const int nSize)
{
    int i;
    for (i = 0; i < nSize; i++)
        {
            naOutFreqHist[naInputList[i]]++;
        }
}
