#include "alphaTree.hpp"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <limits>
#include <sstream>

#include "edge_hierarchy.hpp"
#include "incidence_matrix.hpp"
#include "timer.hpp"

/*!
 * \brief This function computes the bridge edges of the tree.
 *
 * \param[in] myEntry  The entry timestamp of the alpha-edge.
 * \param[in] myExit   The exit timestamp of the alpha-edge.
 *
 * \return A vector of size numAlphaEdges, where each element is the alpha-edge
 *         id of the bridge-edge of the corresponding alpha-edge.
 *
 * \details This function computes the bridge edges of the tree.Ï
 *1. It starts at the entry timestamp of the current alpha-edge
 *	and jumps to the exit timestamp of the current alpha-edge.
 *2. It goes through all the edges that have entry timestamp greater
 *	than the entry timestamp of the current alpha-edge, and exit
 *	timestamp less than the exit timestamp of the current alpha-edge.
 *3. For each such edge, it marks the bridge-edge of this edge as the
 *	current alpha-edge.
 */
std::vector<int> computeBridgeEdges(std::vector<indexedTimeStamp_t> &idxTS,
                                    std::vector<int> &alphaEulerEntry,
                                    std::vector<int> &alphaEulerExit) 
{
    int numAlphaEdges = idxTS.size() / 2;
    // create a vector with numAlphaEdges
    std::vector<int> bridgeEdges(numAlphaEdges, -1);

    for (int alphaId = 0; alphaId < numAlphaEdges; alphaId++)
    {
        // get the entry and exit timestamp of the current alpha-edge
        int myEntry = alphaEulerEntry[alphaId];
        // start at the entry timestamp + 1
        int myExit = alphaEulerExit[alphaId];
        int timeStamp = myEntry + 1;

        while (timeStamp < myExit)
        {
            // get the alpha-edge id of the neighbour alpha-edge
            int neighbourIdx = idxTS[timeStamp].second / 2;

            // mark the bridge-edge as the current alpha-edgeÏ
            bridgeEdges[neighbourIdx] = alphaId;
            // jump to next children of the current alpha-edge
            timeStamp = alphaEulerExit[neighbourIdx] + 1;
        }
    }

    return bridgeEdges;
}

const float FLOAT_EPSILON = std::numeric_limits<float>::epsilon();

int printEnabled = 0;

int branch2AlphaEdge(int branchId)
{
    if (branchId == 0)
        return -1;
    return (branchId - 1) / 2;
}
/**
 * @brief Find the split edges of a tree.
 * @param[in] m_incMatMST The incidence matrix of the MST.
 * @param[out] splitEdgeList The list of split edges.
 */
std::vector<int> findAlphaEdges(incidenceMatrix_t &m_incMatMST)
{
    int m_npts = m_incMatMST.m_npts;
    // list of split edges
    std::vector<int> splitEdgeList;
    // initialize outsets and other variables
    std::vector<std::vector<int>> outSet(2 * m_npts - 1);
    std::vector<int> nChildrenProcessed(m_npts - 1, 0);
    std::vector<std::vector<int>> childrenList(m_npts - 1,
                                               std::vector<int>(2, -1));
    std::vector<int> isEdgeProcessed(m_npts - 1, 0);

    // iterate through every vertex and
    // see if the vertex is incident to
    // an edge
    for (int vtxId = 0; vtxId < m_npts; vtxId++)
    {
        // if the vertex is incident to an edge
        // increment the number of processed children
        int parent = m_incMatMST.maxIncidentEdgeId(vtxId);
        nChildrenProcessed[parent]++;
    }
    // number of split edges
    int nSplitEdges = 0;
    // number of leaf edges
    int nLeafEdges = 0;
    // vector of split edges
    std::vector<int> isSplitEdge(m_npts - 1, 0);
    // iterate through every edge and
    // check if the edge is a split edge
    for (int edgeId = 0; edgeId < m_npts - 1; edgeId++)
    {
        // if the number of children processed on the edge
        // is 0 then the edge is a split edge
        if (nChildrenProcessed[edgeId] == 0)
        {
            nSplitEdges++;
            isSplitEdge[edgeId] = 1;
            splitEdgeList.push_back(edgeId);
        }
        // if the number of children processed on the edge
        // is 2 then the edge is a leaf edge
        if (nChildrenProcessed[edgeId] == 2)
            nLeafEdges++;
    }

    return splitEdgeList;
}

std::vector<indexedTimeStamp_t>
createAlphaIdxTimeStamp(std::vector<int> &alphaEdges,
                        std::vector<int> &mstEulerEntry,
                        std::vector<int> &mstEulerExit)
{
    int numAlphaEdges = alphaEdges.size();
    std::vector<indexedTimeStamp_t> idxTS(2 * numAlphaEdges);
    idxTS.resize(2 * numAlphaEdges);

    for (int i = 0; i < numAlphaEdges; i++)
    {
        int iam = alphaEdges[i];
        idxTS[2 * i].first = mstEulerEntry[iam];
        idxTS[2 * i].second = 2 * i;
        idxTS[2 * i + 1].first = mstEulerExit[iam];
        idxTS[2 * i + 1].second = 2 * i + 1;
    }

    // sorting the indexedTimeStamp array
    std::sort(idxTS.begin(), idxTS.end());
    return idxTS;
}

std::pair<std::vector<int>, std::vector<int>>
computeAlphaEulerEntryExit(std::vector<indexedTimeStamp_t> &idxTS)
{
    int numAlphaEdges = idxTS.size() / 2;
    std::vector<int> alphaEulerEntry(numAlphaEdges);
    std::vector<int> alphaEulerExit(numAlphaEdges);

    for (int i = 0; i < 2 * numAlphaEdges; i++)
    {
        int timeStamp = idxTS[i].first;
        int iam = idxTS[i].second / 2;
        if (idxTS[i].second % 2 == 0)
        {
            alphaEulerEntry[iam] = i;
            // printf("Entering %d at %d\n", iam, i);
        }
        else
        {
            alphaEulerExit[iam] = i;
            // printf("Exiting %d at %d\n", iam, i);
        }
    }

    return std::make_pair(alphaEulerEntry, alphaEulerExit);
}

alphaTree_t::alphaTree_t(incidenceMatrix_t &_incMatMST,
                         std::vector<wtEdge_t> &_wtSortedMST)
    : m_incMatMST(_incMatMST), m_wtSortedMST(_wtSortedMST)
{
    m_npts = m_incMatMST.m_npts;
    // do the eulerTour
    eulerTour();
    // get the alpha edges
    alphaEdges = findAlphaEdges(m_incMatMST);
    numAlphaEdges = alphaEdges.size();

    // perform the eulerTour on Alpha-MST
    // sets idxTS, alphaEulerEntry and alphaEulerExit
    mytimer_t tAlphaEulerTour(std::string("Alpha-euler tour"));
    mytimer_t tAlphaTree(std::string("Alpha-tree construction"));

    tAlphaEulerTour.start();
    // create the indexedTimeStamp array
    idxTS = createAlphaIdxTimeStamp(alphaEdges, mstEulerEntry, mstEulerExit);
    auto alphaEulerEntryExit = computeAlphaEulerEntryExit(idxTS);
    alphaEulerEntry = alphaEulerEntryExit.first;
    alphaEulerExit = alphaEulerEntryExit.second;

    tAlphaEulerTour.end();
    // compute bridgeEdge
    alphaBridgeEdges =
        ::computeBridgeEdges(idxTS, alphaEulerEntry, alphaEulerExit);
    // construct Alpha tree
    tAlphaTree.start();
    alphaParents = constructAlphaTree();
    tAlphaTree.end();

    tAlphaEulerTour.print();
    tAlphaTree.print();
}

int alphaTree_t::isLeftOfEdge(int queryEdge, int rootEdge)
{
    if (mstEulerEntry[queryEdge] > mstEulerEntry[rootEdge] &&
        mstEulerExit[queryEdge] < mstEulerExit[rootEdge])
        return 1;

    return 0;
}

/**
 * \brief Construct the tree of split edges.
 *
 * This function takes in a vector of alpha edges, and returns a vector of
 * parents. The vector of alpha edges isn't modified.
 *
 * \param alphaEdges a vector of alpha edges.
 *
 * \return std::vector<int> list of parents \todo This comment is useless. Most
 * of it is from the header file, and the rest is a line-by-line description of
 * every line of the function.
 */

std::vector<int> alphaTree_t::constructAlphaTree()
{

    int numSplitEdges = alphaEdges.size();

    std::vector<int> parent(numSplitEdges, 0);
    std::vector<int> isProcessed(numSplitEdges, 0);
    std::vector<int> myBranch(numSplitEdges, 0);
    std::vector<int> maxVertexMyBranch(2 * m_npts, m_npts);
    std::vector<int> maxVertexAlphaId(2 * m_npts, -1);

    // First, initialise the data structures:

    isProcessed[0] = 1;
    int numNotProcessed = numSplitEdges - 1;
    int height = 0;

    parent[0] = -1;

    // Then, create the tree in two passes: first, assign all the children to the
    // nodes, and then assign all the parents. This algorithm is generic and works
    // for the alpha tree and for the voronoi tree.

    while (numNotProcessed > 0)
    {
        for (int edge = 1; edge < numSplitEdges; edge++)
        {
            if (!isProcessed[edge])
            {
                // First, figure out who is my parent
                int iam = alphaEdges[edge];

                myBranch[edge] = 2 * parent[edge];

                if (isLeftOfEdge(iam, alphaEdges[parent[edge]]))
                    myBranch[edge]++;

                if (maxVertexMyBranch[myBranch[edge]] > iam)
                {
                    maxVertexMyBranch[myBranch[edge]] = iam;
                    maxVertexAlphaId[myBranch[edge]] = edge;
                }
                // maxVertexMyBranch[myBranch[edge]] =
                // std::min(maxVertexMyBranch[myBranch[edge]],iam ); printf(" %d : \t {
                // %d } --> \t %d \n ", iam, myBranch[edge],
                // maxVertexMyBranch[myBranch[edge]]);
            }
        }

        for (int edge = 0; edge < numSplitEdges; edge++)
        {
            if (!isProcessed[edge])
            {
                int iam = alphaEdges[edge];
                if (maxVertexMyBranch[myBranch[edge]] == iam)
                {
                    isProcessed[edge] = 1;
                    numNotProcessed--;
                    // std::cout<<iam<<"	\t -> \t "<< parent[edge] << "\n";
                    if (printEnabled)
                        std::cout << iam << "	\t -> \t " << alphaEdges[parent[edge]]
                                  << "\n";
                }
                else
                {
                    // parent[edge] = maxVertexMyBranch[myBranch[edge]];
                    parent[edge] = maxVertexAlphaId[myBranch[edge]];
                }
            }
        }

        height++;
    }

    std::cout << "Height of the split edge tree is \t:" << height << "\n";
    return parent;
}

std::vector<int> alphaTree_t::constructEdgeTree()
{
    int printEnabled = 0;
    // Assumption no edges have been processed

    std::vector<branchEdge_t> branchEdge = computeBranchEdge();
    // initialize the member variable
    m_sortedBranchEdgeArr = branchEdge;
    // mark the parents
    std::vector<int> globalParents(m_npts - 1, -1);
    mytimer_t tMarkingParents("Marking Parents");
    tMarkingParents.start();
    for (int edgeId = 1; edgeId < m_npts - 1; edgeId++)
    {
        // we are on same branch
        if (branchEdge[edgeId - 1].branch == branchEdge[edgeId].branch)
        {
            globalParents[branchEdge[edgeId].edgeIdx] =
                branchEdge[edgeId - 1].edgeIdx;
            if (printEnabled)
                printf(" %d -> %d \n", branchEdge[edgeId].edgeIdx,
                       branchEdge[edgeId - 1].edgeIdx);
        }
        else // Iam start of the branch
        {
            // globalParents[branchEdge[edgeId].edgeIdx] =
            // alphaEdges[branchEdge[edgeId].branch/2];
            globalParents[branchEdge[edgeId].edgeIdx] =
                alphaEdges[branch2AlphaEdge(branchEdge[edgeId].branch)];

            if (printEnabled)
                printf(" %d -> %d \n", branchEdge[edgeId].edgeIdx,
                       globalParents[branchEdge[edgeId].edgeIdx]);
        }
    }
    tMarkingParents.end();
    tMarkingParents.print();

    return globalParents;
}

int alphaTree_t::findAlphaParent(int edgeId)
{

    if (edgeId < alphaEdges[0])
        return -1;
    indexedTimeStamp_t idxEntry(mstEulerEntry[edgeId], 0);
    indexedTimeStamp_t idxExit(mstEulerExit[edgeId], 0);
    int myBridge;
    int myEntry;
    int myExit;
    // int firstTimeStamp = idxTS[0].first;
    int lastTimeStamp = idxTS[idxTS.size() - 1].first;
    if (mstEulerEntry[edgeId] > lastTimeStamp)
    {
        myBridge = -1;
        myEntry = idxTS.size();
        myExit = idxTS.size();
    }
    else
    {
        myEntry = std::distance(
            idxTS.begin(), std::lower_bound(idxTS.begin(), idxTS.end(), idxEntry));
        myExit = std::distance(
            idxTS.begin(), std::lower_bound(idxTS.begin(), idxTS.end(), idxExit));

        // quick return
        assert(myEntry >= 0);
        if (!(myEntry < idxTS.size()))
        {
            printf("edgeId%d, (myEntry %d myExit %d) idxTS.size() %d, lastTS %d, "
                   "myTs =%d \n",
                   edgeId, myEntry, myExit, idxTS.size(),
                   idxTS[idxTS.size() - 1].first, mstEulerEntry[edgeId]);
        }
        assert(myEntry < idxTS.size()); // this is failing
        if (edgeId == alphaEdges[idxTS[myEntry].second / 2])
        {
            return alphaParents[idxTS[myEntry].second / 2];
        }

        if (printEnabled)
            printf("%d :\t mstEntry %d mstExit %d newEntry %d newExit %d: \
\n \t\t%d  %d\n",
                   edgeId, mstEulerEntry[edgeId], mstEulerExit[edgeId], myEntry,
                   myExit, alphaEdges[idxTS[myEntry].second / 2],
                   idxTS[myEntry].second);

        if (myEntry == 0)
            myBridge = -1;
        else if (idxTS[myEntry].second % 2 == 0) // starting point
            myBridge = alphaBridgeEdges[idxTS[myEntry].second / 2];
        else
            myBridge = idxTS[myEntry].second / 2;
    }

    int minDescendent = m_npts;
    int minDescIdx = -1;
    int maxAncestor = -1;
    int maxAncIdx = -1;

    int startTime, endTime, minDesc;

    if (myBridge == -1)
    {
        startTime = 0;
        endTime = 2 * numAlphaEdges;
    }
    else
    {
        startTime = alphaEulerEntry[myBridge] + 1;
        endTime = alphaEulerExit[myBridge];

        int bridgeId = alphaEdges[myBridge];
        if (bridgeId > edgeId and bridgeId < minDescendent)
        {
            minDescendent = bridgeId;
            minDescIdx = myBridge;
        }

        if (bridgeId < edgeId and bridgeId > maxAncestor)
        {
            maxAncestor = bridgeId;
            maxAncIdx = myBridge;
        }
    }

    int ts = startTime;
    if (printEnabled)
        printf("%d :\t bridge %d startTime %d endTime %d max %d min %d\n ", edgeId,
               myBridge, startTime, endTime, maxAncestor, minDescendent);
    while (ts < endTime)
    {
        int neighbourIdx = idxTS[ts].second / 2;
        int neighbourEdgeId = alphaEdges[neighbourIdx];
        if (printEnabled)
            printf("%d : \t checking %d\n", edgeId, neighbourEdgeId);
        if (neighbourEdgeId > edgeId and neighbourEdgeId < minDescendent)
        {
            minDescendent = neighbourEdgeId;
            minDescIdx = neighbourIdx;
        }

        if (neighbourEdgeId < edgeId and neighbourEdgeId > maxAncestor)
        {
            maxAncestor = neighbourEdgeId;
            maxAncIdx = neighbourIdx;
        }

        ts = alphaEulerExit[neighbourIdx] + 1;
    }
    if (printEnabled)
        printf(":\t min %d max %d \n ", minDescendent, maxAncestor);

    if (minDescIdx == -1) // leaf node
        return maxAncIdx;

    //
    int myParent = minDescendent;
    int parentIdx = minDescIdx;

    while (myParent > edgeId)
    {
        parentIdx = alphaParents[parentIdx];
        myParent = alphaEdges[parentIdx];
    }

    return parentIdx;
}

std::vector<branchEdge_t> alphaTree_t::computeBranchEdge()
{
    if (m_sortedBranchEdgeArr.size() != 0)
        return m_sortedBranchEdgeArr;
    std::vector<branchEdge_t> branchEdge(m_npts - 1);
    // for each edge find its branch
    mytimer_t tFindAlphaParent("Finding alpha parents");
    tFindAlphaParent.start();
    for (int edgeId = 0; edgeId < m_npts - 1; edgeId++)
    // for(int edgeId =0; edgeId< 10; edgeId++)
    {
        int alParent = findAlphaParent(edgeId);
        if (alParent != -1 and printEnabled)
            printf(" alphaParent(%d)  :\t %d (%d) \n", edgeId, alphaEdges[alParent],
                   alParent);
        // TODO: -1 case
        if (alParent == -1)
            branchEdge[edgeId] = branchEdge_t(0, edgeId);
        else
        {
            // int mstParent = alphaParents[alParent];
            int mstParent = alphaEdges[alParent];
            int myBranch = 2 * alParent + 1;
            myBranch += isLeftOfEdge(edgeId, mstParent);
            branchEdge[edgeId] = branchEdge_t(myBranch, edgeId);
        }
    }
    tFindAlphaParent.end();
    tFindAlphaParent.print();
    // exit(0);
    // TODO: branchEdgeComparator
    mytimer_t tSortBranchEdge("Branch Edge-sorting");
    tSortBranchEdge.start();
    std::sort(branchEdge.begin(), branchEdge.end(), branchEdgeComparator());
    tSortBranchEdge.end();
    tSortBranchEdge.print();

    m_sortedBranchEdgeArr = branchEdge;
    return branchEdge;
}

int alphaTree_t::alphaParent2Branch(int edgeId, int alphaParentIdx)
{
    if (alphaParentIdx < 0)
        return 0;
    int myBranch = 2 * alphaParentIdx + 1;
    if (isLeftOfEdge(edgeId, alphaEdges[alphaParentIdx]))
        myBranch++;

    return myBranch;
}

int alphaTree_t::branchParent(int branchId)
{
    if (branchId == 0)
        return -1;
    int alphaEdgeId = branch2AlphaEdge(branchId);
    int edgeId = alphaEdges[alphaEdgeId];
    // if edge 0 is an alpha edge then its branch id is 0
    if (edgeId == 0)
        return 0;
    return alphaParent2Branch(edgeId, alphaParents[alphaEdgeId]);
}

#include "stability.hpp"

std::vector<int> alphaTree_t::computeFlatClustering(int minClusterSize)
{
    int verbose = 1;
    std::vector<int> flatClustering;
    m_numBranches = 2 * numAlphaEdges + 1;
    // std::vector<int> brHeads(m_numBranches, -1);
    brIHeads.resize(m_numBranches, -1);
    brDescendents.resize(m_numBranches, 0);

    brStabilityScores.resize(m_numBranches, 0.0);
    brSHatScores.resize(m_numBranches, 0.0);
    brDelta.resize(m_numBranches, 1);
    brLength.resize(m_numBranches, 0);

    // std::vector<float> chainStabilityContribution(m_numBranches, 0.0);
    std::vector<double> chainStabilityContribution(m_numBranches, 0.0);
    std::vector<stability_t> SscoreBr(
        m_numBranches); // stability score of each branch
    std::vector<stability_t> HscoreBr(
        m_numBranches); // s-hat score of each branch
    std::vector<stability_t> CScontri(
        m_numBranches); // chain stability contribution of each branch
    // first compute branch edge
    std::vector<branchEdge_t> branchEdge = computeBranchEdge();
    std::vector<int> isLeafBranch(m_numBranches, 0);
    // computing chain lengths
    // why is it starting from 1?
    for (int iEdge = 0; iEdge < m_npts - 1; iEdge++)
    {
        // branchEdge is not sorted so iEdge is not equal to
        // edgeIdx of MST
        int myBranch = branchEdge[iEdge].branch;
        int edgeId = branchEdge[iEdge].edgeIdx;
        int myAlParent = branch2AlphaEdge(myBranch);

        if (m_incMatMST.isLeafEdge(edgeId))
        {

            brLength[myBranch] += 2;
            isLeafBranch[myBranch] = 1;
            if (verbose)
                printf("Stb %f: adding %f\n  ", chainStabilityContribution[myBranch],
                       2.0 / m_wtSortedMST[edgeId].second);
            chainStabilityContribution[myBranch] +=
                2.0 / m_wtSortedMST[edgeId].second;
            CScontri[myBranch] = CScontri[myBranch] + std::make_pair(edgeId, 2);
        }
        else if (!m_incMatMST.isAlphaEdge(edgeId))
        {
            brLength[myBranch] += 1;
            if (verbose)
                printf("Stb %f: adding %f\n  ", chainStabilityContribution[myBranch],
                       1.0 / m_wtSortedMST[edgeId].second);
            chainStabilityContribution[myBranch] +=
                1.0 / m_wtSortedMST[edgeId].second;
            CScontri[myBranch] = CScontri[myBranch] + std::make_pair(edgeId, 1);
        }

        // Top of the chain
        if (iEdge == 0 || branchEdge[iEdge].branch != branchEdge[iEdge - 1].branch)
        {
            brIHeads[myBranch] = iEdge;
        }
    }

    // compute brDescendents
    std::vector<int> numBranchChildrens(m_numBranches, 0);
    std::vector<int> numValidChildrens(m_numBranches, 0);
    std::vector<int> numChildBranchProcessed(m_numBranches, 0);
    std::vector<int> qBranchProcessed(m_numBranches, 0);
    for (int i = 0; i < m_numBranches; i++)
    {
        if (branchParent(i) != -1)
        {
            numBranchChildrens[branchParent(i)]++;
        }
    }
    // assert(0);

    // bottom-up traversal
    int numBranch2Process = m_numBranches;

    while (numBranch2Process)
    {
        for (int i = 0; i < m_numBranches; i++)
        {
            if (qBranchProcessed[i] == 0 and
                numBranchChildrens[i] == numChildBranchProcessed[i])
            {
                if (verbose)
                    std::cout << i << " -> Parent =" << branchParent(i) << "\n";
                if (branchParent(i) != -1)
                {

                    // update descendents
                    int totalDescendents = brDescendents[i] + brLength[i];
                    brDescendents[branchParent(i)] += (totalDescendents);

                    // so there are three cases
                    // 1. either you are invalid
                    // or 2. your number of descdends are greater than mincluster
                    // or 3. it becomes a valid cluster at this point
                    // the third case is the difficult one
                    if (verbose)
                        printf(" Processing %d: stabScore =%f, chain contri =%g\n", i,
                               brStabilityScores[i], chainStabilityContribution[i]);

                    if (totalDescendents < minClusterSize)
                    {
                        // invalid cluster
                        brStabilityScores[i] = 0;
                        brDelta[i] = 0;
                        brSHatScores[i] = 0;

                        SscoreBr[i].setZero();
                        HscoreBr[i].setZero();
                    }
                    else
                    {
                        // update number of valid childrens
                        numValidChildrens[branchParent(i)]++;

                        // float ancStabContribution =
                        double ancStabContribution =
                            1.0 * totalDescendents /
                            m_wtSortedMST[alphaEdges[branch2AlphaEdge(i)]].second;
                        score_t ancSContri(alphaEdges[branch2AlphaEdge(i)],
                                           totalDescendents);

                        if (verbose)
                            printf("Ancestor = %d: length=%d, totaldescend=%d, contribution= "
                                   "%f\n",
                                   alphaEdges[branch2AlphaEdge(i)], brLength[i],
                                   totalDescendents, ancStabContribution);

                        int chainLength =
                            isLeafBranch[i] ? brLength[i] - 1 : brLength[i] + 1;
                        int branchStartEdge = brIHeads[i];
                        if (brDescendents[i] >= minClusterSize)
                        {

                            int etxFinalEdgeId =
                                branchEdge[branchStartEdge + chainLength - 1].edgeIdx;
                            float etWt = 1.0 / m_wtSortedMST[etxFinalEdgeId].second;
                            chainStabilityContribution[i] += etWt * brDescendents[i];
                            CScontri[i] = CScontri[i] +
                                          std::make_pair(etxFinalEdgeId, brDescendents[i]);
                        }
                        else
                        {
                            // print edgeId etc
                            int eta = minClusterSize - brDescendents[i];
                            if (isLeafBranch[i])
                                eta -= 2;
                            if (eta < 0)
                                eta = 0;
                            int chainLenX= chainLength - 1;
                            // int chainLenX= chainLength;
                            int etStart = chainLenX - eta;
                            
                            if (isLeafBranch[i])
                            std::cout<< "This new print: " << i << " " << brDescendents[i] << " " << brLength[i] << " " << chainLength << " " << etStart << "\n ";
                            
                            
                            // if(isLeafBranch[i]) etStart++;
                            int etxFinalEdgeId =
                                branchEdge[branchStartEdge + etStart].edgeIdx;
                            float etWt = 1.0 / m_wtSortedMST[etxFinalEdgeId].second;

                            /* Process upto last-but one edge*/ 
                            for (int etx = etStart; etx < chainLenX; etx++)
                            {
                                int edgeId = branchEdge[branchStartEdge + etx].edgeIdx;
                                chainStabilityContribution[i] +=
                                    (etWt - 1.0 / m_wtSortedMST[edgeId].second);
                                score_t sEtWt(std::make_pair(etxFinalEdgeId, 1));
                                CScontri[i] = (CScontri[i] + sEtWt) - std::make_pair(edgeId, 1);
                            }
                            /* Is the last edge is a leaf edge then the */ 
                            if (isLeafBranch[i] and minClusterSize != 1)
                            {
                                int edgeId =
                                    branchEdge[branchStartEdge + (chainLength - 1)].edgeIdx;
                                chainStabilityContribution[i] +=
                                    2 * (etWt - 1.0 / m_wtSortedMST[edgeId].second);
                                score_t sEtWt(std::make_pair(etxFinalEdgeId, 2));
                                CScontri[i] = (CScontri[i] + sEtWt) - std::make_pair(edgeId, 2);
                            }
                            chainStabilityContribution[i] += etWt * brDescendents[i];
                            score_t sEtWt(std::make_pair(etxFinalEdgeId, brDescendents[i]));
                            CScontri[i] = (CScontri[i] + sEtWt);
                            if (verbose)
                                printf("%d : eta=%d, etStart=%d etxFinalEdgeId =%d \n", i, eta,
                                       etStart, etxFinalEdgeId);
                        }

                        brDelta[i] = 1;

                        if (numValidChildrens[i] == 0)
                        {
                            brStabilityScores[i] =
                                chainStabilityContribution[i] - ancStabContribution;
                            SscoreBr[i] = CScontri[i] - ancSContri;
                            if(isLeafBranch[i])
                            {
                                std::cout<<i<< " Stability scors ";
                                SscoreBr[i].print();
                                // convert SscoreBr[i] to float
                                float score = 0;
                                for (auto it = SscoreBr[i].begin(); it != SscoreBr[i].end(); it++)
                                {
                                    score += it->second * 1.0 / m_wtSortedMST[it->first].second;
                                }
                                std::cout<<i<< " Stability scors " << score << " alterCalc "<<brStabilityScores[i]<< 
                                "Error tol="<< std::max(chainStabilityContribution[i], ancStabContribution) *
                                    FLOAT_EPSILON<<" \t:\t";
                                std::cout<<i<< " Chain contri scors ";
                                CScontri[i].print();
                                std::cout<<i<< " Anc contri scors "<<ancSContri.first<<" "<<ancSContri.second<<"\n";
                                // ancSContri.print();
                                std::cout<<"\n";

                                // brStabilityScores[i] += 1.0 / m_wtSortedMST[alphaEdges[branch2AlphaEdge(i)]].second;
                                // SscoreBr[i] = SscoreBr[i] + std::make_pair(alphaEdges[branch2AlphaEdge(i)], 1);
                            }
                            // TODO: check the following
                            if (SscoreBr[i].isZero())
                            {
                                assert(brStabilityScores[i] <=
                                       std::max(chainStabilityContribution[i],
                                                ancStabContribution) *
                                           FLOAT_EPSILON);
                            }
                            if (brStabilityScores[i] <=
                                std::max(chainStabilityContribution[i], ancStabContribution) *
                                    FLOAT_EPSILON)
                            {
                                brStabilityScores[i] = 0;
                                brDelta[i] = 0;
                            }

                            brSHatScores[i] = brStabilityScores[i];
                            HscoreBr[i] = SscoreBr[i];
                            if (verbose)
                                printf(" %d : 0: chain = %f, ancContri =%f, ,%d, %d \n", i,
                                       chainStabilityContribution[i], ancStabContribution,
                                       brDescendents[i], brLength[i]);
                        }
                        else if (numValidChildrens[i] == 1)
                        {
                            auto selfContri =
                                brStabilityScores[i] + chainStabilityContribution[i];
                            auto selfSContri = SscoreBr[i] + CScontri[i];
                            brStabilityScores[i] =
                                brStabilityScores[i] +
                                (chainStabilityContribution[i] - ancStabContribution);
                            SscoreBr[i] = (SscoreBr[i] + CScontri[i]) - ancSContri;
                            if (SscoreBr[i].isZero())
                            {
                                assert(brStabilityScores[i] <=
                                       std::max(selfContri, ancStabContribution) *
                                           FLOAT_EPSILON);
                            }
                            if (brStabilityScores[i] <=
                                std::max(selfContri, ancStabContribution) * FLOAT_EPSILON)
                            {
                                brStabilityScores[i] = 0;
                                brDelta[i] = 0;
                            }
                            if (verbose)
                                printf(" edge=%d, s =%f, hat =%f, diff=%g \n", i, brStabilityScores[i],
                                       brSHatScores[i], brSHatScores[i] - brStabilityScores[i]);
                            if (brSHatScores[i] > brStabilityScores[i])
                            {
                                brDelta[i] = 0; // not a final cluster
                            }
                            else
                            {
                                brSHatScores[i] = brStabilityScores[i];
                                HscoreBr[i] = SscoreBr[i];
                            }
                            if (verbose)
                                printf(" %d : 1: descdn =%f, chain = %f, ancContri =%f \n", i,
                                       brStabilityScores[i], chainStabilityContribution[i],
                                       ancStabContribution);
                        }
                        else /* two valid childrens */
                        {
                            brStabilityScores[i] =
                                chainStabilityContribution[i] - ancStabContribution;
                            SscoreBr[i] = CScontri[i] - ancSContri;
                            if (SscoreBr[i].isZero())
                            {
                                assert(brStabilityScores[i] <=
                                       std::max(chainStabilityContribution[i],
                                                ancStabContribution) *
                                           FLOAT_EPSILON);
                            }
                            if (brStabilityScores[i] <=
                                std::max(chainStabilityContribution[i], ancStabContribution) *
                                    FLOAT_EPSILON)
                            {
                                brStabilityScores[i] = 0;
                                brDelta[i] = 0;
                                // std::cout<<
                            }
                            if (verbose)
                            {
                                printf(" %d : 2: chain = %f, ancContri =%f \n", i,
                                       chainStabilityContribution[i], ancStabContribution);
                                printf("%d : stb %f, shat %f\n ", i, brStabilityScores[i],
                                       brSHatScores[i]);
                            }

                            // if(HscoreBr[i]>SscoreBr[i] )
                            // {

                            // }
                            if (brSHatScores[i] > brStabilityScores[i])
                            {
                                auto tmpScore = HscoreBr[i] - SscoreBr[i];
                                assert(!tmpScore.isZero()); // tmpScore should not be zero
                                brDelta[i] = 0;             // not a final cluster
                            }
                            else
                            {
                                brSHatScores[i] = brStabilityScores[i];
                                HscoreBr[i] = SscoreBr[i];
                            }
                        }

                        brStabilityScores[branchParent(i)] += brStabilityScores[i];
                        brSHatScores[branchParent(i)] += brSHatScores[i];
                        SscoreBr[branchParent(i)] = SscoreBr[branchParent(i)] + SscoreBr[i];
                        HscoreBr[branchParent(i)] = HscoreBr[branchParent(i)] + HscoreBr[i];
                    }
                    if (verbose)
                    {
                        // stabScore =
                        if (SscoreBr[i].getScore() != 0)
                        {
                            printf("  NONZERO total SCORE !!!:%d\t", SscoreBr[i].getScore());
                            SscoreBr[i].print();
                        }
                        else
                        {
                            printf("  ZERO SCORE:\t");
                        }
                        printf(" Final %d: stabScore =%g, chain contri =%g\n", i,
                               brStabilityScores[i], chainStabilityContribution[i]);
                    }

                    numChildBranchProcessed[branchParent(i)]++;
                }
                // else
                // {
                //     brDelta[i]=0;
                // }
                if (verbose)
                    printf(" branch = %d, stab = %f shat =%f, delta=%d \n", i,
                           brStabilityScores[i], brSHatScores[i], brDelta[i]);
                qBranchProcessed[i] = 1;
                numBranch2Process--;
            }
        }
    }

    for (int i = 0; i < m_numBranches; i++)
    {
        // only if brnach-i is a true cluster, delta should be one
        if (branchParent(i) == -1 || numValidChildrens[branchParent(i)] != 2)
            brDelta[i] = 0;

        int totalDescendents = brDescendents[i] + brLength[i];
        if (totalDescendents == minClusterSize)
        {
            printf(" i =%d, branchHead=%d,branchParent = %d,  "
                   "numvalidChildrenParent=%d, parentdelta=%d \n",
                   i, m_sortedBranchEdgeArr[brIHeads[i]].edgeIdx, branchParent(i),
                   numValidChildrens[branchParent(i)], brDelta[branchParent(i)]);
        }
    }

    // To compute Flat clustering for each branch
    std::vector<int> brFlatCluster = branchFlatCluster(brDelta, m_numBranches);
    // To compute the final clusters
    std::vector<int> edgeFlatMap =
        computeFlatMap(m_npts, m_sortedBranchEdgeArr, brFlatCluster, brIHeads);

    // for(int vtx=0; vtx< m_npts; vtx++)
    // {
    //     auto parentEdge = ??
    //     flatClustering[vtx] = edgeFlatMap[parentEdge];

    // }
    // TODO: dummy return ;
    // edgeFlatMap
    // return flatClustering;
    return edgeFlatMap;
}

/**
 * @brief This function is used to map the edges to their corresponding flat
 * edges. The function is called in the function "FlatMap".
 *
 * @param npts The number of points
 * @param sortedBranchEdgeArr The sorted branch edge array
 * @param brFlatCluster The cluster of the flat branch
 * @param brIHeads The index of the head of the branch
 */
std::vector<int>
computeFlatMap(int npts, const std::vector<branchEdge_t> &sortedBranchEdgeArr,
               const std::vector<int> &brFlatCluster,
               const std::vector<int> &brIHeads)
{
    // Initialize the edge flat map
    std::vector<int> edgeFlatMap(npts - 1, -1);
    // The number of edges
    int nedges = npts - 1;
    // For each edge
    for (int edgeId = 0; edgeId < nedges; edgeId++)
    {
        // Get the global edge index
        auto gEdgeId = sortedBranchEdgeArr[edgeId].edgeIdx;
        // Get the branch index
        auto branchId = sortedBranchEdgeArr[edgeId].branch;
        // If the branch is flat
        if (brFlatCluster[branchId] != -1)
        {
            // Set the edge flat map
            edgeFlatMap[gEdgeId] =
                sortedBranchEdgeArr[brIHeads[brFlatCluster[branchId]]].edgeIdx;
        }
        else
        {
            // Set the edge flat map
            edgeFlatMap[gEdgeId] = -1;
        }
    }
    return edgeFlatMap;
}

/**
 * @brief This function returns a vector of flat cluster ids for each branch
 * @param brDelta a vector of branch delta values
 * @param m_numBranches number of branches
 * @return a vector of flat cluster ids for each branch
 */
std::vector<int> alphaTree_t::branchFlatCluster(std::vector<int> &brDelta,
                                                int m_numBranches)
{
    std::vector<int> brFlatCluster(m_numBranches, -1);
    for (int i = 0; i < m_numBranches; i++)
    {
        auto brParent = i;
        while (brParent != -1)
        {

            if (brDelta[brParent] == 1)
                brFlatCluster[i] = brParent;
            brParent = branchParent(brParent);
        }
    }
    return brFlatCluster;
}