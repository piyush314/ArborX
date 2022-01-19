#include "alphaTree.hpp"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <sstream>

#include "edge_hierarchy.hpp"
#include "incidence_matrix.hpp"
#include "timer.hpp"

int printEnabled = 0;

int branch2AlphaEdge(int branchId)
{
	if (branchId == 0)
		return -1;
	return (branchId - 1) / 2;
}
std::vector<int> alphaTree_t::getAlphaEdges()
{
	std::vector<int> splitEdgeList;
	// initialize outsets
	std::vector<std::vector<int>> outSet(2 * m_npts - 1);
	std::vector<int> nChildrenProcessed(m_npts - 1, 0);
	std::vector<std::vector<int>> childrenList(m_npts - 1,
                 std::vector<int>(2, -1));
	std::vector<int> isEdgeProcessed(m_npts - 1, 0);
	// std::vector<int> vertexMaxIncidentEdgeId = m_incMatMST.maxIncidentEdgeId()

	for (int vtxId = 0; vtxId < m_npts; vtxId++)
	{
		int parent = m_incMatMST.maxIncidentEdgeId(vtxId);
		nChildrenProcessed[parent]++;
		// m_numDescendents[parent]++;
	}
	int nSplitEdges = 0;
	int nLeafEdges = 0;
	std::vector<int> isSplitEdge(m_npts - 1, 0);
	for (int edgeId = 0; edgeId < m_npts - 1; edgeId++)
	{
		// if(edgeId==0)

		if (nChildrenProcessed[edgeId] == 0)
		{
			nSplitEdges++;
			isSplitEdge[edgeId] = 1;
			splitEdgeList.push_back(edgeId);
		}

		if (nChildrenProcessed[edgeId] == 2)
			nLeafEdges++;
	}

	return splitEdgeList;
}

void alphaTree_t::alphaEulerTour()
{
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
	alphaEulerEntry.resize(numAlphaEdges);
	alphaEulerExit.resize(numAlphaEdges);

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
}

alphaTree_t::alphaTree_t(incidenceMatrix_t &_incMatMST,
												 std::vector<wtEdge_t> &_wtSortedMST)
		: m_incMatMST(_incMatMST)
		, m_wtSortedMST(_wtSortedMST)
{
	m_npts = m_incMatMST.m_npts;
	// do the eulerTour
	eulerTour();
	// get the alpha edges
	alphaEdges = getAlphaEdges();
	numAlphaEdges = alphaEdges.size();

	// perform the eulerTour on Alpha-MST
	// sets idxTS, alphaEulerEntry and alphaEulerExit
	timer_t tAlphaEulerTour(std::string("Alpha-euler tour"));
	timer_t tAlphaTree(std::string("Alpha-tree construction"));

	tAlphaEulerTour.start();
	alphaEulerTour();
	tAlphaEulerTour.end();
	// compute bridgeEdge
	alphaBridgeEdges = computeBrideEdges();
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

std::vector<int> alphaTree_t::constructAlphaTree()
{

	//  = listOfSplitEdge();
	int numSplitEdges = alphaEdges.size();

	std::vector<int> parent(numSplitEdges, 0);
	std::vector<int> isProcessed(numSplitEdges, 0);
	std::vector<int> myBranch(numSplitEdges, 0);
	std::vector<int> maxVertexMyBranch(2 * m_npts, m_npts);
	std::vector<int> maxVertexAlphaId(2 * m_npts, -1);
	isProcessed[0] = 1;
	int numNotProcessed = numSplitEdges - 1;
	int height = 0;

	parent[0] = -1;
	while (numNotProcessed > 0)
	{
		for (int edge = 1; edge < numSplitEdges; edge++)
		{
			if (!isProcessed[edge])
			{
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
					// std::cout<<iam<<"  \t -> \t "<< parent[edge] << "\n";
					std::cout << iam << "  \t -> \t " << alphaEdges[parent[edge]] << "\n";
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
	timer_t tMarkingParents("Marking Parents");
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

std::vector<int> alphaTree_t::computeBrideEdges()
{
	std::vector<int> bridgeEdges(numAlphaEdges, -1);

	for (int alphaId = 0; alphaId < numAlphaEdges; alphaId++)
	{
		int myEntry = alphaEulerEntry[alphaId];
		int myExit = alphaEulerExit[alphaId];

		int timeStamp = myEntry + 1;

		while (timeStamp < myExit)
		{
			int neighbourIdx = idxTS[timeStamp].second / 2;
			bridgeEdges[neighbourIdx] = alphaId;

			timeStamp = alphaEulerExit[neighbourIdx] + 1;
		}
	}

	return bridgeEdges;
}

int alphaTree_t::findAlphaParent(int edgeId)
{

	if (edgeId < alphaEdges[0])
		return -1;
	indexedTimeStamp_t idxEntry(mstEulerEntry[edgeId], 0);
	indexedTimeStamp_t idxExit(mstEulerExit[edgeId], 0);
	int myEntry = std::distance(
			idxTS.begin(), std::lower_bound(idxTS.begin(), idxTS.end(), idxEntry));
	int myExit = std::distance(
			idxTS.begin(), std::lower_bound(idxTS.begin(), idxTS.end(), idxExit));

	// quick return
	if (edgeId == alphaEdges[idxTS[myEntry].second / 2])
	{
		return alphaParents[idxTS[myEntry].second / 2];
	}

	if (printEnabled)
		printf("%d :\t mstEntry %d mstExit %d newEntry %d newExit %d: \
		\n \t\t%d  %d\n",
					 edgeId, mstEulerEntry[edgeId], mstEulerExit[edgeId], myEntry, myExit,
					 alphaEdges[idxTS[myEntry].second / 2], idxTS[myEntry].second);

	int myBridge;

#if 0
		if(myEntry==0)
				myBridge =-1;
		else if(idxTS[myEntry-1].second %2 ==0 ) //starting point
				myBridge = idxTS[myEntry-1].second/2;
		else 
				myBridge = alphaBridgeEdges[idxTS[myEntry-1].second/2];
#else
	if (myEntry == 0)
		myBridge = -1;
	else if (idxTS[myEntry].second % 2 == 0) // starting point
		myBridge = alphaBridgeEdges[idxTS[myEntry].second / 2];
	else
		myBridge = idxTS[myEntry].second / 2;
#endif
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
    if(m_sortedBranchEdgeArr.size() != 0)
         return m_sortedBranchEdgeArr;
	std::vector<branchEdge_t> branchEdge(m_npts - 1);
	// for each edge find its branch
	timer_t tFindAlphaParent("Finding alpha parents");
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
	timer_t tSortBranchEdge("Branch Edge-sorting");
	tSortBranchEdge.start();
	std::sort(branchEdge.begin(), branchEdge.end(), branchEdgeComparator());
	tSortBranchEdge.end();
	tSortBranchEdge.print();

    m_sortedBranchEdgeArr = branchEdge;
	return branchEdge;
}

int alphaTree_t::alphaParent2Branch(int edgeId, int alphaParentIdx)
{
    if(alphaParentIdx<0) return 0;
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
    if(edgeId==0) return 0;
	return alphaParent2Branch(edgeId, alphaParents[alphaEdgeId]);
}

std::vector<int> alphaTree_t::computeFlatClustering(int minClusterSize)
{

    std::vector<int> flatClustering;
	m_numBranches = 2 * numAlphaEdges + 1;
	// std::vector<int> brHeads(m_numBranches, -1);
    brIHeads.resize(m_numBranches, -1);
	brDescendents.resize(m_numBranches, 0);
    

	brStabilityScores.resize(m_numBranches, 0.0);
	brSHatScores.resize(m_numBranches, 0.0);
	brDelta.resize(m_numBranches, 1);
	brLength.resize(m_numBranches, 0);

	// first compute branch edge
	std::vector<branchEdge_t> branchEdge = computeBranchEdge();
	// computing chain lengths
	// why is it starting from 1?
	for (int iEdge = 0; iEdge < m_npts - 1; iEdge++)
	{
		// branchEdge is not sorted so iEdge is not equal to
		// edgeIdx of MST
		int myBranch = branchEdge[iEdge].branch;
		int edgeId = branchEdge[iEdge].edgeIdx;
		int myAlParent = branch2AlphaEdge(myBranch);

		// TODO: edgeId is leaf
        // if (m_incMatMST.isAlphaEdge(edgeId))
		if (m_incMatMST.isLeafEdge(edgeId))
		{
			// alphaDescendents[myAlParent]+=2;
			// brDescendents[myBranch] += 2;
            std::cout<<edgeId<<" Adding 2 to branch "<<myBranch<< "  len ="<<brLength[myBranch]<< "\n";
			brLength[myBranch] += 2;
		}
		else if (!m_incMatMST.isAlphaEdge(edgeId))
		{
			// alphaDescendents[myAlParent] +=1;
			// brDescendents[myBranch] += 1;
            std::cout<<edgeId<<" Adding 1 to branch "<<myBranch<<"  len ="<<brLength[myBranch]<< "\n";
			brLength[myBranch] += 1;
		}

		/*Adding stability score to branch */
		if (myAlParent != -1)
		{
			// int parent = alphaEdges[myAlParent];
			brStabilityScores[myBranch] += 1.0 / m_wtSortedMST[edgeId].second;
		}

		/*
		 *   branch head calculation
		 *   TODO: is it needed?
		 *  Its needed for mapping flat cluster to original one
		 */

		if (iEdge == 0 || branchEdge[iEdge].branch != branchEdge[iEdge - 1].branch)
        {
            // brHeads[myBranch] = branchEdge[iEdge].edgeIdx;
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
            std::cout<<" i ="<<i<<"  breanch Head "<<getBranchHead(i)<<" "<<getBranchHead(branchParent(i))<<
            "  len "<<brLength[i] <<"\n";
        
			numBranchChildrens[branchParent(i)]++;
        }
        
    }
		

	// bottom-up traversal
	int numBranch2Process = m_numBranches;

	while (numBranch2Process)
	{
		for (int i = 0; i < m_numBranches; i++)
		{
            // std::cout<<i<<"   "<<qBranchProcessed[i]<<" "<<numBranchChildrens[i]<<" "<<numChildBranchProcessed[i]<<"\n";
            
			if (qBranchProcessed[i] == 0 and
					numBranchChildrens[i] == numChildBranchProcessed[i])
			{
                std::cout<<i<<" -> Parent ="<< branchParent(i) <<"\n";
				if (branchParent(i) != -1)
				{
                    std::cout<<i<<"   "<<qBranchProcessed[i]<<" "<<numBranchChildrens[i]<<" "<<numChildBranchProcessed[i]<<"\n";
                    // update descendents
                    int totalDescendents = brDescendents[i] + brLength[i];
					brDescendents[branchParent(i)] += (totalDescendents);
                    std::cout<<" Descend adding "<< i << " ndesc="<<brDescendents[branchParent(i)]<<" to "<< branchParent(i)<<"\n";
                    // so there are three cases 
                    // 1. either you are invalid
                    // or 2. your number of descdends are greater than mincluster
                    // or 3. it becomes a valid cluster at this point 
                    // the third case is the difficult one
					#if 0
					if(totalDescendents <= minClusterSize   )
                    {
                        // invalid cluster 
                        brStabilityScores[i] =0; 
                        // sHat and delta are don't care 
                        brDelta[i] =0; 
                        brSHatScores[i] =0; 
                    }
                    else
                    {
                        // update number of valid childrens 
                        numValidChildrens[branchParent(i)]++;  
                        // let selfContri = chainContri[i] + brDescendents[i]/myLength
                        if(brDescendents[i] < minClusterSize)
                        {
                            // chainContri[i] =?? use different formula; 
                            
                            for( int etx=0; etx< totalDescendents-minClusterSize; etx++)
                            {
                                int branchStartEdge = brIHeads[i];
                                int edgeId = branchEdge[branchStartEdge+etx].edgeIdx;
                                chainContri[i] += 1.0 / m_wtSortedMST[edgeId].second;
                            }
                            int etxFinal = totalDescendents-minClusterSize;
                        }
                        
                        brDelta[i] =1;
                            // very valid cluster 
                        if(numValidChildrens[i]==0)
                        {
                            brStabilityScores[i] = brSelfContri[i] - ancestorContriution;
                            brSHatScores[i] =brStabilityScores[i]; 
                            
                        }
                        else if(numValidChildrens[i]==1)
                        {
                            brStabilityScores[i] = brStabilityScores[i]+ brSelfContri[i] - ancestorContriution;
                            brSHatScores[i] =brStabilityScores[i]; 
                            
                        }
                        else /* two valid childrens */
                        {
                            brStabilityScores[i] =  brSelfContri[i] - ancestorContriution;
                            
                            if(brSHatScores[i]> brStabilityScores[i])
                            {
                                brDelta[i] =0; // not a final cluster 
                            }
                            else 
                            {
                                brSHatScores[i]= brStabilityScores[i];
                            }
                            
                        }

                        brStabilityScores[branchParent(i)] += brStabilityScores[i];
                        brSHatScores[branchParent(i)] += brSHatScores[i];
                        
                        #endif 
                            // very valid cluster 
                        numChildBranchProcessed[branchParent(i)]++;
                    } 
					

					
                    qBranchProcessed[i] = 1;
                    numBranch2Process--;
                    // std::cout<<" Left to process "<<numBranch2Process<<"\n";
            }
            
        }
    }

    //TODO: dummy return ;
    return flatClustering;
}

