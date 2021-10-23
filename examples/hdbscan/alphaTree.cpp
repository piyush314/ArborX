#include <algorithm>
#include <iostream>
#include <cassert> 
#include <sstream>
#include <cstdio> 
#include "edge_hierarchy.hpp"
#include "incidence_matrix.hpp"
#include "alphaTree.hpp"

std::vector<int> alphaTree_t::getAlphaEdges()
{
    std::vector<int> splitEdgeList;
    // initialize outsets 
    std::vector<std::vector<int>> outSet(2*m_npts-1);
    std::vector<int>  nChildrenProcessed(m_npts-1, 0);
    std::vector<std::vector<int>>  childrenList(m_npts-1, std::vector<int>(2,-1)); 
    std::vector<int>  isEdgeProcessed(m_npts-1, 0);
    // std::vector<int> vertexMaxIncidentEdgeId = m_incMatMST.maxIncidentEdgeId()

    for(int vtxId =0; vtxId< m_npts; vtxId++)
    {
        int parent = m_incMatMST.maxIncidentEdgeId(vtxId); 
        for (int k = 0; k < m_incMatMST.numIncidentEdge(vtxId); k++)
            if(m_incMatMST.k_thIncidentEdge(vtxId,k) != parent)
                outSet[m_npts-1 + vtxId].push_back(m_incMatMST.k_thIncidentEdge(vtxId,k));
        
        
        childrenList[parent][nChildrenProcessed[parent]] = m_npts-1 + vtxId;     

        nChildrenProcessed[parent]++;
        // m_numDescendents[parent]++; 
    }
    int nSplitEdges=0;
    int nLeafEdges=0;
    std::vector<int> isSplitEdge(m_npts-1,0);
    for ( int edgeId = 0; edgeId < m_npts-1; edgeId++)
    {
        if(nChildrenProcessed[edgeId]==0)
        {
            nSplitEdges++;
            isSplitEdge[edgeId]=1;
            splitEdgeList.push_back(edgeId);
        }
            
        if(nChildrenProcessed[edgeId]==2)
            nLeafEdges++;
    }

    return splitEdgeList;
}



void alphaTree_t::alphaEulerTour()
{
    idxTS.resize(2*numAlphaEdges);

    for(int i=0; i<numAlphaEdges; i++)
    {
        int iam = alphaEdges[i];
        idxTS[2*i].first = mstEulerEntry[iam];
        idxTS[2*i].second = 2*i;
        idxTS[2*i+1].first = mstEulerExit[iam];
        idxTS[2*i+1].second = 2*i+1; 
    }
    
    //sorting the indexedTimeStamp array 
    std::sort(idxTS.begin(), idxTS.end());
    alphaEulerEntry.resize(numAlphaEdges);
    alphaEulerExit.resize(numAlphaEdges);

    for(int i=0; i<2*numAlphaEdges; i++)
    {
        int timeStamp = idxTS[i].first;
        int iam = idxTS[i].second/2;
        if(idxTS[i].second%2==0)
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

alphaTree_t::alphaTree_t(incidenceMatrix_t& _incMatMST, std::vector<wtEdge_t>& _wtSortedMST) :
        m_incMatMST(_incMatMST),
        m_wtSortedMST(_wtSortedMST)
{
    m_npts = m_incMatMST.m_npts; 
    // do the eulerTour
    eulerTour(); 
    // get the alpha edges 
    alphaEdges = getAlphaEdges();
    numAlphaEdges = alphaEdges.size();
    
    
    // perform the eulerTour on Alpha-MST
    // sets idxTS, alphaEulerEntry and alphaEulerExit
    alphaEulerTour(); 
    // compute bridgeEdge
    alphaBridgeEdges= computeBrideEdges();
    // construct Alpha tree 
    alphaParents = constructAlphaTree();

}



int alphaTree_t::isLeftOfEdge(int queryEdge, int rootEdge)
{
    if( mstEulerEntry[queryEdge] > mstEulerEntry[rootEdge] &&
        mstEulerExit[queryEdge] < mstEulerExit[rootEdge])
    return 1;

    return 0; 
}


// void edgeHierarchy_t::constructEdgeTreeSplitEdges()
std::vector<int> alphaTree_t::constructAlphaTree()
{

    //  = listOfSplitEdge();
    int numSplitEdges = alphaEdges.size();
    
    std::vector<int> parent(numSplitEdges,0);
    std::vector<int> isProcessed(numSplitEdges,0);
    std::vector<int> myBranch(numSplitEdges,0);
    std::vector<int> maxVertexMyBranch(2*m_npts,m_npts);
    std::vector<int> maxVertexAlphaId(2*m_npts,-1);
    isProcessed[0] =1;
    int numNotProcessed = numSplitEdges-1; 
    int height=0;

    parent[0] =-1;
    while(numNotProcessed>0)
    {
        for(int edge=1; edge<numSplitEdges; edge++ )
        {
            if(!isProcessed[edge])
            {
                int iam = alphaEdges[edge];
                
                myBranch[edge] = 2*parent[edge];

                if(isLeftOfEdge(iam, alphaEdges[parent[edge]]))
                    myBranch[edge]++; 

                if(maxVertexMyBranch[myBranch[edge]]>iam)
                {
                    maxVertexMyBranch[myBranch[edge]] = iam;
                    maxVertexAlphaId[myBranch[edge]] = edge;     
                }
                // maxVertexMyBranch[myBranch[edge]] = std::min(maxVertexMyBranch[myBranch[edge]],iam );
                // printf(" %d : \t { %d } --> \t %d \n ", iam, myBranch[edge], maxVertexMyBranch[myBranch[edge]]);
            }
        }

        for(int edge=0; edge<numSplitEdges; edge++ )
        {
            if(!isProcessed[edge])
            {
                int iam = alphaEdges[edge];
                if(maxVertexMyBranch[myBranch[edge]] == iam)
                {
                    isProcessed[edge]=1; 
                    numNotProcessed--;
                    // std::cout<<iam<<"  \t -> \t "<< parent[edge] << "\n";   
                    std::cout<<iam<<"  \t -> \t "<< alphaEdges[parent[edge]] << "\n";
                }
                else 
                {
                    // parent[edge] = maxVertexMyBranch[myBranch[edge]];
                    parent[edge] =maxVertexAlphaId[myBranch[edge]];
                }
            }
        }

        height++; 
    }
    
    std::cout<< "Height of the split edge tree is \t:"<<height << "\n";
    return parent; 
}

// void alphaTree_t::alpha2EdgeTree(std::vector<int>& alphaEdges, std::vector<int>& alphaParent)
std::vector<int> alphaTree_t::constructEdgeTree()
{
    //Assumption no edges have been processed
    std::vector<branchEdge_t> branchEdge(m_npts-1);
    // for each edge find its branch 
    
    for(int edgeId =0; edgeId< m_npts-1; edgeId++) 
    // for(int edgeId =0; edgeId< 10; edgeId++) 
    {
        int alParent = findAlphaParent(edgeId );
        if(alParent!=-1)
        printf(" alphaParent(%d)  :\t %d (%d) \n", edgeId, alphaEdges[alParent], alParent);
        //TODO: -1 case 
        if(alParent==-1)
            branchEdge[edgeId] = branchEdge_t(-1, edgeId);     
        else 
        {
            // int mstParent = alphaParents[alParent];
            int mstParent = alphaEdges[alParent];
            int myBranch = 2*alParent; 
            myBranch += isLeftOfEdge(edgeId, mstParent);
            branchEdge[edgeId] = branchEdge_t(myBranch, edgeId); 
        }
    }
    // exit(0);
    // TODO: branchEdgeComparator 
    std::sort(branchEdge.begin(),branchEdge.end(), branchEdgeComparator());

    // mark the parents 
    std::vector<int> globalParents(m_npts-1,-1);
    for(int edgeId =1; edgeId< m_npts-1; edgeId++) 
    {
        // we are on same branch
        if(branchEdge[edgeId-1].branch == branchEdge[edgeId].branch)
        {
            globalParents[branchEdge[edgeId].edgeIdx] = 
            branchEdge[edgeId-1].edgeIdx;
            printf(" %d -> %d \n",branchEdge[edgeId].edgeIdx,branchEdge[edgeId-1].edgeIdx);
        }
        else // Iam start of the branch
        {
            globalParents[branchEdge[edgeId].edgeIdx] = alphaEdges[branchEdge[edgeId].branch/2];
            printf(" %d -> %d \n",branchEdge[edgeId].edgeIdx,globalParents[branchEdge[edgeId].edgeIdx]);
        }
    }
    return globalParents;
}

std::vector<int> alphaTree_t::computeBrideEdges()
{
    std::vector<int> bridgeEdges(numAlphaEdges,-1);

    for(int alphaId =0; alphaId <numAlphaEdges; alphaId++ )
    {
        int myEntry = alphaEulerEntry[alphaId];
        int myExit = alphaEulerExit[alphaId];

        int timeStamp = myEntry+1;

        while(timeStamp < myExit)
        {
            int neighbourIdx= idxTS[timeStamp].second/2;
            bridgeEdges[neighbourIdx] = alphaId; 

            timeStamp = alphaEulerExit[neighbourIdx]+1; 
        }
    }

    return bridgeEdges;
    
}


int alphaTree_t::findAlphaParent(int edgeId)
{
    indexedTimeStamp_t idxEntry(mstEulerEntry[edgeId], 0); 
    indexedTimeStamp_t idxExit(mstEulerExit[edgeId], 0);
    int myEntry = std::distance(idxTS.begin(), 
                        std::lower_bound(idxTS.begin(), idxTS.end(), idxEntry)); 
    int myExit  = std::distance(idxTS.begin(), 
                        std::lower_bound(idxTS.begin(), idxTS.end(), idxExit));

    
    // quick return
    if(edgeId == alphaEdges[idxTS[myEntry].second/2]  )
    {
        return alphaParents[idxTS[myEntry].second/2];
    }

    printf("%d :\t mstEntry %d mstExit %d newEntry %d newExit %d: \
    \n \t\t%d  %d %d %d\n",
    edgeId, mstEulerEntry[edgeId], mstEulerExit[edgeId],
         myEntry, myExit
         , alphaEdges[idxTS[myEntry].second/2]
         , idxTS[myEntry].second
         , alphaEdges[idxTS[myEntry-1].second/2]
         , alphaBridgeEdges[idxTS[myEntry-1].second/2]   );

    int myBridge; 


#if 0
    if(myEntry==0)
        myBridge =-1;
    else if(idxTS[myEntry-1].second %2 ==0 ) //starting point
        myBridge = idxTS[myEntry-1].second/2;
    else 
        myBridge = alphaBridgeEdges[idxTS[myEntry-1].second/2];
#else 
    if(myEntry==0)
        myBridge =-1;
    else if(idxTS[myEntry].second %2 ==0 ) //starting point
        myBridge = alphaBridgeEdges[idxTS[myEntry].second/2];
    else 
        myBridge = idxTS[myEntry].second/2;
#endif 
    int minDescendent= m_npts;
    // int minDescendent= -1;
    int minDescIdx=-1;
    // int maxAncestor = m_npts;
    int maxAncestor = -1;
    int maxAncIdx=-1;

    int startTime, endTime, minDesc;

    if(myBridge ==-1)
    {
        startTime=0; 
        endTime = 2*numAlphaEdges; 
    }
    else
    {
        startTime=alphaEulerEntry[myBridge]+1; 
        endTime = alphaEulerExit[myBridge];

        int bridgeId = alphaEdges[myBridge];
        if(bridgeId  >edgeId and 
            bridgeId<minDescendent )
        {   
            minDescendent =     bridgeId;
            minDescIdx = myBridge; 
        }

        if(bridgeId<edgeId and 
            bridgeId>maxAncestor)
        {   
            maxAncestor =     bridgeId;
            maxAncIdx = myBridge; 
        }

    } 


    
    
    
    int ts= startTime;
    printf("%d :\t bridge %d startTime %d endTime %d max %d min %d\n ", edgeId, myBridge,
         startTime, endTime, maxAncestor, minDescendent);
    while(ts<endTime)
    {
        int neighbourIdx= idxTS[ts].second/2;
        int neighbourEdgeId = alphaEdges[neighbourIdx];
        printf("%d : \t checking %d\n", edgeId, neighbourEdgeId);
        if(neighbourEdgeId>edgeId and 
            neighbourEdgeId<minDescendent)
        {   
            minDescendent =     neighbourEdgeId;
            minDescIdx = neighbourIdx; 
        }

        if(neighbourEdgeId<edgeId and 
            neighbourEdgeId>maxAncestor )
        {   
            maxAncestor =     neighbourEdgeId;
            maxAncIdx = neighbourIdx; 
        }

        ts = alphaEulerExit[neighbourIdx]+1; 
    }

    printf(":\t min %d max %d \n ", minDescendent, maxAncestor);

    if(minDescIdx==-1) // leaf node
        return maxAncIdx;

    //
    int myParent = minDescendent;
    int parentIdx=minDescIdx;

    while(myParent>edgeId)
    {
        parentIdx = alphaParents[parentIdx];
        myParent = alphaEdges[parentIdx];
    }

    return parentIdx; 
}


#if 0
std::vector<int> alphaTree_t::getAlphaEdges(int maxLevel)
{
    // TODO: do this later std::vector<int>  nChildrenProcessed(m_npts-1, 0);
    std::vector<int> splitEdgeList;
    // initialize outsets 
    std::vector<std::vector<int>> outSet(2*m_npts-1);
    std::vector<int>  nChildrenProcessed(m_npts-1, 0);
    std::vector<std::vector<int>>  childrenList(m_npts-1, std::vector<int>(2,-1)); 
    std::vector<int>  isEdgeProcessed(m_npts-1, 0);
    // initialize outsets 

    for(int vtxId =0; vtxId< m_npts; vtxId++)
    {
        int parent = m_vertexMaxIncidentEdgeId[vtxId];
        for (int k = 0; k < m_incMatMST.numIncidentEdge(vtxId); k++)
            if(m_incMatMST.k_thIncidentEdge(vtxId,k) != parent)
                outSet[m_npts-1 + vtxId].push_back(m_incMatMST.k_thIncidentEdge(vtxId,k));
        
        
        childrenList[parent][nChildrenProcessed[parent]] = m_npts-1 + vtxId;     

        nChildrenProcessed[parent]++;
        // m_numDescendents[parent]++; 
    }
    int nSplitEdges=0;
    int nLeafEdges=0;
    std::vector<int> isSplitEdge(m_npts-1,0);
    for ( int edgeId = 0; edgeId < m_npts-1; edgeId++)
    {
        if(nChildrenProcessed[edgeId]==0)
        {
            nSplitEdges++;
            isSplitEdge[edgeId]=1;
        }
            
        if(nChildrenProcessed[edgeId]==2)
            nLeafEdges++;
    }
        
    // std::cout<<"Number of split edges ="<<nSplitEdges<<" as percentage  "<<100.0*nSplitEdges/(m_npts-1) << "\n";
    int level=0; 
    std::cout<<"| Level \t| #LeafEdges \t| #SplitEdges \t|" << "\n";
    
    int edgesNotProcessed = m_npts-1;
    while(edgesNotProcessed>1)      // no need to process vertex zero
    {
        nLeafEdges=0;
        
        for ( int edgeId = 0; edgeId < m_npts-1; edgeId++)
        {
            /* code */
            if(!isEdgeProcessed[edgeId]  // edge is not processed
                && nChildrenProcessed[edgeId]==2)        // both its childrens are processed
            {
                nLeafEdges++; 
                if(isSplitEdge[edgeId])
                    nSplitEdges--;  
                // std::cout<<"Processing " <<edgeId<<"\n";
                int left = childrenList[edgeId][0];
                int right = childrenList[edgeId][1];
                // main compute step 
                outSet[edgeId] = mergeOutSets(outSet[left], outSet[right]);
                auto parentIt = std::max_element(outSet[edgeId].begin(), outSet[edgeId].end());

                // up
                m_edgeParent[edgeId] = *parentIt;
                // std::cout<<"Parent is "<< m_edgeParent[edgeId]<<"\n";
                int outsetSize = outSet[edgeId].size(); 
                if(parentIt != outSet[edgeId].end()-1)
                    *parentIt =outSet[edgeId][outsetSize-1]; // set to last element 
                outSet[edgeId].resize(outsetSize-1);

                int parent = m_edgeParent[edgeId];
                childrenList[parent][nChildrenProcessed[parent]] = edgeId; 
                nChildrenProcessed[parent]++;
                m_numDescendents[parent] += m_numDescendents[edgeId];
                if (isValidCluster(edgeId))
                    m_numChildCluster[parent] += 1;
                
                isEdgeProcessed[edgeId] =1; 
                // reduce index 
                edgesNotProcessed--;
            }
        }
        if(level%10==0)
        std::cout<<"|"<<level<<"\t|"<< nLeafEdges<< "\t|" << nSplitEdges <<"\t|" << "\n";
        level++; 
        // std::cout<<"edgesNotProcessed" << edgesNotProcessed<<"\n"; 
        // loopCount++;
        // if(loopCount==100) exit(0);
        if(level==maxLevel) break;
    }

    // int nSplitEdges=0;
    // int nLeafEdges=0;
    // std::vector<int> isSplitEdge(m_npts-1,0);
    for ( int edgeId = 0; edgeId < m_npts-1; edgeId++)
    {
        if(nChildrenProcessed[edgeId]==0)
        {
            // nSplitEdges++;
            // isSplitEdge[edgeId]=1;
            splitEdgeList.push_back(edgeId);
            printf("Adding %d\n", edgeId);
        }
    }

    return splitEdgeList;
}
#endif 