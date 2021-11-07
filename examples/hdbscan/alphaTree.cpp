#include <algorithm>
#include <iostream>
#include <cassert> 
#include <sstream>
#include <cstdio> 
#include "edge_hierarchy.hpp"
#include "incidence_matrix.hpp"
#include "alphaTree.hpp"
#include "timer.hpp"

int printEnabled =0;

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
    timer_t tAlphaEulerTour(std::string("Alpha-euler tour"));
    timer_t tAlphaTree(std::string("Alpha-tree construction"));

    tAlphaEulerTour.start();
    alphaEulerTour(); 
    tAlphaEulerTour.end();
    // compute bridgeEdge
    alphaBridgeEdges= computeBrideEdges();
    // construct Alpha tree 
    tAlphaTree.start(); 
    alphaParents = constructAlphaTree();
    tAlphaTree.end();

    tAlphaEulerTour.print();
    tAlphaTree.print();

}



int alphaTree_t::isLeftOfEdge(int queryEdge, int rootEdge)
{
    if( mstEulerEntry[queryEdge] > mstEulerEntry[rootEdge] &&
        mstEulerExit[queryEdge] < mstEulerExit[rootEdge])
    return 1;

    return 0; 
}



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


std::vector<int> alphaTree_t::constructEdgeTree()
{
    int printEnabled=0;
    //Assumption no edges have been processed
    #if 1
    std::vector<branchEdge_t> branchEdge = computeBranchEdge();
    #else 

    std::vector<branchEdge_t> branchEdge(m_npts-1);
    // for each edge find its branch 
    timer_t tFindAlphaParent("Finding alpha parents");
    tFindAlphaParent.start();
    for(int edgeId =0; edgeId< m_npts-1; edgeId++) 
    // for(int edgeId =0; edgeId< 10; edgeId++) 
    {
        int alParent = findAlphaParent(edgeId );
        if(alParent!=-1 and printEnabled)
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
    tFindAlphaParent.end();
    tFindAlphaParent.print();
    // exit(0);
    // TODO: branchEdgeComparator 
    timer_t tSortBranchEdge("Branch Edge-sorting");
    tSortBranchEdge.start();
    std::sort(branchEdge.begin(),branchEdge.end(), branchEdgeComparator());
    tSortBranchEdge.end();
    tSortBranchEdge.print();
    #endif 
    // mark the parents 
    std::vector<int> globalParents(m_npts-1,-1);
    timer_t tMarkingParents("Marking Parents");
    tMarkingParents.start();
    for(int edgeId =1; edgeId< m_npts-1; edgeId++) 
    {
        // we are on same branch
        if(branchEdge[edgeId-1].branch == branchEdge[edgeId].branch)
        {
            globalParents[branchEdge[edgeId].edgeIdx] = 
            branchEdge[edgeId-1].edgeIdx;
            if(printEnabled)
            printf(" %d -> %d \n",branchEdge[edgeId].edgeIdx,branchEdge[edgeId-1].edgeIdx);
        }
        else // Iam start of the branch
        {
            globalParents[branchEdge[edgeId].edgeIdx] = alphaEdges[branchEdge[edgeId].branch/2];
            if(printEnabled)
            printf(" %d -> %d \n",branchEdge[edgeId].edgeIdx,globalParents[branchEdge[edgeId].edgeIdx]);
        }
    }
    tMarkingParents.end();
    tMarkingParents.print();

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
    
    if(edgeId <alphaEdges[0]) return -1;
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
    
    if(printEnabled)
    printf("%d :\t mstEntry %d mstExit %d newEntry %d newExit %d: \
    \n \t\t%d  %d\n",
    edgeId, mstEulerEntry[edgeId], mstEulerExit[edgeId],
         myEntry, myExit
         , alphaEdges[idxTS[myEntry].second/2]
         , idxTS[myEntry].second
           );

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
    int minDescIdx=-1;
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
    if(printEnabled)
    printf("%d :\t bridge %d startTime %d endTime %d max %d min %d\n ", edgeId, myBridge,
         startTime, endTime, maxAncestor, minDescendent);
    while(ts<endTime)
    {
        int neighbourIdx= idxTS[ts].second/2;
        int neighbourEdgeId = alphaEdges[neighbourIdx];
        if(printEnabled)
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
    if(printEnabled)
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


std::vector<branchEdge_t> alphaTree_t::computeBranchEdge()
{
    std::vector<branchEdge_t> branchEdge(m_npts-1);
    // for each edge find its branch 
    timer_t tFindAlphaParent("Finding alpha parents");
    tFindAlphaParent.start();
    for(int edgeId =0; edgeId< m_npts-1; edgeId++) 
    // for(int edgeId =0; edgeId< 10; edgeId++) 
    {
        int alParent = findAlphaParent(edgeId );
        if(alParent!=-1 and printEnabled)
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
    tFindAlphaParent.end();
    tFindAlphaParent.print();
    // exit(0);
    // TODO: branchEdgeComparator 
    timer_t tSortBranchEdge("Branch Edge-sorting");
    tSortBranchEdge.start();
    std::sort(branchEdge.begin(),branchEdge.end(), branchEdgeComparator());
    tSortBranchEdge.end();
    tSortBranchEdge.print();

    return branchEdge;
}
// std::vector<int> alphaTree_t::computeFlatClustering(int minClusterSize)
// {
//     std::vector<int> chainHeads;
//     std::vector<int> chainLengths;
//     std::vector<int> chainDescendents;
//     std::vector<int> chainStabilityScores;
//     std::vector<int> chainDelta;

//     // computing chain heads 


// }


