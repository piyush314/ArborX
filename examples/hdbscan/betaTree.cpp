#include <iostream>
#include <cassert>
#include "recursiveAlphaTree.hpp"
#include "timer.hpp"


int printBeta=0;
int CheckDandrogramAlphaTree =0;

rAlphaTree_t::rAlphaTree_t (eulerTourTree_t& inTour, int level) :  rlevel(level)
{
    // Compute alpha edges 
    mytimer_t tComputeAlphaEdges(std::string("Computing AlphaEdges-")+std::to_string(level));
    tComputeAlphaEdges.start();
    alphaEdges= inTour.computeAlphaEdges();
    invAlphaEdges.resize(inTour.nEdges(), -1);
    for(int idx=0; idx<alphaEdges.size(); idx++)
        invAlphaEdges[alphaEdges[idx]] = idx;
    tComputeAlphaEdges.end();
    tComputeAlphaEdges.print();
    // compute the alpha MST by computing the Euler Tour of alpha-MST
    mytimer_t tAlphaEulerTour(std::string("Alpha-euler tour-")+std::to_string(level));
    tAlphaEulerTour.start();
    eAlphaTour = computeAlphaEulerTour(alphaEdges, inTour);
    tAlphaEulerTour.end();
    tAlphaEulerTour.print();

    #define MAX_ALPHA_RECURSION_LEVEL 2
    if(rlevel<MAX_ALPHA_RECURSION_LEVEL)
    {
        // compute using alpha tree 
        mytimer_t tBetaTree(std::string("BetaTree-")+std::to_string(level));
        tBetaTree.start();
        betaTree = new rAlphaTree_t (eAlphaTour, rlevel+1);
        tBetaTree.end();
        tBetaTree.print();
        // compute dendrogram 
        // alphaParents =computeDandrogramUsingAlphaTree(inTour, *betaTree);
        mytimer_t tDandrogram(std::string("AlphaTreeDandrogram using Betatree-")+std::to_string(level));
        tDandrogram.start();
        alphaParents =computeDandrogramUsingAlphaTree(eAlphaTour, *betaTree);
        tDandrogram.end();
        tDandrogram.print();
        // TODO: remove the following checks 
        if(CheckDandrogramAlphaTree)
        {
            auto alternative = computeDandrogramTopDown(eAlphaTour);
            // check alternative and alphaParents are same
            assert(alternative.size()==alphaParents.size());
            int count =0; 
            for(int i=0;i<alternative.size();i++)
            {
                if(alternative[i]!=alphaParents[i])
                {
                    std::cout<<"Error in location "<< i<<" \t";
                    std::cout<<"alternative[i]"<< alternative[i]<<" \t";
                    std::cout<<"alphaParents[i]"<< alphaParents[i]<<std::endl;

                    // check if i in *betaTree.alphaEdges
                    auto it = std::find(betaTree->alphaEdges.begin(), betaTree->alphaEdges.end(), i);
                    if(it!=betaTree->alphaEdges.end())
                    {
                        std::cout<<"i is in betaTree->alphaEdges"<<std::endl;
                    }
                    count++;
                }
                    // exit(1);
                // }
                // assert(alternative[i]==alphaParents[i]);
            }
            if(count>0)
            {
                std::cout<<"Error in "<<count<<" locations"<<std::endl;
                exit(1);
            }
        }
        
        
    }
    else // Control the recursion 
    {
        // alphaParents = computeDandrogramTopDown(inTour);
        mytimer_t tDandrogram(std::string("AlphaTreeDandrogram Topdown-")+std::to_string(level));
        tDandrogram.start();
        alphaParents= computeDandrogramTopDown(eAlphaTour);
        tDandrogram.end();
        tDandrogram.print();
    }
        
    mytimer_t tLineageShortcut(std::string("LineageShortcut-")+std::to_string(level));
    tLineageShortcut.start();
    alphaAncestors = computeLineageShortcut(alphaParents);
    tLineageShortcut.end();
    tLineageShortcut.print();

}


std::vector<int> computeDandrogramUsingAlphaTree(eulerTourTree_t& inTour, rAlphaTree_t& alphaTree)
{
    int numEdges = inTour.nEdges();
    // compute Alpha tree
    // rAlphaTree_t   alphaTree(inTour, level+1);
    // For each edge in inTour, compute its branch 
    std::vector<branchEdge_t> branchEdge(numEdges);

    mytimer_t tBranchEdge(std::string("BranchEdge-")+std::to_string(alphaTree.rlevel));
    tBranchEdge.start();
    for (int edgeId = 0; edgeId < numEdges; edgeId++)
        branchEdge[edgeId] = branchEdge_t( alphaTree.computeBranch(edgeId,inTour), edgeId);
    tBranchEdge.end();
    tBranchEdge.print();
    mytimer_t tSortBranchEdge(std::string("SortBranchEdge-")+std::to_string(alphaTree.rlevel));
    tSortBranchEdge.start();
    std::sort(branchEdge.begin(), branchEdge.end(), branchEdgeComparator());
    tSortBranchEdge.end();
    tSortBranchEdge.print();
    

    // Compute the Dandrogram from the branch
    std::vector<int> parent(numEdges, -1);
    mytimer_t tMarkParents(std::string("Marking Parents-")+std::to_string(alphaTree.rlevel));
    tMarkParents.start();
    for (int beIdx = 1; beIdx < numEdges; beIdx++)
    {
        int myBranch = branchEdge[beIdx].branch;
        int myEdgeId = branchEdge[beIdx].edgeIdx;
        // if left branchEdge and Iam  on same branch
        if (branchEdge[beIdx - 1].branch == myBranch)
            parent[myEdgeId] = branchEdge[beIdx - 1].edgeIdx;
        else // If not then Iam the head (i.e.) start of the branch
            parent[myEdgeId] = alphaTree.alphaEdges[branch2AlphaEdge(myBranch)];

        
    }
    tMarkParents.end();
    tMarkParents.print();
    
    return parent; 


}

std::vector<int> computeDandrogramTopDown(eulerTourTree_t& inTour);


// computeBranch: find which branch of alphaTree inputEdge \in parentGraph belongs
int rAlphaTree_t::computeBranch(int edgeId, eulerTourTree_t& pMST)
{
    int alParent = findAlphaParent(edgeId, pMST);
    if (alParent == -1)
        return 0;
    
    int mstParent = alphaEdges[alParent];
    int myBranch = 2 * alParent + 1;
    myBranch += pMST.isLeftOfEdge(edgeId, mstParent);
    return myBranch;
    

}


int rAlphaTree_t::findAlphaParent(int edgeId, eulerTourTree_t& pMST)
{

    /*
    *        QUICK return 
    */
    if (edgeId < alphaEdges[0])
    {
        // std::cout<<edgeId<<":\tQuick return \n";
        return -1;
    }
        
    
    /**
     *  Simulating edge insertion into alphaMST
     *  To do so: get the entry and exit time in the alphaEulerTour
     **/
    // indexedTimeStamp_t idxEntry(pMST.entry(edgeId), 0);
    // indexedTimeStamp_t idxExit(pMST.exit(edgeId), 0);
    // int myEntry = std::distance(
    //     idxTS.begin(), std::lower_bound(idxTS.begin(), idxTS.end(), idxEntry));
    // int myExit = std::distance(
    //     idxTS.begin(), std::lower_bound(idxTS.begin(), idxTS.end(), idxExit));
    
    // Check if the edge is in the alphaMST
    if(invAlphaEdges[edgeId]!=-1)
    {
        // std::cout<<edgeId<<":\tEdge is in alphaMST \n";
        return alphaParents[invAlphaEdges[edgeId]];
    }

    /* What are the assumptions here!? alo pMST.idxTs added later not checked*/ 
    // TODO: The following checks if the edgeId is an alpha edge or not
    // if it is, the parent is already computed, and we can return it
    // So using an array that store alphaIdx[i]= eAlphaTour.idxTS[myEntry].second / 2
    // if i is an alpha edge, can solve it in O(1) time
    // commenting for time time being because simInsert has a bug 
    int myBridge;
    if(0)
    {
        auto myEntryExit = eAlphaTour.simInsert(pMST.entry(edgeId), pMST.exit(edgeId));
        int myEntry = myEntryExit.first;
        int myExit = myEntryExit.second;
        if (edgeId == alphaEdges[eAlphaTour.idxTS[myEntry].second / 2])
        {
            std::cout<<edgeId<<":\tEdge is in alphaMST \n";
            return alphaParents[eAlphaTour.idxTS[myEntry].second / 2];
        }
        myBridge = eAlphaTour.simBridge(myEntry);   
    }
    else 
    {
        myBridge = eAlphaTour.computeSimulatedBridge(pMST.entry(edgeId), pMST.exit(edgeId));   
    }

    // int myBridge ;
    // if (myEntry == 0)
    //     myBridge = -1;
    // else if (idxTS[myEntry].second % 2 == 0) // starting point
    //     myBridge = alphaBridgeEdges[idxTS[myEntry].second / 2];
    // else
    //     myBridge = idxTS[myEntry].second / 2;
    
    // int myBridge = eAlphaTour.simBridge(myEntry);
    #if 0
    int minDescendent = pMST.nEdges()+1;
    int minDescIdx = -1;
    int maxAncestor = -1;
    int maxAncIdx = -1;

    // maxAncestor and minDescdent
    

    // int startTime, endTime;
    int startTime = 0;
    int endTime = 2 * numAlphaEdges;
    
    if(myBridge != -1)
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
    // if (printEnabled)
    //     printf("%d :\t bridge %d startTime %d endTime %d max %d min %d\n ", edgeId,
    //         myBridge, startTime, endTime, maxAncestor, minDescendent);
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

    #endif 

    #define MERGE_BRIDGE_PARENT_SEARCH
    #ifdef MERGE_BRIDGE_PARENT_SEARCH
        // If it is a leaf then
        if(alphaEdges[eAlphaTour.maxIdxBridge(myBridge)] < edgeId)
        {
            // std::cout<<edgeId<<":\tLeaf \n";
            return eAlphaTour.maxIdxBridge(myBridge);
        }
            

        return ::alphaTreeBottomUpParentSearch(edgeId, eAlphaTour.maxIdxBridge(myBridge), this->alphaEdges, this->alphaAncestors);
    #else 
    std::pair<int,int> maxAncMinDesc = computeMaxAncestorMinDescendent(myBridge, edgeId);
    auto maxAncIdx = maxAncMinDesc.first;   
    auto minDescIdx = maxAncMinDesc.second;   
    if (minDescIdx == -1) // leaf node
        return maxAncIdx;

    /*
    * The following can be done efficiently using extra storage to store 
    * alphaParent[i], alphaParent^{2}[i], alphaParent^{4}[i], 
    **/
    int myParent = minDescendent;
    int parentIdx = minDescIdx;

    while (myParent > edgeId)
    {
        parentIdx = alphaParents[parentIdx];
        myParent = alphaEdges[parentIdx];
    }

    return parentIdx;
    #endif 
}


int alphaTreeBottomUpParentSearch(int edgeId, int maxDescIdx, 
    std::vector<int>& alphaEdges, 
    std::vector<std::vector<int>>& alphaAncestors)
{
    // alphaAncestors[k] consists Parent^{2^k}[i] for every i
    assert(edgeId < alphaEdges[maxDescIdx]);
    

    #ifdef ALPHA_ANCESTOR_TRANSPOSED
    assert(alphaEdges.size() == alphaAncestors.size());
    auto mylineage= alphaAncestors[maxDescIdx];
    int logH = mylineage.size(); // number of levels
    // you can do a binary search to insert,
    // but for time being its linear search
    if(alphaEdges[mylineage[0]] < edgeId)
        return mylineage[0];
    
    for( aLvl =1; alvl< logH; alvl++)
    {
        if (mylineage[alvl] < 0 || alphaEdges[mylineage[alvl]] < edgeId )
            return alphaTreeBottomUpParentSearch(edgeId, mylineage[alvl-1], alphaEdges, alphaAncestors);
    }

    return alphaTreeBottomUpParentSearch(edgeId, mylineage[logH-1], alphaEdges, alphaAncestors);
    #else 
    assert(alphaEdges.size() == alphaAncestors[0].size());
    int logH =alphaAncestors.size();
    // DEFINE a lambda function mylineage
    //XXX mylineage(x) = alphaAncestors[x][maxDescIdx]
    // auto mylineage = [alphaAncestors, maxDescIdx](int k) { return alphaAncestors[k][maxDescIdx];};
    // if(alphaEdges[mylineage(0)] < edgeId)
    #define mylineage(k) alphaAncestors[k][maxDescIdx]
    if(alphaEdges[mylineage(0)] < edgeId)
    {
        // std::cout<<edgeId<<":\t alphaEdges["<<mylineage(0) <<"] < edgeId\n";
        return mylineage(0);
    }
        
    
    for( int alvl =1; alvl< logH; alvl++)
    {
        if (mylineage(alvl) < 0 || alphaEdges[mylineage(alvl)] < edgeId )
            return alphaTreeBottomUpParentSearch(edgeId, mylineage(alvl-1), alphaEdges, alphaAncestors);
    }

    return alphaTreeBottomUpParentSearch(edgeId, mylineage(logH-1), alphaEdges, alphaAncestors);

    #endif    
}

std::vector<std::vector<int>> computeLineageShortcut(std::vector<int>& parents)
{
    int numEdges = parents.size();
    // for each edge, parent[edge] contains its parents in edge hierarchy 
    std::vector<std::vector<int>> lineage; 
    lineage.push_back(parents);
    int lvl=0;
    int loop=1;
    while(loop)
    {
        std::vector<int> ancestors2lvl(numEdges,-1);
        loop =0; 
        for( int idx=0; idx< numEdges; idx++ )
        {
            if(lineage[lvl][idx] != -1)
            {
                ancestors2lvl[idx] = lineage[lvl][lineage[lvl][idx]];
                loop=1;
            }
                
        }

        lvl++;
        lineage.push_back(ancestors2lvl);

        
    }

    return lineage; 

}



/** 
  It returns the parent of each edge, the
 *  branch that each edge belongs to, the maximum vertex on the branch of each
 *  edge, and the alpha id of the maximum vertex on the branch of each edge.
 *  
 *  The tree is created as follows:
 *  
 *  1. The edge with alpha value zero is the root.
 *  2. The parent of each edge is the maximum vertex on its branch.
 *  3. The branch of each edge is determined by the parent of the edge.
 *  4. The maximum vertex on the branch of each edge is determined by the
 *     vertices of the edges on the branch.
 *  5. The alpha id of the maximum vertex on the branch of each edge is
 *     determined by the alpha id of the edge that has the maximum vertex as its
 *     vertex.
 *  
 *  The tree is created using two passes:
 *  1. First pass: assign all the children of each edge
 *  2. Second pass: assign all the parents of each edge
 *  
 *  The first pass is done as follows:
 *  1. First, all the edges except the root are considered unprocessed.
 *  2. The branch of each edge is determined.
 *  3. The maximum vertex on the branch of each edge is determined.
 *  4. The maximum vertex on the branch of each edge is stored.
 *  5. The alpha id of the maximum vertex on the branch of each edge is stored.
 *  
 *  The second pass is done as follows:
 *  1. First, all the edges except the root are considered unprocessed.
 *  2. The current edge is considered processed if the maximum vertex on its
 *     branch is its vertex.
 *  3. The current edge is considered unprocessed if the maximum vertex on its
 *     branch is not its vertex.
 *  4. The parent of the current edge is the edge with alpha id equal to the
 *     alpha id of the maximum vertex on the branch of the current edge.
 *  
*/

//  Type of Inputs:
//    1. alphaEdges: vector<int>
//    2. isLeftOfEdge: function

//  Type of Outputs:
//    1. parent: vector<int>
//    2. myBranch: vector<int>
//    3. maxVertexMyBranch: vector<int>
//    4. maxVertexAlphaId: vector<int>


// computeDandrogramTopDown
std::vector<int> computeDandrogramTopDown(eulerTourTree_t& inTour)
{
    // std::vector<int> alphaEdges, 
    
    
    // bool (*isLeftOfEdge)(int, int)
    int numSplitEdges = inTour.nEdges();
    int numPoints = numSplitEdges+1;

    std::vector<int> parent(numSplitEdges, 0);
    std::vector<int> isProcessed(numSplitEdges, 0);
    std::vector<int> myBranch(numSplitEdges, 0);
    std::vector<int> maxVertexMyBranch(2 * numPoints, numPoints);   // this was m_npts earlier
    std::vector<int> maxVertexAlphaId(2 * numPoints, -1);

    // First, initialise the data structures:

    isProcessed[0] = 1;
    int numNotProcessed = numSplitEdges - 1;
    int height = 0;
    parent[0] = -1;

    // Then, create the tree in two passes: first, assign all the children to the
    // nodes, and then assign all the parents. 

    while (numNotProcessed > 0)
    {
        for (int edge = 1; edge < numSplitEdges; edge++)
        {
            if (!isProcessed[edge])
            {
                // First, figure out who is my parent
                // int iam = alphaEdges[edge];
                

                myBranch[edge] = 2 * parent[edge];

                if (inTour.isLeftOfEdge(edge, parent[edge]))
                    myBranch[edge]++;

                if (maxVertexMyBranch[myBranch[edge]] > edge)
                {
                    maxVertexMyBranch[myBranch[edge]] = edge;
                    maxVertexAlphaId[myBranch[edge]] = edge;
                }
            }
        }

        for (int edge = 0; edge < numSplitEdges; edge++)
        {
            if (!isProcessed[edge])
            {
                // int edge = alphaEdges[edge];
                if (maxVertexMyBranch[myBranch[edge]] == edge)
                {
                    isProcessed[edge] = 1;
                    numNotProcessed--;
                    if(printBeta)
                        std::cout << edge << "  \t -> \t " << parent[edge] << "\n";
                }
                else
                {
                
                    parent[edge] = maxVertexAlphaId[myBranch[edge]];
                }
            }
        }

        height++;
    }
    std::cout<<" Height of the tree is "<< height <<"\n";
    return parent;
}
