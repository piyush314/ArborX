#pragma once 
#include<vector>
#include "parallel_boruvka.hpp"
#include "incidence_matrix.hpp"

typedef std::pair<int,int> indexedTimeStamp_t; 
struct  branchEdge_t
{
    int branch, edgeIdx;
    branchEdge_t()
    { branch=-1; edgeIdx=-1; }
    branchEdge_t(int a , int b): branch(a), edgeIdx(b)
    {return ;};
}; 

class alphaTree_t
{

public: 
    int numAlphaEdges;
    int m_npts;  
    std::vector<int> alphaEdges;
    std::vector<int> alphaParents; 
    std::vector<int> alphaEulerEntry;
    std::vector<int> alphaEulerExit;
    std::vector<int> alphaBridgeEdges;
    std::vector<indexedTimeStamp_t> idxTS;
    std::vector<int> mstEulerEntry; 
    std::vector<int> mstEulerExit;
    incidenceMatrix_t& m_incMatMST;
    std::vector<wtEdge_t>& m_wtSortedMST;

    

    void alpha2EdgeTree();
    int isLeftOfEdge(int queryEdge, int rootEdge);
    alphaTree_t(incidenceMatrix_t& _incMatMST, std::vector<wtEdge_t>& _wtSortedMST);

    void eulerTour();
    int timeStampDFS(int vtxId, int time, int bridgeEdge);
    std::vector<int> getAlphaEdges();
    std::vector<int> getAlphaEdges(int maxLevel);
    void alphaEulerTour();
    std::vector<int> computeBrideEdges();
    std::vector<int> constructAlphaTree();
    std::vector<int> constructEdgeTree();
    int findAlphaParent(int edgeId);

    // std::vector<int> computeFlatClustering();
    std::vector<int> computeFlatClustering(int minClusterSize);
    std::vector<branchEdge_t> computeBranchEdge();


    int alphaParent2Branch(int edgeId, int alphaParentIdx);
    int branchParent(int branchId);

    // branch related stuffs
    std::vector<branchEdge_t> m_sortedBranchEdgeArr; 
    int m_numBranches; 
	// std::vector<int> brHeads(numBranches, -1);
    std::vector<int> brIHeads;
	std::vector<int> brDescendents;
    std::vector<int> brLength;

	std::vector<float> brStabilityScores;
	std::vector<float> brSHatScores;
	std::vector<int> brDelta;
	

    int numBrDescedents(int branchId)
    {
        return brDescendents[branchId] + brLength[branchId];
    }
    
    int numBranches()
    {
        return m_numBranches;
    }

    int getBranchHead( int branchId)
    {
        return m_sortedBranchEdgeArr[brIHeads[branchId]].edgeIdx;
    }
}; 

struct branchEdgeComparator
{
    bool operator() (const branchEdge_t& a, const branchEdge_t& b )
    {
        if(a.branch == b.branch)
        {
            return a.edgeIdx < b.edgeIdx;
        }
        else
        {
            return a.branch < b.branch; 
        }
    }
};

