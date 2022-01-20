#pragma once 
#include<vector>
#include "parallel_boruvka.hpp"
#include "incidence_matrix.hpp"

std::vector<int> mergeOutSets(std::vector<int>&A, std::vector<int>&B);

class edgeHierarchy_t
{
    public:
    int m_npts;     
    std::vector<int> m_edgeParent;
    std::vector<int> m_vertexMaxIncidentEdgeId;
    incidenceMatrix_t m_incMatMST;
    std::vector<wtEdge_t>& m_wtSortedMST;
    std::vector<int> m_lastVisitedEdge;
    std::vector<int> m_numDescendents;

    std::vector<int> m_lastVisitVertexValue;
    
    // clsuter related information 

    std::vector<int>  m_numChildCluster;
    std::vector<float> m_stabilityScore;
    std::vector<int>  m_parentCluster;
    std::vector<int>  m_flatClusterMap;

    // related to eulerTour
    std::vector<int> m_eulerEntry;
    std::vector<int> m_eulerExit;
    int timeStampDFS(int vertexId, int startTime, int bridgeEdge);
    void constructEdgeTreeSplitEdges();
    std::vector<int> constructAlphaTree(std::vector<int> &alphaEdges);
    void alpha2EdgeTree(std::vector<int>& alphaEdges, std::vector<int>& alphaParent);
    std::vector<int> listOfSplitEdge();
    std::vector<int> listOfSplitEdge(int maxLevel);
    int isLeftOfEdge(int queryEdge, int rootEdge);
    
    int m_minClusterSize=3;
    bool isValidCluster(int edgeIdx)
    {
        return m_numDescendents[edgeIdx]>=m_minClusterSize; 
    }

    bool isTrueCluster(int edgeId)
    {
        if (edgeId==0) return true;
        return isValidCluster(edgeId) &&
            m_numChildCluster[m_edgeParent[edgeId]] == 2 ;
        
    }
    
    

    edgeHierarchy_t(std::vector<wtEdge_t>& wtSortedMST, int minClusterSize);

    void constructEdgeTree();
    void constructEdgeTreeBottomUP();
    // constructing the euler tour 
    void EulerTour();

    std::vector<int>  constructFlatMap(std::vector<int>& delta);
    std::vector<int> computeDelta();


    int visit(int vertexId, int edgeId);
    std::vector<std::vector<int>> eh2vecofvec();
    void writeGroups(std::ofstream& ofileName);
    void writeMSTdot(std::ofstream& ofileName);
    void writeClusterMaps(std::ofstream &ofileName);
};

