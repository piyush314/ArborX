#pragma once 

#include <vector>
#include <limits>
#include <algorithm>
#include "parallel_boruvka.hpp"

class incidenceMatrix_t
{
    int m_npts; 
    std::vector<wtEdge_t>& m_wtMST;
    std::vector<int> m_vEdgeListIndex;
    std::vector<int> m_edgeIdx; 
public: 
    incidenceMatrix_t(std::vector<wtEdge_t>& wtSortedMST);
    int numIncidentEdge(int k);
    int k_thIncidentEdge(int vertexId, int k);
    std::vector<int> maxIncidentEdgeId();
};


