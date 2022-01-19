#include<vector>
#include "parallel_boruvka.hpp"
#include "incidence_matrix.hpp"


incidenceMatrix_t::incidenceMatrix_t(std::vector<wtEdge_t>& wtSortedMST): m_wtMST(wtSortedMST)
{
    m_npts = m_wtMST.size()+1;
    m_vEdgeListIndex.resize(m_npts+1, 0);
    m_edgeIdx.resize(2*m_npts -2); // double the number of edges  

    // std::vector<int> edgeCounts(m_npts+1, 0); // number of edge count

    // get edge count 
    for(auto wte: m_wtMST)
    {
        m_vEdgeListIndex[wte.first.first]++;
        m_vEdgeListIndex[wte.first.second]++;
    }

    // do a exclusive prefix on edge counts 
    // prefixSumExclusive(m_npts+1, m_vEdgeListIndex.data());
    prefixSumExclusive(m_npts, m_vEdgeListIndex.data());


    std::vector<int> runningEdgeCounts(m_npts, 0);
    for(int i=0; i<m_wtMST.size(); i++)
    {
        int v1 = m_wtMST[i].first.first;
        int v2 = m_wtMST[i].first.second;
        // add edge i to edgelistof v1 and v2
        m_edgeIdx[m_vEdgeListIndex[v1] + runningEdgeCounts[v1]++ ] =i;
        m_edgeIdx[m_vEdgeListIndex[v2] + runningEdgeCounts[v2]++ ] =i;
    }
}

int incidenceMatrix_t::numIncidentEdge(int k)
{
    return m_vEdgeListIndex[k+1] -m_vEdgeListIndex[k];
}

int incidenceMatrix_t::k_thIncidentEdge(int vertexId, int k)
{
    return m_edgeIdx[m_vEdgeListIndex[vertexId] + k];
}

// TODO: can be made a function and array can be part of class
std::vector<int> incidenceMatrix_t::maxIncidentEdgeId()
{
    std::vector<int> out(m_npts, -1);
    for(int vtx=0; vtx< m_npts; vtx++)
    {
        out[vtx] = *std::max_element(&m_edgeIdx[m_vEdgeListIndex[vtx]],&m_edgeIdx[m_vEdgeListIndex[vtx+1]]);
    }
    return out; 
}


int incidenceMatrix_t::maxIncidentEdgeId(int vtx)
{
    return  *std::max_element(&m_edgeIdx[m_vEdgeListIndex[vtx]],
            &m_edgeIdx[m_vEdgeListIndex[vtx+1]]);
}


int incidenceMatrix_t::isLeafEdge(int edgeIdx)
{
    int termVertex1 = m_wtMST[edgeIdx].first.first;
    int termVertex2 = m_wtMST[edgeIdx].first.second;

    if(edgeIdx < maxIncidentEdgeId(termVertex1) ||
        edgeIdx < maxIncidentEdgeId(termVertex2) )
        return 0; 
    return 1; 
}


int incidenceMatrix_t::isAlphaEdge(int edgeIdx)
{
    int termVertex1 = m_wtMST[edgeIdx].first.first;
    int termVertex2 = m_wtMST[edgeIdx].first.second;

    if(edgeIdx < maxIncidentEdgeId(termVertex1) &&
        edgeIdx < maxIncidentEdgeId(termVertex2) )
        return 1; 
    return 0; 
}