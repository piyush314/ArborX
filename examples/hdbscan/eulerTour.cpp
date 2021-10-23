#include <algorithm>
#include <iostream>
#include <cassert> 
#include <sstream>
#include <cstdio> 
#include "edge_hierarchy.hpp"
#include "incidence_matrix.hpp"
#include "alphaTree.hpp"

void alphaTree_t::eulerTour()
{
    int time=0; 
    int root=0;
    mstEulerEntry.resize(m_npts-1, 0);
    mstEulerExit.resize(m_npts-1, 0);
    time = timeStampDFS(m_wtSortedMST[root].first.first, time, root);
    mstEulerEntry[root] = time;
    time++;
    time = timeStampDFS(m_wtSortedMST[root].first.second, time, root);
    mstEulerExit[root] = time;
    time++;
    std::cout<<"NumTime steps in euler tour ="<<time<<"\n";
}

int alphaTree_t::timeStampDFS(int vtxId, int time, int bridgeEdge)
{
    for (int k = 0; k < m_incMatMST.numIncidentEdge(vtxId); k++)
    {
        auto kthEdge = m_incMatMST.k_thIncidentEdge(vtxId, k);
        if(kthEdge==bridgeEdge)
            continue; 
        //other end of the edge 
        mstEulerEntry[kthEdge] = time;
        time++; 
        int utx = m_wtSortedMST[kthEdge].first.first + (m_wtSortedMST[kthEdge].first.second - vtxId);
        time = timeStampDFS(utx, time, kthEdge);
        mstEulerExit[kthEdge] = time;
        time++;
    }
    return time; 
}
