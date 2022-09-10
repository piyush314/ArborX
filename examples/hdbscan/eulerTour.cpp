#include <algorithm>
#include <iostream>
#include <cassert> 
#include <sstream>
#include <cstdio> 
#include "edge_hierarchy.hpp"
#include "incidence_matrix.hpp"
#include "alphaTree.hpp"

int timeStampDFS(int vtxId, int time, int bridgeEdge,std::pair<std::vector<int>,std::vector<int>>& eulerEntryExit, 
incidenceMatrix_t& incMatMST,std::vector<wtEdge_t>& wtSortedMST)
{
    for (int k = 0; k < incMatMST.numIncidentEdge(vtxId); k++)
    {
        auto kthEdge = incMatMST.k_thIncidentEdge(vtxId, k);
        if(kthEdge==bridgeEdge)
            continue; 
        //other end of the edge 
        eulerEntryExit.first[kthEdge] = time;
        time++; 
        int utx = wtSortedMST[kthEdge].first.first + (wtSortedMST[kthEdge].first.second - vtxId);
        time = timeStampDFS(utx, time, kthEdge,eulerEntryExit, incMatMST, wtSortedMST);
        eulerEntryExit.second[kthEdge] = time;
        time++;
    }
    return time; 
}

std::pair<std::vector<int>,std::vector<int>> eulerTour(incidenceMatrix_t& incMatMST,std::vector<wtEdge_t>& wtSortedMST)
{
    // std::vector<int>,std::vector<int>
    int time=0; 
    int root=0;
    
    std::vector<int>mstEulerEntry(wtSortedMST.size(),0); // Earlier (m_npts-1, 0);
    std::vector<int>mstEulerExit(wtSortedMST.size(),0); // Earlier (m_npts-1, 0);
    std::pair<std::vector<int>,std::vector<int>> eulerEntryExit(std::make_pair(mstEulerEntry, mstEulerExit));
    time = timeStampDFS(wtSortedMST[root].first.first, time, root, eulerEntryExit, incMatMST, wtSortedMST);
    eulerEntryExit.first[root] = time; // mstEulerEntry[root] = time;
    time++;
    time = timeStampDFS(wtSortedMST[root].first.second, time, root, eulerEntryExit, incMatMST, wtSortedMST);
    eulerEntryExit.second[root]= time; // mstEulerExit[root] = time;
    time++;
    std::cout<<"NumTime steps in euler tour ="<<time<<"\n";


    return eulerEntryExit;

}


#define TEST_BETA

void alphaTree_t::eulerTour()
{
    #ifdef TEST_BETA
    auto eulerEntryExit= ::eulerTour(m_incMatMST,m_wtSortedMST);
    mstEulerEntry = eulerEntryExit.first;
    mstEulerExit = eulerEntryExit.second;
    #else  
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
    #endif 
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
