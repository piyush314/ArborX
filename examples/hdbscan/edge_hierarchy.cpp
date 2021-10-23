#include <algorithm>
#include <iostream>
#include <cassert> 
#include <sstream>
#include <cstdio> 
#include "edge_hierarchy.hpp"
#include "incidence_matrix.hpp"
#include "alphaTree.hpp"

std::vector<int> mergeOutSets(std::vector<int>&A, std::vector<int>&B)
{
    std::vector<int> AB;
    AB.reserve(A.size() + B.size());
    AB.insert(AB.end(), A.begin(), A.end());
    AB.insert(AB.end(), B.begin(), B.end());
    
    return AB;
}

int visitloopCount = 0;



void edgeHierarchy_t::constructEdgeTree()
{

    // can be done in parallel 
    for (int edgeId = m_npts - 2; edgeId > 0; edgeId--)
    {
        
        auto x1 = visit(m_wtSortedMST[edgeId].first.first, edgeId);
        auto x2 = visit(m_wtSortedMST[edgeId].first.second, edgeId);
        m_edgeParent[edgeId] = std::max(x1, x2);
        // Add descendents to parent
        m_numDescendents[m_edgeParent[edgeId]] += m_numDescendents[edgeId];
        // if I am a valid cluster
        if (isValidCluster(edgeId))
            m_numChildCluster[m_edgeParent[edgeId]] += 1;
    }
    std::cout<<"== Total visits are " <<visitloopCount<< "\n";
}



int loopCount =0; 
void edgeHierarchy_t::constructEdgeTreeBottomUP()
{

    // first n-1, corresponds to edges, following n corresponds to vertices 
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
    }
}


edgeHierarchy_t::edgeHierarchy_t(std::vector<wtEdge_t> &wtSortedMST, int minClusterSize) : 
m_wtSortedMST(wtSortedMST), m_incMatMST(wtSortedMST), m_minClusterSize(minClusterSize)
{
    m_npts = m_wtSortedMST.size() + 1;
    m_edgeParent.resize(m_npts - 1, -1);
    m_lastVisitedEdge.resize(m_npts, -1);
    m_vertexMaxIncidentEdgeId = m_incMatMST.maxIncidentEdgeId();
    m_numDescendents.resize(m_npts - 1, 0);
    m_numChildCluster.resize(m_npts - 1, 0);
    m_stabilityScore.resize(m_npts - 1, 0);
    m_parentCluster.resize(m_npts - 1, 0);
    m_flatClusterMap.resize(m_npts,-1);
    m_lastVisitVertexValue.resize(m_npts,m_npts);

    m_eulerEntry.resize(m_npts-1, 0);
    m_eulerExit.resize(m_npts-1, 0);

    auto sHatScore =m_stabilityScore;
    // std::vector<int> delta(m_npts - 1, 1);   
    // .resize(m_npts,-1);

    /* Perform an euler Tour */ 
    #define SPLIT_EDGE_TREE 
    #ifdef SPLIT_EDGE_TREE 
    // EulerTour();
    // std::vector<int> alphaEdges= listOfSplitEdge();
    // constructAlphaTree(alphaEdges);
    
    // construct alpha Tree
    alphaTree_t alphaTree(m_incMatMST, m_wtSortedMST);
    // TODO: get the edge Tree 
    
    m_edgeParent = alphaTree.constructEdgeTree();
    exit(0);
    // for(int edgeId=0; edgeId<m_npts-1;edgeId++)
    // {
    //     std::cout<<" \t | "<<edgeId<<" \t | "<< m_eulerEntry[edgeId] <<" \t | "<< m_eulerExit[edgeId]<<" \t | " << "\n";
    // }
    #else 
    // add num descendents to leaf nodes
    for (int vtx = 0; vtx < m_npts; vtx++)
    {
        int parent = m_vertexMaxIncidentEdgeId[vtx];
        m_numDescendents[parent]++;
        m_lastVisitVertexValue[vtx] = parent; 
    }

    constructEdgeTreeBottomUP();
    #endif 
    auto delta = computeDelta();
    m_flatClusterMap = constructFlatMap(delta);
    
}


int edgeHierarchy_t::visit(int vtxId, int edgeId)
{
    // loopCount++;
    
    
    if (m_lastVisitedEdge[vtxId] == edgeId)
        return -1;
    m_lastVisitedEdge[vtxId] = edgeId;
    
    if (m_lastVisitVertexValue[vtxId] < edgeId)
    {
        return m_lastVisitVertexValue[vtxId];
    }
    else
    {
        int vst_v_e = -1;
        for (int k = 0; k < m_incMatMST.numIncidentEdge(vtxId); k++)
        {
            visitloopCount++;
            auto kthEdge = m_incMatMST.k_thIncidentEdge(vtxId, k);
            // std::cout<<kthEdge<<" "<<k<<" ";
            if (kthEdge == edgeId)
                continue;
            if (kthEdge < edgeId)
                vst_v_e = std::max(vst_v_e, kthEdge);
            else
            {
                // traverse through the edge
                // utx is the other end of the kthedge (one end is vtxId)

                int utx = m_wtSortedMST[kthEdge].first.first + (m_wtSortedMST[kthEdge].first.second - vtxId);
                // std::cout << "Visting " << utx << " " << vtxId << "\n";
                vst_v_e = std::max(vst_v_e,
                                   visit(utx, edgeId));
            }
        }
        m_lastVisitVertexValue[vtxId] = vst_v_e;
        return vst_v_e;
    }
}

std::vector<std::vector<int>> edgeHierarchy_t::eh2vecofvec()
{
    std::vector<std::vector<int>> out(m_npts - 1);

    for (int vtx = 0; vtx < m_npts; vtx++)
    {
        std::cout << "Vertex=" << vtx << "\n";
        int parent = m_vertexMaxIncidentEdgeId[vtx];
        int loopCount = 0;
        while (parent > -1)
        {
            std::cout << "\t Parent=" << parent << "\n";
            out[parent].push_back(vtx);
            parent = m_edgeParent[parent];
        }
    }

    return out;
}

void edgeHierarchy_t::writeGroups(std::ofstream &ofileName)
{
    std::vector<std::vector<int>> groupVec = eh2vecofvec();

    for (auto grp : groupVec)
    {
        for (auto pts : grp)
        {
            ofileName << pts << " ";
        }
        ofileName << "\n";
    }
}

std::string svertex(int i)
{
  std::ostringstream stringStream;
  stringStream << "v_"<<i;
  std::string copyOfStr = stringStream.str();
  return copyOfStr;
}

std::string sedge(int i)
{
  std::ostringstream stringStream;
  stringStream << "E_"<<i;
  std::string copyOfStr = stringStream.str();
  return copyOfStr;
}

void edgeHierarchy_t::writeMSTdot(std::ofstream& ofileName)
{
    ofileName << "digraph Edge_Hierarchy           \{ \n"; 
    ofileName << "\t label = \"Edge Hierarchy\" ;\n";
    ofileName << "\t rankdir = BT;\n";
    ofileName << "\t ranksep=.75; size = \"11,8\";\n";

    ofileName << "graph [bgcolor=\"#41444A\"];\n";;
	ofileName << "edge [color=white];\n";;

    // now add vertices and edges
    for(int i=0; i<m_npts; i++)
    {
        ofileName << svertex(i) <<" ";
        ofileName << "[shape=doublecircle  , style=filled,fillcolor=grey];\n";

        if(i!=m_npts-1)
        {
            ofileName << sedge(i) <<" ";
            std::string colorName;
            if( isValidCluster(i) && m_numChildCluster[i]==2)
            {
                colorName = std::string("#1F9D97");    
            }
            
            else if (isTrueCluster(i) && i != 0)
            {
                /* code */
                colorName = std::string("#EDAA0F"); 
                
            }
            else
            {
                colorName = std::string("#1A98C4");    
            }
            

            ofileName << "[shape=cds  , style=filled,fillcolor=\""<<colorName<<"\"];\n"; 

            
        }

    }

    for(int vtx=0; vtx<m_npts; vtx++)
    {
        int parent = m_vertexMaxIncidentEdgeId[vtx];   
        ofileName << svertex(vtx) <<" ->" << sedge(parent)<<";\n";
        
        if(vtx!=0 && vtx!=m_npts-1)
        {
            int parentEdge = m_edgeParent[vtx];
            ofileName << sedge(vtx) <<" -> "<< sedge(parentEdge)<<";\n";
            // ofileName << "[shape=cds  , style=filled,fillcolor=white];\n";
        }
    }


    ofileName << "} \n"; 

}


void edgeHierarchy_t::writeClusterMaps(std::ofstream &ofileName)
{
    // std::vector<std::vector<int>> groupVec = eh2vecofvec();

    for (auto cmap : m_flatClusterMap)
    {
        
        ofileName << cmap << "\n";
        
        
    }
}