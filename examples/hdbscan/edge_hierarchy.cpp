#include <algorithm>
#include <iostream>
#include <cassert> 
#include <sstream>
#include "edge_hierarchy.hpp"
#include "incidence_matrix.hpp"

int visitCount = 0;

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
    std::cout<<"== Total visits are " <<visitCount<< "\n";
}

std::vector<int>  edgeHierarchy_t::constructFlatMap(std::vector<int>& delta)
{
    std::vector<int>flatClusterMap(m_npts,-1);
    for(int vtxId=0; vtxId<m_npts;vtxId++)
    {
        int parent = m_vertexMaxIncidentEdgeId[vtxId];
        // flatClusterMap[vtxId] =-1;
        int first_true_cluster =-1;
        
        while (parent > 0)
        {
            if(isValidCluster(parent) && first_true_cluster == -1)
                first_true_cluster = parent; 
                
            if (isTrueCluster(parent) && delta[parent]==1)
                flatClusterMap[vtxId] = parent;
                    
            parent = m_edgeParent[parent];
        }

        // Marking noise points     
        // if(flatClusterMap[vtxId] == first_true_cluster)
        //     flatClusterMap[vtxId]=-1;
    }

    return flatClusterMap;
}


std::vector<int> edgeHierarchy_t::computeDelta()
{
    auto sHatScore =m_stabilityScore;
    std::vector<int> delta(m_npts - 1, 1);   

    std::vector<int> clDescends(m_npts - 1, 0);
    
    // Leaf cluster stabilty caclulation
    for (int vtx = 0; vtx < m_npts; vtx++)
    {
        int parent = m_vertexMaxIncidentEdgeId[vtx];
        int lambdaMaxCluster = -1;
        while (parent > -1)
        {
            if (isValidCluster(parent))
            {
                clDescends[parent]++;
                break;
            }
            parent = m_edgeParent[parent];
        }
    }

    std::vector<int> descCount(m_npts,0); 

    // calculate stability
    
    for (int edgeId = m_npts - 2; edgeId > 0; edgeId--)
    {
        // std::cout<<"------ Edge  ="<<edgeId << "\n";
        int parent = edgeId;
        if (isValidCluster(edgeId))
        {
            while (parent > -1)
            {
                if (isTrueCluster(parent))
                {
                    
                    if(0)
                    {
                        int grandParent = m_edgeParent[parent];
                        m_stabilityScore[parent] += m_numDescendents[edgeId] *
                        ( 1.0/ m_wtSortedMST[edgeId].second - 1.0/ m_wtSortedMST[grandParent].second);
                    }
                    else
                    {
                        int grandParent = m_edgeParent[parent];
                        double tmpScore =( 1.0/ m_wtSortedMST[edgeId].second - 1.0/ m_wtSortedMST[grandParent].second);
                        if(m_numChildCluster[edgeId]==2)
                            m_stabilityScore[parent] += m_numDescendents[edgeId] * tmpScore;
                        else 
                            m_stabilityScore[parent] += clDescends[edgeId] *tmpScore;

                        if(m_numChildCluster[edgeId]==2)
                            descCount[parent] += m_numDescendents[edgeId] ;
                        else 
                            descCount[parent] += clDescends[edgeId];
                    }
                        
                    break;
                }
                else
                {
                    parent = m_edgeParent[parent];
                }
            }
        }

        
        if (isTrueCluster(edgeId) && edgeId != 0)
        {
            parent = m_edgeParent[parent];
            assert(descCount[edgeId]==m_numDescendents[edgeId]);
            
            // m_stabilityScore[edgeId] -= m_numDescendents[edgeId] / m_wtSortedMST[parent].second;
            std::cout << edgeId << " \t : stability score =" << m_stabilityScore[edgeId] << "\n";
            if(sHatScore[edgeId]> m_stabilityScore[edgeId])
            {
                delta[edgeId] =0;
            }
            else
            {
                sHatScore[edgeId]= m_stabilityScore[edgeId];
            }
            std::cout<<"sHat =  "<<sHatScore[edgeId]<<" "<<delta[edgeId]<< "\n";
            sHatScore[parent] += sHatScore[edgeId];
        }
    }

    return delta; 
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

    auto sHatScore =m_stabilityScore;
    // std::vector<int> delta(m_npts - 1, 1);   
    // .resize(m_npts,-1);

    // add num descendents to leaf nodes
    for (int vtx = 0; vtx < m_npts; vtx++)
    {
        int parent = m_vertexMaxIncidentEdgeId[vtx];
        m_numDescendents[parent]++;
        m_lastVisitVertexValue[vtx] = parent; 
    }

    constructEdgeTree();

    auto delta = computeDelta();
    m_flatClusterMap = constructFlatMap(delta);
    
}


int edgeHierarchy_t::visit(int vtxId, int edgeId)
{
    // count++;
    visitCount++;
    
    if (m_lastVisitedEdge[vtxId] == edgeId)
        return -1;
    m_lastVisitedEdge[vtxId] = edgeId;
    
    if (m_lastVisitVertexValue[vtxId] < edgeId)
        return m_lastVisitVertexValue[vtxId];
    else
    {
        int vst_v_e = -1;
        for (int k = 0; k < m_incMatMST.numIncidentEdge(vtxId); k++)
        {
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
        int count = 0;
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
    ofileName << "\t rankdir = RL;\n";
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
            if (isTrueCluster(i) && i != 0)
            {
                /* code */
                colorName = std::string("#EDAA0F"); 
                
            }
            else if( isValidCluster(i) && m_numChildCluster[i]==2)
            {
                colorName = std::string("#1F9D97");    
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