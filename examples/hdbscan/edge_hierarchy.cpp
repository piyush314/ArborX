#include <algorithm>
#include <iostream>
#include <cassert> 
#include <sstream>
#include "edge_hierarchy.hpp"
#include "incidence_matrix.hpp"

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

std::vector<int> mergeOutSets(std::vector<int>&A, std::vector<int>&B)
{
    std::vector<int> AB;
    AB.reserve(A.size() + B.size());
    AB.insert(AB.end(), A.begin(), A.end());
    AB.insert(AB.end(), B.begin(), B.end());
    
    return AB;
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

    int edgesNotProcessed = m_npts-1;
    while(edgesNotProcessed>1)      // no need to process vertex zero
    {
        for ( int edgeId = 0; edgeId < m_npts-1; edgeId++)
        {
            /* code */
            if(!isEdgeProcessed[edgeId]  // edge is not processed
                && nChildrenProcessed[edgeId]==2)        // both its childrens are processed
            {
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
        
        // std::cout<<"edgesNotProcessed" << edgesNotProcessed<<"\n"; 
        // loopCount++;
        // if(loopCount==100) exit(0);
    }
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

    std::vector<int> descloopCount(m_npts,0); 

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
                        if(grandParent!=-1)
                        {
                            double tmpScore =( 1.0/ m_wtSortedMST[edgeId].second - 1.0/ m_wtSortedMST[grandParent].second);
                            if(m_numChildCluster[edgeId]==2)
                            {
                                m_stabilityScore[parent] += m_numDescendents[edgeId] * tmpScore;
                                std::cout<<edgeId<<"  "<<parent<<" "<< grandParent << "\n";
                                std::cout<< m_wtSortedMST[edgeId].second<<"  "<<m_wtSortedMST[parent].second<<" "<< m_wtSortedMST[grandParent].second << "\n";
                                std::cout<<1/m_wtSortedMST[edgeId].second<<"  "<<1/m_wtSortedMST[parent].second<<" "<< 1/m_wtSortedMST[grandParent].second << "\n";
                            }
                                
                            else 
                                m_stabilityScore[parent] += clDescends[edgeId] *tmpScore;

                            if(m_numChildCluster[edgeId]==2)
                                descloopCount[parent] += m_numDescendents[edgeId] ;
                            else 
                                descloopCount[parent] += clDescends[edgeId];
                        }
                        
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
            
            // assert(descloopCount[edgeId]==m_numDescendents[edgeId]);
            
            // m_stabilityScore[edgeId] -= m_numDescendents[edgeId] / m_wtSortedMST[parent].second;
            std::cout << edgeId << " \t : stability score =" << m_stabilityScore[edgeId] <<
            "\t sHat ="<<sHatScore[edgeId] << "\n";
            std::cout<<"loopCOunted desc "<<descloopCount[edgeId]<< " atcual descendent " << m_numDescendents[edgeId]<<"\n";
            if(sHatScore[edgeId]> m_stabilityScore[edgeId])
            {
                delta[edgeId] =0;
            }
            else
            {
                sHatScore[edgeId]= m_stabilityScore[edgeId];
            }
            std::cout<<"sHat =  "<<sHatScore[edgeId]<<" "<<delta[edgeId]<< "\n";
            parent = m_edgeParent[parent];
            while(parent >-1)
            {
                if(isTrueCluster(parent))
                {
                    sHatScore[parent] += sHatScore[edgeId];
                    break;         
                }
                else
                    parent = m_edgeParent[parent];
            }
            
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

    constructEdgeTreeBottomUP();

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