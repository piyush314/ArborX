#include <algorithm>
#include <iostream>
#include <cassert> 
#include <sstream>
#include <cstdio> 
#include "edge_hierarchy.hpp"
#include "incidence_matrix.hpp"

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
#include "stability.hpp"
#if 1

std::vector<int> edgeHierarchy_t::computeDelta()
{
    // auto sHatScore =m_stabilityScore;
    sHatScore.resize(m_npts - 1, 0.0);
    std::vector<int> delta(m_npts - 1, 1);   

    std::vector<int> clDescends(m_npts - 1, 0);
    std::vector<stability_t> SscoreBr(m_npts - 1);
    
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

    // calculate stability
    for (int edgeId = m_npts - 2; edgeId > 0; edgeId--)
    {
        if (isValidCluster(edgeId))
        {
            std::cout<<"------ Edge  ="<<edgeId << "  #Desc " << m_numDescendents[edgeId]<<"\n";
        }
    }
    for (int edgeId = m_npts - 2; edgeId > 0; edgeId--)
    {
        // std::cout<<"------ Edge  ="<<edgeId << "  clDesc " << clDescends[edgeId]<<"\n";
        int parent = edgeId;
        if (isValidCluster(edgeId))
        {
            while (parent > -1)
            {
                if (isTrueCluster(parent))
                {
                    
                    int grandParent = m_edgeParent[parent];
                    if(grandParent!=-1)
                    {
                        double tmpScore =( 1.0/ m_wtSortedMST[edgeId].second - 1.0/ m_wtSortedMST[grandParent].second);

                        if(m_numChildCluster[edgeId]==2)
                        {
                            m_stabilityScore[parent] += m_numDescendents[edgeId] * tmpScore;
                            
                        }
                        else 
                        {
                            m_stabilityScore[parent] += clDescends[edgeId] *tmpScore;
                            SscoreBr[parent] = SscoreBr[parent] + std::make_pair(edgeId, clDescends[edgeId]);
                            SscoreBr[parent] = SscoreBr[parent] - std::make_pair(grandParent, clDescends[edgeId]);
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
            if(m_stabilityScore[edgeId]==0)
            {
                std::cout<<"Edge "<<edgeId<<" has zero stability score: \t ";
                SscoreBr[edgeId].print();
                std::cout<<"\n";
            }
            
            if(m_stabilityScore[edgeId]<0)
            {
                std::cout<<"NONZERO STABILITY SCORE edgeId = "<<edgeId<<" "<<m_stabilityScore[edgeId]<<"\n";
                // SscoreBr[edgeId].print();
            }
            assert(m_stabilityScore[edgeId]>=0);
            if(sHatScore[edgeId]>= m_stabilityScore[edgeId])
            {
                delta[edgeId] =0;
            }
            else
            {
                sHatScore[edgeId]= m_stabilityScore[edgeId];
            }
            // if(delta[edgeId])
            //     std::cout<<"sHat["<<edgeId <<"] =  "<<sHatScore[edgeId]<< "\n";
            // std::cout<<"sHat["<<edgeId <<"] =  "<<sHatScore[edgeId]<<" "<<delta[edgeId]<< "\n";
            std::cout<<"sHat["<<edgeId <<"] =  "<<sHatScore[edgeId]<<" "<<m_stabilityScore[edgeId]<<" "<<delta[edgeId]<< "\n";
            parent = m_edgeParent[edgeId];
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

#else 

std::vector<int> edgeHierarchy_t::computeDelta()
{
    int verbose=0;
    // TODO: for checking shat from alpha tree delete it later 
    #if 0
    auto sHatScore =m_stabilityScore;
    #else 
        sHatScore.resize(m_npts - 1, 0.0);
    #endif 
    std::vector<int> delta(m_npts - 1, 1);   

    std::vector<stability_t> SscoreBr(m_npts - 1);

    // Leaf cluster stabilty caclulation
    for (int vtx = 0; vtx < m_npts; vtx++)
    {
        int parent = m_vertexMaxIncidentEdgeId[vtx];
        if (parent > -1 && isValidCluster(parent) )
        {
            m_stabilityScore[parent] += 1.0/ m_wtSortedMST[parent].second; 
            SscoreBr[parent] = SscoreBr[parent] + parent;
            
        }
            
    }

    // calculate stability
    
    for (int edgeId = m_npts - 2; edgeId > 0; edgeId--)
    {
        int parent = m_edgeParent[edgeId];
        if (!isValidCluster(edgeId) and parent >-1 and  isValidCluster(parent))
        {
            m_stabilityScore[parent] += m_numDescendents[edgeId]*1.0/ m_wtSortedMST[parent].second;
            SscoreBr[parent] = SscoreBr[parent] + std::make_pair(parent, m_numDescendents[edgeId]);
        }
        if (isValidCluster(edgeId))
        {
            // leaf valid cluster 
            if(m_numChildCluster[edgeId]==0)
            {
                if (isTrueCluster(edgeId))
                {
                    m_stabilityScore[edgeId] = m_numDescendents[edgeId]*
                    (1.0/ m_wtSortedMST[edgeId].second -  1.0/ m_wtSortedMST[parent].second);
                    SscoreBr[edgeId] = SscoreBr[edgeId] + std::make_pair(edgeId, m_numDescendents[edgeId]);
                    SscoreBr[edgeId] = SscoreBr[edgeId] - std::make_pair(parent, m_numDescendents[edgeId]);
                    
                }
                m_stabilityScore[parent] += m_numDescendents[edgeId]*1.0/ m_wtSortedMST[edgeId].second;
                SscoreBr[parent] = SscoreBr[parent] + std::make_pair(edgeId, m_numDescendents[edgeId]);
                
            }

            if(m_numChildCluster[edgeId]==1)
            {
                if (isTrueCluster(edgeId))
                {
                    m_stabilityScore[edgeId] -= m_numDescendents[edgeId]*
                    ( 1.0/ m_wtSortedMST[parent].second);
                    SscoreBr[edgeId] = SscoreBr[edgeId] - std::make_pair(parent, m_numDescendents[edgeId]);
                    
                    
                }
                m_stabilityScore[parent] += m_stabilityScore[edgeId];
                SscoreBr[parent] = SscoreBr[parent] + SscoreBr[edgeId];
                
            }
            
            if(m_numChildCluster[edgeId]==2)
            {
                if (isTrueCluster(edgeId))
                {
                    m_stabilityScore[edgeId] = m_numDescendents[edgeId]*
                    (1.0/ m_wtSortedMST[edgeId].second -  1.0/ m_wtSortedMST[parent].second);
                    SscoreBr[edgeId] = SscoreBr[edgeId] + std::make_pair(edgeId, m_numDescendents[edgeId]);
                    SscoreBr[edgeId] = SscoreBr[edgeId] - std::make_pair(parent, m_numDescendents[edgeId]);
                    
                }
                m_stabilityScore[parent] += m_numDescendents[edgeId]*1.0/ m_wtSortedMST[edgeId].second;
                SscoreBr[parent] = SscoreBr[parent] + std::make_pair(edgeId, m_numDescendents[edgeId]);
                
            }
        }



        
        if (isTrueCluster(edgeId) && edgeId != 0)
        {
            
            if(SscoreBr[edgeId].getScore() != 0)
            {
                std::cout<<"SscoreBr["<<edgeId <<"] =  "<<SscoreBr[edgeId].getScore()<< "\n";
            }
            
            if(m_stabilityScore[edgeId]<0)
            {
                std::cout<<"NONZERO STABILITY SCORE edgeId = "<<edgeId<<" "<<m_stabilityScore[edgeId]<<"\n";
                SscoreBr[edgeId].print();
            }
            assert(m_stabilityScore[edgeId]>=0);
            if(sHatScore[edgeId]>= m_stabilityScore[edgeId])
            {
                delta[edgeId] =0;
            }
            else
            {
                sHatScore[edgeId]= m_stabilityScore[edgeId];
            }

            if(verbose)
            std::cout<<"sHat["<<edgeId <<"] =  "<<sHatScore[edgeId]<<" "<<m_stabilityScore[edgeId]<<" "<<delta[edgeId]<< "\n";
            auto parent = m_edgeParent[edgeId];
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
#endif
