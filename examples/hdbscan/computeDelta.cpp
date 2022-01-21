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

#if 0

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
                            m_stabilityScore[parent] += clDescends[edgeId] *tmpScore;

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
    int verbose=1;
    auto sHatScore =m_stabilityScore;
    std::vector<int> delta(m_npts - 1, 1);   

    
    // Leaf cluster stabilty caclulation
    for (int vtx = 0; vtx < m_npts; vtx++)
    {
        int parent = m_vertexMaxIncidentEdgeId[vtx];
        if (parent > -1 && isValidCluster(parent) )
        {
            m_stabilityScore[parent] += 1.0/ m_wtSortedMST[parent].second; 
            if(verbose)
            std::cout<<"0. s["<<parent<<"] \t+= 1.0/"<<m_wtSortedMST[parent].second
            <<" \t="<<1.0/ m_wtSortedMST[parent].second<<" \t="<<m_stabilityScore[parent]
             << "\n";
        }
            
    }

    // calculate stability
    
    for (int edgeId = m_npts - 2; edgeId > 0; edgeId--)
    {
        int parent = m_edgeParent[edgeId];
        if (!isValidCluster(edgeId) and parent >-1 and  isValidCluster(parent))
        {
            m_stabilityScore[parent] += m_numDescendents[edgeId]*1.0/ m_wtSortedMST[parent].second;
            if(verbose)
            std::cout<<"1. s["<<parent<<"] \t+=" << m_numDescendents[edgeId] <<".0/"<<m_wtSortedMST[parent].second
            <<" \t="<<m_numDescendents[edgeId]*1.0/ m_wtSortedMST[parent].second<<" \t="<<m_stabilityScore[parent]
             << "\n";
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
                    if(verbose)
                    printf("2. s[%d] \t= %d x (1/%f - 1/%f) \t= %f\n", edgeId,m_numDescendents[edgeId],
                    m_wtSortedMST[edgeId].second, m_wtSortedMST[parent].second, m_stabilityScore[edgeId] );
                }
                m_stabilityScore[parent] += m_numDescendents[edgeId]*1.0/ m_wtSortedMST[edgeId].second;
                if(verbose)
                printf("a. s[%d] \t+= %d x (1/%f ) \t= %f\n", parent,m_numDescendents[edgeId],
                     m_wtSortedMST[edgeId].second, m_stabilityScore[parent] );
            }

            if(m_numChildCluster[edgeId]==1)
            {
                if (isTrueCluster(edgeId))
                {
                    m_stabilityScore[edgeId] -= m_numDescendents[edgeId]*
                    ( 1.0/ m_wtSortedMST[parent].second);
                    if(verbose)
                    printf("3. s[%d] \t-= %d x (1/%f ) \t= %f\n", edgeId,m_numDescendents[edgeId],
                    m_wtSortedMST[parent].second, m_stabilityScore[edgeId] );
                    
                }
                m_stabilityScore[parent] += m_stabilityScore[edgeId];
                if(verbose)
                printf("b. s[%d] \t+= s[%d] \t= %f\n", parent, edgeId,
                      m_stabilityScore[parent] );
            }
            
            if(m_numChildCluster[edgeId]==2)
            {
                if (isTrueCluster(edgeId))
                {
                    m_stabilityScore[edgeId] = m_numDescendents[edgeId]*
                    (1.0/ m_wtSortedMST[edgeId].second -  1.0/ m_wtSortedMST[parent].second);
                    if(verbose)
                    printf("4. s[%d] \t= %d x (1/%f - 1/%f) \t= %f\n", edgeId,m_numDescendents[edgeId],
                    m_wtSortedMST[edgeId].second, m_wtSortedMST[parent].second, m_stabilityScore[edgeId] );
                }
                m_stabilityScore[parent] += m_numDescendents[edgeId]*1.0/ m_wtSortedMST[edgeId].second;
                if(verbose)
                printf("c. s[%d] \t+= %d x (1/%f ) \t= %f\n", parent,m_numDescendents[edgeId],
                     m_wtSortedMST[edgeId].second, m_stabilityScore[parent] );
            }
        }



        
        if (isTrueCluster(edgeId) && edgeId != 0)
        {
            
            
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
            std::cout<<"sHat["<<edgeId <<"] =  "<<sHatScore[edgeId]<<" "<<m_stabilityScore[edgeId]<<" "<<delta[edgeId]<< "\n";
            auto parent = m_edgeParent[edgeId];
            while(parent >-1)
            {
                if(isTrueCluster(parent))
                {
                    sHatScore[parent] += sHatScore[edgeId];
                    // printf("%d: sHat[%d] += %f = %f\n", edgeId, parent, sHatScore[edgeId], sHatScore[parent]);
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
