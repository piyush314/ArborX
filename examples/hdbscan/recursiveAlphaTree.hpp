#pragma once 
#include "alphaTree.hpp"

class eulerTourTree_t
{
    public: 
        std::vector<indexedTimeStamp_t> idxTS; // index with time stamp 
        std::vector<int> entryTime;
        std::vector<int> exitTime;
        std::vector<int> bridgeEdges;
        std::vector<int> maxIncidentOnBridge;
        std::vector<int>& computemaxIncidentOnBridge();
        int maxIncidentOnRoot=-1;
        int maxIdxBridge(int k);
        
        std::vector<int> computeAlphaEdges();
        int nEdges()
        {
            return entryTime.size();
        }
        std::pair<int,int> simInsert(int pEntry, int pExit);
        // constructor using incidence matrix
        eulerTourTree_t(incidenceMatrix_t& incMatMST,std::vector<wtEdge_t>& wtSortedMST);
        eulerTourTree_t() {};

        bool isLeftOfEdge( int queryEdge, int rootEdge);
        int simBridge(int lEntry);
        int entry(int k)
        {
            return entryTime[k];
        }

        int exit(int k)
        {
            return exitTime[k];
        }
        int computeSimulatedBridge(int pEntry, int pExit);
};

eulerTourTree_t computeAlphaEulerTour(std::vector<int> listOfAlphaEdges, eulerTourTree_t& inTour);
int timeStampDFS(int vtxId, int time, int bridgeEdge,std::pair<std::vector<int>,std::vector<int>>& eulerEntryExit, 
incidenceMatrix_t& incMatMST,std::vector<wtEdge_t>& wtSortedMST);
std::pair<std::vector<int>,std::vector<int>> eulerTour(incidenceMatrix_t& incMatMST,std::vector<wtEdge_t>& wtSortedMST);
std::vector<int> computeDandrogramTopDown(eulerTourTree_t& inTour);


std::vector<std::vector<int>> computeLineageShortcut(std::vector<int>& parents);
class rAlphaTree_t;
class rAlphaTree_t
{

public: 
    std::vector<int> alphaEdges;    // list of alpha edges relative to input MST
    std::vector<int> invAlphaEdges;   // isAlphaEdge[i]=-1. if i is not alpha edge, 
                                    // =idx, where i=alphaEdges[idx] if alpha edge
    eulerTourTree_t eAlphaTour;     // MST of alphaTree in euler tour format
     
    rAlphaTree_t* betaTree=NULL;      // alphatree of alpha tree
    std::vector<int> alphaParents;  // computing dandrogram e.g. the alpha tree

    int rlevel =0;  // a variable to control recursion

    // It can be constructed using either incidence matrix or an eulerTour
    // eulerTourTree_t& parentMST;     // stores the parent MST
    rAlphaTree_t() {}; 
    rAlphaTree_t (incidenceMatrix_t &_incMatMST); 


    // It can be constructed using either incidence matrix or an eulerTour
    rAlphaTree_t (eulerTourTree_t& inTour,int level); 

    // computeBranch: find which branch of alphaTree inputEdge \in parentGraph belongs
    int computeBranch(int pEdge, eulerTourTree_t& pMST);
    int findAlphaParent(int edgeId, eulerTourTree_t& pMST);

private:
    int heightAlphaTree;
    int logHAlpha;
    std::vector<std::vector<int>> alphaAncestors;  // alpha poly ancestors log(Height)*numAlphaEdges
    // std::vector<std::vector<int>> computePolyAncestors();
    int alphaTreeBottomUpParentSearch(int edgeId, int maxDescIdx, 
        std::vector<int> alphaEdges, std::vector<std::vector<int>> alphaAncestors);


};


std::vector<int> computeDandrogramUsingAlphaTree(eulerTourTree_t& inTour, rAlphaTree_t& alphaTree);
int branch2AlphaEdge(int branchId);
int alphaTreeBottomUpParentSearch(int edgeId, int maxDescIdx, 
    std::vector<int>& alphaEdges, 
    std::vector<std::vector<int>>& alphaAncestors);