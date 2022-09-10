#include <iostream>
#include <cassert>
#include "alphaTree.hpp"
#include "recursiveAlphaTree.hpp"

std::vector<int>& eulerTourTree_t::computemaxIncidentOnBridge()
{
    if(this->maxIncidentOnBridge.size()) // already computed 
        return this->maxIncidentOnBridge; 
    
    int numEdges = this->nEdges();
    this->maxIncidentOnBridge.resize(numEdges);
    /*  === Initialize maxIncidentOnBridge === */
    for(int idx=0; idx < numEdges; idx++)
        this->maxIncidentOnBridge[idx] = idx; 
    /*  === Compute maxIncidentOnBridge === */
    for(int idx=0; idx<numEdges; idx++)
    {
        int bridge= bridgeEdges[idx];
        assert(bridge>=-1);
        if (bridge == -1)
        {
            maxIncidentOnRoot = std::max(maxIncidentOnRoot, idx);
        }
        else 
            this->maxIncidentOnBridge[bridge] = std::max(this->maxIncidentOnBridge[bridge], idx);
    }

    return this->maxIncidentOnBridge; 
}

int eulerTourTree_t::maxIdxBridge(int k)
{
    if(!this->maxIncidentOnBridge.size())
    {
        std::cout<<" MaxIncidentOnBridge is not computed !!!! \n";
        exit(0);
    } 
    // assert( !(k>=0 and k<this->maxIncidentOnBridge.size()));
    if(k==-1)
        return this->maxIncidentOnRoot;
    else 
        return this->maxIncidentOnBridge[k];
    // return this->maxIncidentOnBridge[k]; 
    
}

std::vector<int> eulerTourTree_t::computeAlphaEdges()
{
    // Iam an alphaEdge if I have an incident edge on both sides
    // if Iam < maxIncidentEdge(Iam) (left side)
    // and Iam < maxIncident(bridge(Iam)) U bridge(Iam)

    /*
    *   STEP-0: compute maxIncident
    **/ 
    int numEdges = this->nEdges();
    #if 0
        std::vector<int> maxIncident(numEdges);
        /*  === Initialize maxIncident === */
        for(int idx=0; idx < numEdges; idx++)
            maxIncident[idx] = idx; 
        /*  === Compute maxIncident === */
        for(int idx=0; idx<numEdges; idx++)
        {
            int bridge= bridgeEdges[idx];
            maxIncident[bridge] = std::max(maxIncident[bridge], idx);
        }
    #endif 
    std::vector<int> maxIncident= computemaxIncidentOnBridge();

    /*
    *   STEP-1: Compute alpha edges
    **/ 
    
    std::vector<int> alphaEdges; 
    for(int idx=0; idx<numEdges; idx++)
    {
        int bridge= bridgeEdges[idx];
        int maxIncidentOnMyBridge = bridge<0? maxIncidentOnRoot: maxIncident[bridge];
        if(idx < maxIncident[idx] and idx< maxIncidentOnMyBridge )
            alphaEdges.push_back(idx);
    }
        

    return alphaEdges; 
}


std::pair<int,int> eulerTourTree_t::simInsert(int pEntry, int pExit)
{
    #warning Deprecated: use computeSimBridge instead
    indexedTimeStamp_t idxEntry(pEntry, 0);
    indexedTimeStamp_t idxExit(pExit, 0);
    int myEntry = std::distance(
        idxTS.begin(), std::lower_bound(idxTS.begin(), idxTS.end(), idxEntry));
    int myExit = std::distance(
        idxTS.begin(), std::lower_bound(idxTS.begin(), idxTS.end(), idxExit));
    
    return std::make_pair(myEntry, myExit);
}

int eulerTourTree_t::simBridge(int lEntry)
{
    #warning Deprecated: use computeSimBridge
    int myBridge ;

    if (lEntry == 0)
        myBridge = -1;
    else if (idxTS[lEntry].second % 2 == 0) // starting point
        #warning check the following correction
        // myBridge = alphaBridgeEdges[idxTS[lEntry].second / 2];
        myBridge = bridgeEdges[idxTS[lEntry].second / 2];
    else
        myBridge = idxTS[lEntry].second / 2;
    
    return myBridge;
}


int eulerTourTree_t::computeSimulatedBridge(int pEntry, int pExit)
{
    // Replaces simInsert and simBridge with a single call
    indexedTimeStamp_t idxEntry(pEntry, 0);
    indexedTimeStamp_t idxExit(pExit, 0);
    int myBridge;
    int myEntry;
    int myExit;
    // int firstTimeStamp = idxTS[0].first;
    int lastTimeStamp = idxTS[idxTS.size() - 1].first;
    if (pEntry > lastTimeStamp)
    {
        myBridge = -1;
        myEntry = idxTS.size();
        myExit = idxTS.size();
        return -1; 
    }
    else
    {
        myEntry = std::distance(
            idxTS.begin(), std::lower_bound(idxTS.begin(), idxTS.end(), idxEntry));
        myExit = std::distance(
            idxTS.begin(), std::lower_bound(idxTS.begin(), idxTS.end(), idxExit));

        assert(myEntry >= 0);
        assert(myEntry < idxTS.size());
        // if (!(myEntry < idxTS.size()))
        // {
        //     printf("edgeId%d, (myEntry %d myExit %d) idxTS.size() %d, lastTS %d, "
        //             "myTs =%d \n",
        //             edgeId, myEntry, myExit, idxTS.size(),
        //             idxTS[idxTS.size() - 1].first, pEntry);
        // }
    // assert(myEntry < idxTS.size()); // this is failing
    /*
    *  The following if the inserted edge is an alpha edge 
    */
    // if (edgeId == alphaEdges[idxTS[myEntry].second / 2])
    // {
    //   return alphaParents[idxTS[myEntry].second / 2];
    // }

    // if (printEnabled)
    //   printf("%d :\t mstEntry %d mstExit %d newEntry %d newExit %d: \
	// 	\n \t\t%d  %d\n",
    //          edgeId, mstEulerEntry[edgeId], mstEulerExit[edgeId], myEntry,
    //          myExit, alphaEdges[idxTS[myEntry].second / 2],
    //          idxTS[myEntry].second);

    if (myEntry == 0)
        myBridge = -1;
    else if (idxTS[myEntry].second % 2 == 0) // starting point
        myBridge = bridgeEdges[idxTS[myEntry].second / 2];
    else
        myBridge = idxTS[myEntry].second / 2;
    return myBridge;
    }
    // indexedTimeStamp_t idxEntry(pEntry, 0);
    // indexedTimeStamp_t idxExit(pExit, 0);
    // int myEntry = std::distance(
    //     idxTS.begin(), std::lower_bound(idxTS.begin(), idxTS.end(), idxEntry));
    // int myExit = std::distance(
    //     idxTS.begin(), std::lower_bound(idxTS.begin(), idxTS.end(), idxExit));
    
    // return std::make_pair(myEntry, myExit);
}


eulerTourTree_t computeAlphaEulerTour(std::vector<int> listOfAlphaEdges, eulerTourTree_t& inMST)
{
    eulerTourTree_t alphaTour;
    alphaTour.idxTS = createAlphaIdxTimeStamp(listOfAlphaEdges, inMST.entryTime, inMST.exitTime);
    auto alphaEulerEntryExit = computeAlphaEulerEntryExit(alphaTour.idxTS);
    alphaTour.entryTime = alphaEulerEntryExit.first;
    alphaTour.exitTime = alphaEulerEntryExit.second;

    // tAlphaEulerTour.end();
    // compute bridgeEdge
    alphaTour.bridgeEdges = computeBridgeEdges(alphaTour.idxTS , alphaTour.entryTime, alphaTour.exitTime);
    alphaTour.computemaxIncidentOnBridge();

    // return 
    return alphaTour; 
  
}

bool eulerTourTree_t::isLeftOfEdge( int queryEdge, int rootEdge)
{
  if (entryTime[queryEdge] > entryTime[rootEdge] &&
      exitTime[queryEdge] < exitTime[rootEdge])
    return 1;

  return 0;
}

eulerTourTree_t::eulerTourTree_t(incidenceMatrix_t& incMatMST,std::vector<wtEdge_t>& wtSortedMST)
{
    // Compute the euler Tour
    int numEdges = wtSortedMST.size();
    auto eulerEntryExit = eulerTour(incMatMST, wtSortedMST);
    this->entryTime = eulerEntryExit.first;
    this->exitTime = eulerEntryExit.second;
    // Compute the indexTimeStamp
    this->idxTS.resize(2*numEdges);
    for (int i = 0; i < numEdges; i++)
    {
        
        idxTS[2 * i].first = this->entry(i);
        idxTS[2 * i].second = 2 * i;
        idxTS[2 * i + 1].first = this->exit(i);
        idxTS[2 * i + 1].second = 2 * i + 1;
    }

    std::sort(idxTS.begin(), idxTS.end());

    this->bridgeEdges = computeBridgeEdges(this->idxTS , this->entryTime, this->exitTime);
    this->computemaxIncidentOnBridge();
}