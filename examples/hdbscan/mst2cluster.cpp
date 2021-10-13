#include <iostream>
#include <vector>
#include <fstream> 
#include <string>
// #include "hdbscan.hpp"
#include "parallel_boruvka.hpp"
#include "edge_hierarchy.hpp"


std::vector<wtEdge_t> readMST(std::ifstream& inFile)
{
    
    std::vector<wtEdge_t > out;
    int count =0; 
    while(1)
    {
        double dsrc;
        double ddest;
        double wt;   
        inFile>>dsrc;
        inFile>>ddest;
        inFile>>wt;
        if(inFile.eof()) break; 
        // std::cout<<dsrc<< " "<<ddest<<" " <<wt << "\n";
        out.push_back( std::make_pair(std::make_pair((int)dsrc, (int)ddest), wt) ); 
        count++;
    }
        
        std::cout<<"Read MST with \t: "<<count << " edges.\n";
        
    return out; 
}

int main(int argc, char* argv[])
{

    // read the points 
    if(argc<2)
    {
        std::cout<<" usage: "<<argv[0]<<" MSTFile minclusterSize(integer, optional)\n";
        return -1;
    }
    int minClusterSize =2; 
    if(argc>2)
    {
        minClusterSize = atoi(argv[2]);
    }
    std::cout<<"Opening "<<argv[1]<<"\n";
    std::ifstream infile(argv[1]);

    
    // creating edge_Hierarchy 
    auto wtMST = readMST(infile);
    std::sort(wtMST.begin(),wtMST.end(),compareWtEdges() );
    edgeHierarchy_t edgeHrr(wtMST,minClusterSize);
    
    if(0)
    {
        std::string outFileName(argv[1]);
        outFileName = outFileName + ".grp";
        std::cout<<"writing groups to: "<<outFileName<<"\n";
        std::ofstream outFile(outFileName);
        edgeHrr.writeGroups(outFile);
    }

    if(1)
    {
        std::string outFileName(argv[1]);
        outFileName = outFileName + ".dot";
        std::cout<<"writing EdgeTree in dotformat to : "<<outFileName<<"\n";
        std::ofstream outFile(outFileName);
        edgeHrr.writeMSTdot(outFile);
    }
// edgeHierarchy_t::writeClusterMaps(std::ofstream &ofileName)

    if(1)
    {
        std::string outFileName(argv[1]);
        outFileName = outFileName + ".map";
        std::cout<<"writing Clusters  to : "<<outFileName<<"\n";
        std::ofstream outFile(outFileName);
        edgeHrr.writeClusterMaps(outFile);
    }
    // perform hdbscan

    // return output
    return 0;
}