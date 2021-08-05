#include <iostream>
#include <vector>
#include <fstream> 
#include <string>
// #include "hdbscan.hpp"
#include "parallel_boruvka.hpp"
#include "edge_hierarchy.hpp"

int getWriteMST()
{
    //TODO: use getnev 
    return 1; 
}
std::vector<std::vector<double> > readFile(std::ifstream& inFile)
{
    char c;
    int npts;
    int dim; 
    inFile>>c;
    inFile>>npts;
    inFile>>dim;
    std::cout<<c<<" number of pts is "<<npts <<" and dimemsion is "<<dim<<"\n";
    std::vector<std::vector<double> > out(npts, std::vector<double>(dim, 0.0) );
    for(int i=0; i<npts; i++)
        for(int j=0; j<dim; j++)
            inFile>>out[i][j];
        
    return out; 
}

int main(int argc, char* argv[])
{

    // read the points 
    if(argc<2)
    {
        std::cout<<" usage: "<<argv[0]<<" inputFile minclusterSize(integer, optional)\n";
        return -1;
    }
    int minClusterSize =2; 
    if(argc>2)
    {
        minClusterSize = atoi(argv[2]);
    }
    std::cout<<"Opening "<<argv[1]<<"\n";
    std::ifstream infile(argv[1]);

    std::vector<std::vector<double> > points = readFile(infile);
    int kpts =3;
    // hdbscan_t hdbscan(points,kpts);
    parallelBoruvka_t pB(points);
    // writing MST back to file 
    
    if(getWriteMST())
    {
        std::string outFileName(argv[1]);
        outFileName = outFileName + ".mst";
        std::cout<<"writing to: "<<outFileName<<"\n";
        std::ofstream outFile(outFileName);
        pB.writeMST(outFile);
    }
    
    // creating edge_Hierarchy 
    auto wtMST = pB.weightedMST();
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
        std::cout<<"writing MST in dotformat to : "<<outFileName<<"\n";
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