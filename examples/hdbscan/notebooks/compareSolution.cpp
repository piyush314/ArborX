#include <iostream>
#include <vector>
#include <map>
#include "boost/python/numpy.hpp"

namespace py = boost::python;
namespace np = boost::python::numpy;


int compareSolution(std::vector<int> A, std::vector<int> B)
{
    if(A.size() != B.size())
    {
        std::cout<<"Array size mismatch" << "\n";
        return -1;
    }

    std::map<int, int> A2B;
    int count =0; 

    for(int i=0; i< A.size(); i++)
    {
        if(A2B.find(A[i]) == A2B.end() )
        {
            A2B[A[i]] =B[i];
        }
        else
        {
            if(A2B[A[i]] !=B[i])
            {
                count++; 
                std::cout<<"Error at "<<i<<"\n";
            }
        }  
    }

    return count;
}


int compareHDBSCAN()
{
    auto wtMST = readMST(infile);
    std::sort(wtMST.begin(),wtMST.end(),compareWtEdges() );
    edgeHierarchy_t edgeHrr(wtMST,minClusterSize);
}