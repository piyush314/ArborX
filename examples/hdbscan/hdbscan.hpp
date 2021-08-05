#pragma once
#include <vector>
#include <algorithm>
#include "parallel_boruvka.hpp"
#include "incidence_matrix.hpp"

typedef std::pair<edge_t, double> wtEdge_t;

class hdbscan_t
{
    int m_npts, m_dim; //
    int m_kpts;        // used for HDBSCAN core distance computation
    std::vector<std::vector<double>> m_points;
    std::vector<wtEdge_t> wtMST;
    incidenceMatrix_t m_inMatMST;

public:
    hdbscan_t(std::vector<std::vector<double>> &pts, int k_pts);
};



hdbscan_t::hdbscan_t(std::vector<std::vector<double>> &pts, int k_pts) : m_points(pts), m_kpts(k_pts)
{
    m_npts = m_points.size();
    m_dim = m_points[0].size();

    parallelBoruvka_t pBv(pts);
    wtMST = pBv.weightedMST();

    // sort the edges 
    std::sort(wtMST.begin(),wtMST.end(),compareWtEdges() );
    // create the incidence matrix 
    incidenceMatrix_t icmat(wtMST);
    // m_inMatMST = icmat;
};