import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn.datasets as data
import subprocess as sub 
import pandas as pd 
import hdbscan


def plotCluster(binaryExe, dataFileName, minClusterSize=2):
    dataPoints = np.loadtxt(dataFileName, dtype="double", comments="#") 
    argString = dataFileName + " %d"%(minClusterSize)
    print(argString)
    cmd  = binaryExe + " " + dataFileName + " %d"%(minClusterSize)
    try:
        print("Running: " + cmd)
        out = sub.check_output(cmd,shell=True)
        
    except sub.CalledProcessError as e:
        print(e.output)
    # sub.run([binaryExe, argString ], capture_output=True)
    mstFile = dataFileName+".mst"
    mst = np.loadtxt(mstFile, dtype="int")
    ### Now plot the output 
    fig, ax = plt.subplots( figsize=(6,6))
    ax.set_aspect('equal')
    # ax.scatter(dataPoints.T[0], dataPoints.T[1], color='b', **plot_kwds)
    
    ## Annotate 
    npts = dataPoints.shape[0]
    # for pts in range(npts):
    #     str = "%d"%(pts)
    #     ax.annotate(str, dataPoints[pts],fontsize=6)
    # xlim = (-1.5, 2.5)
    # ylim = (-1,3)
    # ax.set_xlim(xlim)
    # ax.set_ylim(ylim)
    ax.set_xlabel("X-coordinate")
    ax.set_ylabel("Y-coordinate")
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.rc('font', size=10  )
    plt.rcParams.update({'font.size': 8})
    # for edge in mst:
    for idx in range(len(mst)):
        edge = mst[idx]
        x_values = [dataPoints[edge[0]][0], dataPoints[edge[1]][0]]
        y_values = [dataPoints[edge[0]][1], dataPoints[edge[1]][1]]
        ax.plot(x_values, y_values, linewidth="1")
        # str = r"$e_{%d}$"%(idx)   # takes too long 
        # str = "e%d"%(idx)
        # dp = [np.sum(x_values)/2, .05+np.sum(y_values)/2]
        # ax.annotate(str, dp,fontsize=6, color='teal')
        # print(edge)dataFileName
    clstFile = dataFileName+".map"
    clst = np.loadtxt(clstFile, dtype="int")
    # print(dataPoints.T[0])
    df = pd.DataFrame(dataPoints)
    df.columns = ["X-coord","Y-coord"]
    df["Cluster"] = clst 
    sns.set_context("paper", rc={"font.size":12,"axes.labelsize":10})
    sns.scatterplot(data=df, x="X-coord", y="Y-coord", hue="Cluster", style="Cluster", palette="tab10", ax=ax)
    plt.legend(loc='lower right')
    plt.show()

def plotSktHdbscan(dataFileName, min_cluster_size=5, saveFig=False, writeMST=False):
    clusterer = hdbscan.HDBSCAN(min_cluster_size, gen_min_span_tree=True)
    dataPoints = np.loadtxt(dataFileName, dtype="double", comments="#")
    clusterer.fit(dataPoints)
    clst = clusterer.labels_
        # print(dataPoints.T[0])
    df = pd.DataFrame(dataPoints)
    df.columns = ["X-coord","Y-coord"]
    df["Cluster"] = clst 

    figs, axs = plt.subplots(2,2, figsize=(10,10))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.rc('font', size=10  )
    sns.set_context("paper", rc={"font.size":12,"axes.labelsize":10})

    sns.scatterplot(data=df, x="X-coord", y="Y-coord", hue="Cluster", 
        style="Cluster",  ax=axs[0,0],palette="tab10")
    plt.legend(loc='lower right')
    plt.sca(axs[0,1])
    clusterer.minimum_spanning_tree_.plot(edge_cmap='viridis', 
                                        edge_alpha=0.6, 
                                        node_size=10, 
                                        edge_linewidth=1)
    plt.sca(axs[1,0])
    clusterer.single_linkage_tree_.plot(cmap='viridis', colorbar=True)
    plt.sca(axs[1,1])
    # clusterer.condensed_tree_.plot()
    clusterer.condensed_tree_.plot(select_clusters=True, selection_palette=sns.color_palette())
    if(saveFig):
        print("Saving outputplot to "+outFile)
        outfile = dataFileName+"-skt.pdf"
        plt.savefig(outfile)
    plt.show()
    if(writeMST):
        outFile = dataFileName +"-skt.mst"
        print("Writing mst to "+outFile)
        np.savetxt(outFile, clusterer.minimum_spanning_tree_._mst, 
            fmt='%.6f', delimiter='\t', newline='\n', footer='', comments='# ', encoding=None)
    return clst; 