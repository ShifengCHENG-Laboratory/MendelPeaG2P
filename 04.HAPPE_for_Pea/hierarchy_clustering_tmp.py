## clustering using pheat file
## usage: python __file__ in.newpheat out.newick out.cluster

# in.pheat format
# name ref alt  sample1 sample2 sample3 ...
# chr1A_131 A G 0 9 -9 ...

import pandas as pd
import numpy as np
from scipy.cluster import hierarchy
import sys
from sklearn import cluster
sys.setrecursionlimit(100000)

def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick
cluster_res = [] #[ [],... ]

def clustering(node,  parentdist, leaf_names,threshold):
    if node.is_leaf():
        return ( [ leaf_names[node.id] ] , parentdist - node.dist)
    else:
        left_guys,left_branch_len = clustering(node.get_left(),  node.dist, leaf_names,threshold)
        right_guys,right_branch_len = clustering(node.get_right(),  node.dist, leaf_names,threshold)
        # print(left_branch_len,right_branch_len,left_branch_len+right_branch_len , len(left_guys+right_guys) )
        if left_branch_len+right_branch_len > threshold:
            cluster_res.append(left_guys)
            cluster_res.append(right_guys)
            # print(len(left_guys))
            # print(len(right_guys))
            return ([], parentdist - node.dist)
        else:
            return (left_guys+right_guys, parentdist - node.dist)


def get_linkage_matrix(model):
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count


    linkage_matrix = np.column_stack([model.children_, model.distances_, counts]).astype(float)
    return linkage_matrix

if __name__ == "__main__":
    # read table file
    inf = open(sys.argv[1],"r")
    header = inf.readline().strip().split()
    sample_gt_d = {}
    for sample in header[3:]:
        sample_gt_d[sample]  = []

    for line in inf:
        ls = line.strip().split()
        for index,value in enumerate( ls[3:] ):
            if value == "NA":
                value = 99
            sample_gt_d[header[ index + 3] ] .append(value)

    inf.close()
    # print(sklearn.__version__)
    pd_df = pd.DataFrame(sample_gt_d)
    col_name = pd_df.columns.tolist()
    pd_df = pd_df.T

    hclustering = cluster.AgglomerativeClustering(distance_threshold=0, n_clusters=None).fit(pd_df)
    # print(len(clustering.labels_))
    Z = get_linkage_matrix(hclustering)
    tree = hierarchy.to_tree(Z, rd=False)
    newick_str = getNewick(tree, "", tree.dist, col_name )

    remains_guys,branch_len = clustering(tree, tree.dist, col_name,1 )
    cluster_res.append(remains_guys)
    # remove empty cluster res
    cluster_res = [x for x in cluster_res if len(x) > 0]
    #output cluster res
    outf = open(sys.argv[3],"w")
    outf.write("sample\tcluster\n")
    cluster_num = 0
    for cluster in cluster_res:
        cluster_num += 1
        for sample in cluster:
            outf.write("%s\t%d\n"%(sample,cluster_num))



    outf.close()
    ouf = open(sys.argv[2],"w")
    ouf.write(newick_str)
    ouf.close()
