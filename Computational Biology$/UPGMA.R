#########################################################
#######        COMPUTATIONAL BIOLOGY         ############
#######             HOMEWORK 3               ############
#########################################################
#                                                       #
# Reconstruct the phylogenetic tree of given sequences  #
# using UPGMA with hamming, JC69, K80 distances.        #
#                                                       #
#########################################################
#########################################################

#########################################################
######    Code below this should not be changed   #######
#########################################################

library(ape)

transform_to_phylo = function(sequences, named_edges, edge_lengths, node_description) {
    # Produce a tree of the phylo class from the matrix of edges, vector of lengths and
    # the dataframe  containing node descriptions.
    #    sequences: a list of the original sequences as strings;
    #    named_edges: an Mx2 matrix of pairs of nodes connected by an edge,
    #                 where the M rows are the different edges and the 2 columns
    #                 are the parent node and the child node of an edge as numeric values.
    #    edge_lengths: a vector of length M of the corresponding edge lengths as numeric values.
    #    node_description: a data frame of the node descriptions as defined in
    #                      initialize_node_description and extended by the upgma code.
    
    edges = named_edges
    for (name in rownames(node_description)) {
        index = which(rownames(node_description) == name)
        edges[which(edges == name)] = as.numeric(index)
    }
    edges = matrix(as.numeric(edges), ncol = 2)
    
    edges[which(edges > length(sequences))] = - edges[which(edges > length(sequences))]
    root = setdiff(edges[,1], edges[,2])
    edges[which(edges==root)] = length(sequences) + 1
    
    k = length(sequences) + 2
    for(x in unique(edges[which(edges < 0)])) {
        edges[which(edges==x)] = k
        k = k + 1
    }
    
    tree = list()
    class(tree) = "phylo"
    tree$edge = edges
    tree$edge.length = edge_lengths
    tree$Nnode = as.integer(length(sequences) - 1)
    tree$tip.label = names(sequences)
    
    # Return the tree in the form of the phylo class from ape
    return(tree)
}

plot_tree = function(tree) {
    # Plot the phylogenetic tree with node labels and edge lengths.
    #    tree: an object of class phylo from the ape package

    plot(tree)
    edgelabels(format(tree$edge.length, digits = 2))
}

initialize_node_description = function(sequences) {
    # Initialize the structure that will hold node descriptions.
    # The created structure is a data frame where the rows are node names, and the columns are
    # node_height -- distance from the node to any tip in the ultrametric tree as a numeric value, and
    # node_size -- number of tips that this node is ancestral to as a numeric value.
    #    sequences: a list of the original sequences as strings;

    N = length(sequences)
    node_names = names(sequences)
    node_sizes = rep(1, times = N)
    node_heights = rep(0, times = N)
    node_description = data.frame(node_sizes, node_heights)
    rownames(node_description) = node_names

    # Return a data frame that contains information on currently existing tip nodes.
    # node_description: a dataframe containing information on the currently existing nodes.
    #                   The row names are the names of the currently existing tips, i.e.
    #                   are the same as the names in the sequence list, node_height is
    #                   0 and node_size is 1 as the all the currently existing nodes are tips.
    return(node_description)
}

add_new_node = function(node_description, merging_nodes) {
    # Add new merged node to the node description data frame.
    # The new node is a combination of the nodes supplied in the merging_nodes,
    # e.g. if one needs to merge nodes "bird" and "fish", the new node in the
    # dataframe will be called "bird.fish".
    #    node_description: the dataframe created by initialize_node_description, containing
    #                      current node sizes and node heights as numeric values.
    #    merging_nodes: a vector of two names of the nodes being merged as strings.

    new_node_name = paste(merging_nodes, collapse = ".")
    new_node_row = data.frame(node_sizes = 0, node_heights = 0)
    rownames(new_node_row) = new_node_name
    new_node_description = rbind(node_description, new_node_row)
    
    # Return the node_description dataframe with a row for the new node added, and
    # the new node name.
    #    node_description: the dataframe where the rows are labelled by current nodes and columns
    #                      contain the node heights and sizes as numeric values.
    #    new_node_name: the name of the newly added node as a string, created from names in merging_nodes.
    return(list(node_description = new_node_description, new_node_name = new_node_name))
}

#########################################################
######    Code above this should not be changed   #######
#########################################################

get_hamming_distance = function(sequence1, sequence2) {
    # Compute the Hamming distance between two sequences.
    #    sequence1: first sequence (string)
    #    sequence2: second sequence (string)
    distance=0
    for (i in 1:nchar(sequence1)){
      if (substr(sequence1, i, i)!=substr(sequence2, i, i)){
        distance=distance+1
      }
    }

    # Return the numerical value of the distance
    return(distance)
}

get_JC69_distance = function(sequence1, sequence2) {
    # Compute the JC69 distance between two sequences.
    #    sequence1: first sequence (string)
    #    sequence2: second sequence (string)
    hamming=get_hamming_distance(sequence1,sequence2)
    distance=(-3/4)*log(1-(4/3)*(hamming/nchar(sequence1)))

    # Return the numerical value of the distance
    return(distance)
}

get_K80_distance = function(sequence1, sequence2) {
    # Compute the K80 distance between two sequences.
    #    sequence1: first sequence (string)
    #    sequence2: second sequence (string)
    S=0
    V=0

    for (i in 1:nchar(sequence1)){
      if (substr(sequence1, i, i)!=substr(sequence2, i, i)){
        #transitions
        if (substr(sequence1, i, i)=="A" & substr(sequence2, i, i)=="G"){
          S=S+1
        }
        else if (substr(sequence1, i, i)=="G" & substr(sequence2, i, i)=="A"){
          S=S+1
        }
        else if (substr(sequence1, i, i)=="T" & substr(sequence2, i, i)=="C"){
          S=S+1
        }
        else if (substr(sequence1, i, i)=="C" & substr(sequence2, i, i)=="T"){
          S=S+1
        }
        #transversion
        else if (substr(sequence1, i, i)=="A" & substr(sequence2, i, i)=="C"){
          V=V+1
        }
        else if (substr(sequence1, i, i)=="C" & substr(sequence2, i, i)=="A"){
          V=V+1
        }
        else if (substr(sequence1, i, i)=="T" & substr(sequence2, i, i)=="A"){
          V=V+1
        }
        else if (substr(sequence1, i, i)=="A" & substr(sequence2, i, i)=="T"){
          V=V+1
        }
        else if (substr(sequence1, i, i)=="C" & substr(sequence2, i, i)=="G"){
          V=V+1
        }
        else if (substr(sequence1, i, i)=="G" & substr(sequence2, i, i)=="C"){
          V=V+1
        }
        else if (substr(sequence1, i, i)=="T" & substr(sequence2, i, i)=="G"){
          V=V+1
        }
        else if (substr(sequence1, i, i)=="G" & substr(sequence2, i, i)=="T"){
          V=V+1
        }
      }
    }
    S=S/nchar(sequence1)
    V=V/nchar(sequence1)
    distance=(-0.5)*log(1-2*S-V)-(0.25)*log(1-2*V)

    # Return the numerical value of the distance
    return(distance)
}

compute_initial_distance_matrix = function(sequences, distance_measure) {
    # Compute the initial distance matrix using one of the distance measures.
    # The matrix is of dimension NxN, where N is the number of sequences.
    # The matrix columns and rows should be labelled with tip names, each row and column
    # corresponding to the appropriate sequence.
    # The matrix can be filled completely (i.e symmetric matrix) or only the upper half
    # (as shown in the lecture).
    # The diagonal elements of the matrix should be Inf.
    #    sequences: the sequence alignment in the form of a list of species names and 
    #               the associated genetic sequences as strings
    #    distance_measure: a string indicating whether the 'hamming', 'JC69' or 'K80' 
    #                      distance measure should be used
  
    N <- length(sequences)
    distance_matrix <- matrix(nrow = N, ncol = N)
    diag(distance_matrix) <- Inf
    rownames(distance_matrix) <- names(sequences)
    colnames(distance_matrix) <- names(sequences)

    for (i in 1:length(sequences)) {
      for (j in 1:length(sequences)){
        if (i != j){
          sequence1=sequences[i]
          sequence2=sequences[j]
          if (distance_measure=='hamming'){
            distance_matrix[i,j]=get_hamming_distance(sequence1,sequence2)  
          }
          else if (distance_measure=='JC69'){
            distance_matrix[i,j]=get_JC69_distance(sequence1,sequence2)
          }
          else if (distance_measure=='K80'){
            distance_matrix[i,j]=get_K80_distance(sequence1,sequence2)
          }
        }
      }
    }
        

    # Return the NxN matrix of numeric inter-sequence distances with Inf on the diagonal
    return(distance_matrix)
}

get_merge_node_distance = function(node_description, distance_matrix, merging_nodes, existing_node) {
    # Compute the new distance between the newly created merge node and an existing old node in the tree
    #    node_description: a dataframe containing information on the currently existing nodes
    #    distance_matrix: the matrix of current distances between nodes
    #    merging_nodes: a vector of two node names that are being merged in this step
    #    existing_node: one of the previously existing nodes, not included in the new node
  
    name1=merging_nodes[1]
    name2=merging_nodes[2]
    for (i in 1:nrow(distance_matrix)){
      if (rownames(distance_matrix)[i]==existing_node){
        row=i
      }
    }
    for (j in 1:ncol(distance_matrix)){
      if (colnames(distance_matrix)[j]==name1){
        dist1=distance_matrix[row,j]
        node_size1=node_description[j,1]
      }
      if (colnames(distance_matrix)[j]==name2){
        dist2=distance_matrix[row,j]
        node_size2=node_description[j,1]
      }
    }
    
    new_distance=(dist1*node_size1+dist2*node_size2)/(node_size1+node_size2)
  
    # Returns the numeric distance between the newly created merge node and the existing node
    return(new_distance)
}

update_distance_matrix = function(node_description, distance_matrix, merging_nodes, new_node_name) {
    # Update the distance matrix given that two nodes are being merged.
    #    node_description: a dataframe containing information on the currently existing nodes
    #    distance_matrix: the current distance matrix that needs to be updated
    #    merging_nodes: a vector of two node names that need to be merged in this step
    #    new_node_name: the name with which the merged node should be labelled
    # The resulting matrix should be one column and one row smaller, i.e. if the given distance matrix
    # was MxM, then the updated matrix will be M-1xM-1, where the 2 rows and cols represent the separate
    # nodes undergoing the merge are taken out and a new row and col added that represents the new node.

    name1=merging_nodes[1]
    name2=merging_nodes[2]
    new_node=c()
    for (i in 1:ncol(distance_matrix)){
      if (rownames(distance_matrix)[i]!=name1 & rownames(distance_matrix)[i]!=name2){
        existing_node=rownames(distance_matrix)[i]
        new_node=append(new_node, get_merge_node_distance(node_description, distance_matrix, merging_nodes, existing_node))
      }else{
        new_node=append(new_node, 0)
      }  
    }
    updated_distance_matrix=rbind(distance_matrix,new_node)
    new_node2=append(new_node, Inf, (length(new_node)+1))
    updated_distance_matrix=cbind(updated_distance_matrix,new_node2)
    n=ncol(updated_distance_matrix)
    rownames(updated_distance_matrix)[n] = new_node_name 
    colnames(updated_distance_matrix)[n] = new_node_name
    for (i in 1:ncol(updated_distance_matrix)){
      if (colnames(updated_distance_matrix)[i]==name1){
        index1=i
      }
      if (colnames(updated_distance_matrix)[i]==name2){
        index2=i
      }
    }
    
    if (nrow(updated_distance_matrix)>3){
      updated_distance_matrix = updated_distance_matrix[-c(index1,index2),]
      updated_distance_matrix = updated_distance_matrix[,-c(index1,index2)]
    } else {
      updated_distance_matrix = updated_distance_matrix[-c(index1,index2),]
      updated_distance_matrix = updated_distance_matrix[-c(index1,index2)]
    }
    
    # Returns the updated matrix of numeric cluster distances
    return(updated_distance_matrix)
}


upgma_one_step = function(node_description, distance_matrix, edges, edge_lengths) {
    # Performs one step of the upgma algorithm, i.e. the nodes with the smallest distance are merged, 
    # the node height of the new node is calculated and the distance matrix is newly calculated.
    # Values that are expected to be returned are listed below.
    #    node_description: a dataframe containing information on the currently existing nodes
    #    distance_matrix: the current distance matrix that needs to be updated (LxL)
    #    edges: an Mx2 matrix of pairs of nodes connected by an edge, where the M rows are the different
    #           edges and the 2 columns are the parent node and the child node of an edge.
    #    edge_lengths: a vector of length M of the corresponding edge lengths.
  
    min_dist=Inf
    for (i in 1:nrow(distance_matrix)){
      for (j in i:ncol(distance_matrix)){
        if (distance_matrix[i,j]<min_dist){
          min_dist=distance_matrix[i,j]
          row=i
          col=j
        }
      }
    }
    merging_nodes=c(colnames(distance_matrix)[row],colnames(distance_matrix)[col])
    new_node=add_new_node(node_description, merging_nodes)
    new_node_name=new_node[[2]]
    node_description=new_node[[1]]
    n=nrow(node_description)
    for (i in 1:nrow(node_description)){
      if (rownames(distance_matrix)[col]==rownames(node_description)[i]){
        ind1=i
      }
      if (rownames(distance_matrix)[row]==rownames(node_description)[i]){
        ind2=i
      }
    }
    node_description[n,1]=node_description[ind2,1]+node_description[ind1,1]  # fill node size for new node
    node_description[n,2]=distance_matrix[row,col]/2 #fill in node height of new node
    edge_new1=distance_matrix[row,col]/2-node_description[ind1,2]
    edge_new2=distance_matrix[row,col]/2-node_description[ind2,2]
    edge_lengths=append(edge_lengths, edge_new1)
    edge_lengths=append(edge_lengths, edge_new2)
    for (i in 1:nrow(node_description)){
      if (rownames(distance_matrix)[col]==rownames(node_description)[i]){
        child1_index=i
      }
      if (rownames(distance_matrix)[row]==rownames(node_description)[i]){
        child2_index=i
      } 
    }
    parent_index=rownames(node_description)[n]
    edges=rbind(edges, c(parent_index, child1_index))
    edges=rbind(edges, c(parent_index, child2_index))
    distance_matrix=update_distance_matrix(node_description, distance_matrix, merging_nodes, new_node_name)
  

  
    # Return the updated distance matrix, edge description matrix, edge length vector and 
    # node_description data frame
    #    node_description: data frame containing sizes and heights of all nodes 
    #        (to be updated using add_new_node())
    #    distance_matrix: the updated matrix of numeric distances between nodes (L-1xL-1)
    #    edges: an Mx2 matrix of pairs of nodes connected by an edge, where the M rows are
    #    the different edges and the 2 columns are the parent node and the child node of an edge as numeric values.
    #    edge_lengths: a vector of length M of the corresponding numeric edge lengths.
    return(list(node_description = node_description, distance_matrix = distance_matrix,
                edges = edges, edge_lengths = edge_lengths))
}

build_upgma_tree = function(sequences, distance_measure) {
    # Build the tree from given sequences using the UPGMA algorithm.
    #    sequences: the sequences in the format of a list of species names and the associated genetic 
    #               sequences as strings.
    #    distance_measure: a string indicating whether the 'hamming', 'JC69' or 'K80' distance measure
    #                      should be used
    N <- length(sequences)
    node_description <- initialize_node_description(sequences)
    edges <- matrix(nrow = 0, ncol = 2)
    edge_lengths <- vector(mode = "numeric", length = 0)
    distance_matrix=compute_initial_distance_matrix(sequences, distance_measure)
  
    while (is.matrix(distance_matrix)){
      results=upgma_one_step(node_description, distance_matrix, edges, edge_lengths)
      node_description=results[[1]]
      distance_matrix=results[[2]]
      edges=results[[3]]
      edge_lengths=results[[4]]
    }
  
    # Return the UPGMA tree of sequences in the form of the phylo class from ape
    tree <- transform_to_phylo(sequences, edges, edge_lengths, node_description)
    return(tree)
    
}

test_tree_building = function() {
    sequences <- list(orangutan = "TCACACCTAAGTTATATCTATATATAATCCGGGCCGGG",
                     chimp = "ACAGACTTAAAATATACATATATATAATCCGGGCCGGG",
                     human = "AAAGACTTAAAATATATATATATATAATCCGGGCCGGG",
                     gorilla = "ACACACCCAAAATATATATATATATAATCCGGGCCGGG",
                     unicorn = "ACACACCCAAAATATATACGCGTATAATCCGGGCCGAA")
    distance_measure <- 'hamming'
    distance_measure <- 'JC69'
    distance_measure <- 'K80'

    tree <- build_upgma_tree(sequences, distance_measure)
    plot_tree(tree)
}

test_tree_building()
