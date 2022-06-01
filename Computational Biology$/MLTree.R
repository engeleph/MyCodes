#########################################################
#######        COMPUTATIONAL BIOLOGY         ############
#######             HOMEWORK 4               ############
#########################################################
#                                                       #
# Compute the log likelihood of an alignment on a given #
# tree using the Felsenstein's tree-pruning algorithm   #
#########################################################
#########################################################

#########################################################
######    Code below this should not be changed   #######
######    or deleted.                             #######
#########################################################

library(ape)
library(Matrix)

nucleotides = c("T", "C", "A", "G")

get_number_from_nucleotide = function(nuc) {
    # Transform a nucleotide letter into the appropriate number
    return(which(nucleotides == nuc))
}

transform_to_numbers = function(sequence) {
    # Transform a sequence of nucleotide characters into the appropriate number sequence
    num_sequence = simplify2array(lapply(strsplit(sequence, split=c())[[1]], get_number_from_nucleotide))
    return(num_sequence)
}

get_sequence_at_tip_node = function(node, tree, sequences) {
    # Get the sequence at the appropriate tip, represented by a numeric vector
    if (node > length(tree$tip.label)) {
        return(NULL)
    }
    tip_name = tree$tip.label[node]
    sequence = sequences[[tip_name]]
    return(transform_to_numbers(sequence))
}

get_node_children = function(node, tree) {
    # Return the two child nodes of the given node.
    # In a bifurcating tree, a node will have either two child nodes, if it is an internal node,
    # or zero child nodes if it is a tip node.
    # For an internal node, return the two child node numbers and the branch lengths leading to them,
    # for a tip node return NULL.
    #   node: the node for which child nodes should be returned
    #   tree: the tree encoded in an object of class phylo.

    if (node <= length(tree$tip.label)) {
        return(NULL)
    }

    children = which(tree$edge[,1] == node)
    child1 = tree$edge[children[1],2]
    branch_length1 = tree$edge.length[children[1]]
    child2 = tree$edge[children[2],2]
    branch_length2 = tree$edge.length[children[2]]

    # Return the child nodes and the branch lengths leading to them.
    #   child1: first child of the given node
    #   branch_length1: the length of the branch leading to child1
    #   child2: second child of the given node
    #   branch_length2: the length of the branch leading to child2
    return(list(child1 = child1, branch_length1 = branch_length1,
                child2 = child2, branch_length2 = branch_length2))
}

#########################################################
######    Code above this should not be changed   #######
######    or deleted.                             #######
#########################################################

create_TN93_Q_matrix = function(pi, alpha1, alpha2, beta) {
    # Create the TN93 transition rate matrix Q as specified in the assignment.
    # The nucleotide order should be (T, C, A, G).
    #   pi: the stationary frequencies of nucleotides
    #   alpha1: rate coefficient for the C <-> T transition
    #   alpha2: rate coefficient for the A <-> G transition
    #   beta: rate coefficient for transversions

    Q=matrix(0,nrow=4,ncol=4)
    Q[,1]=pi[1]
    Q[,2]=pi[2]
    Q[,3]=pi[3]
    Q[,4]=pi[4]
    Q[1,2]=Q[1,2]*alpha1
    Q[2,1]=Q[2,1]*alpha1
    Q[3,4]=Q[3,4]*alpha2
    Q[4,3]=Q[4,3]*alpha2
    Q[1:2,3:4]=Q[1:2,3:4]*beta
    Q[3:4,1:2]=Q[3:4,1:2]*beta
    Q[1,1]=0-(Q[1,2]+Q[1,3]+Q[1,4])
    Q[2,2]=0-(Q[2,1]+Q[2,3]+Q[2,4])
    Q[3,3]=0-(Q[3,1]+Q[3,2]+Q[3,4])
    Q[4,4]=0-(Q[4,1]+Q[4,2]+Q[4,3])

    # Return the transition rate matrix
    # Q: 4 by 4 matrix of rates as numeric values.
    return(Q)
}


calculate_likelihood_from_subtree_likelihoods = function(N, Q,
                                                         subtree1_likelihood, subtree1_branch_length,
                                                         subtree2_likelihood, subtree2_branch_length) {
    # Calculate the likelihood of each nucleotide at each site on an internal node, from the two subtree likelihoods.
    #   N: number of sites in the alignment
    #   Q: the substitution rate matrix
    #   subtree1_likelihood: an N by 4 matrix containing the per site per nucleotide likelihoods
    #                        at the first child node
    #   subtree1_branch_length: the length of the branch leading to the first child node
    #   subtree2_likelihood: an N by 4 matrix containing the per site per nucleotide likelihoods
    #                        at the second child node
    #   subtree2_branch_length: the length of the branch leading to the second child node

    P1=expm(subtree1_branch_length*Q)
    P2=expm(subtree2_branch_length*Q)
    likelihood_per_site=matrix(0,N,4)
    for (i in 1:N){
      for (j in 1:4){
        likelihood1=0
        likelihood2=0
        for (k in 1:4){
          likelihood1=likelihood1+P1[j,k]*subtree1_likelihood[i,k]  
          likelihood2=likelihood2+P2[j,k]*subtree2_likelihood[i,k]
        }
        likelihood_per_site[i,j]=likelihood1*likelihood2
      }
    }

    # Return the per site nucleotide likelihoods on the internal node.
    #   likelihood_per_site: an N by 4 matrix representing the nucleotide likelihoods per site as numeric values.
    return(likelihood_per_site)
}

get_likelihood_from_sequence = function(N, sequence) {
    # Get the matrix of likelihoods of each nucleotide at each site at a tip given the sequence at this tip.
    # The matrix should be of size N x 4, where rows represent the sites in the alignment and
    # the columns represent a nucleotide in the order (T, C, A, G).
    # For each site, the likelihood should be 1 in the column representing
    # the nucleotide at this site and 0 in all other columns.
    #   N: number of sites in the alignment
    #   sequence: the sequence at the tip, encoded as a vector of nucleotide numeric values.

    likelihood_per_site=matrix(0,N,4)
    for (i in 1:N){
      if (sequence[i]==1){
        likelihood_per_site[i,1]=1
      }
      else if (sequence[i]==2){
        likelihood_per_site[i,2]=1
      }
      else if (sequence[i]==3){
        likelihood_per_site[i,3]=1
      }
      else if (sequence[i]==4){
        likelihood_per_site[i,4]=1
      }
    }

    # Return the per site nucleotide likelihoods at the tip.
    #   likelihood_per_site: an N by 4 matrix representing the nucleotide likelihoods per site as numeric values.
    return(likelihood_per_site)
}

get_subtree_likelihood = function(node, tree, sequences, Q) {
    # Return the matrix of likelihoods per site per nucleotide for the subtree below the given node.
    # The matrix should be of size N x 4, where rows represent the sites in the alignment and
    # the columns represent a nucleotide in the order (T, C, A, G).
    # If the node is a tip, the likelihood per site should be a vector with a 1 in the column representing
    # the nucleotide at this site and 0 in all other columns.
    # Use the get_sequence_at_tip_node function to get the numeric representation of the sequence at the
    # given tip.
    # If the node is an internal node, the likelihood for each nucleotide should be computed
    # from the likelihoods of each of the child subtrees.
    #   node: the node for which likelihoods should be computed
    #   tree: the tree encoded in an object of class phylo
    #   sequences: a list of aligned sequences (as character strings) at the tips of the tree
    #   Q: the substitution rate matrix

    N = nchar(sequences[[1]])
    

    if (node <= length(tree$tip.label)) {
        # node is a tip: compute the likelihood for each site and each nucleotide on this node
      sequence=get_sequence_at_tip_node(node, tree, sequences)
      likelihood_per_site=get_likelihood_from_sequence(N, sequence)
      
    } else {
        # node is an internal node: compute the likelihood for each site and each nucleotide on this node
      result=get_node_children(node, tree)
      child1=result[[1]]
      subtree1_likelihood=get_subtree_likelihood(child1, tree, sequences, Q)
      branch_length1=result[[2]]
      child2=result[[3]]
      subtree2_likelihood=get_subtree_likelihood(child2, tree, sequences, Q)
      branch_length2=result[[4]]
      likelihood_per_site=calculate_likelihood_from_subtree_likelihoods(N, Q, subtree1_likelihood, branch_length1, subtree2_likelihood, branch_length2)
    }

    # Return the per site nucleotide likelihoods at this node.
    #   likelihood_per_site: an. N by 4 matrix representing the nucleotide likelihoods per site as numeric values.
    return(likelihood_per_site)
}

Felstensteins_pruning_loglikelihood = function(pi, alpha1, alpha2, beta, newick_tree, sequences) {
    # Compute the log likelihood of a sequence alignment on the given tree under the TN93 model
    # using Felsenstein's pruning algorithm.
    #   pi: the stationary frequencies of nucleotides
    #   alpha1: rate coefficient for the C <-> T transition
    #   alpha2: rate coefficient for the A <-> G transition
    #   beta: rate coefficient for transversions
    #   newick_tree: the tree in newick text format
    #   sequences: a list of aligned sequences (as character strings) at the tips of the tree
  
    # Transfrom the tree from text format to an object of the phylo class which represents the tree in R
    tree = read.tree(text = newick_tree)
    # Reorder the tree for easier traversing
    tree = reorder(tree, order = "cladewise");
    
    # Number of sites in the alignment
    N = nchar(sequences[[1]])

    # Initialize the Q matrix
    Q=create_TN93_Q_matrix(pi, alpha1, alpha2, beta)
    
    # Compute the likelihoods per site of the tree starting from the root
    root = length(tree$tip.label) + 1
    likelihood_per_site = get_subtree_likelihood(root, tree, sequences, Q)
    
    # Sum up the log likelihoods of each site
    loglikelihood=0
    for (i in 1:N){
      likelihood=0
      for (j in 1:4){
        likelihood=likelihood+pi[j]*likelihood_per_site[i,j]
      }
      loglikelihood=loglikelihood+log(likelihood)
    }
    
    

    # Return the final log likelihood of the alignment on the given tree, which is a single number 
    # computed based on the per site, per nucleotide likelihoods at the root of the tree.
    # loglikelihood: the full log likelihood of the alignment given the tree
    return(loglikelihood)
}

test_tree_loglikelihood = function() {
    beta = 0.00135
    alpha1 = 0.5970915
    alpha2 = 0.2940435
    pi = c(0.22, 0.26, 0.33, 0.19)
    newick_tree = "(unicorn:15,(orangutan:13,(gorilla:10.25,(human:5.5,chimp:5.5):4.75):2.75):2);"
    sequences = list(orangutan = "ACCCCTCCCCTCATGTGTAC",
                     chimp = "ACCCCTCCCCTCATGTGTAC",
                     human = "ACCCCTCCCCTCATGTGTAC",
                     gorilla = "ACCCCTCCCCTCATGTGTAC",
                     unicorn = "TGCCCTCCCCTCATGTGTAC")
    expected_result = -89.4346738736

    result = Felstensteins_pruning_loglikelihood(pi, alpha1, alpha2, beta, newick_tree, sequences)
    print(sprintf("Log likelihood value: %.*f" , 10, result))

    comparison = all.equal(result, expected_result)
    if (isTRUE(comparison)) {
        print("The result matched the expected value.")
    } else {
        print(comparison)
    }
}

test_tree_loglikelihood()
