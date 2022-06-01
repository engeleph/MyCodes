#####################################################
#######        COMPUTATIONAL BIOLOGY         ########
#######             HOMEWORK 1               ########
#####################################################
#                                                   #
# Implement the pairwise alignment algorithms       #
# Needleman-Wunsch and Smith-Waterman.              #
#                                                   #
#####################################################
#####################################################

# In all functions the following parameters are the same:
# seqA: the first sequence to align
# seqB: the second sequence to align
# score_gap: score for a gap
# score_match: score for a character match
# score_mismatch: score for a character mismatch
# local: (logical) True if alignment is local, False otherwise

init_score_matrix = function(nrow, ncol, local, score_gap) {
  score_matrix=matrix(0,nrow,ncol)                         #creates matrix filled with 0
  if (local==F){
    for (x in 2:nrow){
      score_matrix[x,1]=score_matrix[(x-1),1]+score_gap     # the gap score is added from row to row in the first column
    }
    for (x in 2:ncol){
      score_matrix[1,x]=score_matrix[1,(x-1)]+score_gap   # the gap score is added from row to row in the first column
    }
  }
  return(score_matrix)
}


init_path_matrix = function(nrow, ncol, local) {
    # Initialize the path matrix with empty values (""), except the top row and the leftmost column if global alignment.
    # If global alignment, the top row has "left" on all positions except 1st.
    # Similarly, leftmost column has "up" on all positions except 1st.
    # nrow: (numeric) number of rows in the matrix
    # ncol: (numeric) number of columns in the matrix
    path_matrix=matrix(,nrow,ncol)
    if (local==F){
        path_matrix[2:nrow,1]="up"
        path_matrix[1,2:ncol]="left"
    }

    # Return the initialized empty path matrix
    # path_matrix: (character) nrow by ncol matrix
    return(path_matrix)
}


get_best_score_and_path = function(row, col, nucA, nucB, score_matrix, score_gap, score_match, score_mismatch, local) {
    # Compute the score and the best path for a particular position in the score matrix
    # nucA: (character) nucleotide in sequence A
    # nucB: (character) nucleotide in sequence B
    # row: (numeric) row-wise position in the matrix
    # col: (numeric) column-wise position in the matrix
    # score_matrix: (double) the score_matrix that is being filled out

    if (nucA==nucB){
      s=score_match
    } else{
      s=score_mismatch
    }
    if (local==T){
      score= max(score_matrix[row-1,col-1]+s , score_matrix[row-1,col]+score_gap,  score_matrix[row,col-1]+score_gap, 0)
      if (score==0){
        path=" "
      }
    } 
    if (local==F){
      score= max(score_matrix[row-1,col-1]+s , score_matrix[row-1,col]+score_gap,  score_matrix[row,col-1]+score_gap)
    }
    score_matrix[row,col]=score
    if (score_matrix[row,col]== score_matrix[row-1,col-1]+s){
      path="diag"
    }
    else if (score_matrix[row,col]== score_matrix[row-1,col]+score_gap){
      path="up"
    }
    else if (score_matrix[row,col]== score_matrix[row,col-1]+score_gap){
      path="left"
    }

    # Return the best score for the particular position in the score matrix
    # In the case that there are several equally good paths available, return any one of them.
    # score: (numeric) best score at this position
    # path: (character) path corresponding to the best score, one of ["diag", "up", "left"] in the global case and of ["diag", "up", "left", "-"] in the local case
    return(list("score"=score, "path"=path))
}

fill_matrices = function(seqA, seqB, score_gap, score_match, score_mismatch, local, score_matrix, path_matrix) {
    # Compute the full score and path matrices
    # score_matrix: (numeric)  initial matrix of the scores
    # path_matrix: (character) initial matrix of paths
    for (i in 2:nrow(score_matrix)){
      for (j in 2:ncol(score_matrix)){
        out=get_best_score_and_path(i,j,substr(seqA,i-1,i-1),substr(seqB,j-1,j-1),score_matrix, score_gap, score_match, score_mismatch, local)
        score=out$score
        path=out$path
        score_matrix[i,j]=score
        path_matrix[i,j]=path
      }
    }

    # Return the full score and path matrices
    # score_matrix: (numeric) filled up matrix of the scores
    # path_matrix: (character) filled up matrix of paths
    return(list("score_matrix"=score_matrix, "path_matrix"=path_matrix))
}


get_best_move = function(nucA, nucB, path, row, col) {
    # Compute the aligned characters at the given position in the score matrix and return the new position,
    # i.e. if the path is diagonal both the characters in seqA and seqB should be added,
    # if the path is up or left, there is a gap in one of the sequences.
    # nucA: (character) nucleotide in sequence A
    # nucB: (character) nucleotide in sequence B
    # path: (character) best path pre-computed for the given position
    # row: (numeric) row-wise position in the matrix
    # col: (numeric) column-wise position in the matrix

    if (path=="up"){
      char1=nucA
      char2="-"
      newrow=row-1
      newcol=col
    }
    if (path=="left"){
      char1="-"
      char2=nucB
      newrow=row
      newcol=col-1
    }
    if (path=="diag"){
      char1=nucA
      char2=nucB
      newrow=row-1
      newcol=col-1
    }

    # Return the new row and column and the aligned characters
    # newrow: (numeric) row if gap in seqA, row - 1 otherwise
    # newcol: (numeric) col if gap in seqB, col  1 otherwise
    # char1: (character) '-' if gap in seqA, appropriate character if a match
    # char2: (character) '-' if gap in seqB, appropriate character if a match
    return(list("newrow"=newrow, "newcol"=newcol, "char1"=char1, "char2"=char2))
}


get_best_alignment = function(seqA, seqB, score_matrix, path_matrix, local) {
    # Return the best alignment from the pre-computed score matrix
    # score_matrix: (numeric) filled up matrix of the scores
    # path_matrix: (character) filled up matrix of paths
    seq1=character()
    seq2=character()
    if (local==T){
      score=max(score_matrix)
      vec=which(score_matrix == score, arr.ind = TRUE)[1,]
      row=vec[1]
      col=vec[2]
      while (score_matrix[row,col]>0){
        out=get_best_move(substr(seqA,row-1,row-1), substr(seqB,col-1,col-1), path_matrix[row,col], row, col)
        char1=out$char1
        char2=out$char2
        newrow=out$newrow
        newcol=out$newcol
        seq1=c(seq1,char1)
        seq2=c(seq2, char2)
        row=newrow
        col=newcol
      }
      
    }
    
    else{
      score=score_matrix[nrow(score_matrix),ncol(score_matrix)]
      row=nrow(score_matrix)
      col=ncol(score_matrix)
      while (row>1 | col>1){
        out=get_best_move(substr(seqA,row-1,row-1), substr(seqB,col-1,col-1), path_matrix[row,col], row, col)
        char1=out$char1
        char2=out$char2
        newrow=out$newrow
        newcol=out$newcol
        seq1=c(seq1,char1)
        seq2=c(seq2, char2)
        row=newrow
        col=newcol
      }
    }
    seq1=rev(seq1)
    seq2=rev(seq2)
    alignment=rbind(seq1,seq2)

    # Return the best score and alignment (or one thereof if there are multiple with equal score)
    # score: (numeric) score of the best alignment
    # alignment: (character) the actual alignment in the form of a vector of two strings
    return(list("score"=score, "alignment"=alignment))
}


align = function(seqA, seqB, score_gap, score_match, score_mismatch, local) {
    # Align the two sequences given the scoring scheme
    # For testing purposes, use seqA for the rows and seqB for the columns of the matrices
  
    # Initialize score and path matrices
    score_matrix=init_score_matrix(nchar(seqA)+1, nchar(seqB)+1, local, score_gap)
    path_matrix=init_path_matrix(nchar(seqA)+1, nchar(seqB)+1, local)
  
    # Fill in the matrices with scores and paths using dynamic programming
    output=fill_matrices(seqA, seqB, score_gap, score_match, score_mismatch, local, score_matrix, path_matrix)
    score_matrix=output$score_matrix
    path_matrix=output$path_matrix
  
    # Get the best score and alignment (or one thereof if there are multiple with equal score)
    output=get_best_alignment(seqA, seqB, score_matrix, path_matrix, local)
    score=output$score
    alignment=output$alignment
    # Return the best score and alignment (or one thereof if there are multiple with equal score)
    # Returns the same value types as get_best_alignment
    result=list("score"=score, "alignment"=alignment)
    return(result)
}

test_align = function() {
    seqA = "TCACACTAC"
    seqB = "AGCACAC"
    score_gap = -2
    score_match = +3
    score_mismatch = -1
    local = F
    result = align(seqA, seqB, score_gap, score_match, score_mismatch, local)
    print(result$alignment)
    print(result$score)
}

test_align()
