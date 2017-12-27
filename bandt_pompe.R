suppressMessages(library(igraph))
require(e1071, quietly=TRUE)

# Implementation of the bandt-pope methodology

# Bandt-Pompe transformation
#
# data: the dataset (vector)
# D:    the embedding dimension (size of sliding window)
# tau:  the embedding delay ('step' value)
#
# return The list of symbols from the transformed dataset
bandt_pompe = function(data, D=4, tau=1, by=1)
{
    # the list of symbols to be returned
    symbols = c()
    
    # dicovering the sequences of order n
    for (s in seq(1, length(data)-(D-1)*tau, by=by))
    {
        # the indices for the subsequence
        ind = seq( s, s+(D-1)*tau, by=tau)
        
        # get the sub-timeseries (sliding window)
        sub = data[ind]

        # the current permutation pattern
        pattern = paste(order(sub), collapse='')

        # adding the current pattern to the list of symbols
        symbols = c(symbols, pattern)
    }

    return(symbols)
}

# Bandt-Pompe distribution
#
# Parameters:
# data: the dataset (vector)
# D:    the embedding dimension (size of sliding window)
# tau:  the embedding delay ('step' value)
# numred: numerosity reduction, similar to BOSS algorithm, do not count
#         repetitions of the same symbol
bandt_pompe_distribution = function(data, D=4, tau=1, numred=FALSE, by=1)
{
    # the list of symbols from the BP transformation
    symbols = bandt_pompe(data, D=D, tau=tau, by=by)

    # the distribution of permutations (pi)
    dpi = rep(0, factorial(D))
    
    # to get the index of the permutation pi
    perms = sort(apply(permutations(D), 1, paste, collapse=''))
    
    if (numred == FALSE)
    {
        # dicovering the sequences of order n
        for (i in 1:length(perms))
        {
            # counting the pattern
            dpi[i] = sum(symbols == perms[i])
        }
    }
    else
    {
        # only counts different subsequence of symbols
        oldsymbol = symbols[1]
        ind = which(perms == oldsymbol)
        dpi[ind] = dpi[ind] + 1
        
        for (i in 2:length(symbols))
        {
            if (oldsymbol != symbols[i])
            {
                ind = which(perms == symbols[i])
                dpi[ind] = dpi[ind] + 1
            }
            oldsymbol = symbols[i]
        }
    }
    
    # the returning format
    output = data.frame(patterns = perms, 
                        frequencies = dpi,
                        probabilities = dpi/sum(dpi))
    return(output)
}



# The Band-Pompe Transition (adjacency matrix)
#
# Parameters:
# data: the dataset (vector)
# D:    the embedding dimension (size of sliding window)
# tau:  the embedding delay ('step' value)
# loop: enbale (TRUE) the creation of loopings in the graph
# normalized: returns as percentage
# vect: returns as a vector instead of matrix
bandt_pompe_transition = function(data, D=4, tau=1, 
                                  normalized=TRUE, quiroz=FALSE,
                                  loop=TRUE, vect=FALSE, by=1)
{
    # the number of permutations
    dfact = factorial(D)

    # to get the index of the permutation pi
    perms = sort(apply(e1071::permutations(D), 1, paste, collapse=''))
    
    # the transitions matrix
    M = matrix(0, ncol=dfact, nrow=dfact)
    rownames(M) = perms
    colnames(M) = perms
    
    # the list of symbols from the BP transformation
    symbols = bandt_pompe(data, D=D, tau=tau, by=by)

    # dicovering the sequences of order n
    for (i in 2:length(symbols))
    {
        # the previous permutation pattern
        from = symbols[i-1]

        # the current permutation pattern
        to = symbols[i]
        
        # checking if the creation of loopings are allowed
        if (from != to | loop == TRUE)
        {
            # incrementing the counting for this transition
            M[from, to] = M[from, to] + 1
        }
    }

    if (normalized == TRUE)
    {
        # this is the normalization suggested by Quiroz (2015)
        # similar to a markov chain probabilities
        if (quiroz == TRUE)
        {
            # normalization by row: right stochastic matrix
            M = t(apply(M, 1, function(x)
                            {
                                if(sum(x) > 0)
                                {
                                    x/sum(x)
                                } else {
                                    rep(0, length(x))
                                }
                            }
                            ))
        }
        else
        {
            # doing a general normalization
            M = M/sum(M)
        }
    }

    if (vect == TRUE)
    {
        M = c(M)
    }
    
    return(M)
}

# get the BP transition graph
# - D: embedded dimension
# - tau: embedded delay
# - empty: TRUE to remove the vertices with 0 degree
# - quiroz: normalization proposed by quiroz, similar to a markov-chain
bandt_pompe_transition_graph = function(x, D=4, tau=1, 
                                        empty=TRUE, quiroz=FALSE)
{
    # adjacency matrix from bandt pompe transition
    A = bandt_pompe_transition(x, D=D, tau=tau, quiroz=quiroz)
    
    # the graph
    gA = graph_from_adjacency_matrix(A, mode="directed", weighted=TRUE)
    
    # removing vertices without transition
    if (empty == TRUE)
    {
        gA = delete_vertices(gA, which(igraph::degree(gA) == 0))
    }

    return(gA)
}


