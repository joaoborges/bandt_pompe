
# Metrics to be computed

# Shannon entropy 'H'
#
# probs: vector of probabilities from the Bandt-Pompe distribution
shannon_entropy = function(probs, normalized=FALSE)
{
    # consider only the non-zero values
    p = which(probs > 1e-30)
    
    # considering log base n
    entropy = -sum(probs[p]*log(probs[p]))
    
    # if normalized, the log base is irrelevant
    if(normalized)
    {
        entropy = entropy/log(length(probs))
    }

    return(entropy)
}

# Statistical Complexity
#
# probs: vector of probabilities from the Bandt-Pompe distribution
# entropy: the Shannon entropy for probs
complexity = function(probs, entropy = NULL, normalized=TRUE)
{
    if (is.null(entropy))
    {
        entropy = shannon_entropy(probs, normalized=normalized)
    }

    # the length of the probabilities, 
    N = length(probs)

    # the reference distribution (uniform)
    P_u = rep(1/N, N)

    # the Jensen-shannon divergence
    JS = shannon_entropy( (probs + P_u) / 2  ) - 
        shannon_entropy(probs)/2 - shannon_entropy(P_u)/2

    # the statistical complexity
    aux = (   ((N+1)/N) * log(N + 1) - 2*log(2*N) + log(N)    )
    Q_0 = -2*(1/aux)
    Q = Q_0 * JS
    C = Q*entropy

    # "piggypacking" the JS variable as an 'attr'
    attr(C, "JS") = JS

    return(C)
}

# Fisher Information
# probs: vector of probabilities from the Bandt-Pompe distribution
fisher_information = function(probs)
{
    # the number of probabilities
    N = length(probs)

    # the normalization constant
    if (probs[1] == 1 | probs[N] == 1)
    {
        F0 = 1
    } else {
        F0 = 1/2
    }

    #print((sqrt(probs[2:N]) - sqrt(probs[1:(N-1)]))^2)

    # the fisher information
    f = F0 * sum((sqrt(probs[2:N]) - sqrt(probs[1:(N-1)]))^2)

    return(f)
}

# computes the probability of self-transitions (pst) from a given
# transition graph
pst = function(g)
{
    # self-transitions
    st = E(g)$weight[which_loop(g)]

    # pst
    return(sum(st))
}

