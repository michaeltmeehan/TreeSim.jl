"""
    forward_algorithm(times::Vector{Float64},
                      leaves::Vector{Int64},
                      Nₑ::Float64,
                      bound::Float64) -> Matrix{Float64}

Calculate the forward probabilities of transitioning between different numbers of lineages 
over a series of time intervals using a forward algorithm.

# Arguments
- `times::Vector{Float64}`: A vector of time points.
- `leaves::Vector{Int64}`: A vector indicating the number of leaves at each time point.
- `Nₑ::Float64`: The effective population size.
- `bound::Float64`: The lower bound for the time interval.

# Returns
- `Matrix{Float64}`: A matrix of forward probabilities, where each element represents the probability of transitioning from one number of lineages to another over a given time interval.
"""
function forward_algorithm(times::Vector{Float64},
                           leaves::Vector{Int64},
                           Nₑ::Float64,
                           bound::Float64)::Matrix{Float64}

    n_times = length(times)
    sum_leaves = leaves[end]

    # Initiate forward algorithm
    forward_probs = zeros(sum(leaves), n_times+1)
    forward_probs[leaves[end], end] = 1.0

    # Forward recursion through sampling times
    for k in n_times:-1:2
        dt = times[k] - times[k - 1]
        for n_end in 1:sum_leaves
            for n_start in 1:sum_leaves
                transition_prob = homochronous_probability(n_start, n_end, dt, Nₑ)
                forward_probs[n_end + leaves[k-1], k] += transition_prob * forward_probs[n_start, k+1]
            end
        end
        sum_leaves += leaves[k-1]
    end

    # Bound probabilities
    dt = times[1] - bound

    for n_end in 1:sum_leaves
        for n_start in 1:sum_leaves
            transition_prob = homochronous_probability(n_start, n_end, dt, Nₑ)
            forward_probs[n_end, 1] += transition_prob * forward_probs[n_start, 2]
        end
    end
    return forward_probs
end