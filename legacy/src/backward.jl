"""
    backward_sampler!(lineages::Vector{Int64},
                      forward_probs::Matrix{Float64},
                      times::Vector{Float64},
                      leaves::Vector{Int64},
                      Nₑ::Float64,
                      bound::Float64,
                      bound_size::Int64) -> Float64

Perform backward sampling to generate lineages and calculate the likelihood of the sample.

# Arguments
- `lineages::Vector{Int64}`: A vector to store the number of lineages at each time step.
- `forward_probs::Matrix{Float64}`: A matrix of forward probabilities.
- `times::Vector{Float64}`: A vector of time points.
- `leaves::Vector{Int64}`: A vector indicating the number of leaves at each time point.
- `Nₑ::Float64`: The effective population size.
- `bound::Float64`: The lower bound for the time interval.
- `bound_size::Int64`: The initial number of lineages at the bound.

# Returns
- `Float64`: The likelihood of the sampled lineages.
"""
function backward_sampler!(lineages::Vector{Int64},
                           forward_probs::Matrix{Float64},
                           times::Vector{Float64},
                           leaves::Vector{Int64},
                           Nₑ::Float64,
                           bound::Float64,
                           bound_size::Int64)::Float64
    # Initialise with bound condition
    lineages[1] = bound_size

    # Likelihood of sample
    likelihood = 1.0

    # Total the leaves
    total_leaves = sum(leaves)

    # Initiate backwards sampling
    dt = times[1] - bound

    # Generate random number for sampling
    u = rand()

    # Track cumulative probability for sampling
    sum_prob = 0.0

    # Evaluate smoothed probabilities, sample, and update likelihood
    for n_start in 1:total_leaves
        transition_prob = homochronous_probability(n_start, lineages[1], dt, Nₑ)
        smoothed_prob = (transition_prob * forward_probs[n_start, 2]) /
                            forward_probs[lineages[1], 1]

        sum_prob += smoothed_prob

        if (u < sum_prob)
            lineages[2] = n_start
            likelihood *= smoothed_prob
            break
        end
    end

    # Backwards sampling recursion
    for k in 2:length(times)
        dt = times[k] - times[k - 1]
        # Generate random number for sampling
        u = rand()

        # Track cumulative probability for sampling
        sum_prob = 0.0

        # Evaluate smoothed probabilities, sample, and update likelihood
        for n_start in 1:total_leaves
            # Number of lineages before leaves are added
            n_end = lineages[k] - leaves[k - 1]
            transition_prob = homochronous_probability(n_start, n_end, dt, Nₑ)
            smoothed_prob = (transition_prob * forward_probs[n_start, k + 1]) /
                                forward_probs[lineages[k], k]

            sum_prob += smoothed_prob
            if (u < sum_prob)
                lineages[k + 1] = n_start
                likelihood *= smoothed_prob
                break
            end
        end
    end
    # Return likelihood
    return likelihood
end



backward_sampler!(lineages::Vector{Int64},
                  forward_probs::Matrix{Float64},
                  times::Vector{Float64},
                  leaves::Vector{Int64},
                  Nₑ::Float64,
                  bound::Float64)::Float64 = backward_sampler!(lineages::Vector{Int64},
                                                               forward_probs::Matrix{Float64},
                                                               times::Vector{Float64},
                                                               leaves::Vector{Int64},
                                                               Nₑ::Float64,
                                                               bound::Float64,
                                                               1)::Float64

