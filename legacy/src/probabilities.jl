"""
    homochronous_probability(n_start::Int, 
                             n_end::Int, 
                             dt::Float64, 
                             Nₑ::Float64) -> Float64

Calculate the probability of transitioning from `n_start` lineages to `n_end` lineages over a 
time interval `dt` under a constant population size `Nₑ`.

# Arguments
- `n_start::Int`: The number of lineages at the start of the interval.
- `n_end::Int`: The number of lineages at the end of the interval.
- `dt::Float64`: The time interval over which the transition occurs.
- `Nₑ::Float64`: The effective population size.

# Returns
- `Float64`: The probability of transitioning from `n_start` to `n_end` lineages over the time interval `dt`. Returns 0.0 if the inputs are invalid or if `n_start` is less than `n_end`.
"""
function homochronous_probability(n_start::Int, 
                                  n_end::Int, 
                                  dt::Float64, 
                                  Nₑ::Float64)::Float64

    # Check that inputs are valid
    (n_start <= 0 || n_end <= 0 || n_start < n_end || dt < 0 || Nₑ <= 0) && return 0.0
    (n_start == 1 && n_end == 1) && return 1.0

    if n_end == 1
        # Initialise total probability
        prob_total = 0.0
        for k in 2:n_start
            dk = convert(Float64, k)
            # Initialise probability increment
            prob_increment = 1.0

            for l in 2:n_start
                dl = convert(Float64, l)
                # Calculate coefficients
                if l != k 
                    prob_increment *= (dl * (dl - 1.0)) /
                                                (dl * (dl - 1.0) - dk * (dk - 1.0))
                end
            end

            # Calculate exponential
            prob_increment *= 1.0 - exp(- ((dk * (dk - 1.0)) / (2.0 * Nₑ)) * dt)

            # Update total probability
            prob_total += prob_increment
        end
    else
        # Initialise total probability
        prob_total = 0.0
        for k in n_end:n_start
            dk = convert(Float64, k)

            # Initialise probability increment
            prob_increment = (dk * (dk - 1.0)) / (convert(Float64, n_end) * (convert(Float64, n_end) - 1.0))
            
            for l in n_end:n_start
                dl = convert(Float64, l)
                # Calculate coefficients
                if l != k
                    prob_increment *= (dl * (dl - 1.0)) /
                                                (dl * (dl - 1.0) - dk * (dk - 1.0))
                end
            end

            # Calculate exponential
            prob_increment *= exp(- ((dk * (dk - 1.0)) / (2.0 * Nₑ)) * dt)

            # Update total probability
            prob_total += prob_increment
        end
    end
    return prob_total
end


"""
    bounded_times_likelihood(leaf_times::Vector{Float64},
                             leaves::Vector{Int64},
                             coalescence_times::Vector{Float64},
                             Nₑ::Float64,
                             bound::Float64) -> Float64

Calculate the likelihood of the coalescence times under the bounded coalescent model.

# Arguments
- `leaf_times::Vector{Float64}`: A vector of leaf (sampling) times.
- `leaves::Vector{Int64}`: A vector indicating the number of leaves at each time point.
- `coalescence_times::Vector{Float64}`: A vector of coalescence times.
- `Nₑ::Float64`: The effective population size.
- `bound::Float64`: The lower bound for the time interval.

# Returns
- `Float64`: The likelihood of the coalescence times given the bounded coalescent model.
"""
function bounded_times_likelihood(leaf_times::Vector{Float64},
                                  leaves::Vector{Int64},
                                  coalescence_times::Vector{Float64},
                                  Nₑ::Float64,
                                  bound::Float64)::Float64
    # Counters
    k = length(leaf_times)
    c = length(coalescence_times)
    # Current time
    current_time = leaf_times[k]
    # Current lineages
    current_lineages = leaves[k]
    k -= 1
    # Next event times
    next_leaf_time = k >= 1 ? leaf_times[k] : bound
    next_coalescence_time = coalescence_times[c]
    likelihood = 1.0
    while current_time > coalescence_times[1]
        if next_leaf_time > next_coalescence_time
            dt = current_time - next_leaf_time
            coef = (convert(Float64, current_lineages) * (convert(Float64,current_lineages) - 1.0)) / (2.0 * Nₑ)
            likelihood *= exp(-coef * dt)
            current_lineages += leaves[k]
            current_time = next_leaf_time
            k -=1
            next_leaf_time = k >= 1 ? leaf_times[k] : bound
        else
            dt = current_time - next_coalescence_time
            coef = (convert(Float64, current_lineages) * (convert(Float64, current_lineages) - 1.0)) / (2.0 * Nₑ)
            likelihood *= coef * exp(-coef * dt)
            current_lineages -= 1
            current_time = next_coalescence_time
            c -=1
            if c >= 1 
                next_coalescence_time = coalescence_times[c]
            end
        end
    end
    total_leaves = sum(leaves)
    # forward_probs = zero(total_leaves, length(leaf_times) + 1)
    forward_probs = forward_algorithm(leaf_times, leaves, Nₑ, bound)
    likelihood /= forward_probs[1, 1]
    return likelihood
end