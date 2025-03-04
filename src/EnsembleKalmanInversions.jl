module EnsembleKalmanInversions

using Distributions
using ProgressBars
using Random
using Suppressor
using LinearAlgebra
using EnsembleKalmanProcesses.EnsembleKalmanProcessModule
using EnsembleKalmanProcesses.ParameterDistributionStorage

using ..InverseProblems: n_ensemble, observation_map, forward_map

function lognormal_with_mean_std(mean, std)
    k = std^2 / mean^2 + 1
    μ = log(mean / sqrt(k))
    σ = sqrt(log(k))
    return LogNormal(μ, σ)
end

struct ConstrainedNormal{FT}
    # θ is the original constrained paramter, θ̃ is the unconstrained parameter ~ N(μ, σ)
    # θ = lower_bound + (upper_bound - lower_bound)/（1 + exp(θ̃)）
    μ::FT
    σ::FT
    lower_bound::FT
    upper_bound::FT
end

# Model priors are sometimes constrained; EKI deals with unconstrained, Normal priors.
convert_prior(prior::LogNormal) = Normal(prior.μ, prior.σ)
convert_prior(prior::Normal) = prior
convert_prior(prior::ConstrainedNormal) = Normal(prior.μ, prior.σ)


# Convert parameters to unconstrained for EKI
forward_parameter_transform(::LogNormal, parameter) = log(parameter)
forward_parameter_transform(::Normal, parameter) = parameter
forward_parameter_transform(cn::ConstrainedNormal, parameter) = log((cn.upper_bound - parameter)/(cn.upper_bound - cn.lower_bound))

# Convert parameters from unconstrained (EKI) to constrained
inverse_parameter_transform(::LogNormal, parameter) = exp(parameter)
inverse_parameter_transform(::Normal, parameter) = parameter
inverse_parameter_transform(cn::ConstrainedNormal, parameter) = cn.lower_bound+(cn.upper_bound - cn.lower_bound)/(1 + exp(parameter))

# Convert covariance from unconstrained (EKI) to constrained
inverse_covariance_transform(::Tuple{Vararg{LogNormal}}, parameters, covariance) = Diagonal(exp.(parameters)) * covariance * Diagonal(exp.(parameters))
inverse_covariance_transform(::Tuple{Vararg{Normal}}, parameters, covariance) = covariance
function inverse_covariance_transform(cn::Tuple{Vararg{ConstrainedNormal}}, parameters, covariance)
    upper_bound = [cn[i].upper_bound for i = 1:length(cn)]
    lower_bound = [cn[i].lower_bound for i = 1:length(cn)]
    dT = Diagonal(-(upper_bound - lower_bound) .* exp.(parameters)./(1.0 .+ exp.(parameters)).^2) 
    return dT * covariance * dT'
end

mutable struct EnsembleKalmanInversion{I, P, E, M, O, F, S, D}
    inverse_problem :: I
    parameter_distribution :: P
    ensemble_kalman_process :: E
    mapped_observations :: M
    noise_covariance :: O
    inverting_forward_map :: F
    iteration :: Int
    iteration_summaries :: S
    dropped_ensemble_members :: D
end

Base.show(io::IO, eki::EnsembleKalmanInversion) = 
    print(io, "EnsembleKalmanInversion", '\n',
              "├── inverse_problem: ", typeof(eki.inverse_problem).name.wrapper, '\n',
              "├── parameter_distribution: ", typeof(eki.parameter_distribution).name.wrapper, '\n',
              "├── ensemble_kalman_process: ", typeof(eki.ensemble_kalman_process), '\n',
              "├── mapped_observations: ", typeof(eki.mapped_observations), '\n',
              "├── noise_covariance: ", typeof(eki.noise_covariance), '\n',
              "├── inverting_forward_map: ", typeof(eki.inverting_forward_map).name.wrapper, '\n',
              "├── iteration: $(eki.iteration)", '\n',
              "├── iteration_summaries: $(eki.iteration_summaries)", '\n',
              "└── dropped_ensemble_members: $(eki.dropped_ensemble_members)", '\n'
              )

construct_noise_covariance(noise_covariance::AbstractMatrix, y) = noise_covariance

function construct_noise_covariance(noise_covariance::Number, y)
    # Independent noise for synthetic observations
    n_obs = length(y)
    return noise_covariance * Matrix(I, n_obs, n_obs)
end

"""
    EnsembleKalmanInversion(inverse_problem; noise_covariance = 1e-2)

Return an object that interfaces with EnsembleKalmanProcesses.jl and uses Ensemble
Kalman Inversion to iteratively "solve" the inverse problem

```math
y = G(θ) + η,
```

for the parameters ``θ``, where ``y`` is a "normalized" vector of observations,
``G(θ)`` is a forward map that predicts the observations, and ``η ∼ N(0, Γ_y)`` is zero-mean
random noise with covariance matrix ``Γ_y`` representing uncertainty in the observations.

By "solve", we mean that the iteration finds  ``θ`` that minimizes the distance
between ``y`` and ``G(θ)``.

(For more details on the Ensemble Kalman Inversion algorithm refer to the
[EnsembleKalmanProcesses.jl Documentation](https://clima.github.io/EnsembleKalmanProcesses.jl/stable/ensemble_kalman_inversion/).)

The "forward map output" `G` can have many interpretations. The specific statistics that `G`
computes have to be selected for each use case to provide a concise summary of the complex
model solution that contains the values that we would most like to match to the corresponding
truth values `y`. For example, in the ocean surface boundary layer context, this summary 
could be, for example, a vector of concatenated `u`, `v`, `b`, `e` profiles at all or some time steps
of the CATKE solution.

Arguments
=========
- `inverse_problem :: InverseProblem`: an inverse problem representing the comparison between
                                       synthetic observations generated by Oceananigans and model
                                       predictions also generated by Oceananigans.

- `noise_covariance` (`AbstractMatrix` or `Number`): normalized covariance representing observational
                                                     uncertainty. If `noise_covariance isa Number` then
                                                     it's converted to an identity matrix scaled by
                                                     `noise_covariance`.
"""
function EnsembleKalmanInversion(inverse_problem; noise_covariance = 1e-2)
    free_parameters = inverse_problem.free_parameters
    original_priors = free_parameters.priors

    transformed_priors = [Parameterized(convert_prior(prior)) for prior in original_priors]
    no_constraints = [[no_constraint()] for _ in transformed_priors]
    parameter_distribution = ParameterDistribution(transformed_priors, no_constraints, collect(string.(free_parameters.names)))

    # Seed for pseudo-random number generator for reproducibility
    initial_ensemble = construct_initial_ensemble(parameter_distribution, n_ensemble(inverse_problem); rng_seed = Random.seed!(41))

    # Build EKP-friendly observations "y" and the covariance matrix of observational uncertainty "Γy"
    y = dropdims(observation_map(inverse_problem), dims=2) # length(forward_map_output) column vector
    Γy = construct_noise_covariance(noise_covariance, y)

    # The closure G(θ) maps (N_params, ensemble_size) array to (length(forward_map_output), ensemble_size)
    function G(θ) 
        batch_size = size(θ, 2)
        inverted_parameters = [inverse_parameter_transform.(values(original_priors), θ[:, i]) for i in 1:batch_size]
        return forward_map(inverse_problem, inverted_parameters)
    end

    ensemble_kalman_process = EnsembleKalmanProcess(initial_ensemble, y, Γy, Inversion())

    return EnsembleKalmanInversion(inverse_problem, parameter_distribution, ensemble_kalman_process, y, Γy, G, 0, [], Set())
end

"""
    UnscentedKalmanInversion(inverse_problem, prior_mean, prior_cov;
                             noise_covariance = 1e-2, α_reg = 1, update_freq = 0)

Return an object that interfaces with EnsembleKalmanProcesses.jl and uses Unscented
Kalman Inversion to iteratively "solve" the inverse problem

```math
y = G(θ) + η,
```

for the parameters ``θ``, where ``y`` is a "normalized" vector of observations,
``G(θ)`` is a forward map that predicts the observations, and ``η ∼ N(0, Γ_y)`` is zero-mean
random noise with covariance matrix ``Γ_y`` representing uncertainty in the observations.

By "solve", we mean that the iteration finds  ``θ`` that minimizes the distance
between ``y`` and ``G(θ)``.

(For more details on the Unscented Kalman Inversion algorithm refer to the
[EnsembleKalmanProcesses.jl Documentation](https://clima.github.io/EnsembleKalmanProcesses.jl/stable/unscented_kalman_inversion/).)
    
Arguments
=========

- `inverse_problem :: InverseProblem`: an inverse problem representing the comparison between 
                                       synthetic observations generated by Oceananigans and model
                                       predictions also generated by Oceananigans.
- `prior_mean :: Vectorß{Float64}`: prior mean
- `prior_cov :: Array{Float64, 2}`: prior covariance
- `noise_covariance :: Float64`: observation error covariance
- `α_reg :: Float64`: regularization parameter toward the prior mean (0 < `α_reg` ≤ 1);
                      default `α_reg=1` implies no regularization
- `update_freq::IT`: set to 0 when the inverse problem is not identifiable (default), namely the
                     inverse problem has multiple solutions, the covariance matrix will represent
                     only the sensitivity of the parameters, instead of posterior covariance information;
                     set to 1 (or anything > 0) when the inverse problem is identifiable, and 
                     the covariance matrix will converge to a good approximation of the 
                     posterior covariance with an uninformative prior
"""
function UnscentedKalmanInversion(inverse_problem, prior_mean, prior_cov;
                                  noise_covariance = 1e-2, α_reg = 1, update_freq = 0)
    free_parameters = inverse_problem.free_parameters
    original_priors = free_parameters.priors

    transformed_priors = [Parameterized(convert_prior(prior)) for prior in original_priors]
    no_constraints = [[no_constraint()] for _ in transformed_priors]
    parameter_distribution = ParameterDistribution(transformed_priors, no_constraints, collect(string.(free_parameters.names)))

    # Build EKP-friendly observations "y" and the covariance matrix of observational uncertainty "Γy"
    y = dropdims(observation_map(inverse_problem), dims=2) # length(forward_map_output) column vector
    Γy = construct_noise_covariance(noise_covariance, y)

    # The closure G(θ) maps (N_params, ensemble_size) array to (length(forward_map_output), ensemble_size)
    function G(θ) 
        batch_size = size(θ, 2)
        inverted_parameters = [inverse_parameter_transform.(values(original_priors), θ[:, i]) for i in 1:batch_size]
        return forward_map(inverse_problem, inverted_parameters)
    end

    ensemble_kalman_process = EnsembleKalmanProcess(y, Γy, Unscented(prior_mean, prior_cov, α_reg, update_freq))

    return EnsembleKalmanInversion(inverse_problem, parameter_distribution, ensemble_kalman_process, y, Γy, G, 0, [], Set())
end

# return:
# mean : N_iter × Nθ mean matrix
# cov  : N_iter vector of Nθ × Nθ covariance matrix
# std  : N_iter × Nθ standard deviation matrix
# err  : N_iter error array
function UnscentedKalmanInversionPostprocess(eki)

    original_priors = eki.inverse_problem.free_parameters.priors
    θ_mean_raw = hcat(eki.ensemble_kalman_process.process.u_mean...)
    θθ_cov_raw = eki.ensemble_kalman_process.process.uu_cov
    
    θ_mean = similar(θ_mean_raw)
    θθ_cov= similar(θθ_cov_raw)
    θθ_std_arr = similar(θ_mean_raw)
    for i = 1:size(θ_mean, 2)
        θ_mean[:, i] = inverse_parameter_transform.(values(original_priors), θ_mean_raw[:, i])
        θθ_cov[i] = inverse_covariance_transform(values(original_priors), θ_mean_raw[:, i], θθ_cov_raw[i])
    
        for j in 1:size(θ_mean, 1)
            θθ_std_arr[j, i] = sqrt(θθ_cov[i][j, j])
        end
    end

    return θ_mean, θθ_cov, θθ_std_arr, eki.ensemble_kalman_process.err
end

struct IterationSummary{P, E}
    parameters :: P
    mean_square_errors :: E
end

function IterationSummary(parameters, forward_map, observations)
    N_observations, N_ensemble = size(forward_map)

    mean_square_errors = [
        mapreduce((x, y) -> (x - y)^2, +, observations, view(forward_map, m, :)) / N_observations
        for m in 1:N_ensemble
    ]

    return IterationSummary(parameters, mean_square_errors)
end

function drop_ensemble_member!(eki, member)
    parameter_ensemble = eki.iteration_summary[end].parameters
    new_parameter_ensemble = vcat(parameter_ensemble[1:member-1], parameter_ensemble[member+1:end])
    new_ensemble_kalman_process = EnsembleKalmanProcess(new_parameter_ensemble, eki.mapped_observations, eki.noise_covariance, Inversion())

    push!(eki.dropped_ensemble_members, memeber)
    eki.ensemble_kalman_process = new_ensemble_kalman_process

    return nothing
end

"""
    iterate!(eki::EnsembleKalmanInversion; iterations=1)

Iterate the ensemble Kalman inversion problem `eki` forward by `iterations`.
"""
function iterate!(eki::EnsembleKalmanInversion; iterations = 1)

    first_iteration = eki.iteration + 1
    final_iteration = eki.iteration + 1 + iterations

    @suppress begin
        for iter = ProgressBar(first_iteration:final_iteration)
            θ = get_u_final(eki.ensemble_kalman_process) # (N_params, ensemble_size) array
            G = eki.inverting_forward_map(θ) # (len(G), ensemble_size)

            # Save the parameter values and mean square error between forward map
            # and observations at the current iteration
            summary = IterationSummary(θ, G, eki.mapped_observations)

            @show summary.mean_square_errors
            
            eki.iteration = iter
            push!(eki.iteration_summaries, summary)

            update_ensemble!(eki.ensemble_kalman_process, G)
        end
    end

    #=
    # All unconstrained
    params = mean(get_u_final(ensemble_kalman_process), dims=2) # ensemble mean

    # losses = [norm([G([mean(get_u(ensemble_kalman_process, i), dims=2)...])...]  - y) for i in 1:n_iterations] # particle loss
    mean_vars = [diag(cov(get_u(ensemble_kalman_process, i), dims=2)) for i in 1:iterations] # ensemble variance for each parameter

    original_priors = eki.inverse_problem.free_parameters.priors
    params = inverse_parameter_transform.(values(original_priors), params)

    # return params, mean_vars, mean_us_constrained
    mean_us_eki = [[mean(get_u(ensemble_kalman_process, i), dims=2)...] for i in 1:iterations]
    mean_us_constrained = [inverse_parameter_transform.(values(original_priors), mean_u) for mean_u in mean_us_eki]
    return [params...], mean_vars, mean_us_constrained
    =#

    return nothing
end

end # module