using Random, LinearAlgebra, Statistics, Distributions, StatsBase

data = vcat(
    rand(Normal(-1.0, 0.5), 500),
    rand(Normal( 2.0, 0.5), 1000)
)

hist = append!(Histogram(-2:0.1:4), data)

using Plots
plot(
    normalize(hist, mode=:density),
    st = :steps, label = "Data",
    title = "Data"
)

function fit_function(p::NamedTuple{(:a, :mu, :sigma)}, x::Real)
    p.a[1] * pdf(Normal(p.mu[1], p.sigma), x) +
    p.a[2] * pdf(Normal(p.mu[2], p.sigma), x)
end

true_par_values = (a = [500, 1000], mu = (-1.0, 2.0), sigma = 0.5)
plot(
    normalize(hist, mode=:density),
    st = :steps, label = "Data",
    title = "Data and True Statistical Model"
)
plot!(
    -4:0.01:4, x -> fit_function(true_par_values, x),
    label = "Truth"
)

using BAT, IntervalSets

likelihood = let h = hist, f = fit_function
    # Histogram counts for each bin as an array:
    observed_counts = h.weights

    # Histogram binning:
    bin_edges = h.edges[1]
    bin_edges_left = bin_edges[1:end-1]
    bin_edges_right = bin_edges[2:end]
    bin_widths = bin_edges_right - bin_edges_left
    bin_centers = (bin_edges_right + bin_edges_left) / 2

    params -> begin
        # Log-likelihood for a single bin:
        function bin_log_likelihood(i)
            # Simple mid-point rule integration of fit function `f` over bin:
            expected_counts = bin_widths[i] * f(params, bin_centers[i])
            logpdf(Poisson(expected_counts), observed_counts[i])
        end

        # Sum log-likelihood over bins:
        idxs = eachindex(observed_counts)
        ll_value = bin_log_likelihood(idxs[1])
        for i in idxs[2:end]
            ll_value += bin_log_likelihood(i)
        end

        # Wrap `ll_value` in `LogDVal` so BAT knows it's a log density-value.
        return LogDVal(ll_value)
    end
end

likelihood(true_par_values)

using ValueShapes

prior = NamedTupleDist(
    a = [Weibull(1.1, 5000), Weibull(1.1, 5000)],
    mu = [-2.0..0.0, 1.0..3.0],
    sigma = Weibull(1.2, 2)
)

posterior = PosteriorDensity(likelihood, prior)
samples = bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^7, nchains = 4)).result

plot(
    samples, :(mu[1]),
    mean = true, std = true, globalmode = true, marginalmode = true,
    nbins = 50, title = "Marginalized Distribution for mu[1]"
)