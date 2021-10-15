module Observations

using Oceananigans
using Oceananigans: short_show, fields
using Oceananigans.Grids: AbstractGrid
using Oceananigans.Fields
using Oceananigans.Utils: SpecifiedTimes
using JLD2

import Oceananigans.Fields: set!

abstract type AbstractObservation end

include("normalization.jl")

"""
    OneDimensionalTimeSeries{F, G, T, P, M} <: AbstractObservation

A time series of horizontally-averaged observational or LES data
gridded as Oceananigans fields.
"""
struct OneDimensionalTimeSeries{F, G, T, P, M, N} <: AbstractObservation
    field_time_serieses :: F
                   grid :: G
                  times :: T
                   path :: P
               metadata :: M
          normalization :: N
end

obs_str(ts::OneDimensionalTimeSeries) = "OneDimensionalTimeSeries of $(keys(ts.field_time_serieses)) on $(short_show(ts.grid))"

tupleit(t) = try
    Tuple(t)
catch
    tuple(t)
end

const not_metadata_names = ("serialized", "timeseries")

read_group(group::JLD2.Group) = NamedTuple(Symbol(subgroup) => read_group(group[subgroup]) for subgroup in keys(group))
read_group(group) = group

get_grid(field_time_serieses) = field_time_serieses[1].grid

function OneDimensionalTimeSeries(path; field_names, normalize=IdentityNormalization, times=nothing)
    field_names = tupleit(field_names)
    field_time_serieses = NamedTuple(name => FieldTimeSeries(path, string(name); times) for name in field_names)
    grid = get_grid(field_time_serieses)
    times = first(field_time_serieses).times

    # validate_data(fields, grid, times) # might be a good idea to validate the data...
    file = jldopen(path)
    metadata = NamedTuple(Symbol(group) => read_group(file[group]) for group in filter(n -> n ∉ not_metadata_names, keys(file)))
    close(file)

    normalization = Dict(name => normalize(field_time_serieses[name]) for name in keys(field_time_serieses))

    return OneDimensionalTimeSeries(field_time_serieses, grid, times, path, metadata, normalization)
end

observation_times(observation) = observation.times
observation_times(obs::Vector) = observation_times(first(obs))

#####
##### set! for simulation models and observations
#####

default_initial_condition(ts, name) = 0

function set!(model, ts::OneDimensionalTimeSeries, index=1)
    # Set initial condition
    for name in keys(fields(model))

        model_field = fields(model)[name]

        if name in keys(ts.field_time_serieses)
            ts_field = ts.field_time_serieses[name][index]
            set!(model_field, ts_field)
        else
            set!(model_field, 0) #default_initial_condition(ts, Val(name)))
        end
    end

    return nothing
end

"""
    column_ensemble_interior(observations::Vector{<:OneDimensionalTimeSeries}, field_name, time_indices::Vector, N_ens)

Returns an `N_cases × N_ens × Nz` array of the interior of a field `field_name` defined on a 
`OneDimensionalEnsembleGrid` of size `N_cases × N_ens × Nz`, given a list of `OneDimensionalTimeSeries` objects
containing the `N_cases` single-column fields at time index in `time_index`.
"""
function column_ensemble_interior(observations::Vector{<:OneDimensionalTimeSeries}, field_name, time_index, ensemble_size)
    batch = [ts.field_time_serieses[field_name][time_index] for ts in observations]

    for ts in observations
        
    end

    batch = cat(batch..., dims = 2) # (n_batch, n_z)
    return cat([batch for i = 1:ensemble_size]..., dims = 1) # (ensemble_size, n_batch, n_z)
end

function set!(model, ts::Vector{<:OneDimensionalTimeSeries}, index=1)
    # Set initial condition
    for name in keys(fields(model))

        model_field = fields(model)[name]

        if name in keys(ts.field_time_serieses)
            ts_field = ts.field_time_serieses[name][index]
            set!(model_field, ts_field)
        else
            set!(model_field, 0) #default_initial_condition(ts, Val(name)))
        end
    end

    return nothing
end

## should maybe include interpolation functionality
function set!(field::AbstractField, data::AbstractArray)

    arch = architecture(field)

    # Reshape `data` to the size of `field`'s interior
    reshaped_data = arch_array(arch, reshape(data, size(field)))

    # Sets the interior of field `field` to values of `data`
    field .= reshaped_data
end

struct FieldTimeSeriesCollector{G, D, F, T}
    grid :: G
    field_time_serieses :: D
    collected_fields :: F
    times :: T
end

"""
    FieldTimeSeriesCollector(fields)

Returns a `FieldTimeSeriesCollector` for `fields` of `simulation`.

`fields` is a `NamedTuple` of `AbstractField`s that are to be collected.
"""
function FieldTimeSeriesCollector(collected_fields, times; architecture=CPU())

    grid = collected_fields[1].grid
    field_time_serieses = Dict{Symbol, Any}()

    for name in keys(collected_fields)
        field = collected_fields[name]
        LX, LY, LZ = location(field)
        field_time_series = FieldTimeSeries{LX, LY, LZ}(architecture, field.grid, times)
        field_time_serieses[name] = field_time_series
    end

    # Convert to NamedTuple
    field_time_serieses = NamedTuple(name => field_time_serieses[name] for name in keys(collected_fields))

    return FieldTimeSeriesCollector(grid, field_time_serieses, collected_fields, times)
end

function (collector::FieldTimeSeriesCollector)(simulation)
    for field in collector.collected_fields
        compute!(field)
    end

    current_time = simulation.model.clock.time
    time_index = findfirst(t -> t >= current_time, collector.times)

    for name in keys(collector.collected_fields)
        field_time_series = collector.field_time_serieses[name]
        set!(field_time_series[time_index], collector.collected_fields[name])
    end

    return nothing
end

function initialize_simulation!(simulation, ts, time_series_collector, time_index=1)
    set!(simulation.model, ts, time_index) 

    times = observation_times(ts)

    initial_time = ts.times[time_index]
    simulation.model.clock.time = initial_time
    simulation.model.clock.iteration = 0

    # Zero out time series data
    for time_series in time_series_collector.field_time_serieses
        time_series.data .= 0
    end

    simulation.callbacks[:data_collector] = Callback(time_series_collector, SpecifiedTimes(times...))

    simulation.stop_time = times[end]

    return nothing
end

summarize_metadata(::Nothing) = ""
summarize_metadata(metadata) = keys(metadata)

Base.show(io::IO, ts::OneDimensionalTimeSeries) =
    print(io, "OneDimensionalTimeSeries with fields $(propertynames(ts.field_time_serieses))", '\n',
              "├── times: $(ts.times)", '\n',    
              "├── grid: $(short_show(ts.grid))", '\n',
              "├── path: \"$(ts.path)\"", '\n',
              "├── metadata: ", summarize_metadata(ts.metadata), '\n',
              "└── normalization: $(short_show(ts.normalization))")

end # module