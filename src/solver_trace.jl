mutable struct EnergyColumn{T<:AbstractFloat} <: TraceColumn
    E::T
    fmt::FormatExpr
    header::String
    eV::Bool
end

function EnergyColumn(E::T, eV::Bool=true) where T
    header = "Energy"
    fmt = if eV
        FormatExpr("{1:+10.5f} Ha {2:+10.5f} eV")
    else
        FormatExpr("{1:+10.5f} Ha")
    end
    EnergyColumn(E, fmt,
                 rpad(header, length(eV ? format(fmt, E, E) : format(fmt, E))), eV)
end

(e::EnergyColumn)(::Integer) =
    e.eV ? (e.E,27.211e.E) : (e.E,)

mutable struct RelaxationColumn{T<:AbstractFloat} <: TraceColumn
    ω::T
    fmt::FormatExpr
    header::String
end

function RelaxationColumn(ω::T) where T
    header = "1-ω"
    fmt = FormatExpr("{1:.2f}×10{2:3s}")
    RelaxationColumn(ω, fmt,
                     rpad(header, length(format(fmt, base_exp(abs(1.0 - ω))...))))
end

function (r::RelaxationColumn{T})(::Integer) where T
    ωb,ωe = base_exp(abs(1.0 - r.ω))
    (ωb,to_superscript(ωe))
end
