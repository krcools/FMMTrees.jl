export reconstruct, storedentries

struct LowRankMatrix{T}
    A::Matrix{T}
    B::Matrix{T}
end

Base.:*(A::LowRankMatrix, x::Vector) = (A*(B*x))

struct LowRankBlock{T}
    matrix::LowRankMatrix{T}
    τ::Vector{Int}
    σ::Vector{Int}
end

struct HMatrix{T}
    blocks::Vector{LowRankBlock{T}}
    I::Vector{Int}
    J::Vector{Int}
end

Base.eltype(A::HMatrix{T}) where {T} = T
Base.size(A::HMatrix) = (length(A.I), length(A.J))

function Base.:*(A::HMatrix, x)
    y = zeros(length(A.I))
    for block in A.blocks
        y[block.τ] .+= block.matrix * x[block.σ]
    end
    return y
end


function storedentries(A::HMatrix)
    sz = 0
    for block in A.blocks
        sz += length(block.matrix.A) + length(block.matrix.B)
    end
    return sz
end


function reconstruct(A::HMatrix)
    M = zeros(eltype(A), size(A)...)
    for block in A.blocks
        M[block.τ,block.σ] = block.matrix.A * block.matrix.B
    end
    return M
end
