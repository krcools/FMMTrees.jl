export reconstruct, storedentries

struct LowRankMatrix{T}
    A::Matrix{T}
    B::Matrix{T}
end

Base.:*(M::LowRankMatrix, x::Vector) = (M.A*(M.B*x))

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



storedentries(lrm::LowRankMatrix) = length(lrm.A) + length(lrm.B)
storedentries(lrb::LowRankBlock) = storedentries(lrb.matrix)

function storedentries(A::HMatrix)
    sz = 0
    for block in A.blocks
        sz += storedentries(block)
    end
    return sz
end

function reconstruct(A::HMatrix)
    M = zeros(eltype(A), size(A)...)
    for block in A.blocks
        @assert norm(M[block.τ,block.σ]) == 0
        M[block.τ,block.σ] = block.matrix.A * block.matrix.B
    end
    return M
end
