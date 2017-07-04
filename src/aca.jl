# μ is a lazy array. After all, the whole point of doing ACA is to
# avoid explicit computation of all entries of the matrix under consideration.
function adaptive_cross_approximation(μ,ϵ=sqrt(eps(eltype(μ))))

    T = eltype(μ)
    α = zeros(T)
    a = zeros()

    I = Set(indices(μ,1))
    J = Set(indices(μ,2))

    l, j = 0, first(J)
    while true
        l += 1
        i = indmax(abs(M[k,j] - R[k,j]) for k in 1:m)
    end

end
