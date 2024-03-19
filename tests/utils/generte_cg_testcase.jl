### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ d45aa6a5-9cfa-4016-9375-c55343307996
using Random

# ╔═╡ c499c797-a2e1-4cb0-bff5-0d8c47bed734
using LinearAlgebra

# ╔═╡ ec0a14bd-21f3-4e88-8609-0c9c251a7bb1
using IterativeSolvers

# ╔═╡ 22373066-a4c4-4d12-94fb-b6acad925cb6
md"
Julia/Pluto notebook to generate and validate test cases for the implementation of conjugate gradients.
"

# ╔═╡ 045688d1-c90f-47b3-a663-35b5ce702925
n = 3

# ╔═╡ 43b1423a-e540-11ee-00de-4bfda149915c
begin
	Random.seed!(0)
	R = rand(Float32, (n,n))
end

# ╔═╡ 76666054-9ded-4838-b150-494afeb1cedf
A = R' * R + 2*n*I

# ╔═╡ d6ead66e-d307-4f75-8fa2-287c250ebad7
begin
	Random.seed!(0)
	b = rand(Float32, 3)
end

# ╔═╡ 81a575a4-a491-404f-9fd3-7e114e3df602
x = inv(A)*b

# ╔═╡ 1cb479aa-9f94-4720-bd55-17c5ded594c3
[abs(x1 - x2) < 0.00001 for (x1, x2) in zip(b, A*x)]

# ╔═╡ 232c9445-2c71-439e-b61e-368c35e07637
x2 = cg(A, b)

# ╔═╡ cf454a74-2d97-4956-b9ec-a104e963bc8a
[abs(x1 - x2) < 0.00001 for (x1, x2) in zip(b, A*x2)]

# ╔═╡ 299bd1aa-3450-4f37-8dd8-f74753cd1ae0
"""
    conjugate_gradient!(A, b, x)

Return the solution to `A * x = b` using the conjugate gradient method.
"""
function conjugate_gradient!(
    A::AbstractMatrix, b::AbstractVector, x::AbstractVector; tol=eps(eltype(b))
)
    # Initialize residual vector
    residual = b - A * x
    # Initialize search direction vector
    search_direction = residual
    # Compute initial squared residual norm
	norm(x) = sqrt(sum(x.^2))
    old_resid_norm = norm(residual)

    # Iterate until convergence
    while old_resid_norm > tol
        A_search_direction = A * search_direction
        step_size = old_resid_norm^2 / (search_direction' * A_search_direction)
        # Update solution
        @. x = x + step_size * search_direction
        # Update residual
        @. residual = residual - step_size * A_search_direction
        new_resid_norm = norm(residual)
        
        # Update search direction vector
        @. search_direction = residual + 
            (new_resid_norm / old_resid_norm)^2 * search_direction
        # Update squared residual norm for next iteration
        old_resid_norm = new_resid_norm
    end
    return x
end

# ╔═╡ a8460106-1302-4574-8c3d-1efed8f27623
x3 = conjugate_gradient!(A, b, [0.0;0.0;0.0])

# ╔═╡ b30da3e1-30fc-43e1-b949-47cd1f75eed5
[abs(x1 - x2) < 0.0001 for (x1, x2) in zip(b, A*x3)]

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
IterativeSolvers = "~0.9.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.2"
manifest_format = "2.0"
project_hash = "336d853d46b26fbc38a3aac7a62e3fe045151dd0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.0+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "59545b0a2b27208b0650df0a46b8e3019f85055b"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.4"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"
"""

# ╔═╡ Cell order:
# ╟─22373066-a4c4-4d12-94fb-b6acad925cb6
# ╠═d45aa6a5-9cfa-4016-9375-c55343307996
# ╠═c499c797-a2e1-4cb0-bff5-0d8c47bed734
# ╠═ec0a14bd-21f3-4e88-8609-0c9c251a7bb1
# ╠═045688d1-c90f-47b3-a663-35b5ce702925
# ╠═43b1423a-e540-11ee-00de-4bfda149915c
# ╠═76666054-9ded-4838-b150-494afeb1cedf
# ╠═d6ead66e-d307-4f75-8fa2-287c250ebad7
# ╠═81a575a4-a491-404f-9fd3-7e114e3df602
# ╠═1cb479aa-9f94-4720-bd55-17c5ded594c3
# ╠═232c9445-2c71-439e-b61e-368c35e07637
# ╠═cf454a74-2d97-4956-b9ec-a104e963bc8a
# ╠═299bd1aa-3450-4f37-8dd8-f74753cd1ae0
# ╠═a8460106-1302-4574-8c3d-1efed8f27623
# ╠═b30da3e1-30fc-43e1-b949-47cd1f75eed5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
