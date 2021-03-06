### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 449faf07-3214-47b6-bb02-824cf900bc07
begin
	using PlutoUI # for the @bind macro
	using Plots # plotting front end
end

# ╔═╡ 31bc44b6-ade1-4277-87dc-02379091d2c2
html"""
<div style="
position: absolute;
width: calc(100% - 30px);
border: 50vw solid #282936;
border-top: 100px solid #282936;
border-bottom: none;
box-sizing: content-box;
left: calc(-50vw + 15px);
top: -500px;
height: 200px;
pointer-events: none;
"></div>

<div style="
height: 200px;
width: 100%;
background: #282936;
color: #fff;
padding-top: 0px;
">
<span style="
font-family: Vollkorn, serif;
font-weight: 700;
font-feature-settings: 'lnum', 'pnum';
"> <p style="
font-size: 1.5rem;
opacity: .8;
"><em>Section x.x</em></p>
<p style="text-align: center; font-size: 2rem;">
<em> Timestepping </em>
</p>

<p style="
font-size: 1.5rem;
text-align: center;
opacity: .8;
"><em>Finite-difference equations</em></p>
<div style="display: flex; justify-content: center;">
</div>
</div>

<style>
body {
overflow-x: hidden;
}
</style>"""

# ╔═╡ b44f6459-6816-4747-9823-4a1575a77578
md"""
# Notebook packages
"""

# ╔═╡ 439520d7-7453-43f6-92ef-91524e598d12
md"""
Using `Pluto`'s package management system for reproducible results.
"""

# ╔═╡ 2dd2d4e3-14a0-43f1-b207-969b98ee5346
PlutoUI.TableOfContents(aside=true)

# ╔═╡ 64d42218-59a4-49b0-b260-f07ee0d0d5bd
md"""
# Abstract
"""

# ╔═╡ 89f1748e-593f-473f-907d-b3a64073e197
md"""
## Science and math concepts

- Finite difference equations
- Their solution
- Time stepping
"""

# ╔═╡ f878dcc5-06db-4dd4-8df8-c4666cf8f605
md"""
## Julia and numerical concepts

- Array comprehension
- Plotting, using `Plots` package
- Functions
"""

# ╔═╡ f5e86f2a-076d-11ec-1080-c341c92d3fd1
md"""
# Physical laws
"""

# ╔═╡ 26f6893e-17ab-48e9-88d0-98a79ec0a1b0
md"""
## Change is what matters
"""

# ╔═╡ 40b5474e-a8eb-4f6e-bcba-d49f714c83f2
md"""
A physical law for a quantity (such as an electric field, temperature, quantum mechanical wavefunction) is almost always (or always?) formulated in terms of how that quantity *changes*. The *rate of change* of that quantity is the thing that appears in the equation that formulates the law.

For instance, some quantity ``x`` may evolve in time ``t`` according to
```math
x(t+\Delta t) = x(t) + f(x(t),t)\Delta t,
```
where ``f(x(t),t)`` is the rate of change of ``x``. It can depend on the value of ``x`` at the time ``t`` as well as explicitly on the time. By the way, this type of equation is called a **finite difference equation**.
"""

# ╔═╡ ecb886b8-cf77-4d82-865b-d2af63f22a2e
md"""
## Time stepping

Note that this equation itself is not enough to find ``x`` at all times, which is what we want to do. The equation tells us what is involved in order to *update* the value of ``x``, which means we need a value to start with, say at ``t=0``.

If we know ``x(0)``, then the equation tells us how to find the value at some future time ``\Delta t``. We can then use ``x(\Delta t)`` to find ``x(2\Delta t)``, and so on to the final time we're interested in. The starting value ``x(0)`` (or whatever the starting time is) is called the **initial value**, and this is referred to as an **initial value problem**.

But while we know what determines how ``x`` evolves (namely, ``f``), the law does not specify the technicalities of *how* to do the updating or time stepping. It can sometimes be more complicated that it appears...
"""

# ╔═╡ 42986972-b891-46ae-8462-93a31924e6fe
md"""
### ``f`` is constant in time
"""

# ╔═╡ 2d13c6f7-4212-4bde-8e3d-1e46bc7992a6
md"""
Suppose ``f`` does not depend on time or ``x``. Then it's constant over any interval ``\Delta t``. This is a case where the updating is easy. Whenever ``f`` is needed you just use that constant value.
"""

# ╔═╡ 1dd0d637-4419-4475-9f36-4da7b8107e9a
md"""
### ``f`` depends on time
"""

# ╔═╡ c357aaef-725a-40b2-aadf-2bd5a26ab997
md"""
Now, there's a question: If ``f(t)`` changes over the time interval ``\Delta t``, at what time should we evaluate ``f(t)`` at?

A couple of options come to mind:
1. Take ``f`` to be constant and equal to that at the *start* of the time increment, or
2. Take ``f`` to be constant and equal to that in the *middle* of the time increment, or
3. Take ``f`` to be constant and equal to the *average* speed during the time increment.
"""

# ╔═╡ 8e1fafdb-bdf0-4649-aff7-7720280f466a
md"""
### ``f`` depends on ``x(t)`` and maybe also time
"""

# ╔═╡ 2e7cbd7c-c650-4088-ac92-3490f2db8e2c
md"""
This is almost the same situation as the last one, but now ``f`` depends on the actual solution ``x(t)`` and perhaps also has explicit dependence on ``t``. A couple of examples are ``f = x^2`` or ``f = t\, \sin(x)``.
"""

# ╔═╡ 2b1e6d9b-abed-4633-a52e-a84514a6c27a
md"""
# Computer experiments
"""

# ╔═╡ 07c5d7e1-574c-4f3b-b879-6b2842d862c1
md"""
## ``f`` is constant in time
"""

# ╔═╡ 6e215ec5-dafb-4eb0-a159-b9028d6de08d
md"""
Let's first illustrate the calculation of the solution ``x(t)`` for the simplest case of a constant ``f``.
"""

# ╔═╡ c87acb83-84e6-484a-bed7-319bdb3d0feb
md"""
We will split time up into a number of intervals, yielding an array of time values from the intial to final times. This is often called the **time grid**. We'll set it up using **array comprehensions**, which we've seen before.
"""

# ╔═╡ 263c3ed1-18b4-45eb-9de9-48b4f4f5c417
begin
	Δt = 0.2
	nmax = 10  # final time = Δt*nmax
	t = [n*Δt for n in 0:nmax]
end

# ╔═╡ a1f373d4-4c8c-448b-b336-82288f712f7c
md"""
> 🤔 Question: How do you reference or call the initial time? What about the final time? (careful) How would you make the initial time to be 5?
>
> Answer:
"""

# ╔═╡ fcba2327-2306-4fc3-bb70-14e94e5d55c5
md"""
Define a function to calculate the change in ``x`` over a time interval. Allow the value of ``f`` to be passed, but use a default value.
"""

# ╔═╡ 809e639a-8492-4b0a-bc98-2226c2404478
function Δx(Δt, f=10.0)
	Δt*f
end

# ╔═╡ 8d015906-914e-4ea5-935f-9cde818000af
md"""
> ☡ Coding pointer: The output of the function definition says you've defined a function with 2 methods. This is often called function **overloading**. It's one function that can be executed 2 different ways: with or without passing the function ``f`` in this case.
"""

# ╔═╡ d6c69e40-e383-4588-8130-09ac2a2a69a9
begin
	x = similar(t)
	x[1] = 0.0 # initial value
	for i = 2:nmax+1
		x[i] = x[i-1] + Δx(Δt)
	end
end

# ╔═╡ 1fbc2b57-aa8a-47e8-8886-603e15474861
md"""
> ☡ Coding pointer: The `nmax+1` can be tricky. Best to use `length(x)`; that way you don't have to remember how long the array is, or scroll back to an earlier point in the code to look for the variable name for it. We will do that from now on.
>
>Also, the `similar` command makes an array of the same size, but the values it contains are *uninitialized*. You need to make any initial value settings.
"""

# ╔═╡ 6fd61abc-6bf3-4df7-b6bc-8adfafc461e0
md"""
We now have ``x`` at each ``t`` value and can plot.
"""

# ╔═╡ 8d22438a-17ce-4bb6-814f-eecbd371fb71
begin
	plot(t,x,label="x(t)")
	plot!(title="x vs time",xlabel="time [some units]",ylabel="x [some units]")
end

# ╔═╡ 97216c01-ac6f-439e-8bf0-8ade70e08e03
md"""
## ``f`` depends on time
"""

# ╔═╡ 1edc180d-c2fe-4e30-8650-de2952b97a14
md"""
We mentioned three options -- and there are more -- for handling this situation. Again, the problem is that as we step over ``\Delta t``, ``f`` is changing. So, where should we evaluate ``f``? At what time?
"""

# ╔═╡ 002fab49-b0c0-4973-b676-473daaf2c698
md"""
Pick some ``f`` that depends on time for illustration.
"""

# ╔═╡ 0656317b-3470-4d97-8b5d-9eddaab3ee8c
function f(t)
	2*t
end

# ╔═╡ dd4a6136-75a7-460d-9bbc-97a374e23a4f
function Δx_start(t₀, Δt, f::Function)
	f(t₀)*Δt
end

# ╔═╡ 44762363-cc3e-422a-8044-670bc04b1690
function Δx_mid(t₀, Δt, f::Function)
	f(t₀ + Δt/2)*Δt
end

# ╔═╡ 9059e6d7-ad27-4b06-b740-a91b9bc9f349
function Δx_ave(t₀, Δt, f::Function)
	f_ave = (f(t₀) + f(t₀ + Δt))/2
	f_ave*Δt
end

# ╔═╡ f28a6182-2310-4444-a90e-402afd923bf9
md"""
> 🤔 Question: Instead of this usual average, what other sort of average could you consider?
>
> Answer:
"""

# ╔═╡ 453ad11f-7bed-49f7-8255-6b4e7a86ae57
md"""
Let's now calculate the ``x`` values for these three cases. We'll take the initial condition to be again ``x(0) = 0``.
"""

# ╔═╡ cc9737b7-33ee-48c4-9cbd-dde3f9b28afc
begin
	x_start = similar(t)
	x_mid = similar(t)
	x_ave = similar(t)
	x_start[1] = 0.0  # initial value
	x_mid[1] = 0.0  # initial value
	x_ave[1] = 0.0  # initial value
end

# ╔═╡ e59e1646-d949-456d-8e28-a4f9e1f55a8a
for i in 2:length(t)
	x_start[i] = x_start[i-1] + Δx_start(t[i-1], Δt, f)
	x_mid[i] = x_mid[i-1] + Δx_mid(t[i-1], Δt, f)
	x_ave[i] = x_ave[i-1] + Δx_ave(t[i-1], Δt, f)
end

# ╔═╡ 421f19c1-38eb-4cff-b8d4-7a6ae4f438d0
md"""
As a quick check of the results, show the arrays and make sure nothing diverges, for instance. 
"""

# ╔═╡ a87dd44d-89ca-41cc-b077-d444ea37c48d
x_start

# ╔═╡ becfff3a-b130-44e8-8f9b-f63d693fdc0d
x_mid

# ╔═╡ 20e80196-d302-44ad-af58-c2c6211bcc95
x_ave

# ╔═╡ 202922e8-9d1a-4a23-afdc-f91deceb899e
md"""
> 🤔 Question: Do you notice anything interesting here, even without plotting and going any further?
>
> Answer:
"""

# ╔═╡ d98ad456-1f29-4047-8724-78a8cb35187a
md"""
Let's now make some plots of ``x`` as a function of time.
"""

# ╔═╡ bafa445a-f3c8-4e83-bfa7-2b907bd33319
begin
	plot(t,x_start,label="start")
	plot!(t,x_mid,label="mid")
	plot!(t,x_ave,label="ave")
	plot!(title="3 choices for evaluating f",xlabel="time [some units]",ylabel="x [some units]")
end

# ╔═╡ 61f10f8b-9f83-45ae-815f-22a86180dcfe
md"""
It looks like two of the curves coincide: the ones for `x_mid` and `x_ave`. This is in fact the case, which can be verified by commenting out the `plot` command for one of them. Be sure to try it and convince yourself!

*And it's definitely the case that the choice of where to evaluate ``f(t)`` affects the solution!*
"""

# ╔═╡ ec533632-19dc-47e6-8491-3da3c2a3fe8a
md"""
> 🤔 Question: How does the choice of Δt change the results? Make is smaller and larger than the above value of 0.2. But reset it to 0.2 when you're done for the next section.
>
> Answer:
"""

# ╔═╡ c416bf8c-8096-4ee8-bbfe-c5287f52d732
md"""
### Exact result
"""

# ╔═╡ c39566f3-f184-4cba-9dfb-390a9a5a8448
md"""
We know the exact result for ``x(t)``. It's
```math
x(t) = t^2 .
```
We can compare this exact result with the choices above.
"""

# ╔═╡ d844d41b-35ee-4391-a36e-cffc5820cf30
begin
	x_texact = similar(t)
	x_texact[1] = 0.0
	for i = 2:length(x_texact)
		x_texact[i] = t[i]^2
	end
end

# ╔═╡ 934a56c2-001f-4d7b-9f22-44b67c1b0481
begin
	plot(t,x_start,label="start")
	plot!(t,x_mid,label="mid")
	#plot!(t,x_ave,label="ave")
	plot!(t,x_texact,label="exact")
	plot!(title="With exact result",xlabel="time [some units]",ylabel="x [some units]")
end

# ╔═╡ 6a0e0537-db4a-46b2-b5a1-c9afbf98a028
md"""
The exact result matches the midpoint or average evaluation of ``f``. When we let ``Δt`` approach zero, all results converge. *But otherwise, we see that evaluating at the beginning of the interval is not so good.*
"""

# ╔═╡ 6a167ed2-45b4-4a09-a7e1-e24478c6ed38
md"""
## ``f`` depends on ``x(t)`` and maybe also time
"""

# ╔═╡ a48c3f1a-012c-4de3-b710-6c24bfb025ce
md"""
Consider another illustrative example of the function ``f``, in this case ``f = x^2``.
"""

# ╔═╡ 02f6c327-048a-42ad-8d51-116fb2ba42fe
md"""
> ☡ Coding note: This is Pluto notebook issue... we can't redefine a quantity, in this case the function ``f``. That would destroy the feature that makes Pluto unique! Namely, it's reactivity. So we'll call it ``g`` for this case.
"""

# ╔═╡ 858520cf-cc81-4fa4-afa1-b412bc43e062
function g(x)
	x^2
end

# ╔═╡ 78a73608-de13-438c-822f-d2eb4b4a0a5b
md"""
We'll have the same Pluto issue when defining the functions to increment ``x``, so we need to change the name a bit. And while we don't need the time at the start of the interval here, we'll include it for generality and ease of modification later.
"""

# ╔═╡ 5785a365-738f-41f9-907a-133ca9cddc67
function Δx_fstart(t₀, x₀, Δt, f::Function)
	f(x₀)*Δt
end

# ╔═╡ 07651aee-71aa-46c1-9e70-636db648b0fd
md"""
> ☡ Coding pointer: Notice that within this updating function we can still use ``f`` as a variable name. We don't need to change it to ``g``. This is because ``f`` is **local** to `Δx_gstart`.
"""

# ╔═╡ 75b4d476-9a25-41c3-97ad-c61b7bba7027
function Δx_fmid(t₀, x₀, Δt, f::Function)
	x_mid = x₀ + f(x₀)*Δt/2
	f(x_mid)*Δt
end

# ╔═╡ aa27d3ae-12da-45f1-af6b-3fdfed216835
function Δx_fave(t₀, x₀, Δt, f::Function)
	x_end = x₀ + f(x₀)*Δt
	f_ave = (f(x_end) + f(x₀))/2
	f_ave*Δt
end

# ╔═╡ 6afbfaed-9e8e-4b20-96d1-15b3832ceb78
md"""
Let's now calculate the ``x`` values for these three cases. We'll take the initial condition to be this time ``x(0) = 0.4``.
"""

# ╔═╡ a9bdbc6d-2cfd-47dd-840c-f495130601c1
begin
	x_fstart = similar(t)
	x_fmid = similar(t)
	x_fave = similar(t)
	x_fstart[1] = 0.4  # initial value
	x_fmid[1] = 0.4  # initial value
	x_fave[1] = 0.4  # initial value
end

# ╔═╡ 695669d8-f368-4156-84df-42fd40653758
for i in 2:length(t)
	x_fstart[i] = x_fstart[i-1] + Δx_fstart(t[i-1], x_fstart[i-1], Δt, g)
	x_fmid[i] = x_fmid[i-1] + Δx_fmid(t[i-1], x_fmid[i-1], Δt, g)
	x_fave[i] = x_fave[i-1] + Δx_fave(t[i-1], x_fave[i-1], Δt, g)
end

# ╔═╡ 61b28060-5ee9-4699-a9d5-0628537f14d9
x_fstart

# ╔═╡ 5e79a3d1-d403-437c-8da3-2c0216c8432e
x_fmid

# ╔═╡ 0de24a1d-dbb1-4944-97ed-97cba4428a35
x_fave

# ╔═╡ 9e59294b-40fe-4c21-8496-bd85680f15b6
md"""
Interesting... now *all three choices give different results*, but the middle and average cases are very close.
"""

# ╔═╡ 3861a74a-bb31-439d-a359-233333b5c536
begin
	plot(t,x_fstart,label="start")
	plot!(t,x_fmid,label="mid")
	plot!(t,x_fave,label="ave")
	plot!(title="3 choices for evaluating f",xlabel="time [some units]",ylabel="x [some units]")
end

# ╔═╡ 8d3c3463-c9bd-414d-a0ca-ef6d214155ed
md"""
> 🤔 Question: How does the choice of Δt change the results? Make is smaller and larger than the above value of 0.2, but reset to 0.2 when done.
>
> Answer:
"""

# ╔═╡ bf1a5cb5-4073-4d7a-989e-4b8f6bd97415
md"""
### Exact result
"""

# ╔═╡ 780b02c4-5fc5-400b-9aac-09345c2769a0
md"""
Now, it turns out we know the exact result for ``x(t)`` here. It's
```math
x(t) = \frac{1}{x(0)^{-1} - t} .
```
We can compare this exact result with the choices above.
"""

# ╔═╡ dbf96f31-84d4-45c0-9e37-17a277372de8
begin
	x_exact = similar(t)
	x_exact[1] = 0.4
	for i = 2:length(t)
		x_exact[i] = 1/( (1/x_exact[1]) - t[i] )
	end
end

# ╔═╡ 629aae68-24f2-4cbd-99f3-1d964c723247
begin
	plot(t,x_fstart,label="start")
	plot!(t,x_fmid,label="mid")
	#plot!(t,x_fave,label="ave")
	plot!(t,x_exact,label="exact")
	plot!(title="With exact result",xlabel="time [some units]",ylabel="x [some units]")
end

# ╔═╡ 968338b0-4536-45e5-ac12-573b939b9393
md"""
The exact result is much closer to the middle or average choice, but still shows some divergence at later times. You can experiment with decreasing Δt, and you'll find that the middle and average results converge to the exact result. The starting choice is better, but not as good as the middle or average choices.
"""

# ╔═╡ acda2d19-e948-4b19-9da6-63b5a16c2e08
md"""
# Conclusions
"""

# ╔═╡ ed702970-dde8-494f-8725-2b89bfb4d292
md"""
- A finite difference equation can be solved numerically.
- For a given time step ``\Delta t``, evaluating ``f`` at the middle of the time interval or using an average over the time interval is the most accurate.
- But, as ``\Delta t`` goes to zero (i.e., becomes very small), all of the methods for evaluating ``f`` during a timestep approach the same result. This will lead to the notion of a **differential equation** as opposed to a finite difference one.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Plots = "~1.21.2"
PlutoUI = "~0.7.9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c3598e525718abcc440f69cc6d5f60dda0a1b61e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.6+5"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "e2f47f6d8337369411569fd45ae5753ca10394c6"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.0+6"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "9995eb3977fbf67b86d0a0a0508e83017ded03f2"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.14.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "727e463cfebd0c7b999bbf3e9e7e16f254b94193"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.34.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[DataAPI]]
git-tree-sha1 = "ee400abb2298bd13bfc3df1c412ed228061a2385"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.7.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "92d8f9f208637e8d2d28c664051a00569c01493d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.1.5+1"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "LibVPX_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "3cc57ad0a213808473eafef4845a74766242e05f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.3.1+4"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "35895cf184ceaab11fd778b4590144034a167a2f"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.1+14"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "cbd58c9deb1d304f5a245a0b7eb841a2560cfec6"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.1+5"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "0c603255764a1fa0b61752d2bec14cfbd18f7fe8"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+1"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "182da592436e287758ded5be6e32c406de3a2e47"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.58.1"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d59e8320c2747553788e4fc42231489cc602fa50"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.58.1+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "44e3b40da000eab4ccb1aecdc4801c040026aeb5"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.13"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a4b12a1bd2ebade87891ab7e36fdbce582301a92"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.6"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[LibVPX_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "12ee7e23fa4d18361e7c2cde8f8337d4c3101bc7"
uuid = "dd192d2f-8180-539f-9fb4-cc70b1dcf69a"
version = "1.10.0+0"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "0fb723cd8c45858c22169b2e42269e53271a6df7"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.7"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "2ca267b08821e86c5ef4376cffed98a46c2cb205"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.1"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "438d35d2d95ae2c5e8780b330592b6de8494e779"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.3"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "9ff1c70190c1c30aebca35dc489f7411b256cd23"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.13"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "9e1a400fb1f27b4146fe35dc1a22de6c793b8f20"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.21.2"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "d4491becdc53580c6dadb0f6249f90caae888554"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.4.0"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3240808c6d463ac46f1c1cd7638375cd22abbccb"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.12"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8cbbc098554648c84f79a463c9ff0fd277144b6c"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.10"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "1700b86ad59348c0f9f68ddc95117071f947072d"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.1"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "d0c690d37c73aeb5ca063056283fde5585a41710"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.5.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "acc685bcf777b2202a904cdcb49ad34c2fa1880c"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.14.0+4"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7a5780a0d9c6864184b3a2eeeb833a0c871f00ab"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "0.1.6+4"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d713c1ce4deac133e3334ee12f4adff07f81778f"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2020.7.14+2"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "487da2f8f2f0c8ee0e83f39d13037d6bbf0a45ab"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.0.0+3"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─31bc44b6-ade1-4277-87dc-02379091d2c2
# ╟─b44f6459-6816-4747-9823-4a1575a77578
# ╟─439520d7-7453-43f6-92ef-91524e598d12
# ╠═449faf07-3214-47b6-bb02-824cf900bc07
# ╠═2dd2d4e3-14a0-43f1-b207-969b98ee5346
# ╟─64d42218-59a4-49b0-b260-f07ee0d0d5bd
# ╟─89f1748e-593f-473f-907d-b3a64073e197
# ╟─f878dcc5-06db-4dd4-8df8-c4666cf8f605
# ╟─f5e86f2a-076d-11ec-1080-c341c92d3fd1
# ╟─26f6893e-17ab-48e9-88d0-98a79ec0a1b0
# ╟─40b5474e-a8eb-4f6e-bcba-d49f714c83f2
# ╟─ecb886b8-cf77-4d82-865b-d2af63f22a2e
# ╟─42986972-b891-46ae-8462-93a31924e6fe
# ╟─2d13c6f7-4212-4bde-8e3d-1e46bc7992a6
# ╟─1dd0d637-4419-4475-9f36-4da7b8107e9a
# ╟─c357aaef-725a-40b2-aadf-2bd5a26ab997
# ╟─8e1fafdb-bdf0-4649-aff7-7720280f466a
# ╟─2e7cbd7c-c650-4088-ac92-3490f2db8e2c
# ╟─2b1e6d9b-abed-4633-a52e-a84514a6c27a
# ╟─07c5d7e1-574c-4f3b-b879-6b2842d862c1
# ╟─6e215ec5-dafb-4eb0-a159-b9028d6de08d
# ╟─c87acb83-84e6-484a-bed7-319bdb3d0feb
# ╠═263c3ed1-18b4-45eb-9de9-48b4f4f5c417
# ╟─a1f373d4-4c8c-448b-b336-82288f712f7c
# ╟─fcba2327-2306-4fc3-bb70-14e94e5d55c5
# ╠═809e639a-8492-4b0a-bc98-2226c2404478
# ╟─8d015906-914e-4ea5-935f-9cde818000af
# ╠═d6c69e40-e383-4588-8130-09ac2a2a69a9
# ╟─1fbc2b57-aa8a-47e8-8886-603e15474861
# ╟─6fd61abc-6bf3-4df7-b6bc-8adfafc461e0
# ╠═8d22438a-17ce-4bb6-814f-eecbd371fb71
# ╟─97216c01-ac6f-439e-8bf0-8ade70e08e03
# ╟─1edc180d-c2fe-4e30-8650-de2952b97a14
# ╟─002fab49-b0c0-4973-b676-473daaf2c698
# ╠═0656317b-3470-4d97-8b5d-9eddaab3ee8c
# ╠═dd4a6136-75a7-460d-9bbc-97a374e23a4f
# ╠═44762363-cc3e-422a-8044-670bc04b1690
# ╠═9059e6d7-ad27-4b06-b740-a91b9bc9f349
# ╟─f28a6182-2310-4444-a90e-402afd923bf9
# ╟─453ad11f-7bed-49f7-8255-6b4e7a86ae57
# ╠═cc9737b7-33ee-48c4-9cbd-dde3f9b28afc
# ╠═e59e1646-d949-456d-8e28-a4f9e1f55a8a
# ╟─421f19c1-38eb-4cff-b8d4-7a6ae4f438d0
# ╠═a87dd44d-89ca-41cc-b077-d444ea37c48d
# ╠═becfff3a-b130-44e8-8f9b-f63d693fdc0d
# ╠═20e80196-d302-44ad-af58-c2c6211bcc95
# ╟─202922e8-9d1a-4a23-afdc-f91deceb899e
# ╟─d98ad456-1f29-4047-8724-78a8cb35187a
# ╠═bafa445a-f3c8-4e83-bfa7-2b907bd33319
# ╟─61f10f8b-9f83-45ae-815f-22a86180dcfe
# ╟─ec533632-19dc-47e6-8491-3da3c2a3fe8a
# ╟─c416bf8c-8096-4ee8-bbfe-c5287f52d732
# ╟─c39566f3-f184-4cba-9dfb-390a9a5a8448
# ╠═d844d41b-35ee-4391-a36e-cffc5820cf30
# ╠═934a56c2-001f-4d7b-9f22-44b67c1b0481
# ╟─6a0e0537-db4a-46b2-b5a1-c9afbf98a028
# ╟─6a167ed2-45b4-4a09-a7e1-e24478c6ed38
# ╟─a48c3f1a-012c-4de3-b710-6c24bfb025ce
# ╟─02f6c327-048a-42ad-8d51-116fb2ba42fe
# ╠═858520cf-cc81-4fa4-afa1-b412bc43e062
# ╟─78a73608-de13-438c-822f-d2eb4b4a0a5b
# ╠═5785a365-738f-41f9-907a-133ca9cddc67
# ╟─07651aee-71aa-46c1-9e70-636db648b0fd
# ╠═75b4d476-9a25-41c3-97ad-c61b7bba7027
# ╠═aa27d3ae-12da-45f1-af6b-3fdfed216835
# ╟─6afbfaed-9e8e-4b20-96d1-15b3832ceb78
# ╠═a9bdbc6d-2cfd-47dd-840c-f495130601c1
# ╠═695669d8-f368-4156-84df-42fd40653758
# ╠═61b28060-5ee9-4699-a9d5-0628537f14d9
# ╠═5e79a3d1-d403-437c-8da3-2c0216c8432e
# ╠═0de24a1d-dbb1-4944-97ed-97cba4428a35
# ╟─9e59294b-40fe-4c21-8496-bd85680f15b6
# ╠═3861a74a-bb31-439d-a359-233333b5c536
# ╟─8d3c3463-c9bd-414d-a0ca-ef6d214155ed
# ╟─bf1a5cb5-4073-4d7a-989e-4b8f6bd97415
# ╟─780b02c4-5fc5-400b-9aac-09345c2769a0
# ╠═dbf96f31-84d4-45c0-9e37-17a277372de8
# ╠═629aae68-24f2-4cbd-99f3-1d964c723247
# ╟─968338b0-4536-45e5-ac12-573b939b9393
# ╟─acda2d19-e948-4b19-9da6-63b5a16c2e08
# ╟─ed702970-dde8-494f-8725-2b89bfb4d292
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
