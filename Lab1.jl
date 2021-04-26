### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 83728874-7aff-11eb-361e-27c8a01017c8
md"""
$$\def\com#1{\color{blue}{\textrm{#1}}}$$
"""

# ╔═╡ a1eabf29-f453-445d-9e15-d6809b570b36
md"""

# Lab1: Intro to Julia

In this lab you will gain familiarity with using Julia and the Pluto notebook environment. 
"""

# ╔═╡ 56a74762-59ea-4aa5-a361-64a40b91e3fe
md"""
# Submitting lab assignments
- __Rename your notebook to__ `Lab1_userID.ipynb` before submitting it to the link provided on Blackboard/Slack. For example, if your login name is `365joe`, make the notebook called `Lab1_joe.ipynb`. __If you do not use this naming convention, your lab may be lost__.
- Upload your completed lab to the [365 lab dropbox](https://www.dropbox.com/request/7DxEF6F3HYLl7KEyAYjt)
- __Due:__ 1pm Wednesday 10th March.

$\com{Grading comments will appear like this}$
---
"""

# ╔═╡ 26cc84f7-b12c-46a0-9352-c20c200eabdd
md"""
# Part 1: Julia and Pluto notebooks

First, we want you to do some basic calculations in Julia, and a get better feel for Pluto. _Note:_ you will need to load any packages you need from within this notebook.

- Go to [JuliaAcademy](https://juliaacademy.com/) and sign up for a free account. 
- Watch the course "Introduction to Julia (for programmers)". This will take about 1 hour in short video tutorials, stopping before the last 20 minute video (do this later  if you are interested). Take a minute in between each video to answer the __review questions__.  
- Watch the [Pluto JuliaCon talk](https://www.youtube.com/watch?v=IAF8DjrQSSk) to get a nice overview and a deeper idea of how Pluto works.
"""

# ╔═╡ 6c6c3a6e-7c97-11eb-1e5f-1152bf523a15
md"""
# Part 2: Julia basics
"""

# ╔═╡ 70f484a1-94b7-473b-8fb8-906f8a67ca88
md"""
1) Answer each question using text in a separate cell below.

```julia
a) v=1:4; u = collect(v) # how can we make a column vector [1 2 ... 1000000]^T ?
b) w = transpose(v); w = v' # how can we transpose a matrix?
c) x = rand(3,3)+im*rand(3,3); transpose(x); x' # what is the difference between transpose and '?
d) M = rand(5,5); A=M[2:4,1:3] # what does this command do? 
e) b = collect(1:12); B=reshape(b,3,4) # what does the second argument of reshape stand for? If you are unsure, use ?
```
"""

# ╔═╡ e29ce398-38ef-41e7-9474-e49f1bf7a7d1
md"""
2) Enter each command below in a separate cell, evaluate it, and then answer the questions below
```julia
A = [1 1; -2 6]
B = [3 -5;7 2]
C = [1 2 3; 2 4 6; 10 11 12]
1+A; 1.+A; 1 .+A
one(A)+A
A+B
A-B
A*B
A.*B
@. A*B
A/B
A./B
A^2
A.^2
A.^B
3 .^A; 3.0 .^A
A+C
```
"""

# ╔═╡ 8b303f0b-08b4-45fa-841e-aafe3b55bdb8
md"""
a) What does the second `2` in `2:2:6` mean?

b) What is the difference between `*` and `.*`? How about `^` and `.^`?

c) What does each element of `A.^B` give in `julia`?

d) Could you exponentiate `A` to `B`, i.e. `A^B`, or just element-wise? If "no", then what is the error message? If "yes", what is the result?
"""

# ╔═╡ 0be2d2fa-aae9-4d7e-8747-75e6526f09a5
md"""
# Part 3: Functions
"""

# ╔═╡ 413965f8-2cea-4a0d-885e-c0608a2b60fe
md"""
We have seen in lectures that a function can be defined as simply as `g(x)=x^2`. For more complicated functions, the general way to define a function in Julia is, e.g:

```julia 
function f(x)
    a = cos(x)
    b = 3
    return a,b
end
```
and in this case, would be called for `x=3` with the command `a,b=f(3)`, assigning two numbers to variables `a`, `b` _outside the function_ (but `a,b,x` inside the function are all _private_ to the function). 
"""

# ╔═╡ a789061c-7c98-11eb-053e-99f8f04bc41c
md"""
## To vectorize or not to vectorize?

A common strategy in technical computing is to __vectorize__ your code. This refers to the process of expressing your problem in terms of vectors or matrices and (hopefully) sparse matrices acting on them. This allows high level languages like Matlab to call efficient compiled C-libraries for linear algebra. This is usuall essential because `for` loops are slow in high-level languages. 

Let's find out how what julia has to say about this by investigating a simple example in detail. Pluto gives a simple time for each cell. Better benchmarking will require the package `BenchmarkTools`.
"""

# ╔═╡ ec5a80f4-7c98-11eb-3c75-d1e998ad408a
md"""
Let's say we want to calculate the sum of the reciprocals of the first ten thousand cubes, i.e. 

$$s = 1^{-3}+2^{-3}+3^{-3}+\dots +10000^{-3},$$

an approximation to $\zeta(3)$, where the [Riemann-zeta function](https://en.wikipedia.org/wiki/Riemann_zeta_function)

$$\boxed{\quad\zeta(z)=\sum_{k=1}^\infty\frac{1}{k^z}\quad}$$

is an important function in several branches of Physics ($\zeta(3)$ enters the formula for Bose-Einstein condensation at a critical temperature).
"""

# ╔═╡ 2588aa11-7065-4f78-9dff-8f4284c96c78
md"""
## $\zeta(3)$
Determine $s$ using `julia` in three different ways.

(a) Use a [for loop](https://docs.julialang.org/en/v1/manual/control-flow/#man-loops-1) in a cell (not a function). You will need to initialize the sum outside the loop (`s=0.0`), and then to modify `s` inside the loop as discussed in Lecture 2. Notice also that the statement 
>`s=0.0`
given directly to julia defines a global variable: __variables defined outside functions (i.e. in the REPL) are global__. This also means they can often lead to slower code (harder for the compiler to reason about). Hence an important lesson in julia is: 
>if speed matters use functions. 
This is also good practice because __it makes your code more reusable.__
Pay attention to the time taken to execute the code.

(b) Put your for loop inside a [function](https://docs.julialang.org/en/v1/manual/functions/#) `mysum1()`, and again not the execution time. Pay attention to whether you are seeing compile time, or pure run time: the first time you call a function, the compiler does some extra work writing a fast function. The true timing we are interested in is for any call _after_ the first one where it is defined (redefinitions count as definitions too!). We refer to the execution time without compile time as the __run time__.

(c) Use only a single-line __vectorised__ `julia` expression to calculate $s$ (no function wrapper). If you haven't done this before, the strategy is to define a vector `n=[1 2 3 ... 10000]` (you can do this with a range object `1:10000`, even though it looks a bit different to a vector, it still behaves like one in most circumstances). Then perform local (dot) operations to achieve the result. You should be able to get this down to one very short expression that easily fits on one line. Note the run time.

(d) Wrap your single-line expression in a function `mysum2()` and note the run time. 

(e) Use the package [BenchmarkTools](https://github.com/JuliaCI/BenchmarkTools.jl), and the macro `@btime` to compare more accurate timings for (b),(c),(d) (`@btime code` averages `code` over multiple calls, so the garbage collection is triggered. This is a more realistic time). 

- What can you say about whether to use functions or write code directly in the cell (or REPL)?
- What can you conclude about the speed of vectorizing (writing in terms of array operations), as compared with writing a `for` loop? Would you recommend vectorizing code in julia? Justify your answer.
"""

# ╔═╡ Cell order:
# ╟─83728874-7aff-11eb-361e-27c8a01017c8
# ╟─a1eabf29-f453-445d-9e15-d6809b570b36
# ╟─56a74762-59ea-4aa5-a361-64a40b91e3fe
# ╟─26cc84f7-b12c-46a0-9352-c20c200eabdd
# ╟─6c6c3a6e-7c97-11eb-1e5f-1152bf523a15
# ╟─70f484a1-94b7-473b-8fb8-906f8a67ca88
# ╟─e29ce398-38ef-41e7-9474-e49f1bf7a7d1
# ╟─8b303f0b-08b4-45fa-841e-aafe3b55bdb8
# ╟─0be2d2fa-aae9-4d7e-8747-75e6526f09a5
# ╟─413965f8-2cea-4a0d-885e-c0608a2b60fe
# ╟─a789061c-7c98-11eb-053e-99f8f04bc41c
# ╟─ec5a80f4-7c98-11eb-3c75-d1e998ad408a
# ╟─2588aa11-7065-4f78-9dff-8f4284c96c78
