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

# ╔═╡ d4c40347-5cb2-42fb-bf2b-66c3c057753e
begin
	using PlutoUI # for the @bind macro
	using Images # https://juliaimages.org/latest/pkgs/ for included packages
	using ImageShow # inline display of images
	using ImageView # interactive display of images
	using FileIO 
	using ImageIO
end

# ╔═╡ ac612840-8770-44d8-bbc1-1f842b6fef88
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
<em> Images as Arrays </em>
</p>

<p style="
font-size: 1.5rem;
text-align: center;
opacity: .8;
"><em>An introduction to both</em></p>
<div style="display: flex; justify-content: center;">
</div>
</div>

<style>
body {
overflow-x: hidden;
}
</style>"""

# ╔═╡ 9535ac29-e781-45de-92de-c29ddd50bd77
md"""
# Notebook packages
"""

# ╔═╡ 3afc2a0d-98fa-4110-869a-9a67e837f674
md"""
Using `Pluto`'s package management system for reproducible results.
"""

# ╔═╡ 8d027df7-3dc5-4345-a740-d6f4130ac090
PlutoUI.TableOfContents(aside=true)

# ╔═╡ 5f944d93-c0a1-4918-88ff-7d13c2194e49
md"""
# Abstract
"""

# ╔═╡ 3149a71a-ef3e-44d0-b783-628fcc20858d
md"""
## Science and math concepts

- What is an image?
- Some image processing
"""

# ╔═╡ b234f6a2-54b4-4282-bbff-f62e578edc2f
md"""
## Julia and numerical concepts

- Arrays
- Array indexing
- Julia `function`s
- Broadcasting
"""

# ╔═╡ 9c8b026e-f47f-4c37-8e3f-5ba506297cbc
md"""
# Images
"""

# ╔═╡ 90af1c4e-c709-49c7-bcbf-4150e6a01873
md"""
What is an image on a computer? It's a collection of many (sometimes millions) **pixels** ("picture elements"), in which each one has a single color. The pixels are arranged in a two-dimensional rectangular grid, and the color is specified by indicating how much red, green, blue is in the color. The collection of all these little pixels comprise the image you see.
"""


# ╔═╡ a39d3d02-8ae4-44d6-8915-6530c1153ac9
md"""
## Downloading and viewing an image
"""

# ╔═╡ 6158f257-9f6a-461e-b06d-d36544e99142
md"""
Let's download an image from the internet. Specify its URL and name and assign it to a variable for ease of reference. This type of variable is called a **string**.
"""

# ╔═╡ 491be4c4-4779-4e1d-a102-7d080b04216c
url = "https://raw.githubusercontent.com/jamesamiller/Computational-Thinking/main/Images/bear_small.png"

# ╔═╡ c664cb0b-0cf6-4f67-b00d-1b654b8086f8
md"""
We can now download the image to our local machine. We can also specify a path where it will be kept. If none is specified, then the image will go to some temporary folder (with possibly a strange name, which doesn't matter).
"""

# ╔═╡ 4afaec5d-fbaf-499b-ae58-8d815cdd62b9
bear_path = download(url,"/Users/miller/photo.png")

# ╔═╡ 796d7a6b-ddbe-4bff-ad36-1192dc437734
md"""
Using the package `FileIO` that we loaded at the start of the notebook, we can `load` the data and show the actual image. This is Bear, a 15 year old Shih Tzu, who looks remarkably spry here.
"""

# ╔═╡ 366ae591-39db-4590-b7fd-ee67d1a637f5
bear = load(bear_path)

# ╔═╡ 5035a8d0-8797-497d-b4ee-fe67b4c5b98d
md"""
> 💪 Quick exercise: Load an image on your computer (if running Pluto locally) or another one from the internet.
"""

# ╔═╡ 7503af9b-1a25-45d2-a787-a6c7fbb594fc
md"""
## Inspecting the image
"""

# ╔═╡ 98b2b5fd-c6ed-45f4-b0c8-d51885703fed
md"""
Perhaps the most basic piece of information is how big the image is. That is, how many pixels is it made of. The `size` command is good here. It returns a pair of numbers, the first is the height and the second is the width.
"""

# ╔═╡ 51ac7599-dba3-4444-9fc4-73cab88bca74
bear_size = size(bear)

# ╔═╡ fe1e939d-5219-472d-a5ab-dfa8c5bf0581
bear_height = bear_size[1]

# ╔═╡ 3c43e753-dff3-4396-91ed-946ee778ccf9
bear_width = bear_size[2]

# ╔═╡ ad4b40aa-305e-474a-8e6a-302e675a3a6b
md"""
Congratulations. You've now worked with an **array**. That's what the pair of numbers was, and the square brackets `[  ]` are how you pick out an element of that array. Here, the array was one-dimensional, or just a single row or column of elements.
"""

# ╔═╡ 8499b752-2a04-4039-801e-11233caf7159
md"""
Suppose now we want to examine a portion of the image in more detail. We'll need a way to specify that portion.

The image is a grid of pixels, and this grid can be thought of as a **two-dimensional array**. Now, in order to refer to a pixel, we need to give a pair of positive integers for the coordinates of that single pixel, one for its horizontal location and one for its vertial. Specifying coordinates like this is called **indexing**. Think of the index of a book, which tells you *on which page* an idea is discussed.

As with one-dimensional arrays, the indices are specified inside square brackets `[  ]`.

"""

# ╔═╡ 9edd8313-84fd-4c32-9a8d-e2c32c669d74
md"""
A natural question is where we start the counting for these indices. That is, where is the origin and does the first number refer to a height or width position on the image?

A nice way to see this is to use the interactive command `imshow` from the package `ImageView`, loaded at the start. This opens the image in a new window and you can mouse around to see the coordinates of the pointer.
"""

# ╔═╡ e345a926-eb92-4f3d-b74c-a0ac1d82b202
imshow(bear)

# ╔═╡ 91a57c6f-31c6-4fed-9d03-07c1ee01c093
md"""
> 🤔 Question: A pixel, or element of the array, is specified by the row and column position in the array. In the pair of numbers ``[x, y]``, does ``x`` specify the row or column? What about ``y``? Where does the counting start, at what corner of the image?
>
> Answer: ``x`` refers to the row of the array, of which there are 648. ``y`` refers to the column, of which there are 864. Counting starts in the upper left corner.
"""

# ╔═╡ 3f5128c9-da58-468b-9f4c-cdcdee0d1c60
md"""
Let's pick out a pixel from the toy. We can easily select one we want using the interactive view.
"""

# ╔═╡ b704ad67-9949-4e10-b4dc-9358ffca38be
pixel1 = bear[480,150]

# ╔═╡ 21fd8e6b-9dc0-45d8-b5c7-0e55dc0db5e3
md"""
### Pluto Interactivity with `Slider`
"""

# ╔═╡ 02b97d35-9b1f-4398-b8b6-3711322a2967
md"""
We can assign row and column positions to variables, and can even do pairs of assignments at once.
"""

# ╔═╡ b34898b9-f8aa-4793-83fd-0a8a8dab11c3
begin
	row_i, col_i = 1, 10
	bear[row_i, col_i]
end

# ╔═╡ 10594bf0-b10f-4c2c-9016-dcd1ddba272c
md"""
Suppose we want to vary these values? Doing it manually like this can be time consuming, but there's another way, afforded by the interactivity of Pluto notebooks.

It's to use a `Slider`, which is part of the `PlutoUI` package. We "bind" the slider value to a variable, which can be used later.

Let's show the pixel for the chosen row and column. Notice how the pixel changes as we change their values. This is a consequence of Pluto's **reactivity**: other cells that reference these variables will be immediately updated too. This is especially valuable for learning a new language.
"""
	

# ╔═╡ 82309cf6-7f11-41f2-bfed-b09b61c9fdac
@bind row_slide Slider(1:bear_height, show_value=true)

# ╔═╡ 3ada3cfe-ada4-4e74-84ba-c3052d4e316e
@bind col_slide Slider(1:bear_width, show_value=true)

# ╔═╡ bbc37da2-10a6-435e-8223-922326e034a7
bear[row_slide, col_slide]

# ╔═╡ 23c6acc1-cd56-4b5e-b3a0-1f14b0e3cefd
md"""
## Parts of an image
"""

# ╔═╡ 70441526-acd3-41a1-b114-64c5549eb9fa
md"""
We saw that we can use the row number and column number to index a _single pixel_ of an image. Next, we will use a *range of numbers* to index _multiple rows or columns_ at once, returning a **subarray**.

The `:` symbol in an expression like `x:y` means all the numbers between, and including, `x` and `y`. The `:` symbol alone means to take all allowed values.

Let's take a horizontal slice of the image between rows 200 and 300, two different ways.
"""

# ╔═╡ e2e915fd-6140-4969-9456-2d970e3525e4
bear[200:300, 1:bear_width]

# ╔═╡ f0182ff5-7b69-412d-8b07-60472cf876ed
bear[200:300, :]

# ╔═╡ 019cf2a3-7a0e-4c99-a369-70fa0f681b7f
md"""
Here's a single row of pixels across the image, and just Bear's head.
"""

# ╔═╡ 75e8aec3-27fa-4a5a-90e5-afe6782cdcad
bear[250,:]

# ╔═╡ f2c1336b-a404-45e3-bc32-1a011b6b3663
bear[150:350,300:550]

# ╔═╡ 38557f02-2084-45f8-b926-82d8526ff222
md"""
We can use the `Slider` functionality to make this less manual as well. Here, we'll use a variant called `RangeSlider`, and you can see immediately what it does. *Again, note how convenient Pluto's reactivity is!*
"""

# ╔═╡ 22c3abef-48cf-4bb5-bbda-8c8e31a99278
@bind row_range RangeSlider(1:bear_height)

# ╔═╡ 3431b58a-9ca1-464f-a5a7-05812bcba258
@bind col_range RangeSlider(1:bear_width)

# ╔═╡ 739089a8-5fbf-44b5-be7e-b75aed35e517
bear[row_range, col_range]

# ╔═╡ 47339aed-d351-4956-93d2-0085522b248c
md"""
## Modifying an image
"""

# ╔═╡ 519589ad-0ce3-4d45-815a-d79f5776bd38
md"""
Now that we have access to image data, we can start to *process* that data to extract information and/or modify it in some way.

We might want to detect what type of objects are in the image, say a tree or house or person. To achieve a high-level goal like this, we will need to perform mid-level operations, such as detecting edges that separate different objects based on their color. And, in turn, to carry that out we will need to do low-level operations like comparing colors of neighboring pixels and somehow deciding if they are different.

"""

# ╔═╡ 9751caf8-5b01-4b89-9862-a7496a3b9537
md"""
Let's look more at colors. We can  use indexing to modify a pixel's color. To do so, we need a way to specify a new color.

Color turns out to be a complicated concept, having to do with the interaction of the physical properties of light with the physiological mechanisms and mental processes by which we detect it!

We use a standard method of representing colors in the computer as an **RGB triple**, which is a set of three numbers (r, g, b), giving the amount of red, of green and of blue in a colour, respectively. These are numbers between 0 (none) and 1 (full). 
"""

# ╔═╡ 9b13d289-85bf-4c42-98a1-5a0ce4b6cf90
RGB(1.0,0.0,0.0)

# ╔═╡ 016d2145-db62-4f52-91b3-26580ffa7834
RGB(10.1,0.5,2.0)

# ╔═╡ 033033dc-5ca5-45d4-a6e0-cf702f19b599
md"""
> 💪 Exercise: Use Sliders to create a color.
"""

# ╔═╡ 047dbbf5-8cfc-4029-b734-0defec6f826d
md"""
> ✅ Here is one solution. The `0.1` is called the step size or increment. Allowed `Slider` values are thus 0, 0.1, 0.2, ..., 1.0.
"""

# ╔═╡ 2080e340-62d9-40b4-b1ce-6c919920bced
@bind red_amount Slider(0:0.1:1, show_value=true)

# ╔═╡ 29f93e8d-b8ea-4c8c-8e05-51567cac3f9b
@bind green_amount Slider(0:0.1:1, show_value=true)

# ╔═╡ 47125e13-c22b-4a7a-923f-b06d7cb2bfa3
@bind blue_amount Slider(0:0.1:1, show_value=true)

# ╔═╡ ca4f2907-af4e-4f69-8d35-13a7c79a6d68
RGB(red_amount, green_amount, blue_amount)

# ╔═╡ 95247f0e-85d0-4e1c-af5b-0576efa5829c
md"""
> 💪 Exercise: Let's *invert* a color. Inverting means ``(r, g, b) \rightarrow (1-r, 1-g, 1-b)``. Specifically, write a function to do this. You can call it whatever you want, but perhaps `invert` is a logical choice.
>
>We haven't talked much about `function`s yet, but you can guess what they are: ways of encapsulating some code, which eases its reuse later on for other things. A template follows.
>
> Hint: The way to extract the red, green, or blue color value from `color` is to the commands, `red`, `green`, or `blue`. It's nice when things are the obvious choice...

```julia
# Function template
function my_function_name(color::AbstractRGB)
	# some code
	return # the inverted color
end
```

"""

# ╔═╡ 10b6525c-179d-4b97-a72f-df32e662da99
md"""
> ✅ Answer: Here's an answer. Let's define a color and invert it to test. 
"""

# ╔═╡ ceabd170-9040-4149-a288-8b52e5942105
function invert(color::AbstractRGB)
	inverted_color = RGB(1 - red(color), 1 - green(color), 1 - blue(color))
	return inverted_color
end

# ╔═╡ 8a510354-320f-47c5-8984-55a0e38e732a
testc = RGB(0.2,0,0.7)

# ╔═╡ a033d602-61d7-4b7b-a429-f5ed9f4d8184
invert(testc)

# ╔═╡ bb398e63-1d34-4032-b1dc-7ea37d97a1dc
md"""
We can invert the color of an image pixel too.
"""

# ╔═╡ e3535710-963d-4bab-a27a-763636b0380c
bear[300,450]

# ╔═╡ 2b01bc69-a933-4af0-bb26-bf9db35417b4
invert(bear[300,450])

# ╔═╡ a8bcee58-fc68-4a60-b20b-608fba4b222f
md"""
> 💪 Exercise: This one is more challenging. The function we just defined works on a single pixel, but suppose we want to invert an entire image? Try your hand at defining a function for this.
"""

# ╔═╡ 2493f995-f59e-486a-8d54-dadeb7381129
# Function template
function my_name(image)
	# some code
	return # the inverted image
end

# ╔═╡ 76a74a55-c5ee-4b50-a5bf-038e18136a2e
md"""
> ✅ Answer: Here's an answer. And let's test it.
"""

# ╔═╡ cfc430f2-e5b3-4f3a-9493-c46fb5a3ad6a
function iminvert(image)
	rows = size(image)[1]  # number of rows, or the height in pixels
	cols = size(image)[2]  # number of columns, or the width in pixels
	inverted_image = similar(image)  # allocate the inverted image array
	for i = 1:rows  # loop over the rows
		for j = 1:cols  # loop over the columns for each row
			inverted_image[i,j] = invert(image[i,j])
		end
	end
	return inverted_image
end


# ╔═╡ a6310891-5c66-4a0a-88e3-91db94eae846
inverted_bear = iminvert(bear)

# ╔═╡ 5fff41c0-282d-4999-ada0-3a37e733f472
md"""
> ☡ Programming note: Notice that I used the function we already built `invert` in this function. This is good programming practice. We already have a tool that works, so it will reduce bugs if we simply reuse that.
"""

# ╔═╡ 152e4834-ae82-47a9-b1d1-2bc75bea8590
md"""
Now suppose we want to reduce the size of an image. One way is to **sample** the original image at every ``n`` pixels. This will surely reduce the size, but it will have the obvious effect on the resolution of the image.
"""

# ╔═╡ 3db7f041-ed09-4704-8a6f-79d3e9f259fc
smaller_bear = bear[1:8:end, 1:8:end]

# ╔═╡ 8ecb0f33-2909-409b-b9bb-b3794b069661
md"""
That looks pretty bad, as expected. (And notice the new Julia statement `end` here.)

We can do better by using the function `imresize` from the package `ImageTransformations` that is included in `Images`. This does a better job, and uses an algorithm called "SVD." Unless you dig into the details of how this works, `imresize` should be treated as a "black box." This means code or an algorithm that may be unknown or unclear to us. There is nothing wrong in using properly vetted ones, so be sure they are from a trusted, reputable source! 
"""

# ╔═╡ d1376074-ff05-49ed-80e7-7fe6fb7e792e
better_smaller_bear = imresize(bear, ratio=1/4)

# ╔═╡ 5c2fa757-efd8-4691-b12f-7b007192360c
size(better_smaller_bear)

# ╔═╡ 641840b3-ad20-4810-b0dd-bb64799e7ad1
md"""
Let's now save this image. Your save location will be different!
"""

# ╔═╡ dc7e92a9-1008-427e-81d0-1f126f3db49f
save("/Users/miller/better_smaller_bear.png", better_smaller_bear)

# ╔═╡ 32239004-b079-4dc3-9d94-8daa146887dd
md"""
# Arrays
"""

# ╔═╡ 0e59fd0d-44bb-408c-bd73-35bdbf161470
md"""
## Basics
"""

# ╔═╡ 58088872-d0f9-42b0-94db-3d1e749990b8
md"""
An image is a concrete example of a fundamental concept in computer science, namely an **array**. 

Just as an image is a rectangular grid, where each grid cell contains a single color,
an array is a rectangular grid for storing data. Data is stored and retrieved using indexing, just as in the image examples: each cell in the grid can store a single "piece of data" of a given type.

An array can be one-dimensional, like the strip of pixels above, two-dimensional, three-dimensional, and so on. The dimension tells us the number of indices that we need to specify a unique location in the grid. The array object also needs to know the length of the data in each dimension.

One-dimensional arrays are often called **vectors** and two-dimensional arrays are **matrices**. 

An array is an example of a **data structure**, i.e. a way of arranging data such that we can access it. A key theme in computer science is that of designing different data structures that represent data in different ways.

Conceptually, we can think of an array as a block of data that has a position or location in space. This can be a useful way to arrange data if, for example, we want to represent the fact that values in nearby locations in array are somehow near to one another.

Images are a good example of this: neighbouring pixels often represent different pieces of the same object, for example the rug or floor, or Philip himself, in the photo. We thus expect neighbouring pixels to be of a similar color. On the other hand, if they are not, this is also useful information, since that may correspond to the edge of an object.

"""

# ╔═╡ ad078b68-17fc-4b60-8941-833c5ff5e24d
md"""
## Constructing and manipulations
"""

# ╔═╡ 43c0494e-48a0-4af5-b2a9-4fd9aa7dee99
md"""
One-dimensional arrays, or vectors, are written using square brackets `[ ]`. In the example that follows, click on the sideways triangle before the vector. You'll see the results displayed as a column. Vectors in Julia are actually represented as columns, which is really only important when doing linear algebra.

Vectors (and also matrices) can be composed of objects of different **data types**: integers, floating point numbers, colors, strings, etc.
"""

# ╔═╡ 47fc9219-3231-407f-be8b-cc017733d3f1
vec1 = [1, 4.0, "this is a string"]  # can mix data types

# ╔═╡ dbc99eef-6491-469c-876c-b7472445e933
typeof(vec1)

# ╔═╡ 94a81a04-2d69-46ce-917a-2ba7113ac1d0
vec_pic = [RGB(1,0,0), RGB(0,1,0), RGB(0,0,1)]

# ╔═╡ fd49f09e-15eb-4dcd-b1c8-fcafe656351f
typeof(vec_pic)

# ╔═╡ 31bf42be-aab0-47ac-ac37-c6280bd4c895
md"""
Vectors are one dimensional.
"""

# ╔═╡ 6ae658de-66a4-451c-8dfa-5793986e34be
vec2 = [2.0 5.0 9.0]

# ╔═╡ ad18bf51-cee4-4f1c-99f8-8b89f47eb333
typeof(vec2)

# ╔═╡ 279663dc-abe6-4eca-903f-327c54702458
md"""
Notice what happened when we used spaces to separate elements instead of commas. Now the construction is a matrix and it's two dimensional: 1 row and 3 columns in this case.

We can generalize this using a ";" to make a matrix with more than one row.
"""

# ╔═╡ 26cd57ed-c96b-4a46-a49b-64d951915d67
mat1 = [RGB(1,0,0) RGB(0,1,0) RGB(0,0,1);  # notice the new line character ;
        RGB(0,0,1) RGB(0,1,0) RGB(1,0,0)]

# ╔═╡ 631e9e6d-dedc-4a9c-89a5-5bf449b42574
md"""
It's clear that if we want to create an array with more than a few elements, it will be very tedious to do so by hand like this.
Rather, we want to *automate* the process of creating an array by following some pattern, for example to create a whole palette of colors!

Let's start with all the possible colors interpolating between black, `RGB(0, 0, 0)`, and red, `RGB(1, 0, 0)`.  Since only one of the values is changing, we can represent this as a vector.

A Julia method to do this is an **array comprehension**. Again we use square brackets  to create an array, but now we use a variable that runs over a given range values.
"""

# ╔═╡ e3c22035-a453-4ddc-b107-a4f7e4e16671
[RGB(x, 0, 0) for x = 0:0.1:1]  # 0.1 is the step size or increment here

# ╔═╡ 37cc67d7-d868-4f74-8b46-0550328f918e
md"""
The syntax `0:0.1:1` is a **range**. The first and last numbers are the start and end values, and the middle number is the size of the step. The `=` sign after the `for` statement can be replaced by `in`, but I usually prefer the equal sign.
"""

# ╔═╡ c0a4eb82-7035-438e-b0ce-d6d81517e83a
md"""
Comprehensions generalize to matrices.
"""

# ╔═╡ 6dca0b3f-4f29-4a42-bc6f-81182bb4c9f5
[RGB(i, j, 0) for i = 0:0.1:1, j in 0:0.1:1]

# ╔═╡ 8395676a-4b83-492f-9acf-25518d14912e
md"""
Vectors and matrices can be joined in a very natural way:
"""

# ╔═╡ 85855eff-c04d-4ec0-9d03-e1ce5a307e43
[bear, inverted_bear]

# ╔═╡ 975d5350-5f2f-4943-8e1c-b587905d0679
[bear inverted_bear]

# ╔═╡ bdc46829-e183-41e7-ad86-cdb603d21267
[bear inverted_bear;
 inverted_bear bear]

# ╔═╡ b229214e-63e9-46e5-a4cd-f2f8c2b4f6d4
md"""
## Modifying an image, part 2
"""

# ╔═╡ 118deed1-4415-4747-94e0-521563a6beeb
md"""
Lastly, let's look at **redacting** an image, or erasing certain parts of it.

How about erasing Bear's eyes and replacing them with red? This is a **destructive** process: If we use the original image, it will be lost (unless it wasn't saved to disk yet). So let's copy the original into a temporary image and mess with the temporary image.
"""

# ╔═╡ e98a4172-5294-4d99-9736-6a6117c49ec0
begin
	temp = copy(bear)
	for row_index = 200:250  # loop over rows
		for col_index = 350:500  # loop over columns
			temp[row_index,col_index] = RGB(1.0,0.0,0.0)
		end
	end
	temp
end

# ╔═╡ 431122e9-1967-450e-bc53-7705e962006e
md"""
Or how about this?
"""

# ╔═╡ 3fc2afe3-aae5-4013-9b21-edc0de70a571
begin
	temp1 = copy(bear)
	temp1[200:250, 350:500] .= RGB(1.0, 0.0, 0.0)
	temp1
end

# ╔═╡ c3df7685-1792-48d1-b73e-3b1b9b7d6106
md"""
> ☡ Programming note: This is shorter, the thus less prone to errors in typing and syntax. It's an example of **broadcasting**. I "broadcast" the equal sign over a whole range of pixels. This is very convenient, allows for more concise code, and is quicker to implement. Try it. Retype these two approaches and see which one you want to use.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
ImageIO = "82e4d734-157c-48bb-816b-45c225c6df19"
ImageShow = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
ImageView = "86fae568-95e7-573e-a6b2-d8a6b900c9ef"
Images = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
FileIO = "~1.14.0"
ImageIO = "~0.6.2"
ImageShow = "~0.3.4"
ImageView = "~0.11.0"
Images = "~0.25.2"
PlutoUI = "~0.7.38"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0"
manifest_format = "2.0"

[[deps.ATK_jll]]
deps = ["Artifacts", "Glib_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "58c36d8a1beeb12d63921bcfaa674baf30a1140e"
uuid = "7b86fcea-f67b-53e1-809c-8f1719c154e8"
version = "2.36.1+0"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "cf6875678085aed97f52bfc493baaebeb6d40bcb"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.5"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9950387274246d08af38f6eef8cb5480862a435f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.14.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "1e315e3f4b0b7ce40feded39c73049692126cf53"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.3"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "75479b7df4167267d75294d14b58244695beb2ac"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.14.2"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "3f1f500312161f1ae067abe07d13b40f78f32e07"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.8"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "b153278a25dd42c65abbf4e62344f9d22e59191b"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.43.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "681ea870b918e7cff7111da58791d7f718067a19"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.2"

[[deps.CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "cc1a8e22627f33c789ab60b36a9132ac050bbf75"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.12"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "97f1325c10bd02b1cc1882e9c2bf6407ba630ace"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.12.16+3"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "cbdf14d1e8c7c8aacbe8b19862e0179fd08321c2"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "505876577b5481e50d089c1c68899dfb6faebc62"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.6"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "9267e5f50b0e12fdfd5a2455534345c4cf2c7f7a"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.14.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GTK3_jll]]
deps = ["ATK_jll", "Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Libepoxy_jll", "Pango_jll", "Pkg", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXcomposite_jll", "Xorg_libXcursor_jll", "Xorg_libXdamage_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "Xorg_libXrender_jll", "at_spi2_atk_jll", "gdk_pixbuf_jll", "iso_codes_jll", "xkbcommon_jll"]
git-tree-sha1 = "b080a592525632d287aee4637a62682576b7f5e4"
uuid = "77ec8976-b24b-556a-a1bf-49a033a670a6"
version = "3.24.31+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78e2c69783c9753a91cdae88a8d432be85a2ab5e"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "1c5a84319923bea76fa145d49e93aa4394c73fc2"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "57c021de207e234108a6f1454003120a1bf350c4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.6.0"

[[deps.Gtk]]
deps = ["Cairo", "Cairo_jll", "Dates", "GTK3_jll", "Glib_jll", "Graphics", "JLLWrappers", "Libdl", "Librsvg_jll", "Pkg", "Reexport", "Serialization", "Test", "Xorg_xkeyboard_config_jll", "adwaita_icon_theme_jll", "gdk_pixbuf_jll", "hicolor_icon_theme_jll"]
git-tree-sha1 = "eb7302ec0e1690787f396bc882a7c681d430bfdb"
uuid = "4c0ca9eb-093a-5379-98c5-f87ac0bbbf44"
version = "1.2.1"

[[deps.GtkObservables]]
deps = ["Cairo", "Colors", "Dates", "FixedPointNumbers", "Graphics", "Gtk", "IntervalSets", "LinearAlgebra", "Observables", "Reexport", "RoundingIntegers"]
git-tree-sha1 = "cf87f031fee932b90023ea37207c7a1de8caee6f"
uuid = "8710efd8-4ad6-11eb-33ea-2d5ceb25a41c"
version = "1.2.3"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "c54b581a83008dc7f292e205f4c409ab5caa0f04"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.10"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageContrastAdjustment]]
deps = ["ImageCore", "ImageTransformations", "Parameters"]
git-tree-sha1 = "0d75cafa80cf22026cea21a8e6cf965295003edc"
uuid = "f332f351-ec65-5f6a-b3d1-319c6670881a"
version = "0.3.10"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "9a5c62f231e5bba35695a20988fc7cd6de7eeb5a"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.3"

[[deps.ImageDistances]]
deps = ["Distances", "ImageCore", "ImageMorphology", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "7a20463713d239a19cbad3f6991e404aca876bda"
uuid = "51556ac3-7006-55f5-8cb3-34580c88182d"
version = "0.2.15"

[[deps.ImageFiltering]]
deps = ["CatIndices", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageBase", "ImageCore", "LinearAlgebra", "OffsetArrays", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "15bd05c1c0d5dbb32a9a3d7e0ad2d50dd6167189"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.7.1"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "539682309e12265fbe75de8d83560c307af975bd"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.2"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils", "Libdl", "Pkg", "Random"]
git-tree-sha1 = "5bc1cb62e0c5f1005868358db0692c994c3a13c6"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.2.1"

[[deps.ImageMagick_jll]]
deps = ["Artifacts", "Ghostscript_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f025b79883f361fa1bd80ad132773161d231fd9f"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.12+2"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "36cbaebed194b292590cba2593da27b34763804a"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.8"

[[deps.ImageMorphology]]
deps = ["ImageCore", "LinearAlgebra", "Requires", "TiledIteration"]
git-tree-sha1 = "e7c68ab3df4a75511ba33fc5d8d9098007b579a8"
uuid = "787d08f9-d448-5407-9aad-5290dd7ab264"
version = "0.3.2"

[[deps.ImageQualityIndexes]]
deps = ["ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "OffsetArrays", "Statistics"]
git-tree-sha1 = "1d2d73b14198d10f7f12bf7f8481fd4b3ff5cd61"
uuid = "2996bd0c-7a13-11e9-2da2-2f5ce47296a9"
version = "0.3.0"

[[deps.ImageSegmentation]]
deps = ["Clustering", "DataStructures", "Distances", "Graphs", "ImageCore", "ImageFiltering", "ImageMorphology", "LinearAlgebra", "MetaGraphs", "RegionTrees", "SimpleWeightedGraphs", "StaticArrays", "Statistics"]
git-tree-sha1 = "36832067ea220818d105d718527d6ed02385bf22"
uuid = "80713f31-8817-5129-9cf8-209ff8fb23e1"
version = "1.7.0"

[[deps.ImageShow]]
deps = ["Base64", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "25f7784b067f699ae4e4cb820465c174f7022972"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.4"

[[deps.ImageTransformations]]
deps = ["AxisAlgorithms", "ColorVectorSpace", "CoordinateTransformations", "ImageBase", "ImageCore", "Interpolations", "OffsetArrays", "Rotations", "StaticArrays"]
git-tree-sha1 = "42fe8de1fe1f80dab37a39d391b6301f7aeaa7b8"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.9.4"

[[deps.ImageView]]
deps = ["AxisArrays", "Cairo", "Compat", "Graphics", "Gtk", "GtkObservables", "ImageBase", "ImageCore", "ImageMetadata", "RoundingIntegers", "StatsBase"]
git-tree-sha1 = "d3aeb8f4c6b212cad334d2c650c3a88b836f2716"
uuid = "86fae568-95e7-573e-a6b2-d8a6b900c9ef"
version = "0.11.0"

[[deps.Images]]
deps = ["Base64", "FileIO", "Graphics", "ImageAxes", "ImageBase", "ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "ImageIO", "ImageMagick", "ImageMetadata", "ImageMorphology", "ImageQualityIndexes", "ImageSegmentation", "ImageShow", "ImageTransformations", "IndirectArrays", "IntegralArrays", "Random", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "TiledIteration"]
git-tree-sha1 = "03d1301b7ec885b266c0f816f338368c6c0b81bd"
uuid = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
version = "0.25.2"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.IntegralArrays]]
deps = ["ColorTypes", "FixedPointNumbers", "IntervalSets"]
git-tree-sha1 = "509075560b9fce23fdb3ccb4cc97935f11a43aa0"
uuid = "1d092043-8f09-5a30-832f-7509e371ab51"
version = "0.1.4"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b7bc05649af456efc75d178846f47006c2c4c3c7"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.6"

[[deps.IntervalSets]]
deps = ["Dates", "Statistics"]
git-tree-sha1 = "eb381d885e30ef859068fce929371a8a5d06a914"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.6.1"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "336cc738f03e069ef2cac55a104eb823455dca75"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.4"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "81b9477b49402b47fbe7f7ae0b252077f53e4a08"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.22"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "a77b273f1ddec645d1b7c4fd5fb98c8f90ad10a5"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libepoxy_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "18b65a0eff6b58546bec18065e73f8a04e83758d"
uuid = "42c93a91-0102-5b3f-8f9d-e41de60ac950"
version = "1.5.8+1"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Librsvg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pango_jll", "Pkg", "gdk_pixbuf_jll"]
git-tree-sha1 = "25d5e6b4eb3558613ace1c67d6a871420bfca527"
uuid = "925c91fb-5dd6-59dd-8e8c-345e74382d89"
version = "2.52.4+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "76c987446e8d555677f064aaac1145c4c17662f8"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.14"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "e595b205efd49508358f7dc670a940c790204629"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.0.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.MetaGraphs]]
deps = ["Graphs", "JLD2", "Random"]
git-tree-sha1 = "2af69ff3c024d13bde52b34a2a7d6887d4e7b438"
uuid = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
version = "0.7.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "b34e3bc3ca7c94914418637cb10cc4d1d80d877d"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.3"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "ded92de95031d4a8c61dfb6ba9adb6f1d8016ddd"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.10"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore"]
git-tree-sha1 = "18efc06f6ec36a8b801b23f076e3c6ac7c3bf153"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Observables]]
git-tree-sha1 = "dfd8d34871bc3ad08cd16026c1828e271d554db9"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.1"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "aee446d0b3d5764e35289762f6a18e8ea041a592"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.11.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "e925a64b8585aa9f4e3047b8d2cdc3f0e79fd4e4"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.16"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a121dfbba67c94a5bec9dde613c3d0cbcf3a12b"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.3+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "1285416549ccfcdf0c50d4997a94331e88d68413"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "670e559e5c8e191ded66fa9ea89c97f10376bb4c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.38"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.Quaternions]]
deps = ["DualNumbers", "LinearAlgebra", "Random"]
git-tree-sha1 = "b327e4db3f2202a4efafe7569fcbe409106a1f75"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.5.6"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegionTrees]]
deps = ["IterTools", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "4618ed0da7a251c7f92e869ae1a19c74a7d2a7f9"
uuid = "dee08c22-ab7f-5625-9660-a9af2021b33f"
version = "0.3.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "3177100077c68060d63dd71aec209373c3ec339b"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.3.1"

[[deps.RoundingIntegers]]
git-tree-sha1 = "99acd97f396ea71a5be06ba6de5c9defe188a778"
uuid = "d5f540fe-1c90-5db3-b776-2e2f362d9394"
version = "1.1.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays", "Test"]
git-tree-sha1 = "a6f404cc44d3d3b28c793ec0eb59af709d827e4e"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.2.1"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "cd56bf18ed715e8b09f06ef8c6b781e6cdc49911"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c82aaa13b44ea00134f8c9c89819477bd3986ecd"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.3.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "f90022b44b7bf97952756a6b6737d1a0024a3233"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.5.5"

[[deps.TiledIteration]]
deps = ["OffsetArrays"]
git-tree-sha1 = "5683455224ba92ef59db72d10690690f4a8dc297"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.3.1"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcomposite_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll"]
git-tree-sha1 = "7c688ca9c957837539bbe1c53629bb871025e423"
uuid = "3c9796d7-64a0-5134-86ad-79f8eb684845"
version = "0.4.5+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdamage_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll"]
git-tree-sha1 = "fe4ffb2024ba3eddc862c6e1d70e2b070cd1c2bf"
uuid = "0aeada51-83db-5f97-b67e-184615cfc6f6"
version = "1.1.5+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libXtst_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "Xorg_libXi_jll"]
git-tree-sha1 = "0c0a60851f44add2a64069ddf213e941c30ed93c"
uuid = "b6f176f1-7aea-5357-ad67-1d3e565ea1c6"
version = "1.2.3+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.adwaita_icon_theme_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "hicolor_icon_theme_jll"]
git-tree-sha1 = "37c9a36ccb876e02876c8a654f1b2e8c1b443a78"
uuid = "b437f822-2cd6-5e08-a15c-8bac984d38ee"
version = "3.33.92+5"

[[deps.at_spi2_atk_jll]]
deps = ["ATK_jll", "Artifacts", "JLLWrappers", "Libdl", "Pkg", "XML2_jll", "Xorg_libX11_jll", "at_spi2_core_jll"]
git-tree-sha1 = "f16ae690aca4761f33d2cb338ee9899e541f5eae"
uuid = "de012916-1e3f-58c2-8f29-df3ef51d412d"
version = "2.34.1+4"

[[deps.at_spi2_core_jll]]
deps = ["Artifacts", "Dbus_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXtst_jll"]
git-tree-sha1 = "d2d540cd145f2b2933614649c029d222fe125188"
uuid = "0fc3237b-ac94-5853-b45c-d43d59a06200"
version = "2.34.0+4"

[[deps.gdk_pixbuf_jll]]
deps = ["Artifacts", "Glib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Xorg_libX11_jll", "libpng_jll"]
git-tree-sha1 = "c23323cd30d60941f8c68419a70905d9bdd92808"
uuid = "da03df04-f53b-5353-a52f-6a8b0620ced0"
version = "2.42.6+1"

[[deps.hicolor_icon_theme_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b458a6f6fc2b1a8ca74ed63852e4eaf43fb9f5ea"
uuid = "059c91fe-1bad-52ad-bddd-f7b78713c282"
version = "0.17.0+3"

[[deps.iso_codes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "5ee24c3ae30e006117ec2da5ea50f2ce457c019a"
uuid = "bf975903-5238-5d20-8243-bc370bc1e7e5"
version = "4.3.0+4"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78736dab31ae7a53540a6b752efc61f77b304c5b"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.8.6+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─ac612840-8770-44d8-bbc1-1f842b6fef88
# ╟─9535ac29-e781-45de-92de-c29ddd50bd77
# ╟─3afc2a0d-98fa-4110-869a-9a67e837f674
# ╠═d4c40347-5cb2-42fb-bf2b-66c3c057753e
# ╠═8d027df7-3dc5-4345-a740-d6f4130ac090
# ╟─5f944d93-c0a1-4918-88ff-7d13c2194e49
# ╟─3149a71a-ef3e-44d0-b783-628fcc20858d
# ╟─b234f6a2-54b4-4282-bbff-f62e578edc2f
# ╟─9c8b026e-f47f-4c37-8e3f-5ba506297cbc
# ╟─90af1c4e-c709-49c7-bcbf-4150e6a01873
# ╟─a39d3d02-8ae4-44d6-8915-6530c1153ac9
# ╟─6158f257-9f6a-461e-b06d-d36544e99142
# ╠═491be4c4-4779-4e1d-a102-7d080b04216c
# ╟─c664cb0b-0cf6-4f67-b00d-1b654b8086f8
# ╠═4afaec5d-fbaf-499b-ae58-8d815cdd62b9
# ╟─796d7a6b-ddbe-4bff-ad36-1192dc437734
# ╠═366ae591-39db-4590-b7fd-ee67d1a637f5
# ╟─5035a8d0-8797-497d-b4ee-fe67b4c5b98d
# ╟─7503af9b-1a25-45d2-a787-a6c7fbb594fc
# ╟─98b2b5fd-c6ed-45f4-b0c8-d51885703fed
# ╠═51ac7599-dba3-4444-9fc4-73cab88bca74
# ╠═fe1e939d-5219-472d-a5ab-dfa8c5bf0581
# ╠═3c43e753-dff3-4396-91ed-946ee778ccf9
# ╟─ad4b40aa-305e-474a-8e6a-302e675a3a6b
# ╟─8499b752-2a04-4039-801e-11233caf7159
# ╟─9edd8313-84fd-4c32-9a8d-e2c32c669d74
# ╠═e345a926-eb92-4f3d-b74c-a0ac1d82b202
# ╟─91a57c6f-31c6-4fed-9d03-07c1ee01c093
# ╟─3f5128c9-da58-468b-9f4c-cdcdee0d1c60
# ╠═b704ad67-9949-4e10-b4dc-9358ffca38be
# ╟─21fd8e6b-9dc0-45d8-b5c7-0e55dc0db5e3
# ╟─02b97d35-9b1f-4398-b8b6-3711322a2967
# ╠═b34898b9-f8aa-4793-83fd-0a8a8dab11c3
# ╟─10594bf0-b10f-4c2c-9016-dcd1ddba272c
# ╠═82309cf6-7f11-41f2-bfed-b09b61c9fdac
# ╠═3ada3cfe-ada4-4e74-84ba-c3052d4e316e
# ╠═bbc37da2-10a6-435e-8223-922326e034a7
# ╟─23c6acc1-cd56-4b5e-b3a0-1f14b0e3cefd
# ╟─70441526-acd3-41a1-b114-64c5549eb9fa
# ╠═e2e915fd-6140-4969-9456-2d970e3525e4
# ╠═f0182ff5-7b69-412d-8b07-60472cf876ed
# ╟─019cf2a3-7a0e-4c99-a369-70fa0f681b7f
# ╠═75e8aec3-27fa-4a5a-90e5-afe6782cdcad
# ╠═f2c1336b-a404-45e3-bc32-1a011b6b3663
# ╟─38557f02-2084-45f8-b926-82d8526ff222
# ╠═22c3abef-48cf-4bb5-bbda-8c8e31a99278
# ╠═3431b58a-9ca1-464f-a5a7-05812bcba258
# ╠═739089a8-5fbf-44b5-be7e-b75aed35e517
# ╟─47339aed-d351-4956-93d2-0085522b248c
# ╟─519589ad-0ce3-4d45-815a-d79f5776bd38
# ╟─9751caf8-5b01-4b89-9862-a7496a3b9537
# ╠═9b13d289-85bf-4c42-98a1-5a0ce4b6cf90
# ╠═016d2145-db62-4f52-91b3-26580ffa7834
# ╟─033033dc-5ca5-45d4-a6e0-cf702f19b599
# ╟─047dbbf5-8cfc-4029-b734-0defec6f826d
# ╠═2080e340-62d9-40b4-b1ce-6c919920bced
# ╠═29f93e8d-b8ea-4c8c-8e05-51567cac3f9b
# ╠═47125e13-c22b-4a7a-923f-b06d7cb2bfa3
# ╠═ca4f2907-af4e-4f69-8d35-13a7c79a6d68
# ╟─95247f0e-85d0-4e1c-af5b-0576efa5829c
# ╟─10b6525c-179d-4b97-a72f-df32e662da99
# ╠═ceabd170-9040-4149-a288-8b52e5942105
# ╠═8a510354-320f-47c5-8984-55a0e38e732a
# ╠═a033d602-61d7-4b7b-a429-f5ed9f4d8184
# ╟─bb398e63-1d34-4032-b1dc-7ea37d97a1dc
# ╠═e3535710-963d-4bab-a27a-763636b0380c
# ╠═2b01bc69-a933-4af0-bb26-bf9db35417b4
# ╟─a8bcee58-fc68-4a60-b20b-608fba4b222f
# ╠═2493f995-f59e-486a-8d54-dadeb7381129
# ╟─76a74a55-c5ee-4b50-a5bf-038e18136a2e
# ╠═cfc430f2-e5b3-4f3a-9493-c46fb5a3ad6a
# ╠═a6310891-5c66-4a0a-88e3-91db94eae846
# ╟─5fff41c0-282d-4999-ada0-3a37e733f472
# ╟─152e4834-ae82-47a9-b1d1-2bc75bea8590
# ╠═3db7f041-ed09-4704-8a6f-79d3e9f259fc
# ╟─8ecb0f33-2909-409b-b9bb-b3794b069661
# ╠═d1376074-ff05-49ed-80e7-7fe6fb7e792e
# ╠═5c2fa757-efd8-4691-b12f-7b007192360c
# ╟─641840b3-ad20-4810-b0dd-bb64799e7ad1
# ╠═dc7e92a9-1008-427e-81d0-1f126f3db49f
# ╟─32239004-b079-4dc3-9d94-8daa146887dd
# ╟─0e59fd0d-44bb-408c-bd73-35bdbf161470
# ╟─58088872-d0f9-42b0-94db-3d1e749990b8
# ╟─ad078b68-17fc-4b60-8941-833c5ff5e24d
# ╟─43c0494e-48a0-4af5-b2a9-4fd9aa7dee99
# ╠═47fc9219-3231-407f-be8b-cc017733d3f1
# ╠═dbc99eef-6491-469c-876c-b7472445e933
# ╠═94a81a04-2d69-46ce-917a-2ba7113ac1d0
# ╠═fd49f09e-15eb-4dcd-b1c8-fcafe656351f
# ╟─31bf42be-aab0-47ac-ac37-c6280bd4c895
# ╠═6ae658de-66a4-451c-8dfa-5793986e34be
# ╠═ad18bf51-cee4-4f1c-99f8-8b89f47eb333
# ╟─279663dc-abe6-4eca-903f-327c54702458
# ╠═26cd57ed-c96b-4a46-a49b-64d951915d67
# ╟─631e9e6d-dedc-4a9c-89a5-5bf449b42574
# ╠═e3c22035-a453-4ddc-b107-a4f7e4e16671
# ╟─37cc67d7-d868-4f74-8b46-0550328f918e
# ╟─c0a4eb82-7035-438e-b0ce-d6d81517e83a
# ╠═6dca0b3f-4f29-4a42-bc6f-81182bb4c9f5
# ╟─8395676a-4b83-492f-9acf-25518d14912e
# ╠═85855eff-c04d-4ec0-9d03-e1ce5a307e43
# ╠═975d5350-5f2f-4943-8e1c-b587905d0679
# ╠═bdc46829-e183-41e7-ad86-cdb603d21267
# ╟─b229214e-63e9-46e5-a4cd-f2f8c2b4f6d4
# ╟─118deed1-4415-4747-94e0-521563a6beeb
# ╠═e98a4172-5294-4d99-9736-6a6117c49ec0
# ╟─431122e9-1967-450e-bc53-7705e962006e
# ╠═3fc2afe3-aae5-4013-9b21-edc0de70a571
# ╟─c3df7685-1792-48d1-b73e-3b1b9b7d6106
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
