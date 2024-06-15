# Raster analysis in Julia in 45 minutes

This tutorial is based on the xarray in 45 minutes 
https://tutorial.xarray.dev/overview/xarray-in-45-min.html

In this lesson, we discuss cover the basics of Raster analysis data structures. By the
end of the lesson, we will be able to:

- Understand the basic data structures in Julia
- Inspect `DimArray` and `DimStack` objects.
- Read and write netCDF files using Rasters.jl.
- Understand that there are many packages that build on top of xarray


We'll start by reviewing the various components of the Xarray data model, represented here visually:

<img src="https://docs.xarray.dev/en/stable/_images/dataset-diagram.png" align="center" width="50%">

As example data, we use the data in the xarray-data repository.

Here we'll use `air temperature` from the [National Center for Environmental Prediction](https://www.weather.gov/ncep/). Xarray objects have convenient HTML representations to give an overview of what we're working with:

Before we can start we need to activate the current environment.

```julia
using Pkg
Pkg.activate(@__DIR__)
```

```@example raster
using Rasters
using NCDatasets

path = "data/air_temperature.nc"
# First we download the data locally if needed.
if !isfile(path)
    download("https://github.com/pydata/xarray-data/blob/master/air_temperature.nc", "air_temperature.nc")
end

# Now we can open the data as a RasterStack
ds = RasterStack(path)
```



## What's in a DimStack? 

*Many DimArrays!* 

DimStacks are dictionary-like containers of "DimArray"s. They are a mapping of
variable name to DimArray.
The RasterStack is a special case of a DimStack with some geospatial information.
The DimArray s in the DimStack can share dimensions butdon't have all the same dimensionality.
If layers share a dimension name, this dimension is only stored once for the whole DimStack.


```@example raster
# pull out "air" dataarray with dictionary syntax
ds["air"]
```

You can save some typing by using the "attribute" or "dot" notation and using tab completion. 



```@example raster
# pull out dataarray using dot notation
ds.air
```

## What's in a DimArray? 

*data + (a lot of) metadata*



### Name (optional)


```@example raster
da = ds.air

Rasters.name(da)
```

### Named dimensions 

`dims(da)` correspond to the axes of your data. 

In this case we have 2 spatial dimensions (`X` and `Y`) and one temporal dimension (`Ti`).


```@example raster
dims(da)
```

You can also extract a single dimension.


```@example raster
# extracting coordinate variables
dims(da, X)
```

It is useful to think of the values in the dimensions as axis
"labels" such as "tick labels" in a figure. These are coordinate locations on a
grid at which you have data.


### Arbitrary attributes 

`metadata` is a metadata object that can hold arbitrary attributes that describe the underlying data. 


```@example raster
metadata(da)
```


```@example raster
# assign your own attributes!
metadata(da)["myattrs"] = "Mine"
```

### Underlying data 


A DimensionalData data structures wrap underlying simpler array-like data structures. These arrays have to fit into the Julia Array interface but can be either in memory arrays or DiskArray.jl arrays for lazy access or  part of Xarray is quite extensible allowing for distributed array, GPU arrays, sparse arrays, arrays with units etc. We'll  briefly look at this later in this tutorial.

To access the underlying data use the `parent` function:


```@example raster
parent(da)
```


```@example raster
# what is the type of the underlying data
typeof(parent(da))
```

We can change the underlying data type by using the lazy keyword for opening the data.
This is especially helpful for very large data or data that is hosted online where we would not want to download the whole dataset before starting the analysis.



```@example raster
dsl = RasterStack(path, lazy=true)
```


```@example raster
dal = dsl.air
```


```@example raster
typeof(parent(dal))
```

### Review

DimensionalData provides two main data structures:

1. `DimArrays` that wrap underlying data containers and contain associated metadata in the 
1. `DimStacks` that are dictionary-like containers of DataArrays

`DimArrays` contain underlying arrays and associated metadata:
1. Name
2. Dimension names
3. Lookup values
4. Metadata

---

## Why named dimensions and labeled dimensions? 

Metadata provides context and provides code that is more legible. This reduces the likelihood of errors from typos and makes analysis more intuitive and fun!

### Analysis without named dimensions:


```@example raster
# plot the first timestep
lon = ds.air.dims[1].val.data  # Vector
lat = ds.air.dims[2].val.data  # Vector
temp = parent(da)  # vector
```


```@example raster
using GLMakie
heatmap(lon, lat, temp[1, :, :])
```


```@example raster
using Statistics
mean(temp, dims=3)# On what dimensions did we apply the reduction? I can't tell by looking at this line.
```

### Analysis with DimensionalData

How readable is this code?



```@example raster
plot(ds.air[Ti=1])
```

Use dimension names instead of axis numbers



```@example raster
plot((mean(ds.air, dims=Ti)[Ti=1]))
```

---

## Extracting data or "indexing" 

DimensionalData supports

- label-based indexing using `Selector`s
- position-based indexing using `Integer`

See the [Documentation about Selectors](https://rafaqz.github.io/DimensionalData.jl/v0.27.0/selectors) for more.

### Label-based indexing

DimensionalData implements label based indexing where you can use the name of the dimension and also the labels for the entries in the dimension.



```@example raster
# here's what our dataset looks like
ds
```


```@example raster
# We can extract the Time dimension
dims(ds, Ti)
```


```@example raster
# pull out data for all of 2013-May
ds[Ti=Where(x->yearmonth(x) == (2013, 5))]
```


```@example raster
# demonstrate slicing, extract all time slices between to given dates
ds[Ti=Date(2013,5,1)..Date(2013,8,1)]
```


```@example raster
# demonstrate "nearest" indexing
ds[X=Near(240.2)]
```


```@example raster
# "nearest indexing at multiple points"
ds[X=Near([240.125, 234]), Y=Near([40.1, 50.1])]
```

These selectors can be mixed for different dimensions. So that we could have a `Where` selection for time and a nearest neighbor selection in space.

### Position-based indexing

This is similar to usual array indexing `array[1, 2, 3]` but with the power of named
dimensions!



```@example raster
# pull out time index 0, lat index 2, and lon index 3
ds.air[Ti=1, Y=2, X=3]  #  much better than ds.air[3, 2, 1]
```


```@example raster
# demonstrate slicing
ds.air[X=1:10]
```

---

## Concepts for computation

Consider calculating the *mean air temperature per unit surface area* for this dataset. Because latitude and longitude correspond to spherical coordinates for Earth's surface, each 2.5x2.5 degree grid cell actually has a different surface area as you move away from the equator! This is because *latitudinal length* is fixed ($ \delta Lat = R \delta \phi  $), but *longitudinal length varies with latitude* ($ \delta Lon = R \delta \lambda \cos(\phi) $)

So the [area element for lat-lon coordinates](https://en.wikipedia.org/wiki/Spherical_coordinate_system#Integration_and_differentiation_in_spherical_coordinates) is


$$ \delta A = R^2 \delta\phi \, \delta\lambda \cos(\phi) $$

where $\phi$ is latitude, $\delta \phi$ is the spacing of the points in latitude, $\delta \lambda$ is the spacing of the points in longitude, and $R$ is Earth's radius. (In this formula, $\phi$ and $\lambda$ are measured in radians)


```@example raster
lon
```


```@example raster
# Earth's average radius in meters
R = 6.371e6

# Coordinate spacing for this dataset is 2.5 x 2.5 degrees
dϕ = deg2rad(2.5)
dλ = deg2rad(2.5)

dlat = fill(R * dϕ, dims(ds,X))
```


```@example raster
dlonval = R .* dλ .* cos.(deg2rad.(dims(ds.air, Y)))
dlon = DimArray(reshape(dlonval, (1, length(dlonval))), (X,dims(ds, Y)))
```


```@example raster
cell_area = dlat .* dlon
```

You can apply functions like `cos` and `deg2rad` elementwise on array types by using broadcasting in Julia.

### Broadcasting: expanding data

In DimensionalData the broadcast does not automatically expand the data therefore we have to reshape the underlying data into a row vector to get a 2D array from the two 1D vectors.


```@example raster
cell_area = dlon .* dlat
cell_area
```

---

### Alignment: putting data on the same grid

When broadcasting arithmetic operations DimensionalData automatically "aligns" i.e. puts the
data on the same grid. In this case `cell_area` and `ds.air` are at the same
lat, lon points we end up with a result with the same shape (25x53):



```@example raster
ds.air[Ti=1] ./ cell_area
```

DimensionalData only compares the dimension name and would currently broadcast two arrays with the same size along each dimension together.


```@example raster
# make a copy of cell_area
# then add 1e-5 degrees to latitude
cell_area_bad = similar(cell_area)
```


```@example raster
set(cell_area_bad, set(dims(cell_area_bad, X), 1:53))
```


```@example raster
cell_area_bad .* ds.air[Ti=1]
```

This results in the same values but the dimensions are different. 

---

## High level computation 

(`groupby`, `resample`, `rolling`, `coarsen`, `weighted`)

Xarray has some very useful high level objects that let you do common
computations:

1. `groupby` :
   [Bin data in to groups and reduce](https://docs.xarray.dev/en/stable/groupby.html)
1. `resample` :
   [Groupby specialized for time axes. Either downsample or upsample your data.](https://docs.xarray.dev/en/stable/user-guide/time-series.html#resampling-and-grouped-operations)
1. `rolling` :
   [Operate on rolling windows of your data e.g. running mean](https://docs.xarray.dev/en/stable/user-guide/computation.html#rolling-window-operations)
1. `coarsen` :
   [Downsample your data](https://docs.xarray.dev/en/stable/user-guide/computation.html#coarsen-large-arrays)
1. `weighted` :
   [Weight your data before reducing](https://docs.xarray.dev/en/stable/user-guide/computation.html#weighted-array-reductions)


Below we quickly demonstrate these patterns. See the user guide links above and [the tutorial](https://tutorial.xarray.dev/intermediate/01-high-level-computation-patterns.html) for more.

### groupby



```@example raster
# here's ds
ds
```


```@example raster
groups = groupby(ds, Ti=>seasons())
```


```@example raster
# make a seasonal mean
seasonal_mean = mean.(groups, dims=Ti)
seasonal_mean
```


```@example raster
#TODO: Plot the mean map for every season
```

### resample



```@example raster
# Reduce the time dimension to monthly means 
monthlymeans = dropdims.(mean.(groupby(ds.air, Ti=>yearmonth), dims=Ti), dims=Ti)
```

### weighted



```@example raster
# weight by cell_area and take mean over (time, lon)
#ds.weighted(cell_area).mean(["lon", "time"]).air.plot(y="lat");
weightedmean = dropdims(mean(ds.air, dims=(X, Ti)), dims=(X,Ti))
weightedmean
```


```@example raster
y = dims(weightedmean, Y)
@show y
```


```@example raster
nx = X(reverse(Float32[75.0, 72.5, 70.0, 67.5, 65.0, 62.5, 60.0, 57.5, 55.0, 52.5, 50.0, 47.5, 45.0, 42.5, 40.0, 37.5, 35.0, 32.5, 30.0, 27.5, 25.0, 22.5, 20.0, 17.5, 15.0]))
ny = Y(reverse(Float32[75.0, 72.5, 70.0, 67.5, 65.0, 62.5, 60.0, 57.5, 55.0, 52.5, 50.0, 47.5, 45.0, 42.5, 40.0, 37.5, 35.0, 32.5, 30.0, 27.5, 25.0, 22.5, 20.0, 17.5, 15.0]))

narr = DimArray(rand(25,25), (ny, nx))
fig, ax, pl = plot(narr)
```


```@example raster
nsingle = DimArray(rand(25), ny)
```


```@example raster
fig, ax, pl = plot(nsingle)
```


```@example raster
ax.finallimits
```

---

## Visualization

(`.plot`)


We have seen very simple plots earlier. Xarray also lets you easily visualize
3D and 4D datasets by presenting multiple facets (or panels or subplots) showing
variations across rows and/or columns.


```@example raster
# facet the seasonal_mean
seasonal_mean.air.plot(col="season", col_wrap=2);
```


```@example raster
# contours
seasonal_mean.air.plot.contour(col="season", levels=20, add_colorbar=True);
```


```@example raster
# line plots too? wut
seasonal_mean.air.mean("lon").plot.line(hue="season", y="lat");
```

For more see the [user guide](https://docs.xarray.dev/en/stable/plotting.html), the [gallery](https://docs.xarray.dev/en/stable/examples/visualization_gallery.html), and [the tutorial material](https://tutorial.xarray.dev/fundamentals/04.0_plotting.html).

---

## Reading and writing files

Xarray supports many disk formats. Below is a small example using netCDF. For
more see the [documentation](https://docs.xarray.dev/en/stable/user-guide/io.html)



```@example raster
# write to netCDF
ds.to_netcdf("my-example-dataset.nc")
```

!!! note
    To avoid the `SerializationWarning` you can assign a _FillValue for any NaNs in 'air' array by adding the keyword argument encoding=dict(air={_FillValue=-9999})


```@example raster
# read from disk
fromdisk = xr.open_dataset("my-example-dataset.nc")
fromdisk
```


```@example raster
# check that the two are identical
ds.identical(fromdisk)
```

!!! tip
    A common use case to read datasets that are a collection of many netCDF
    files. See the [documentation](https://docs.xarray.dev/en/stable/user-guide/io.html#reading-multi-file-datasets) for how
    to handle that.


Finally to read other file formats, you might find yourself reading in the data using a different library and then creating a DataArray([docs](https://docs.xarray.dev/en/stable/user-guide/data-structures.html#creating-a-dataarray), [tutorial](https://tutorial.xarray.dev/fundamentals/01.1_creating_data_structures.html)) from scratch. For example, you might use `h5py` to open an HDF5 file and then create a Dataset from that.
For MATLAB files you might use `scipy.io.loadmat` or `h5py` depending on the version of MATLAB file you're opening and then construct a Dataset.

---

## The scientific python ecosystem

Xarray ties in to the larger scientific python ecosystem and in turn many
packages build on top of xarray. A long list of such packages is here:
<https://docs.xarray.dev/en/stable/related-projects.html>.

Now we will demonstrate some cool features.


### Pandas: tabular data structures

You can easily [convert](https://docs.xarray.dev/en/stable/pandas.html) between xarray and [pandas](https://pandas.pydata.org/) structures. This allows you to conveniently use the extensive pandas 
ecosystem of packages (like [seaborn](https://seaborn.pydata.org/)) for your work.



```@example raster
# convert to pandas dataframe
df = ds.isel(time=slice(10)).to_dataframe()
df
```


```@example raster
# convert dataframe to xarray
df.to_xarray()
```

### Alternative array types

This notebook has focused on Numpy arrays. Xarray can wrap [other array](https://docs.xarray.dev/en/stable/user-guide/duckarrays.html) types! For example:

<img src="https://docs.dask.org/en/stable/_images/dask_horizontal.svg" width="20%"> [distributed parallel arrays](https://docs.dask.org/en/latest/array.html) & [Xarray user guide on Dask](https://docs.xarray.dev/en/stable/user-guide/dask.html)

<img src="https://raw.githubusercontent.com/pydata/sparse/master/docs/logo.svg" width="15%"> **pydata/sparse** : [sparse arrays](https://sparse.pydata.org)

<img src="https://raw.githubusercontent.com/cupy/cupy.dev/master/images/cupy_logo.png" width="22%"> [GPU arrays](https://cupy.dev) & [cupy-xarray](https://cupy-xarray.readthedocs.io/)

<img src="https://pint.readthedocs.io/en/stable/_static/logo-full.jpg" width="10%"> **pint** : [unit-aware arrays](https://pint.readthedocs.io) & [pint-xarray](https://github.com/xarray-contrib/pint-xarray)


### Dask

Dask cuts up NumPy arrays into blocks and parallelizes your analysis code across
these blocks

<img src="https://raw.githubusercontent.com/dask/dask/main/docs/source/images/dask-array.svg" style="width:45%">



```@example raster
# demonstrate dask dataset
dasky = xr.tutorial.open_dataset(
    "air_temperature",
    chunks={"time": 10},  # 10 time steps in each block
)

dasky.air
```

All computations with dask-backed xarray objects are lazy, allowing you to build
up a complicated chain of analysis steps quickly



```@example raster
# demonstrate lazy mean
dasky.air.mean("lat")
```

To get concrete values, call `.compute` or `.load`



```@example raster
# "compute" the mean
dasky.air.mean("lat").compute()
```

### HoloViz

Quickly generate interactive plots from your data!

The [`hvplot` package](https://hvplot.holoviz.org/user_guide/Gridded_Data.html) attaches itself to all
xarray objects under the `.hvplot` namespace. So instead of using `.plot` use `.hvplot`


```@example raster
import hvplot.xarray

ds.air.hvplot(groupby="time", clim=(270, 300), widget_location='bottom')
```

```{note}
The time slider will only work if you're executing the notebook, rather than viewing the website
```

### cf_xarray 

[cf_xarray](https://cf-xarray.readthedocs.io/) is a project that tries to
let you make use of other CF attributes that xarray ignores. It attaches itself
to all xarray objects under the `.cf` namespace.

Where xarray allows you to specify dimension names for analysis, `cf_xarray`
lets you specify logical names like `"latitude"` or `"longitude"` instead as
long as the appropriate CF attributes are set.

For example, the `"longitude"` dimension in different files might be labelled as: (lon, LON, long, x…), but cf_xarray let's you always refer to the logical name `"longitude"` in your code:


```@example raster
import cf_xarray
```


```@example raster
# describe cf attributes in dataset
ds.air.cf
```

The following `mean` operation will work with any dataset that has appropriate
attributes set that allow detection of the "latitude" variable (e.g.
`units: "degress_north"` or `standard_name: "latitude"`)



```@example raster
# demonstrate equivalent of .mean("lat")
ds.air.cf.mean("latitude")
```


```@example raster
# demonstrate indexing
ds.air.cf.sel(longitude=242.5, method="nearest")
```

### Other cool packages

- [xgcm](https://xgcm.readthedocs.io/) : grid-aware operations with xarray
  objects
- [xrft](https://xrft.readthedocs.io/) : fourier transforms with xarray
- [xclim](https://xclim.readthedocs.io/) : calculating climate indices with
  xarray objects
- [intake-xarray](https://intake-xarray.readthedocs.io/) : forget about file
  paths
- [rioxarray](https://corteva.github.io/rioxarray/stable/index.html) : raster
  files and xarray
- [xesmf](https://xesmf.readthedocs.io/) : regrid using ESMF
- [MetPy](https://unidata.github.io/MetPy/latest/index.html) : tools for working
  with weather data

Check the Xarray [Ecosystem](https://docs.xarray.dev/en/stable/ecosystem.html) page and [this tutorial](https://tutorial.xarray.dev/intermediate/xarray_ecosystem.html) for even more packages and demonstrations.

## Next

1. Read the [tutorial](https://tutorial.xarray.dev) material and [user guide](https://docs.xarray.dev/en/stable/user-guide/index.html)
1. See the description of [common terms](https://docs.xarray.dev/en/stable/terminology.html) used in the xarray documentation: 
1. Answers to common questions on "how to do X" with Xarray are [here](https://docs.xarray.dev/en/stable/howdoi.html)
1. Ryan Abernathey has a book on data analysis with a [chapter on Xarray](https://earth-env-data-science.github.io/lectures/xarray/xarray_intro.html)
1. [Project Pythia](https://projectpythia.org/) has [foundational](https://foundations.projectpythia.org/landing-page.html) and more [advanced](https://cookbooks.projectpythia.org/) material on Xarray. Pythia also aggregates other [Python learning resources](https://projectpythia.org/resource-gallery.html).
1. The [Xarray Github Discussions](https://github.com/pydata/xarray/discussions) and [Pangeo Discourse](https://discourse.pangeo.io/) are good places to ask questions.
1. Tell your friends! Tweet!


## Welcome!

DimensionalData and the whole Julia data analyiss ecosystem is an open-source project and gladly welcomes all kinds of contributions. This could include reporting bugs, discussing new enhancements, contributing code, helping answer user questions, contributing documentation (even small edits like fixing spelling mistakes or rewording to make the text clearer). Welcome!
