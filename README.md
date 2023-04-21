# Probabilistic and causal modelling of genome-scale data

This is a set of lecture notes on machine learning methods for analyzing genome-scale data, taught as part of the [BINF301 course](https://www.uib.no/en/course/BINF301).

The notes are written in markdown and converted to a static site using [Docsy][], a [Hugo theme module][] for technical documentation sites.

The public URL for the site is: [https://tmichoel.github.io/genome-scale-modelling/](https://tmichoel.github.io/genome-scale-modelling/)

## Running the website locally

You can find detailed theme instructions in the [Docsy user guide][].

Building and running the site locally requires a recent `extended` version of [Hugo](https://gohugo.io).
You can find out more about how to install Hugo for your environment in our
[Getting started](https://www.docsy.dev/docs/getting-started/#prerequisites-and-installation) guide.

Once you've made your working copy of the site repo, from the repo root folder, run:

```
hugo server
```

### Troubleshooting

As you run the website locally, you may run into the following error:

```
➜ hugo server

INFO 2021/01/21 21:07:55 Using config file: 
Building sites … INFO 2021/01/21 21:07:55 syncing static files to /
Built in 288 ms
Error: Error building site: TOCSS: failed to transform "scss/main.scss" (text/x-scss): resource "scss/scss/main.scss_9fadf33d895a46083cdd64396b57ef68" not found in file cache
```

This error occurs if you have not installed the extended version of Hugo.
See this [section](https://www.docsy.dev/docs/get-started/docsy-as-module/installation-prerequisites/#install-hugo) of the user guide for instructions on how to install Hugo.

Or you may encounter the following error:

```
➜ hugo server

Error: failed to download modules: binary with name "go" not found
```

This error occurs if you have not installed the `go` programming language on your system.
See this [section](https://www.docsy.dev/docs/get-started/docsy-as-module/installation-prerequisites/#install-go-language) of the user guide for instructions on how to install `go`.


[alternate dashboard]: https://app.netlify.com/sites/goldydocs/deploys
[deploys]: https://app.netlify.com/sites/docsy-example/deploys
[Docsy user guide]: https://docsy.dev/docs
[Docsy]: https://github.com/google/docsy
[example.docsy.dev]: https://example.docsy.dev
[Hugo theme module]: https://gohugo.io/hugo-modules/use-modules/#use-a-module-for-a-theme
[Netlify]: https://netlify.com

## Code

The "code" subfolder contains example notebooks to illustrate the course material.

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> Genome-scale ML code

To (locally) reproduce this project, do the following:

1. Download this code base. Notice that raw data are typically not included in the git-history and may need to be downloaded independently.
2. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "Genome-scale ML code"
```
which auto-activate the project and enable local path handling from DrWatson.

The main code files are [Pluto notebooks](https://plutojl.org/) and are located in the [notebooks folder](https://github.com/tmichoel/genome-scale-modelling/tree/main/code/notebooks). To run the reactive notebooks, follow the instruction at the [bottom of the Pluto homepage](https://plutojl.org/).
