# ephys-jl

Julia versions of some of the ephys analysis routines. Faster than the python versions, but
stripped-down also. This is beta code and is not fully tested yet.

Note: as of 7 Dec 2023, the config file can be read in pure julia (no pycall)

To run LSPSAnalysis.jl fron a cold start:
--------

  Run julia.

  julia --threads 12 (or some number)

Make sure a filename is defined:
  for example: fn = "/Volumes/Pegasus_002/ManisLab_Data3/Kasten_Michael/VGAT_DCNmap/2018.08.31_000/slice_000/cell_000/VGAT_5mspulses_5Hz_pt045_to_pt1_000"

  include("src/LSPSAnalysis.jl")
  LSPSAnalysis.LSPS_read_and_plot(filename=fn)

That should generate a graph after a bit.

To run IVAnalysis.jl:

include("src/IVAnalysis.jl")
fniv = "/Volumes/Pegasus_002/ManisLab_Data3/Kasten_Michael/VGAT_DCNmap/2018.03.20_000/slice_000/cell_000/CCIV_1nA_max_000"






Cold:
Start Julia
 use the shell to go to the ephys-jl directory (cd ephys-jl)
 then load the package manager: 
] pkg
pkg >  activate .
pkg >  instantiate
pkg >  update
pkg > build ( or rebuild)
pkg > add .


Now build the working environment:
add HDF5_jll
add HDF5
just do add HDF5_jll (yes, that is right)

try : "add Acq4Reader"
add packages as needed

instantiate
activate

using Acq4Reader

Acq4Reader.read_hdf5(filename--- .ma)


======================
Building a compiled version:

function julia_main()::Cint
  # do something based on ARGS?
  return 0 # if things finished successfully
end

this could be useful:

(Acq4Reader) pkg> add AddPackage
   Resolving package versions...

using Pkg
Pkg.generate("myapp")
using PackageCompiler
create_app("myapp", "myappcompiled")

Run with: 
myappcompiled/bin/myapp args...


to use packages:
] PKG
pkg> add ephys





OLD (Don't do this!):
    Do the "brew install hdf5" separately, and maybe use "reinstall" to force update.
    be sure to install from a native M1 terminal, not the rosetta terminal (2023)
    Set Julia path in .zshrc:
    export JULIA_HDF5_PATH=/opt/homebrew/Cellar/hdf5/1.12.2_2
    restart julia to be sure it can be added in the package manager:
    ] pkg > add HDF5



Check that the correct version of Julia is run. May need to modify .zshrc (or .bash_profile or .bash_rc).
V1.7 is the first version that runs natively on M1 hardware.
The manifest is set up for Julia 1.8.5

using Pkg
Pkg.activate(".")
(Pkg.instantiate(), Pkg.upgrade_manifest())activate

or ] (pkg)
activate
instantiate

Do the includes (main Julia prompt):

include("src/Acq4Reader.jl")
First time around, this might take a while
