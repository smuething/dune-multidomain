MultiDomain
===========

This is the 2.0.0 release of [MultiDomain][1], a set of multi-domain
extensions for the PDE solver toolbox [PDELab][2] that is part of the [DUNE][3]
project. See the file README for further details including license information.

MultiDomain 2.0.0 is a minor release, mainly focusing on bug fixes and keeping
current with changes to the underlying DUNE modules and PDELab.

If you need help, please ask via [GitHub][1]) If you find bugs, you can also
submit them to the [bugtracker][4]. Even better, if you have managed to fix a
problem, open a [pull request][5] to get your patch merged into the library.


Changes
-------

### MultiDomainGrid 2.0

* First "official" release.

* Support for building with CMake.

* Compatibility with DUNE 2.3.x, PDELab 2.0.x and MultiDomainGrid 2.3.x.

* Some minor bugfixes.

* Release history

  * `2.0.0` Final release
    * Minor bugfixes

  * `2.0.0-rc1` Initial release candidate


Caveats
-------

The following list is a non-exhaustive overview of possible problems you might
encounter with the current release.


### General

* Compile times can be really long for non-trivial problems. Some developers
  have had good success with using the clang compiler instead of GCC during
  development and bug-testing to reduce compile times.

* After MultiDomain 2.0, the minimum compiler requirement of MultiDomain will
  be increased to GCC 4.7. Please be aware of this change in minimum
  requirements.


Links
-----

[1]: http://github.com/smuething/dune-multidomain
[2]: http://dune-project.org/pdelab/
[3]: http://dune-project.org
[4]: https://github.com/smuething/dune-multidomain/issues
[5]: https://github.com/smuething/dune-multidomain/pulls
