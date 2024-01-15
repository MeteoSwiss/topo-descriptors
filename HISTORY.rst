=======
History
=======

0.3.0 (2024-01-15)
------------------

* Move from DataArray to Dataset for DEM object to allow transferring global attributes.
* Add units as variable attributes.
* Output slope in units of [degree] instead of [m / pixel].
* Fix bug in slope calculation.
* Remove parallelization of scales with multiprocessing for valley and ridge

0.2.1 (2022-10-19)
------------------

* Fix bug in the scale to pixel conversion in case of WGS84 grids.

0.2.0 (2021-06-12)
------------------

* Add Sx descriptor.
* Add STD descriptor.

0.1.2 (2021-05-14)
------------------

* First working release on PyPI.
