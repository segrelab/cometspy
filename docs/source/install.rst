.. Installation documentation page, written by Helen Scott, last
   updated on 2023-08-21.
   You can adapt this file completely to your liking, but it should at
   least contain the root `toctree` directive.

.. The toctree directive controls what links are in the "Navigation" bar
   at the top of each page. The "maxdepth" argument determines how many
   levels deep the links are shown. The "caption" argument determines the
   title of the "Navigation" bar.
Installation
============

How to Install:
----------------------------
pip it:

pip install cometspy

Warning! The current version will work only with Pandas version 1.5.0 or higer. 
Pandas changed the variable line_terminator to lineterminator in the 'csv' modules.
More at: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_csv.html
"Changed in version 1.5.0: Previously was line_terminator, changed for consistency with read_csv and the standard library ‘csv’ module."

If you get an error message about line_terminator, update your pandas package. 


