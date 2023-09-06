# RISM grids

Some jupyter-notebooks to visualise RISM data and align RISM grids to molecules.  Also includes classes to load and manipulate RISM data generated using AMBER.

- `RISM_grids_usage.ipynb`: Loading and manipulating RISM data using the classes in `RISM_grids.py`.
- `RISM_3D_plots.ipynb`: Examples of plotting RISM data in 3D, including interactive plots.

The `html` files contain working examples of the jupyter-notebooks, including the interactive plots in `RISM_3D_plots.ipynb`.  Open them via the links below (uses https://htmlpreview.github.io/):
- [RISM_grids_usage.html](https://htmlpreview.github.io/?https://github.com/jc32173/RISM_grids/blob/master/RISM_grids_usage.html)
- [RISM_3D_plots.html](https://htmlpreview.github.io/?https://github.com/jc32173/RISM_grids/blob/master/RISM_3D_plots.html)

## Examples:
(From `RISM_3D_plots.ipynb`)

- #### Align molecule with 2D slices through 3D RISM data.
  Code will generate VMD commands to align the view of the molecule with the corresponding 2D RISM slice.
<img src="examples/3D_RISM_VMD_plot_modified.png" width="220">

- #### Interactive plot showing the RISM distribution around a molecule.
  The slider moves the 2D surface over the molecule to show how the RISM grid reflects the 3D shape of the molecule.  See: 
[interactive_fig_x](https://htmlpreview.github.io/?https://github.com/jc32173/RISM_grids/blob/master/interactive_fig_x.html), 
[interactive_fig_y](https://htmlpreview.github.io/?https://github.com/jc32173/RISM_grids/blob/master/interactive_fig_y.html), 
[interactive_fig_z](https://htmlpreview.github.io/?https://github.com/jc32173/RISM_grids/blob/master/interactive_fig_z.html).
![RISM_plot](examples/example_slider_plots.png)
