import RISM_grids as rg
mol3d = rg.Grid3D.FromDxFile('mobley_7754849_conf0_guv.O.1.dx')
mol2d = mol3d.slice('i', 34)
mol2d.vmd_cmd()
mol2d.plt()
