# Classes to hold RISM data and read from .dx file.

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

__version__ = '2'

# TODO:
# ADD TITLE TO ReadFromDx so that molecule ID, type of RISM (e.g. guv), atom used etc can be recorded and also passed to slices.
# Take diagonal slices?
# Interpolate between points, especially for diagonal slices?

# Reads dx file from 3D RISM, assigns values to i,j,k grid points and plots heatmap of slice through the box.

class Grid3D():

    def __init__(self, gridsize, gridstep, origin, grid, name=None):
        self.gridsize = gridsize
        self.gridstep = gridstep
        self.origin = origin
        self.grid = grid
        self.name = name

    def __call__(self):
        print('RISM data on 3D grid')

    # Use a classmethod to read dx file and return a Grid3D object:
    @classmethod
    def FromDxFile(cls, dx_filename):
        print('Reading Dx file: {}'.format(dx_filename))
        # Parse file header:
        line = open(dx_filename, 'r').readlines()
	# SHOULD CLOSE THE FILE NOW AFTER READLINES??
        n = 0
        while True:
            #if line[n][:8] == 'object 3':
            if line[n].startswith('object 3'):
                data_start = n + 1
                break
            elif line[n].startswith('object 1 class gridpositions counts') and \
                 line[n+1].startswith('origin') and \
                 line[n+2].startswith('delta') and \
                 line[n+3].startswith('delta') and \
                 line[n+4].startswith('delta'):

                # Check grid is orthogonal:
#               if :
#                   sys.exit('ERROR: Grid is not orthogonal, check "delta" values')
                # THIS SHOULD NOT BE HARD CODED:
                gridstep = np.array([0.5, 0.5, 0.5])

                n_x, n_y, n_z = [int(n_grid) for n_grid in line[n].split()[-3:]]
                gridsize = np.array([n_x, n_y, n_z])
                origin = np.array([float(ori) for ori in line[n+1].split()[-3:]])
                #self.gridstep = [float(line[n+2].split[-3]), float(line[n+3].split[-2]), float(line[n+4].split[-1])]
                n += 5
            else:
                n += 1

        # Read in data:
        grid = np.zeros(gridsize)
        for line_n, line in enumerate(line[data_start:]):
            if line[:6] == 'object':
                break
            #elif len(line.split()) != 3:
            #    sys.exit('ERROR')
            for col_n, val in enumerate(line.split()):
                i =  (line_n*3 + col_n)//(n_z*n_y)
                j = ((line_n*3 + col_n) %(n_z*n_y)) // n_z
                k = ((line_n*3 + col_n) %(n_z*n_y))  % n_z
                grid[i,j,k] = float(val)

        return cls(gridsize, gridstep, origin, grid)

    # Could use miller indices to get diagonal slices?
    # or give normal vector?
    def slice(self, axis, level):
        if axis == 'x':
            #return self.grid[level,:,:]
            return Grid2D(self.gridsize[[1, 2]], 
                          self.gridstep[[1, 2]], 
                          np.array([level*self.gridstep[0], 0, 0]) + self.origin, 
                          self.grid[level,:,:],
                          ('y', 'z'), 
                          [1, 0, 0])
        elif axis == 'y':
            #return self.grid[:,level,:]
            return Grid2D(self.gridsize[[0, 2]], 
                          self.gridstep[[0, 2]], 
                          np.array([0, level*self.gridstep[1], 0]) + self.origin, 
                          self.grid[:,level,:], 
                          ('x', 'z'), 
                          [0, 1, 0])
        elif axis == 'z':
            #return self.grid[:,:,level]
            return Grid2D(self.gridsize[[0, 1]], 
                          self.gridstep[[0, 1]], 
                          np.array([0, 0, level*self.gridstep[2]]) + self.origin, 
                          self.grid[:,:,level], 
                          ('x', 'y'), 
                          [0, 0, 1])
        else:
            raise ValueError('axis must be "x", "y" or "z"')

#    def diag_slice(norm, ori = [0, 0, 0]):
#        plane = []
#        Get diagonal elements
#        np.diagonal()
#        return Grid2D(self.gridsize[[0, 1]],
#                          self.gridstep[[0, 1]],
#                          np.array([0, 0, level*self.gridstep[2]]) + self.origin,
#                          self.grid[:,:,level],
#                          ('x', 'y'),
#                          [0, 0, 1])


    # Reduce the grid size by 1/n:
    def coarsen(self, n):
        coarse_grid = np.zeros(self.gridsize//n)
        print(np.shape(coarse_grid))
        coarse_gridsize = self.gridsize//n
        # SHOULD CENTRE THE RANGE THAT IS KEPT IN CASE GRIDSIZE NOT EXACTLY DIVISIBLE BY n
        for i in range(0, self.gridsize[0]):
            for j in range(0, self.gridsize[1]):
                for k in range(0, self.gridsize[2]):
                    coarse_grid[i,j,k] = self.grid[i*n,j*n,k*n]
        return Grid3D(coarse_gridsize, self.gridstep, self.origin, coarse_grid)
	# Also update origin?
        #pass

    # FOR 2D OR 3D:
    # Reduce the grid size by 1/n:
    #def coarsen(self, n):
    #    coarse_grid = np.zeros(self.gridsize//n)
    #    print(np.shape(coarse_grid))
    #    self.gridsize = self.gridsize//n
    #    # SHOULD CENTRE THE RANGE THAT IS KEPT IN CASE GRIDSIZE NOT EXACTLY DIVISIBLE BY n
    #    for d in len(self.gridsize):
    #        coarse_grid????
    #    for i in range(0, self.gridsize[0]):
    #        for j in range(0, self.gridsize[1]):
    #            for k in range(0, self.gridsize[2]):
    #                coarse_grid[i,j,k] = self.grid[i*n,j*n,k*n]
    #    self.grid = coarse_grid
    #    # Also update origin?
    #    #pass

    def prune(self, r, nullvalue=-1):
        #for d in len(self.gridsize):
        pruned_grid = np.zeros(self.gridsize)
        for i in range(0, self.gridsize[0]):
            for j in range(0, self.gridsize[1]):
        #    for i in range(0, self.gridsize[0]):
                for k in range(0, self.gridsize[2]):
                    if 0.0 not in ((self.grid[max(0, i-r):i+r+1,:,:])[:,max(0, j-r):j+r+1,:])[:,:,max(0, k-r):k+r+1]:
                        #self.grid[i,j,k] = nullvalue #np.NaN #self.grid[i,j,k]
                        # Cannot set existing values to nan, must be values in new zeros array:
                        pruned_grid[i,j,k] = np.nan
                    else:
                        pruned_grid[i,j,k] = self.grid[i,j,k]
        return Grid3D(self.gridsize, self.gridstep, self.origin, pruned_grid)

    def vmd_box(self, vmd_filename=None):

        ori = self.origin
        incr_grid = self.gridstep
        n_x, n_y, n_z = self.gridsize

        # save file for vmd to plot box and slice
        
        # box limits:
        lims = [[ori[0], ori[0] + incr_grid[0]*n_x], [ori[1], ori[1] + incr_grid[1]*n_y], [ori[2], ori[2] + incr_grid[2]*n_z]]
        
        # Draw box:
        vmd_box = 'proc draw_rism_box {} {\n' \
                  'set sel [atomselect top all]\ndraw color yellow\n'
        for x in lims[0]:
           for y in lims[1]:
               vmd_box += 'draw line "%.3f %.3f %.3f" "%.3f %.3f %.3f" style solid width 1\n' % (x, y, lims[2][0], x, y, lims[2][1])
           for z in lims[2]:
               vmd_box += 'draw line "%.3f %.3f %.3f" "%.3f %.3f %.3f" style solid width 1\n' % (x, lims[1][0], z, x, lims[1][1], z)
        for y in lims[1]:
           for z in lims[2]:
               vmd_box += 'draw line "%.3f %.3f %.3f" "%.3f %.3f %.3f" style solid width 1\n' % (lims[0][0], y, z, lims[0][1], y, z)
        vmd_box += '}\n'

        if vmd_filename is None:
            print(vmd_box)
        else:
            vmd_file = open(vmd_filename, 'w')
            vmd_file.write(vmd_box)
            vmd_file.close()


class Grid2D():

    def __init__(self, gridsize, gridstep, origin, grid, gridaxes=('x', 'y'), normvec=[0, 0, 1]):
        self.gridsize = gridsize
        self.gridstep = gridstep
        self.origin = origin
        self.grid = grid
        self.gridaxes = gridaxes
        self.norm = {'x', 'y', 'z'} - set(self.gridaxes)
        # Convert set to string:
        self.norm = self.norm.pop()
#        self.normvec = np.array(normvec)
#        self.axes = []

    def __call__(self):
        print('RISM data on 2D grid')

    # THIS NEEDS FIXING:
    def slice(self, axis, level):
        # NEED TO ALSO OUTPUT THE ORIGIN
        if axis == 'x':
            return Grid1D(self.gridsize[1], self.gridstep, self.grid[level,:])
        elif axis == 'y':
            return Grid1D(self.gridsize[0], self.gridstep, self.grid[:,level])
        else:
            raise ValueError('axis must be "x" or "y"')

    # Plot heatmap:
    def plt(self, **kwargs):
        # Pass **kwargs to heatmap to control plot:
        sns.heatmap(self.grid, **kwargs)
        # Note that for a heatmap, the first array index gives 
        # the row number (position along the y axis), not the x 
        # coordinate, so this is the ylabel:
        plt.xlabel(self.gridaxes[1])
        plt.ylabel(self.gridaxes[0])
        plt.gca().set_aspect('equal')
        plt.show()

    # Generate functions for VMD:
    # Useful VMD commands reference: 
    # https://www.ks.uiuc.edu/Research/vmd/vmd_help.html
    def vmd_cmds(self, vmd_filename=None):

        # Record all vmd functions:
        vmd_fns = ''

        # Draw edges of slice:
        # --------------------

        ind = {'x' : 0, 'y' : 1, 'z' : 2}

        off_corns = []
        # Have to make a copy, otherwise all list elements just point to self.origin
        oris = [self.origin.copy()]
        ori = self.origin.copy()
        for i2d, ax in enumerate(self.gridaxes):
            i3d = ind[ax]
            off_corn = self.origin.copy()
            off_corn[i3d] += self.gridstep[i2d]*self.gridsize[i2d]
            off_corns.append(off_corn)
            ori[i3d] += self.gridstep[i2d]*self.gridsize[i2d]
        oris.append(ori)

        vmd_edges = 'proc draw_rism_edges {} {\n' \
                    'set sel [atomselect top all]\n' \
                    'draw color yellow\n'
        
        for ori in oris:
            for off_corn in off_corns:
                vmd_edges += 'draw line "%.3f %.3f %.3f" "%.3f %.3f %.3f" style solid width 1\n' % (ori[0], ori[1], ori[2], off_corn[0], off_corn[1], off_corn[2])

        vmd_edges += '}\n'

        vmd_fns += vmd_edges
        
        # Draw slice as solid plane:
        # --------------------------

        # Use Trinorm to define normals to edges parallel to the surface 
        # to prevent the surface diasppearing when looked from the side:
        #https://www.ks.uiuc.edu/Research/vmd/vmd-1.7.1/ug/node171.html
        
        #y_coord = ori[1] + 34*incr_grid[1]
        vmd_slice = 'proc draw_rism_area {} {\n' \
                    'set sel [atomselect top all]\n' \
                    'draw color yellow\n'

        # Draw slice plane:
        for ori in oris:
            vmd_slice += 'draw triangle "%.3f %.3f %.3f" "%.3f %.3f %.3f" "%.3f %.3f %.3f"\n' % (
                ori[0]-0.25, ori[1], ori[2], 
                off_corns[0][0]-0.25, off_corns[0][1], off_corns[0][2], 
                off_corns[1][0]-0.25, off_corns[1][1], off_corns[1][2])

        for ori in oris:
            vmd_slice += 'draw triangle "%.3f %.3f %.3f" "%.3f %.3f %.3f" "%.3f %.3f %.3f"\n' % (
                ori[0]+0.25, ori[1], ori[2], 
                off_corns[0][0]+0.25, off_corns[0][1], off_corns[0][2], 
                off_corns[1][0]+0.25, off_corns[1][1], off_corns[1][2])

        vmd_slice += '}\n'

        vmd_fns += vmd_slice

        # Clip plane:
        # -----------

        # https://www.ks.uiuc.edu/Training/Tutorials/vmd-imgmv/imgmv/tutorial-html/node2.html

#        clip_norm = "1 0 0"
#        clip_cen = "1 0 0"
        ind = {'x' : 0, 'y' : 1, 'z' : 2}
        plane_ind = ind[self.norm]
        clip_norm = [0, 0, 0]
        clip_cen = [0, 0, 0]
        clip_norm[plane_ind] = 1
        clip_cen[plane_ind] = self.origin[plane_ind]

        # Functions to add clipplanes:
        # Clip below slice:
        clip_cmd =  'proc clipplane_above {{mol 0}} {\n' \
                    'mol clipplane status 0 0 $mol 2\n'
        clip_cmd += 'mol clipplane normal 0 0 $mol "{}"\n'.format(' '.join([str(i) for i in clip_norm]))
        clip_cmd += 'mol clipplane center 0 0 $mol "{}"\n'.format(' '.join([str(j) for j in clip_cen]))
        clip_cmd += '}\n'
        vmd_fns += clip_cmd

        # Clip above slice:
        clip_norm[plane_ind] = -1
        clip_cmd =  'proc clipplane_below {{mol 0}} {\n' \
                    'mol clipplane status 0 0 $mol 2\n'
        clip_cmd += 'mol clipplane normal 0 0 $mol "{}"\n'.format(' '.join([str(i) for i in clip_norm]))
        clip_cmd += 'mol clipplane center 0 0 $mol "{}"\n'.format(' '.join([str(j) for j in clip_cen]))
        clip_cmd += '}\n'
#                   '}\n'.format(' '.join([str(i) for i in clip_norm]), ' '.join([str(j) for j in clip_cen]))
        vmd_fns += clip_cmd

        # NEED TO SORT OUT THESE FUNCTIONS:
        clip_norm[plane_ind] = 1
        clip_cen[plane_ind] = self.origin[plane_ind] - 1
        clip_cmd = 'proc vmd_clip_up {{mol 0}} {\n' \
                   'mol clipplane status 1 0 0 2\n'
        clip_cmd += 'mol clipplane normal 1 0 0 "{}"\n'.format(' '.join([str(i) for i in clip_norm]))
        clip_cmd += 'mol clipplane center 1 0 0 "{}"\n'.format(' '.join([str(j) for j in clip_cen]))
        clip_cmd += '}\n'
        vmd_fns += clip_cmd

        clip_norm[plane_ind] = -1
        clip_cen[plane_ind] = self.origin[plane_ind] + 1
        clip_cmd = 'proc vmd_clip_down {{mol 0}} {\n' \
                   'mol clipplane status 2 0 0 2\n'
        clip_cmd += 'mol clipplane normal 2 0 0 "{}"\n'.format(' '.join([str(i) for i in clip_norm]))
        clip_cmd += 'mol clipplane center 2 0 0 "{}"\n'.format(' '.join([str(j) for j in clip_cen]))
        clip_cmd += '}\n'
        vmd_fns += clip_cmd

        # Unset all clipping planes:
        rm_clip = 'proc rm_clipplanes {{mol 0}} {\n' \
                  'mol clipplane status 0 0 $mol 0\n' \
                  'mol clipplane status 1 0 $mol 0\n' \
                  'mol clipplane status 2 0 $mol 0\n' \
                  '}\n'
        vmd_fns += rm_clip

        # Rotation function:
        # ------------------

        if self.norm == 'x':
            rot = 'proc rot_view_x {} {\n' \
                  'display resetview\n' \
                  'rotate y by -90\n' \
                  'rotate z by 180\n' \
                  '}\n'
        elif self.norm == 'y':
            rot = 'proc rot_view_y {} {\n' \
                  'display resetview\n' \
                  'rotate x by -90\n' \
                  'rotate z by -90\n' \
                  '}\n'
        elif self.norm == 'z':
            rot = 'proc rot_view_z {} {\n' \
                  'display resetview\n' \
                  'rotate z by -90\n' \
                  '}\n'
        vmd_fns += rot

        # Print or write list of functions:

        if vmd_filename is None:
            print(vmd_fns)
        else:
            vmd_file = open(vmd_filename, 'w')
            vmd_file.write(vmd_fns)
            vmd_file.close()


class Grid1D():

    def __init__(self, gridsize, gridstep, grid, gridaxis='x', name=None):
        self.gridsize = gridsize
        self.gridsize = gridstep
        self.origin = origin
        self.grid = grid
        self.gridaxis = gridaxis
        self.name = name

    def __call__(self):
        print('RISM data on 1D grid')

    def plt(self):
        plt.plot(list(range(self.gridsize)), self.grid)
        plt.xlabel('Distance along {} (grid points)'.format(self.gridaxis))
        plt.ylabel('RISM')
        plt.show()
