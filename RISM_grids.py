# Classes to hold RISM data and read from .dx file.
# Reads dx file from 3D RISM, assigns values to i,j,k grid points and 
# plots heat map of slices through the box.

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

__version__ = '3'

coord_inds = {'x' : 0, 'y' : 1, 'z' : 2}

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
        n = 0
        while True:
            if line[n].startswith('object 3'):
                data_start = n + 1
                break
            elif line[n].startswith('object 1 class gridpositions counts') and \
                 line[n+1].startswith('origin') and \
                 line[n+2].startswith('delta') and \
                 line[n+3].startswith('delta') and \
                 line[n+4].startswith('delta'):

                # Check grid is orthogonal:
                if not np.all(np.array([float(line[n+2].split()[-2]), float(line[n+2].split()[-1]),
                                        float(line[n+3].split()[-3]), float(line[n+3].split()[-1]),
                                        float(line[n+4].split()[-3]), float(line[n+4].split()[-2])]) == 0):
                    sys.exit('ERROR: Grid is not orthogonal, check "delta" values')

                n_x, n_y, n_z = [int(n_grid) for n_grid in line[n].split()[-3:]]
                gridsize = np.array([n_x, n_y, n_z])
                origin = np.array([float(ori) for ori in line[n+1].split()[-3:]])
                gridstep = np.array([float(line[n+2].split()[-3]), 
                                     float(line[n+3].split()[-2]), 
                                     float(line[n+4].split()[-1])])
                n += 5
            else:
                n += 1

        # Read in data:
        grid = np.zeros(gridsize)
        for line_n, line in enumerate(line[data_start:]):
            if line[:6] == 'object':
                break
            for col_n, val in enumerate(line.split()):
                # Maybe better to read all values into a 1D array and then
                # use numpy reshape.
                i =  (line_n*3 + col_n)//(n_z*n_y)
                j = ((line_n*3 + col_n) %(n_z*n_y)) // n_z
                k = ((line_n*3 + col_n) %(n_z*n_y))  % n_z
                grid[i,j,k] = float(val)

        return cls(gridsize, gridstep, origin, grid, name=dx_filename)

    # Take 2D slice through RISM box:
    def slice(self, axis, level):
        if axis == 'x':
            return Grid2D(self.gridsize[[1, 2]], 
                          self.gridstep[[1, 2]], 
                          self.grid[level,:,:],
                          np.array([level*self.gridstep[0], 0, 0]) + self.origin, 
                          ('y', 'z'), 
                          level,
                          name=self.name+", "+axis+" = "+str(level))
        elif axis == 'y':
            return Grid2D(self.gridsize[[2, 0]], 
                          self.gridstep[[2, 0]],
                          # For y, have to take transpose to switch order of x and z axes
                          # to match permutation of axis labels: 
                          self.grid[:,level,:].T, 
                          np.array([0, level*self.gridstep[1], 0]) + self.origin, 
                          ('z', 'x'), 
                          level,
                          name=self.name+", "+axis+" = "+str(level))
        elif axis == 'z':
            return Grid2D(self.gridsize[[0, 1]], 
                          self.gridstep[[0, 1]], 
                          self.grid[:,:,level], 
                          np.array([0, 0, level*self.gridstep[2]]) + self.origin, 
                          ('x', 'y'), 
                          level,
                          name=self.name+", "+axis+" = "+str(level))
        else:
            raise ValueError('axis must be "x", "y" or "z"')

    # Reduce the grid size by 1/n:
    # Should probably centre the range.
    def coarsen(self, n):
        coarse_grid = np.zeros(self.gridsize//n)
        coarse_gridsize = self.gridsize//n
        for i in range(0, self.gridsize[0]):
            for j in range(0, self.gridsize[1]):
                for k in range(0, self.gridsize[2]):
                    coarse_grid[i,j,k] = self.grid[i*n,j*n,k*n]
        return Grid3D(coarse_gridsize, self.gridstep, self.origin, coarse_grid)

    # Only keep data within r of the molecule:
    def prune(self, r, nullvalue=-1):
        pruned_grid = np.zeros(self.gridsize)
        for i in range(0, self.gridsize[0]):
            for j in range(0, self.gridsize[1]):
                for k in range(0, self.gridsize[2]):
                    if 0.0 not in ((self.grid[max(0, i-r):i+r+1,:,:])[:,max(0, j-r):j+r+1,:])[:,:,max(0, k-r):k+r+1]:
                        #self.grid[i,j,k] = nullvalue #np.NaN #self.grid[i,j,k]
                        # Cannot set existing values to nan, must be values in new zeros array:
                        pruned_grid[i,j,k] = np.nan
                    else:
                        pruned_grid[i,j,k] = self.grid[i,j,k]
        return Grid3D(self.gridsize, self.gridstep, self.origin, pruned_grid)

    # Print function to plot VMD box edges:
    def vmd_box(self, vmd_filename=None):
        ori = self.origin
        incr_grid = self.gridstep
        n_x, n_y, n_z = self.gridsize

        # Box limits:
        lims = [[ori[0], ori[0] + incr_grid[0]*n_x], [ori[1], ori[1] + incr_grid[1]*n_y], [ori[2], ori[2] + incr_grid[2]*n_z]]
        
        # Draw box function:
        vmd_box = 'proc draw_rism_box {} {\n' \
                  'set sel [atomselect top all]\ndraw color blue\n'
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

    def __init__(self, gridsize, gridstep, grid, origin=[0, 0, 0], gridaxes=('x', 'y'), norm_intercept=0, name=None):
        self.gridsize = gridsize
        self.gridstep = gridstep
        self.grid = grid
        self.origin = origin
        self.gridaxes = gridaxes
        self.norm = {'x', 'y', 'z'} - set(self.gridaxes)
        # Convert set to string:
        self.norm = self.norm.pop()
        self.norm_intercept = norm_intercept
        self.name = name

    def __call__(self):
        print('RISM data on 2D grid:')

    # Get a 1D profile:
    def slice(self, axis, level):
        parallel_axis = set(self.gridaxes) - set(axis)
        # Convert set to string:
        parallel_axis = parallel_axis.pop()
        new_origin = self.origin.copy()
        if axis == self.gridaxes[0]:
            new_origin[coord_inds[axis]] += level*self.gridstep[0]
            return Grid1D(self.gridsize[1], 
                          self.gridstep[1], 
                          self.grid[level,:], 
                          new_origin, 
                          parallel_axis, 
                          name=self.name+", "+axis+" = "+str(level))
        elif axis == self.gridaxes[1]:
            new_origin[coord_inds[axis]] += level*self.gridstep[1]
            return Grid1D(self.gridsize[0], 
                          self.gridstep[0], 
                          self.grid[:,level], 
                          new_origin, 
                          parallel_axis, 
                          name=self.name+", "+axis+" = "+str(level))
        else:
            raise ValueError('axis must be "{}" or "{}"'.format(self.gridaxes[0], self.gridaxes[1]))

    # Plot heatmap:
    def plt(self, ax=None, **kwargs):
        show = False
        # If axes not passed to function, create some:
        if ax is None:
            show = True
            ax = plt.subplot(1, 1, 1)
        # Pass **kwargs to heatmap to control plot:
        sns.heatmap(self.grid, ax=ax, **kwargs)
        # Note that for a heatmap, the first array index gives 
        # the row number (position along the y axis), not the x 
        # coordinate, so this is the ylabel:
        ax.set_xlabel(str(self.gridaxes[1])+" (grid points)")
        ax.set_ylabel(str(self.gridaxes[0])+" (grid points)")
        ax.set_aspect('equal')
        if show:
            plt.show()
        else:
            return ax

    # Generate functions for VMD:
    # Useful VMD commands reference: 
    # https://www.ks.uiuc.edu/Research/vmd/vmd_help.html
    def vmd_cmds(self, vmd_filename=None):

        # Record all vmd functions:
        vmd_fns = ''

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
                  'rotate x by 90\n' \
                  '}\n'
        elif self.norm == 'z':
            rot = 'proc rot_view_z {} {\n' \
                  'display resetview\n' \
                  'rotate z by -90\n' \
                  '}\n'
        vmd_fns += rot

        # Draw slice:
        # -----------

        off_corners = []
        # Have to make a copy, otherwise all list elements just point to self.origin
        oris = [self.origin.copy()]
        ori = self.origin.copy()
        for i2d, ax in enumerate(self.gridaxes):
            i3d = coord_inds[ax]
            off_corner = self.origin.copy()
            off_corner[i3d] += self.gridstep[i2d]*self.gridsize[i2d]
            off_corners.append(off_corner)
            ori[i3d] += self.gridstep[i2d]*self.gridsize[i2d]
        oris.append(ori)

        # Draw edges:

        vmd_edges = 'proc draw_rism_edges {} {\n' \
                    'set sel [atomselect top all]\n' \
                    'draw color blue\n'
        
        for ori in oris:
            for off_corner in off_corners:
                vmd_edges += 'draw line "%.3f %.3f %.3f" "%.3f %.3f %.3f" style solid width 1\n' % (ori[0], ori[1], ori[2], off_corner[0], off_corner[1], off_corner[2])

        vmd_edges += '}\n'

        vmd_fns += vmd_edges
        
        # Draw solid plane:

        #y_coord = ori[1] + 34*incr_grid[1]
        vmd_slice = 'proc draw_rism_area {} {\n' \
                    'set sel [atomselect top all]\n' \
                    'draw color blue\n'

        for ori in oris:
            vmd_slice += 'draw triangle "%.3f %.3f %.3f" "%.3f %.3f %.3f" "%.3f %.3f %.3f"\n' % (
                ori[0], ori[1], ori[2], 
                off_corners[0][0], off_corners[0][1], off_corners[0][2], 
                off_corners[1][0], off_corners[1][1], off_corners[1][2])

        # Tried using trinorm to define normals to edges parallel to the surface 
        # to prevent the surface diasppearing when looked from the side:
        #https://www.ks.uiuc.edu/Research/vmd/vmd-1.7.1/ug/node171.html
#        for i in [0, 1]:
#            ori = oris[i]
#            ori_op = oris[(i+1)%2]
#            vmd_slice += 'draw trinorm "%.3f %.3f %.3f" "%.3f %.3f %.3f" "%.3f %.3f %.3f" "%.3f %.3f %.3f" "%.3f %.3f %.3f" "%.3f %.3f %.3f"\n' % (
#                ori[0], ori[1], ori[2], 
#                off_corners[0][0], off_corners[0][1], off_corners[0][2], 
#                off_corners[1][0], off_corners[1][1], off_corners[1][2],
#                2*ori[0]-off_corners[0][0]-off_corners[1][0], 2*ori[1]-off_corners[0][1]-off_corners[1][1], 2*ori[2]-off_corners[0][2]-off_corners[1][2],
#                2*off_corners[i][0]-ori[0]-ori_op[0], 2*off_corners[i][1]-ori[1]-ori_op[1], 2*off_corners[i][2]-ori[2]-ori_op[2],
#                2*off_corners[i][0]-ori[0]-ori_op[0], 2*off_corners[i][1]-ori[1]-ori_op[1], 2*off_corners[i][2]-ori[2]-ori_op[2])

        vmd_slice += '}\n'

        vmd_fns += vmd_slice

        # Clip plane:
        # -----------

        # https://www.ks.uiuc.edu/Training/Tutorials/vmd-imgmv/imgmv/tutorial-html/node2.html

        plane_ind = coord_inds[self.norm]
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
        vmd_fns += clip_cmd

        # Clip 2 grid points below slice:
        clip_norm[plane_ind] = 1
        clip_cen[plane_ind] = self.origin[plane_ind] - self.gridsize[plane_ind]*2
        clip_cmd = 'proc clipplane_above_margin {{mol 0}} {\n' \
                    'mol clipplane status 1 0 $mol 2\n'
        clip_cmd += 'mol clipplane normal 1 0 $mol "{}"\n'.format(' '.join([str(i) for i in clip_norm]))
        clip_cmd += 'mol clipplane center 1 0 $mol "{}"\n'.format(' '.join([str(j) for j in clip_cen]))
        clip_cmd += '}\n'
        vmd_fns += clip_cmd

        # Clip 2 grid points above slice:
        clip_norm[plane_ind] = -1
        clip_cen[plane_ind] = self.origin[plane_ind] + self.gridsize[plane_ind]*2
        clip_cmd = 'proc clipplane_below_margin {{mol 0}} {\n' \
                    'mol clipplane status 2 0 $mol 2\n'
        clip_cmd += 'mol clipplane normal 2 0 $mol "{}"\n'.format(' '.join([str(i) for i in clip_norm]))
        clip_cmd += 'mol clipplane center 2 0 $mol "{}"\n'.format(' '.join([str(j) for j in clip_cen]))
        clip_cmd += '}\n'
        vmd_fns += clip_cmd

        # Unset all clipping planes:
        rm_clip = 'proc rm_clipplanes {{mol 0}} {\n' \
                  'mol clipplane status 0 0 $mol 0\n' \
                  'mol clipplane status 1 0 $mol 0\n' \
                  'mol clipplane status 2 0 $mol 0\n' \
                  '}\n'
        vmd_fns += rm_clip

        # Print or write list of functions:

        if vmd_filename is None:
            print(vmd_fns)
        else:
            vmd_file = open(vmd_filename, 'w')
            vmd_file.write(vmd_fns)
            vmd_file.close()


class Grid1D():

    def __init__(self, gridsize, gridstep, grid, origin=None, gridaxis='', name=None):
        self.gridsize = gridsize
        self.gridstep = gridstep
        self.grid = grid
        self.origin = origin
        self.gridaxis = gridaxis
        self.name = name

    def __call__(self):
        print('RISM data on 1D grid')

    # Plot 1D profile:
    def plt(self, ax=None):
        # If axes not passed to function, create some:
        show = False
        if ax is None:
            show = True
            ax = plt.subplot(1, 1, 1)
        ax.plot(list(range(self.gridsize)), self.grid)
        ax.set_xlabel('Distance{}(grid points)'.format(" along "+str(self.gridaxis)+" "))
        ax.set_ylabel('RISM value')
        if show:
            plt.show()
        else:
            return ax

    # Draw line to corresponding to 1D profile in VMD:
    def vmd_line(self, vmd_filename=None):
        if self.origin is None:
            sys.exit('ERROR: No origin given')
        ori = self.origin.copy()
        line_end = self.origin.copy()
        line_end[coord_inds[self.gridaxis]] += self.gridstep*self.gridsize
        vmd_line = 'proc draw_rism_line {} {\n' \
                   'set sel [atomselect top all]\n' \
                   'draw color blue\n' \
                   'draw line "%.3f %.3f %.3f" "%.3f %.3f %.3f" style solid width 1\n' % (ori[0], ori[1], ori[2], line_end[0], line_end[1], line_end[2])
        if vmd_filename is None:
            print(vmd_line)
        else:
            vmd_file = open(vmd_filename, 'w')
            vmd_file.write(vmd_line)
            vmd_file.close()
