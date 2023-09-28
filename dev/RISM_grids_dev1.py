import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

# Reads dx file from 3D RISM, assigns values to i,j,k grid points and plots heatmap of slice through the box.

class Grid3D():

    def __init__(self, gridsize, origin, grid):
        self.gridsize = gridsize
        self.origin = origin
        self.grid = grid

    def __call__(self):
        print('Data on 3D grid')

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

        return cls(gridsize, origin, grid)

    # Could use miller indices to get diagonal slices?
    def slice(self, axis, level):
        if axis == 'i':
            #return self.grid[level,:,:]
            return Grid2D(self.gridsize[[1, 2]], np.array([level*self.gridsize[0], 0, 0]) + self.origin, self.grid[level,:,:])
        elif axis == 'j':
            #return self.grid[:,level,:]
            return Grid2D(self.gridsize[[0, 2]], np.array([0, level*self.gridsize[1], 0]) + self.origin, self.grid[:,level,:])
        elif axis == 'k':
            #return self.grid[:,:,level]
            return Grid2D(self.gridsize[[0, 1]], np.array([0, 0, level*self.gridsize[2]]) + self.origin, self.grid[:,:,level])
        else:
            raise ValueError('axis must be "i" "j" or "k"')

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
        return Grid3D(coarse_gridsize, self.origin, coarse_grid)
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
        return Grid3D(self.gridsize, self.origin, pruned_grid)

class Grid2D():

    def __init__(self, gridsize, origin, grid):
        self.gridsize = gridsize
        self.origin = origin
        self.grid = grid

    def __call__(self):
        print('RISM on 2D grid')

    def slice(self, axis, level):
        # NEED TO ALSO OUTPUT THE ORIGIN
        if axis == 'x':
            return Grid1D(self.gridsize[1], self.grid[level,:])
        elif axis == 'y':
            return Grid1D(self.gridsize[0], self.grid[:,level])
        else:
            raise ValueError('axis must be "x" or "y"')

    def plt(self):
        sns.heatmap(self.grid)
        plt.gca().set_aspect('equal')
        plt.show()

class Grid1D():

    def __init__(self, gridsize, grid):
        self.gridsize = gridsize
        self.grid = grid

    def __call__(self):
        print('Data on 1D grid')

    def plt(self):
        plt.plot(list(range(self.gridsize)), self.grid)
        plt.xlabel('Distance (grid points)')
        plt.ylabel('h')
        plt.show()

#    # Remove values outside a certain radius from the central void left by the molecule:
#    def prune(self, r):
#        for d in len(self.gridsize):
#            
#        for x in self.gridsize[0]:
#            np.argwhere(self.grid == 0.0)
#        self.grid

#
#print('k with most zero values:')
#
#from collections import Counter
#
#occurence_count = Counter([g[2] for g in void])
#print(occurence_count.most_common(1))
#
#from statistics import mode
#print(mode([g[2] for g in void]))
#
## Output slice for certain k:
#print(rism_slice)
#
#sns.heatmap(rism_slice)
#plt.axes().set_aspect('equal')
#plt.show()
#
## save file for vmd to plot box and slice
#
## box limits:
#lims = [[ori[0], ori[0] + incr_grid[0]*n_x], [ori[1], ori[1] + incr_grid[1]*n_y], [ori[2], ori[2] + incr_grid[2]*n_z]]
#
#vmd_out = open('vmd_lines.vmd', 'w')
#
## Draw box:
#vmd_out.write('proc vmd_draw_box {mol} {\n')
#vmd_out.write('set sel [atomselect top all]\ndraw color yellow\n')
#for x in lims[0]:
#   for y in lims[1]:
#       vmd_out.write('draw line "%.3f %.3f %.3f" "%.3f %.3f %.3f" style solid width 1\n' % (x, y, lims[2][0], x, y, lims[2][1]))
#   for z in lims[2]:
#       vmd_out.write('draw line "%.3f %.3f %.3f" "%.3f %.3f %.3f" style solid width 1\n' % (x, lims[1][0], z, x, lims[1][1], z))
#for y in lims[1]:
#   for z in lims[2]:
#       vmd_out.write('draw line "%.3f %.3f %.3f" "%.3f %.3f %.3f" style solid width 1\n' % (lims[0][0], y, z, lims[0][1], y, z))
#vmd_out.write('}\n')
#
## Draw slice:
## Use Trinorm to define normals to edges parallel to the surface to prevent the surface diasppearing when looked from the side:
##https://www.ks.uiuc.edu/Research/vmd/vmd-1.7.1/ug/node171.html
#y_coord = ori[1] + 42*incr_grid[1]
#vmd_out.write('proc vmd_draw_slice {mol} {\n')
#vmd_out.write('set sel [atomselect top all]\ndraw color yellow\n')
#vmd_out.write('draw triangle "%.3f %.3f %.3f" "%.3f %.3f %.3f" "%.3f %.3f %.3f"\n' % (lims[0][0], y_coord, lims[2][0], lims[0][1], y_coord, lims[2][0], lims[0][0], y_coord, lims[2][1]))
#vmd_out.write('draw triangle "%.3f %.3f %.3f" "%.3f %.3f %.3f" "%.3f %.3f %.3f"\n' % (lims[0][1], y_coord, lims[2][1], lims[0][0], y_coord, lims[2][1], lims[0][1], y_coord, lims[2][0]))
#vmd_out.write('}\n')
#
