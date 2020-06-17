"""
Create fgmax_grid.txt input file
"""

from clawpack.geoclaw import fgmax_tools

def make_fgmax_grid1():
    fg = fgmax_tools.FGmaxGrid()
    fg.fgno = 1
    fg.point_style = 2       # will specify a 2d grid of points
    fg.x1 = 203.515
    fg.x2 = 203.545
    fg.y1 = 20.885
    fg.y2 = 20.905
    fg.dx = 1./(3600.)    # 1 arcssecond grid
    fg.tstart_max = 7.9*3600.     # when to start monitoring max values
    fg.tend_max = 1.e10       # when to stop monitoring max values
    fg.dt_check = 30.         # target time (sec) increment between updating 
                               # max values
    fg.min_level_check = 6    # which levels to monitor max on
    fg.arrival_tol = 1.e-2    # tolerance for flagging arrival

    # fg.input_file_name = 'fgmax1.txt'
    fid = open('fgmax1.txt','w')
    fg.write_to_fgmax_data(fid)


if __name__ == "__main__":
    make_fgmax_grid1()


