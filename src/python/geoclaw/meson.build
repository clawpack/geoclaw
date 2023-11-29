package = 'clawpack.geoclaw'
pkg_dir = join_paths(package.split('.'))

python_sources = {
  '': [
    '__init__.py',
    'data.py',
    'dtopotools.py',
    'etopotools.py',
    'fgmax_tools.py',
    'fgout_tools.py',
    'geoplot.py',
    'kmltools.py',
    'marching_front.py',
    'most2geoclaw.py',
    'nonuniform_grid_tools.py',
    'okada.py',
    'plotfg.py',
    'resolution.py',
    'test.py',
    'topotools.py',
    'units.py',
    'util.py',
  ],
  'datatools': [
    '__init__.py',
    'fixdata.py',
    'iotools.py',
  ],
  'multilayer': [
    '__init__.py',
    'plot.py',
  ],
  'surge': [
    '__init__.py',
    'plot.py',
    'quality.py',
    'storm.py',
  ],
}

foreach subpkg, subsrcs: python_sources
  sources = []
  foreach src: subsrcs
    sources += subpkg / src
  endforeach
  py.install_sources(
    sources,
    subdir: pkg_dir / subpkg,
  )
endforeach
