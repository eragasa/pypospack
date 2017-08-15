import pypospack.io.vasp as vasp

# an array of symbols
symbols = ['Mg','O']
xc = 'GGA'
potcar = vasp.Potcar(symbols=symbols,xc='GGA')
potcar.write('POTCAR')

symbols = ['Ni']
xc = 'GGA'
potcar = vasp.Potcar(symbols=symbols,xc=xc)
potcar.write('POTCAR')

