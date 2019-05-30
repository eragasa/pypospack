import pypospack.io.vasp as vasp

    
if __name__ == '__main__':
    symbols = ['Mg','O']
    filename = 'POTCAR'

    # let's try to read
    potcar = vasp.Potcar(symbols=symbols)
    potcar.write('POTCAR')

    # now let's try to write
    potcar = vasp.Potcar(symbols=symbols)
    potcar.write('POTCAR')
    potcar.read('POTCAR')
    print('symbols',potcar.symbols)
    print('encut_min',potcar.encut_min)
    print('encut_max',potcar.encut_max)
    print('lexch',potcar.xc)
    print('models',potcar.models)
