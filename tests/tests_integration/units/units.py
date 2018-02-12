
def convert(src_x,dst_unit):
    joules_per_erg = 1e-7
    joules_per_eV = 1.602176620e-19

    if src_x[1] == 'erg':
        x = src_x[0] * joules_per_erg
    elif src_x[1] == 'eV':
        x = src_x[0] * joules_per_eV

    if dst_unit == 'eV':
        dst_x = x / joules_per_eV
    elif dst_unit == 'erg':
        dst_x = x / joules_per_erg

    m_per_cm = 1e-2
    m_per_Ang = 1e-10
    if src_x[1] == 'cm':
        x = src_x[0] * m_per_cm
    
    if dst_unit == 'Ang':
        dst_x = x / m_per_Ang
    
    return [dst_x,dst_unit]

x0 = [0.35099e-12,'erg']
xf = convert(x0,'eV')
print(x0,xf)

a = [3.25238e-8,'cm']
print('a',a,convert(a,'Ang'))

r0 = [0.71727*a[0],'cm']
print('r0',r0,convert(r0,'Ang'))

alpha_0 = 8.7660/a[0]
alpha_f = 1/convert([1/alpha_0,'cm'],'Ang')[0]
print('alpha',alpha_0,alpha_f)
