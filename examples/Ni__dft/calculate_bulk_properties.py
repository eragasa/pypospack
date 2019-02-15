class PypospackTask(object):
    
structures = {}
structures['Ni_fcc_cubic'] = [['Ni'],'fcc', 3.508,'cubic']
structures['Ni_bcc_cubic'] = [['Ni'],'bcc', 3.508,'cubic']
structures['Ni_hcp_cubic'] = [['Ni'],'hcp', 3.508,'cubic']
structures['Ni_dia_cubic'] = [['Ni'],'diamond', 3.508,'cubic']
structures['Ni_sc_cubic']  = [['Ni'],'sc', 3.508,'cubic']

for k,v in structures.items():
    print(k,v)

