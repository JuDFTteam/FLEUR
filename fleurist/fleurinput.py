def read_inpxml(file):
    input={
        "jspins":1,
        "noco":False,
        "nkpt":0,
        "basis":0
    }
    from lxml import etree
    with open(file,'r') as inpxmlfile:
        tree=etree.parse(inpxmlfile)
        tree.xinclude()

    try:
        input["jspins"]=int(tree.find(".//magnetism").attrib["jspins"])
        input["noco"]=tree.find(".//magnetism").attrib["l_noco"]=='t' or tree.find(".//magnetism").attrib["l_noco"]=='T'
    except:
        pass

    #find number of k-points
    kpointset=tree.find(".//kPointListSelection").attrib['listName']
    klists=tree.findall('.//kPointList')
    for kl in klists:
        if kl.attrib['name']==kpointset:
            input['nkpt']=int(kl.attrib['count'])



    #Read bravais matrix
    import re
    import numpy as np
    mat=tree.find(".//bravaisMatrix")
    matrix=np.zeros((3,3),np.float)
    for i in range(3):
        s=mat[i].text
        m=re.search(" *([\.\d]+) +([\.\d]+) +([\.\d]+)",s)
        matrix[i,:]=[float(m.groups()[0]),float(m.groups()[1]),float(m.groups()[2])]
    #invert matrix and calculate eigenvalues
    eigvals=np.linalg.eigvals(np.linalg.inv(matrix))

    
    #approximate matrix size
    kmax=float(tree.find(".//cutoffs").attrib["Kmax"])
    input['basis']=int(kmax**3/eigvals[0]/eigvals[1]/eigvals[2]/250)  #250 for 2Pi factor
    if input['noco']: input['basis']=2*input['basis']

    return input