def Tetrahedron(element):
    with open('{}_tetrahedron.xyz'.format(element),'w') as f:
        print("""4

{}  -1.97198893029705      0.52869350703018     -0.01999480282656
{}  -1.88609335015628      2.89785563439104      0.21239667365365
{}  -0.64863258875857      1.82514619388775     -1.51740948337694
{}  -3.02558513078810      1.92320466469102     -1.63849238745014""".format(element,element,element,element),file=f)
    return 
