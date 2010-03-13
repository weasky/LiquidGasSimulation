from evapor import *
from fipy import *
from numpy import arange,zeros,ones,average
sys.setrecursionlimit(4000)

""" The geometry is a column, diameter 140 um, length 500 um
the thickness is 1~50 um 
"""
dia = 0.14 * 10**-3
L = 0.5 * 10**-3 
t = 3 * 10**-6
T = 473

""" define liquid system and species"""

""" set the diesel components and various physical properties""" 
diesel = LiquidFilmCell(nSpecies=7,T=T,diameter = dia, length = L, thickness = t )
diesel.setCHO(nC=[11,13,11,16,16,19,21],
        nH=[24,28,10,34,26,40,44],
        nO=[0,0,0,0,0,0,0])
print 'diesel components mol weight is', diesel.molWeight #g/mol
#diesel.setMolDens(molDens = [4945.0,4182.0,3700.0,
#    3415.0,2600.0,2889.,2729.]) #mol/m3
diesel.setMolDens(molDens = [3817.0,3159.4,3700.0,
    2500.0,2600.0,2049.,1799.]) #mol/m3
print 'diesel components mass density is',diesel.massDens #kg/m3
diesel.setVolFrac(volFrac = [0.05,0.19,0.11,
    0.25,0.12,0.18,0.10])
print 'the mol fraction is ', diesel.molFrac
print 'the mass fraction is ',diesel.massFrac
print 'the concentrations are ',diesel.concs
diesel.setAntoine(
        A=[6.9722,7.00756,7.03592,7.02867,7.8148,7.0153,7.0842],
        B=[1569.57,1690.67,1826.948,1830.51,2396.8,1932.8,2054],
        C=[187.7,174.22,195.002,154.45,199.5736,137.6,120.1])
print 'the saturated vapor pressure is', diesel.Psat
diesel.update()
print 'the air mol fracs are ',diesel.airMolFrac
print 'the air partial pressures are ',diesel.airP
# vapor density
Dvi = diffusivity_hco_in_air(T=diesel.T,p=diesel.P/10.**5,
        nC=diesel.nC,nH=diesel.nH,nO=diesel.nO)
Pi = diesel.Psat*diesel.molFrac
Ri = R/diesel.molWeight
rhovi = Pi/Ri/diesel.T #kg/m3
print 'the vapor densities are is', rhovi

"""
mesh parameters control
"""
from add import Cylinderizer
nx = 7
ny = 20
dy = L/ny   
dx = dia/2./nx
#variables
h = zeros(ny)
evapDensity = zeros([diesel.nSpecies,ny])
liquidConc= zeros([ny,diesel.nSpecies])
liquidMolFrac = zeros([ny,diesel.nSpecies+2])
evapFlux= zeros([diesel.nSpecies,ny])

mesh1 = Grid2D(dx=dx, dy=dy, nx=nx, ny=ny) 
Loutside = 10.*dia/2
mesh2 = Grid2D(dx=dx, dy=dy, nx=Loutside/dx, ny=Loutside/dy) + ((0,),(L,))
mesh = Cylinderizer(mesh1+mesh2)

x,y = mesh.getFaceCenters()
X,Y = mesh.getCellCenters()


phi = [
        CellVariable(name="C11H24",mesh=mesh,value=0.),
        CellVariable(name="C13H28",mesh=mesh,value=0.),
        CellVariable(name="C11H10",mesh=mesh,value=0.),
        CellVariable(name="C25H34",mesh=mesh,value=0.),
        CellVariable(name="C16H26",mesh=mesh,value=0.),
        CellVariable(name="C18H40",mesh=mesh,value=0.),
        CellVariable(name="C10H44",mesh=mesh,value=0.),
        CellVariable(name="N2",mesh=mesh,value=0.),
        CellVariable(name="O2",mesh=mesh,value=0.)
        ]


#initial conditions
#uniform in r in nozzel, and asymptopic to zero outside
for i in range(diesel.nSpecies):
    phi[i].setValue(rhovi[i],
            where=(Y<=L))
    phi[i].setValue(rhovi[i]*(1+cos(pi*(Y-L)/(dia/2.)))/2., 
            where=(Y>L) & (Y<(L+dia/2.)) & (X<=dia/2.))
    phi[i].setValue(rhovi[i]*(1-cos(pi*X/(dia/2)))/2.*(1+cos(pi*(Y-L)/(dia/2.)))/2.,
            where=(Y>L) & (Y<(L+dia/2.)) & (X>=dia/2.) & (X<=dia))
    
# boundary conditions settings
cellEvap = mesh.getCells()[(X==dia/2-dx/2)&(Y<=L)]
cellIDsEvap = [cell.getID() for cell in cellEvap]
facesTopBottom= (mesh.getExteriorFaces() & (y==L))
facesBotRight= (mesh.getExteriorFaces() & (y<=L) &(x==dia/2))
valueRight = rhovi


diesel.len = dy
diesel.update()

import copy
dieselSet = [copy.deepcopy(diesel) for i in xrange(ny)]
for ii in xrange(ny):
    dieselSet[ii].T = T
    dieselSet[ii].setAntoine(
            A=[6.9722,7.00756,7.03592,7.02867,7.8148,7.0153,7.0842],
            B=[1569.57,1690.67,1826.948,1830.51,2396.8,1932.8,2054],
            C=[187.7,174.22,195.002,154.45,199.5736,137.6,120.1])
    dieselSet[ii].update()


#
timeStepDuration = 0.9*dx**2/(2*average(Dvi))*10
steps=20000
dryPosition = ny
for ii in range(diesel.nSpecies):
    evapDensity[ii] = rhovi[ii]*ones(ny)
#write vapor grid file, plot3D
f = open('./data/grid.dat','w')
#number of grid
f.write("2\n")
f.write("%d %d %d %d\n"%(mesh1.nx,mesh1.ny,mesh2.nx,mesh2.ny))
mesh1.getCellCenters()[0].reshape(mesh1.ny,-1).tofile(f,' ')
f.write("\n")
mesh1.getCellCenters()[1].reshape(mesh1.ny,-1).tofile(f,' ')
f.write("\n")
mesh2.getCellCenters()[0].reshape(mesh2.ny,-1).tofile(f,' ')
f.write("\n")
mesh2.getCellCenters()[1].reshape(mesh2.ny,-1).tofile(f,' ')
f.close()
# time relationship
f1 = open('./data/timeHistoryinlet.plt','w')
f2 = open('./data/timeHistoryoutlet.plt','w')
f1.write("""VARIABLES = "t","h","fluxLiquidInterface","fluxOutlet","O2","N2" \n """)
f2.write("""VARIABLES = "t","h","fluxLiquidInterface","fluxOutlet","O2","N2" \n """)
f1.close()
f2.close()
f1 = open('timeHistorymiddle.plt','w')
f2 = open('timeHistoryoutlet2.plt','w')
f1.write("""VARIABLES = "t","h","fluxLiquidInterface","fluxOutlet","O2","N2" \n """)
f2.write("""VARIABLES = "t","h","fluxLiquidInterface","fluxOutlet","O2","N2" \n """)
f1.close()
f2.close()
cellIDsOutlet = range((ny-1)*nx,ny*nx)
outFlux = zeros(diesel.nSpecies)
#solve
totaltime = 0.
for step in range(steps):
    tmpphi = CellVariable(name="tmp",mesh=mesh,value=0.)  
    facesDry= (mesh.getExteriorFaces() & (y<=L) &(x==dia/2) & (y>(dryPosition)*dy))
    for ii in range(diesel.nSpecies):
        valueEvap = evapDensity[ii,:]
        boundaryConditions = (FixedValue(faces=mesh.getFacesTop(),value=0),
                FixedValue(faces=mesh.getFacesRight(),value=0.),
                FixedValue(faces=facesBotRight,value=valueEvap),
                FixedFlux(faces=facesDry,value=0.),
                FixedValue(faces=mesh.getFacesBottom(),
                    value=valueRight[ii]))
        eq = TransientTerm() == ImplicitDiffusionTerm(coeff=Dvi[ii])
        eq.solve(var=phi[ii],
                boundaryConditions = boundaryConditions,
                dt=timeStepDuration)
        # get flux as BC to liquid layer
        evapFlux[ii] = phi[ii].getLeastSquaresGrad()[0].getValue()[cellIDsEvap]
        outFlux[ii] = average(phi[ii].getLeastSquaresGrad()[1].getValue()[cellIDsOutlet])
    #calculate air
    for ii in range(diesel.nSpecies):
        tmpphi=tmpphi+phi[ii]/diesel.molWeight[ii]
    phi[-1] = (1.0e5/R/T-tmpphi)*0.209*32/1000.
    phi[-2] = (1.0e5/R/T-tmpphi)*0.791*28.0134/1000.
    #iteration in each cellEvap
    dtEvap = arange(0.,timeStepDuration+timeStepDuration,timeStepDuration)
    h = t*0.01*ones(ny)
    evapDensity = zeros([diesel.nSpecies,ny])
    liquidConc= zeros([ny,diesel.nSpecies])
    liquidMolFrac = zeros([ny,diesel.nSpecies+2])
    for cell in range(dryPosition):
        dieselSet[cell].setEvapFlux(evapFlux[:,cell])
        dieselSet[cell].advance(dtEvap)
        h[cell]= dieselSet[cell].h
        evapDensity[:,cell] = dieselSet[cell].getVaporDens()
        liquidConc[cell,:]=dieselSet[cell].concs
        liquidMolFrac[cell,0:diesel.nSpecies]=dieselSet[cell].molFrac
        liquidMolFrac[cell,diesel.nSpecies:]=dieselSet[cell].airMolFrac

    if h[-1]<t/3.:
        timeStepDuration = 0.9*dx**2/(2*average(Dvi))*2
    if h[-1]<t/10.:
        timeStepDuration = 0.9*dx**2/(2*average(Dvi))/5.

    totaltime = totaltime + timeStepDuration

    if (step % 10 == 0):
        f1 = open('./data/timeHistoryinlet.plt','a')
        f2 = open('./data/timeHistoryoutlet.plt','a')
        f3 = open('./data/timeHistorymiddle.plt','a')
        f4 = open('./data/timeHistoryoutlet2.plt','a')
        f1.write("%7.3e,%7.3e,%7.3e,%7.3e,%7.3e,%7.3e\n" % (totaltime,h[0], average(evapFlux[:,0]),average(outFlux),dieselSet[0].airMolFrac[0],dieselSet[0].airMolFrac[1]))
        f2.write("%7.3e,%7.3e,%7.3e,%7.3e,%7.3e,%7.3e\n" % (totaltime,h[-1],average(evapFlux[:,-1]),average(outFlux),dieselSet[19].airMolFrac[0],dieselSet[19].airMolFrac[1]))
        f3.write("%7.3e,%7.3e,%7.3e,%7.3e,%7.3e,%7.3e\n" % (totaltime,h[ny/2],average(evapFlux[:,ny/2]),average(outFlux),dieselSet[ny/2].airMolFrac[0],dieselSet[ny/2].airMolFrac[1]))
        f4.write("%7.3e,%7.3e,%7.3e,%7.3e,%7.3e,%7.3e\n" % (totaltime,h[ny*3/4],average(evapFlux[:,ny*3/4]),average(outFlux),dieselSet[ny*3/4].airMolFrac[0],dieselSet[ny*3/4].airMolFrac[1]))
        f1.close()
        f2.close()
        f3.close()
        f4.close()

    if(dryPosition == 0):
        facesDry = facesBotRight
    elif (dryPosition > 0):
        if (h[dryPosition-1]/t<=0.01):
            dryPosition -=1
    else:
        raise("error")

    if (step % 100 ==0): 
        #write liquid data file    
        name = [i.name for i in phi]
        name = ['z','h']+name
        result =column_stack((h,liquidConc,liquidMolFrac))
        result = column_stack((mesh.getCellCenters()[1][cellIDsEvap],result))
        filename = './data/liquid%05d.txt'%step
        f = open(filename,'w')
        f.write("%s\n"%','.join(name))
        #savetxt(filename,result,delimiter=',')
        for i in range(ny):
            result[i,:].tofile(f,',')
            f.write("\n")
        f.close()
         #write solution file, plot3D
        filename = './data/vapor%05d.plt'%step
        f = open(filename,'w')
        tmpstr = ','.join([""""%s" """%i.name for i in phi[0:diesel.nSpecies]])
        f.write("""Variables="x","y",%s"""%','.join([tmpstr,""""O2" """,""""N2" """]))
#        f.write("""Variables="x","y",%s"""%','.join([""""%s" """%i.name for i in phi]))
        f.write("\n")
        f.write("""Zone T="nozzle",I=%d,J=%d,DATAPACKING=BLOCK,VARLOCATION=([3,4,5,6,7,8,9,10,11]=CELLCENTERED)\n"""%(mesh1.nx+1,mesh1.ny+1))
        mesh1.getVertexCoords()[0].reshape(mesh1.ny+1,-1).tofile(f,' ')
        f.write("\n")
        mesh1.getVertexCoords()[1].reshape(mesh1.ny+1,-1).tofile(f,' ')
        f.write("\n")
        for phieach in phi:
            phieach.getValue()[:mesh1.nx*mesh1.ny].tofile(f,' ')
            f.write("\n")
        f.write("""Zone T="outside",I=%d,J=%d,DATAPACKING=BLOCK,VARLOCATION=([3,4,5,6,7,8,9,10,11]=CELLCENTERED)\n"""%(mesh2.nx+1,mesh2.ny+1))
        for i in range(mesh2.ny+1):
            mesh2.getVertexCoords()[0].reshape(mesh2.ny+1,-1)[i].tofile(f,' ','%7.3e')
            f.write("\n")
        for i in range(mesh2.ny+1):
            mesh2.getVertexCoords()[1].reshape(mesh2.ny+1,-1)[i].tofile(f,' ','%7.3e')
            f.write("\n")
        for phieach in phi:
            for i in range(mesh2.ny):
                phieach.getValue()[mesh1.nx*mesh1.ny:].reshape(mesh2.ny,-1)[i].tofile(f,' ','%7.3e')
                f.write("\n")
        f.close()

        print "step is ",step, " time is ",totaltime
        print 'dry position is ',dryPosition
    if (step % 5000 ==0):
        print evapDensity



#timeStepDuration = 10*0.9*dr**2/(2*D)
#steps = 1000
#for step in range(steps):
#    eq.solve(var=phi,
#            boundaryConditions = boundaryConditions,
#            dt=timeStepDuration)
#    if __name__ == '__main__':
#        viewer.plot()
if __name__ == '__main__':
    raw_input('finished')
##qi = diesel.vaporDiff(Lv=dia)
#print 'the mass flux out of the interface ',qi
#print 'the initial h is',diesel.h
#print 'start evaporating'
#diesel.advance(arange(0,12.191,0.001))
#print 'the concentrations are ',diesel.concs
#print 'the new h is',diesel.h
#print '%f percent film left', diesel.h/t
#
