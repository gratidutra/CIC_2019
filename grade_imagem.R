data(plethodon) 
plethodon$land
Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
ref<-mshape(Y.gpa$coords)
plotRefToTarget(ref,Y.gpa$coords[,,39])
plotRefToTarget(ref,Y.gpa$coords[,,39],mag=2,outline=plethodon$outline)   #magnify by 2X
plotRefToTarget(ref,Y.gpa$coords[,,39],method="vector",mag=3)
plotRefToTarget(ref,Y.gpa$coords[,,39],method="points",outline=plethodon$outline)
plotRefToTarget(ref,Y.gpa$coords[,,39],gridPars=gridPar(pt.bg = "green", pt.size = 1),
                method="vector",mag=3)
## tri x tri
z.gpa<-gpagen(t_t)
fer<-mshape(z.gpa$coords)
plotRefToTarget(fer,z.gpa$coords[,,89])
plotRefToTarget(fer,z.gpa$coords[,,89],mag=5)
plotRefToTarget(fer,z.gpa$coords[,,89],method="vector",mag=3)
plotRefToTarget(fer,z.gpa$coords[,,89],method="points")
plotRefToTarget(fer,z.gpa$coords[,,89],gridPars=gridPar(pt.bg = "green", 
                pt.size = 1),method="vector",mag=3)

