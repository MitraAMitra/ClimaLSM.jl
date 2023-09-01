# Domain setup
global zmin
global zmax
global nelements
global h_stem
global h_leaf
println("depth of soil domain is  "*string(zmin))
println("dx of soil domain is  "*string(round((zmax-zmin)./nelements,digits=3)))
land_domain =
    LSMSingleColumnDomain(; zlim = (zmin, zmax), nelements = nelements)

compartment_midpoints = [h_stem / 2, h_stem + h_leaf / 2];
compartment_surfaces = [zmax, h_stem, h_stem + h_leaf];
